from __future__ import annotations

import hashlib
import json
import os
import re
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Mapping, Sequence

import numpy as np

from lineage_charisma import condensed as cd
from lineage_charisma.distances import cosine_distance_matrix
from lineage_charisma.embed import EmbeddingResult, chunk_spans, l2_normalize, pool_passages
from lineage_charisma.h2 import decision, share_of_raw
from lineage_charisma.io_utils import atomic_write_text
from lineage_charisma.masking import PILOT_NAME_RANKS, Masker, build_vocabulary
from lineage_charisma.wiki import Section, clean_text, split_sections

# Every rule below is fixed in docs/preregistration_amendment_extensions.md (Amendment 5), designed after the 2026-11-11 lock.
AMENDMENT = "docs/preregistration_amendment_extensions.md"
STATUS = "secondary (Amendment 5, designed after the 2026-11-11 lock)"
MIN_TOKENS = 50            # the registered stub threshold, applied to each section group
MIN_GROUP_SPECIES = 50
COST_STOP_USD = 10.0


def amendment_is_committed(root: Path, rel: str = AMENDMENT) -> tuple[bool, str]:
    """The amendment must be tracked with no uncommitted changes before any statistic or request (section 4)."""
    if not (Path(root) / rel).exists():
        return False, f"{rel} does not exist"
    tracked = subprocess.run(["git", "ls-files", "--error-unmatch", rel], cwd=root, capture_output=True).returncode == 0
    clean = subprocess.run(["git", "diff", "--quiet", "HEAD", "--", rel], cwd=root, capture_output=True).returncode == 0
    if not (tracked and clean):
        return False, f"{rel} must be committed with no uncommitted changes before any extension analysis runs"
    return True, ""


def scale_masker(cfg) -> Masker:
    """The masker scripts/05b_build_corpus.py builds for this scale (pilot or full-scale vocabulary)."""
    from lineage_charisma.io_utils import read_csv

    m = cfg.get("masking") or {}
    vocab = build_vocabulary(read_csv(cfg.mdd_species_path()), read_csv(cfg.mdd_synonyms_path()), extra_terms=m.get("extra_taxonomy_terms", []),
                             keep_terms=m.get("keep_terms", []), keep_when_followed_by=m.get("keep_when_followed_by", {}),
                             head_noun_rule=m.get("head_noun_rule", "pilot"), epithets_anywhere=m.get("epithets_anywhere", True),
                             name_ranks=m.get("name_ranks", PILOT_NAME_RANKS))
    return Masker(vocab, m.get("token", "[TAXON]"))


# ---------------------------------------------------------------- registered tests on a new text matrix (section 3)

@dataclass(frozen=True)
class Registered:
    """The fixed matrices every extension is tested against at one scale."""
    order: list[str]           # H1 species
    phylo: np.ndarray          # on `order`
    h2_species: list[str]      # subset of `order`, same order
    ecology: np.ndarray        # on `h2_species`
    length: np.ndarray         # on `h2_species`

    def restrict(self, keep: Sequence[str]) -> "Registered":
        """The same matrices on the species in `keep` (used when a representation has to drop a species)."""
        keep = set(keep)
        order, h2 = [s for s in self.order if s in keep], [s for s in self.h2_species if s in keep]
        return Registered(order, cd.subset(self.phylo, self.order, order), h2, cd.subset(self.ecology, self.h2_species, h2),
                          cd.subset(self.length, self.h2_species, h2))


def tables_dir(cfg) -> Path:
    p = cfg.reports / "tables" / "extensions" / (cfg.output_subdir or "")
    p.mkdir(parents=True, exist_ok=True)
    return p


def work_dir(cfg) -> Path:
    p = cfg.processed / "extensions"
    p.mkdir(parents=True, exist_ok=True)
    return p


def load_registered(cfg) -> Registered:
    from lineage_charisma.h2 import load_h2_matrices
    from lineage_charisma.distances import load_matrix, phylo_matrix_path

    h2, species, order = load_h2_matrices(cfg)
    phylo = load_matrix(phylo_matrix_path(cfg.matrices_dir()), cfg.species_order_path())[0]
    return Registered(order, phylo, species, h2["ecology"], h2["length"])


def h1_test(text: np.ndarray, phylo: np.ndarray, permutations: int, seed: int, method: str = "spearman", workers: int | None = None) -> dict:
    res = cd.mantel(cd.prepare(text, method, "text"), cd.prepare(phylo, method, "phylogeny"), permutations, seed, workers=workers)
    return {"test": "H1", "r": res.r, "p": res.p, "n_species": res.n, "n_pairs": res.n_pairs, "permutations": permutations}


def registered_tests(text: np.ndarray, reg: Registered, permutations: int, seed: int, alpha: float, method: str = "spearman",
                     workers: int | None = None) -> list[dict]:
    """H1 on the H1 species; raw and partial (ecology + length) Mantel on the H2 species with rules (a) and (b). Phylogeny moves."""
    h1 = h1_test(text, reg.phylo, permutations, seed, method, workers)
    h1["outcome"] = "positive" if h1["r"] > 0 and h1["p"] < alpha else "not positive"
    text_h2, phylo_h2 = cd.subset(text, reg.order, reg.h2_species), cd.subset(reg.phylo, reg.order, reg.h2_species)
    P, X = cd.prepare(text_h2, method, "text"), cd.prepare(phylo_h2, method, "phylogeny")
    C = [cd.prepare(reg.ecology, method, "ecology"), cd.prepare(reg.length, method, "length")]
    null = cd.null_correlations(X, [P, *C], permutations, seed, workers=workers)
    raw = cd.mantel(P, X, seed=seed, null=null[:, 0])
    part = cd.partial_mantel(X, P, C, seed=seed, null_corr=null)
    rules = decision(raw.r, part.r, part.p, alpha)
    h2 = {"test": "H2 partial (ecology + length)", "raw_r": raw.r, "raw_p": raw.p, "r": part.r, "p": part.p, "share_of_raw": share_of_raw(raw.r, part.r),
          "n_species": part.n, "n_pairs": part.n_pairs, "permutations": permutations,
          "rule_a": "survives" if rules["survives"] else "fails",
          "rule_b": "undefined" if rules["confounded"] is None else ("confounded" if rules["confounded"] else "not confounded")}
    return [h1, h2]


def matrix_agreement(matrices: Mapping[str, np.ndarray], method: str = "spearman") -> list[dict]:
    """Descriptive Spearman correlation between every pair of text-distance matrices on the same species (no p-values)."""
    prepared = {k: cd.prepare(v, method, k) for k, v in matrices.items()}
    names = list(prepared)
    return [{"a": a, "b": b, "r": cd.correlation(prepared[a], prepared[b])} for i, a in enumerate(names) for b in names[i + 1:]]


# ---------------------------------------------------------------- B. TF-IDF baseline

TFIDF_SETTINGS = dict(lowercase=True, token_pattern=r"(?u)\b[^\W\d_]{2,}\b", stop_words="english", ngram_range=(1, 1),
                      min_df=5, max_df=1.0, sublinear_tf=True, smooth_idf=True, use_idf=True, norm="l2")


@dataclass(frozen=True)
class TfidfResult:
    distance: np.ndarray       # on the kept documents
    kept: list[int]            # positions of the input texts with a non-zero vector
    dropped: list[int]
    vocabulary_size: int
    median_terms: float


def tfidf_distance(texts: Sequence[str], mask_token: str = "[TAXON]", **overrides) -> TfidfResult:
    from sklearn.feature_extraction.text import TfidfVectorizer

    docs = [t.replace(mask_token, " ") for t in texts]
    vec = TfidfVectorizer(**{**TFIDF_SETTINGS, **overrides})
    X = vec.fit_transform(docs).tocsr()
    nnz = np.diff(X.indptr)
    kept, dropped = np.flatnonzero(nnz > 0).tolist(), np.flatnonzero(nnz == 0).tolist()
    X = X[kept]
    # rows are L2-normalized by the vectorizer, so the sparse Gram matrix is the cosine similarity
    dist = 1.0 - np.clip((X @ X.T).toarray(), -1.0, 1.0)
    dist = (dist + dist.T) / 2.0
    np.fill_diagonal(dist, 0.0)
    return TfidfResult(np.clip(dist, 0.0, 2.0), kept, dropped, len(vec.vocabulary_), float(np.median(nnz[kept])) if kept else 0.0)


# ---------------------------------------------------------------- C. section groups

GROUPS = ("appearance", "behavior", "distribution_habitat", "diet", "other")
LEAD = "lead"                  # the sixth, reference row: the lead alone
REPORTED = (*GROUPS, LEAD)
EXCLUDED = "excluded"

EXACT_HEADINGS: dict[str, tuple[str, ...]] = {
    EXCLUDED: ("references", "external links", "see also", "further reading", "notes", "literature cited", "sources", "bibliography", "footnotes",
               "citations", "gallery", "notes and references", "references and notes", "cited literature", "works cited"),
    "appearance": ("description", "physical description", "characteristics", "physical characteristics", "appearance", "morphology", "anatomy",
                   "anatomy and morphology", "anatomy and physiology", "size", "description and anatomy", "identification", "coloration", "pelage",
                   "dentition", "skull"),
    "behavior": ("behavior", "ecology", "biology", "behavior and ecology", "ecology and behavior", "biology and ecology", "biology and behavior",
                 "behavior and biology", "life history", "natural history", "habits", "reproduction", "breeding", "mating", "reproduction and development",
                 "reproduction and life cycle", "mating and reproduction", "social behavior", "social structure", "communication", "vocalizations",
                 "echolocation", "hibernation", "locomotion", "senses", "predators", "predation", "lifespan"),
    "distribution_habitat": ("distribution and habitat", "distribution", "habitat", "range and habitat", "habitat and distribution", "range",
                             "population and distribution", "habitat and ecology", "distribution and ecology", "range and status",
                             "distribution and status"),
    "diet": ("diet", "feeding", "food habits", "hunting", "foraging", "diet and foraging", "food and foraging", "hunting and diet", "diet and hunting",
             "feeding ecology", "diet and feeding"),
    "other": ("taxonomy", "etymology", "taxonomy and etymology", "taxonomy and evolution", "taxonomy and phylogeny", "taxonomy and naming", "evolution",
              "phylogeny", "fossil record", "classification", "names", "subspecies", "hybrids", "genetics", "history", "discovery", "extinction",
              "conservation", "conservation status", "status", "threats", "threats and conservation", "status and conservation", "population status",
              "relationship with humans", "relationships with humans", "interaction with humans", "interactions with humans", "human interaction",
              "human interactions", "in culture", "in popular culture", "cultural significance", "in captivity", "parasites", "physiology",
              "attacks on humans", "cultural references"),
}
# checked in this order, so `excluded` wins a tie at the same word
KEYWORDS: dict[str, tuple[str, ...]] = {
    EXCLUDED: ("references", "notes", "sources", "bibliography", "footnotes", "citations", "gallery", "literature", "links"),
    "appearance": ("description", "appearance", "characteristics", "morphology", "anatomy", "size", "coloration", "pelage", "fur", "dentition", "skull",
                   "measurements", "physical", "identification"),
    "behavior": ("behavior", "ecology", "biology", "habits", "social", "communication", "reproduction", "breeding", "mating", "echolocation",
                 "vocalizations", "vocalization", "hibernation", "locomotion", "activity", "life", "lifecycle", "lifespan", "predators", "predation", "senses",
                 "development"),
    "distribution_habitat": ("distribution", "habitat", "habitats", "range", "occurrence", "geographic"),
    "diet": ("diet", "feeding", "food", "foraging", "hunting", "prey", "nutrition"),
}
_EXACT = {h: g for g, hs in EXACT_HEADINGS.items() for h in hs}
_KEYWORD = {}
for _g, _ws in KEYWORDS.items():
    for _w in _ws:
        _KEYWORD.setdefault(_w, _g)


def normalize_heading(heading: str) -> str:
    h = heading.strip().lower().replace("&", "and")
    h = re.sub(r"\bbehaviour", "behavior", h)
    h = re.sub(r"\bcolouration\b", "coloration", h)
    return re.sub(r"\s+", " ", h)


def classify_heading(heading: str) -> tuple[str | None, str]:
    """(group, rule) from the exact table or the keyword rule; (None, "none") when neither applies."""
    h = normalize_heading(heading)
    if h in _EXACT:
        return _EXACT[h], "exact"
    for word in re.findall(r"[a-z]+", h):
        if word in _KEYWORD:
            return _KEYWORD[word], "keyword"
    return None, "none"


def assign_groups(sections: Sequence[Section]) -> list[tuple[str, str]]:
    """(group, rule) for every section after the lead; unmatched subsections inherit from their nearest enclosing section."""
    out: list[tuple[str, str]] = []
    stack: list[tuple[int, str]] = []          # (level, group) of enclosing sections
    for sec in sections[1:]:
        while stack and stack[-1][0] >= sec.level:
            stack.pop()
        group, rule = classify_heading(sec.heading)
        if group is None:
            group, rule = (stack[-1][1], "inherited") if stack and sec.level > 2 else ("other", "default")
        out.append((group, rule))
        stack.append((sec.level, group))
    return out


def group_texts(extract: str) -> dict[str, str]:
    """Cleaned text of each section group in article order, plus the cleaned lead under `lead`; excluded sections are left out."""
    sections = split_sections(extract)
    parts: dict[str, list[str]] = {g: [] for g in REPORTED}
    lead = clean_text(sections[0].body)
    if lead:
        parts[LEAD].append(lead)
    for sec, (group, _) in zip(sections[1:], assign_groups(sections)):
        if group == EXCLUDED:
            continue
        body = clean_text(sec.body)
        if body:
            parts[group].append(body)
    return {g: "\n\n".join(p) for g, p in parts.items() if p}


# ---------------------------------------------------------------- A. commercial embeddings

@dataclass(frozen=True)
class Provider:
    name: str
    model: str
    url: str
    key_env: str
    limit: int                 # passage limit in the provider's tokens
    max_batch: int             # texts per request
    max_batch_tokens: int      # tokens per request
    usd_per_mtok: float        # list price, confirmed before any request (amendment 5, section A)

    def body(self, inputs: list[str]) -> dict:
        return {"model": self.model, "input": inputs, "encoding_format": "float"}


PROVIDERS = {
    "openai": Provider("openai", "text-embedding-3-large", "https://api.openai.com/v1/embeddings", "OPENAI_API_KEY", 8191, 96, 100_000, 0.13),
}   # OpenAI only: Voyage was dropped before any request (amendment 5, Change 1)
PESSIMISTIC = {"retry_factor": 2.0}


class Tokenizer:
    """Encode and decode in one provider's tokens."""

    def __init__(self, encode: Callable[[str], list[int]], decode: Callable[[list[int]], str], name: str):
        self.encode, self.decode, self.name = encode, decode, name


def provider_tokenizer(provider: Provider) -> tuple[Tokenizer, int]:
    """(tokenizer, passage limit): OpenAI's cl100k_base at the model's input limit."""
    import tiktoken

    enc = tiktoken.get_encoding("cl100k_base")
    return Tokenizer(enc.encode, enc.decode, "cl100k_base"), provider.limit


def split_passages(texts: Sequence[str], tokenizer: Tokenizer, limit: int, overlap: int = 32) -> tuple[list[str], list[int], list[int]]:
    """(passages, owner document, token count per passage) under the registered chunk rule at the provider's limit."""
    passages, owner, counts = [], [], []
    for doc, text in enumerate(texts):
        ids = tokenizer.encode(text)
        spans = chunk_spans(len(ids), limit, overlap)
        if not spans:
            raise ValueError(f"document {doc} has no tokens")
        for start, end in spans:
            passages.append(text if len(spans) == 1 else tokenizer.decode(ids[start:end]))
            owner.append(doc)
            counts.append(end - start)
    return passages, owner, counts


class ApiError(RuntimeError):
    pass


def _default_post(url: str, headers: dict, body: dict) -> tuple[int, dict]:
    import requests

    resp = requests.post(url, headers=headers, json=body, timeout=120)
    try:
        payload = resp.json()
    except ValueError:
        payload = {"error": resp.text[:500]}
    return resp.status_code, payload


class ApiEmbedder:
    """Embeds passages through one provider, one cached JSON file per passage. The key is read from the environment at the first
    request and never stored on the object, printed, or written."""

    def __init__(self, provider: Provider, cache_dir: Path, post: Callable[[str, dict, dict], tuple[int, dict]] = _default_post,
                 sleep: Callable[[float], None] = time.sleep, max_retries: int = 5, backoff: float = 2.0):
        self.provider, self.cache_dir = provider, Path(cache_dir) / provider.name
        self.post, self.sleep, self.max_retries, self.backoff = post, sleep, max_retries, backoff
        self.requests_sent, self.tokens_billed = 0, 0

    def cache_path(self, passage: str) -> Path:
        digest = hashlib.sha256(f"{self.provider.name}|{self.provider.model}|{passage}".encode()).hexdigest()
        return self.cache_dir / digest[:2] / f"{digest}.json"

    def _headers(self) -> dict:
        key = os.environ.get(self.provider.key_env)
        if not key:
            raise ApiError(f"{self.provider.key_env} is not set")
        return {"Authorization": f"Bearer {key}", "Content-Type": "application/json"}

    def _request(self, inputs: list[str]) -> dict:
        for attempt in range(self.max_retries + 1):
            status, payload = self.post(self.provider.url, self._headers(), self.provider.body(inputs))
            self.requests_sent += 1
            if status == 200:
                return payload
            if status != 429 and status < 500:
                detail = re.sub(r"sk-[\w*-]+", "[key redacted]", str(payload.get("error", payload)))   # providers echo a masked key
                raise ApiError(f"{self.provider.name} returned HTTP {status}: {detail[:300]}")
            if attempt < self.max_retries:
                self.sleep(self.backoff * 2 ** attempt)
        raise ApiError(f"{self.provider.name}: no valid response after {self.max_retries} retries (last HTTP {status})")

    def embed(self, passages: Sequence[str], token_counts: Sequence[int]) -> np.ndarray:
        todo = [i for i, p in enumerate(passages) if not self.cache_path(p).exists()]
        batch: list[int] = []
        for i in [*todo, None]:
            full = i is not None and batch and (len(batch) >= self.provider.max_batch or
                                               sum(token_counts[j] for j in batch) + token_counts[i] > self.provider.max_batch_tokens)
            if batch and (i is None or full):
                self._store([passages[j] for j in batch])
                batch = []
            if i is not None:
                batch.append(i)
        vectors = [json.loads(self.cache_path(p).read_text())["embedding"] for p in passages]
        return np.asarray(vectors, dtype=np.float64)

    def _store(self, inputs: list[str]) -> None:
        payload = self._request(inputs)
        data = sorted(payload.get("data", []), key=lambda d: d["index"])
        if len(data) != len(inputs):
            raise ApiError(f"{self.provider.name}: {len(data)} embeddings returned for {len(inputs)} inputs")
        usage = payload.get("usage", {})
        self.tokens_billed += int(usage.get("total_tokens", usage.get("prompt_tokens", 0)))
        for text, item in zip(inputs, data):
            path = self.cache_path(text)
            path.parent.mkdir(parents=True, exist_ok=True)
            atomic_write_text(path, json.dumps({"provider": self.provider.name, "model_requested": self.provider.model,
                                                "model": payload.get("model"), "usage": usage, "n_inputs": len(inputs),
                                                "embedding": item["embedding"]}))


def embed_with_provider(texts: Sequence[str], embedder: ApiEmbedder, tokenizer: Tokenizer, limit: int, overlap: int = 32) -> EmbeddingResult:
    passages, owner, counts = split_passages(texts, tokenizer, limit, overlap)
    vectors = l2_normalize(embedder.embed(passages, counts))
    pooled = pool_passages(vectors, owner, len(texts))
    n_passages = np.bincount(owner, minlength=len(texts)).tolist()
    totals = [len(tokenizer.encode(t)) for t in texts]
    return EmbeddingResult(pooled.astype(np.float32), n_passages, totals, totals)


def cost_estimate(tokens: Mapping[str, int], provider: Provider) -> list[dict]:
    """Per scale: list-price cost and a pessimistic bound (every request sent twice). The token counts are exact for cl100k_base."""
    rows = []
    for scale, n in tokens.items():
        factor = PESSIMISTIC["retry_factor"]
        rows.append({"provider": provider.name, "model": provider.model, "scale": scale, "tokens": int(n),
                     "usd": n / 1e6 * provider.usd_per_mtok, "usd_pessimistic": n * factor / 1e6 * provider.usd_per_mtok})
    return rows


def text_distance(embeddings: np.ndarray) -> np.ndarray:
    return cosine_distance_matrix(embeddings)
