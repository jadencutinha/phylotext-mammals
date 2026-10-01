from __future__ import annotations

from pathlib import Path
from typing import Callable, Mapping, Sequence

import pandas as pd

from lineage_charisma.io_utils import slugify_species
from lineage_charisma.masking import Masker, text_column

CORPUS_COLUMNS = ["species_id", "family", "text", "n_tokens", "is_stub"]

TokenCounter = Callable[[Sequence[str]], list[int]]


def load_texts(species: Sequence[str], clean_dir: Path) -> dict[str, str]:
    clean_dir = Path(clean_dir)
    texts, missing = {}, []
    for name in species:
        path = clean_dir / f"{slugify_species(name)}.txt"
        if path.exists():
            texts[name] = path.read_text(encoding="utf-8").strip()
        else:
            missing.append(name)
    if missing:
        raise FileNotFoundError(f"{len(missing)} species have no cleaned text in {clean_dir}, e.g. {missing[:5]}")
    return texts


def load_token_counter(tokenizer_name: str) -> TokenCounter:
    from tokenizers import Tokenizer

    tok = Tokenizer.from_pretrained(tokenizer_name)
    tok.no_truncation()
    tok.no_padding()

    def count(texts: Sequence[str]) -> list[int]:
        return [len(enc.ids) for enc in tok.encode_batch(list(texts), add_special_tokens=False)]

    return count


def build_corpus(
    matched: pd.DataFrame,
    texts: Mapping[str, str],
    count_tokens: TokenCounter,
    min_tokens: int,
    masker: Masker | None = None,
    mask_levels: Sequence[str] = (),
) -> pd.DataFrame:
    """One row per species. `text`, `n_tokens` and `is_stub` always describe the unmasked text, so the
    stub set is the same at every mask level; each masked level adds `text_<level>` and `n_masked_<level>`."""
    if matched["species"].duplicated().any():
        dupes = matched.loc[matched["species"].duplicated(), "species"].tolist()
        raise ValueError(f"duplicate species in matched table: {dupes[:5]}")
    missing = [s for s in matched["species"] if s not in texts]
    if missing:
        raise KeyError(f"{len(missing)} matched species have no text, e.g. {missing[:5]}")
    corpus = pd.DataFrame({
        "species_id": matched["species"].astype(str).to_numpy(),
        "family": matched["family"].astype(str).to_numpy(),
    })
    corpus["text"] = [texts[s].strip() for s in corpus["species_id"]]
    corpus["n_tokens"] = pd.array(count_tokens(corpus["text"].tolist()), dtype="int64")
    corpus["is_stub"] = corpus["n_tokens"] < int(min_tokens)
    for level in mask_levels:
        if level == "none":
            continue
        if masker is None:
            raise ValueError(f"mask level {level!r} requested without a masker")
        corpus[text_column(level)] = [masker.mask(t, level) for t in corpus["text"]]
        corpus[f"n_masked_{level}"] = pd.array([masker.count(t) for t in corpus[text_column(level)]], dtype="int64")
    return corpus.sort_values("species_id").reset_index(drop=True)


def stub_table(corpus: pd.DataFrame, extra: pd.DataFrame | None = None) -> pd.DataFrame:
    stubs = corpus.loc[corpus["is_stub"], ["species_id", "family", "n_tokens", "text"]].sort_values("n_tokens")
    if extra is not None:
        stubs = stubs.merge(extra, on="species_id", how="left")
    return stubs.reset_index(drop=True)


def length_summary(corpus: pd.DataFrame, limits: Sequence[int] = (384, 512)) -> dict[str, float | int]:
    n = corpus["n_tokens"]
    out: dict[str, float | int] = {
        "n_species": int(len(corpus)),
        "n_stub": int(corpus["is_stub"].sum()),
        "min": int(n.min()),
        "q1": int(n.quantile(0.25)),
        "median": int(n.median()),
        "mean": round(float(n.mean()), 1),
        "q3": int(n.quantile(0.75)),
        "max": int(n.max()),
    }
    for limit in limits:
        out[f"over_{limit}"] = int((n > limit - 2).sum())
    return out


def write_corpus(corpus: pd.DataFrame, path: Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp")
    extra = [c for c in corpus.columns if c not in CORPUS_COLUMNS]
    corpus[CORPUS_COLUMNS + extra].to_parquet(tmp, index=False)
    tmp.replace(path)
    return path


def read_corpus(path: Path) -> pd.DataFrame:
    corpus = pd.read_parquet(path)
    absent = [c for c in CORPUS_COLUMNS if c not in corpus.columns]
    if absent:
        raise ValueError(f"{path} is missing columns {absent}")
    if corpus["species_id"].duplicated().any():
        raise ValueError(f"{path} has duplicate species_id values")
    return corpus
