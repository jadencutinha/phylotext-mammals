from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Protocol, Sequence

import numpy as np

from lineage_charisma.io_utils import atomic_write_text, is_cached, read_json, sha256, write_json

RULES = ("chunk", "truncate")


def model_slug(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")


def detect_device(preference: str | None = "auto") -> str:
    if preference and preference != "auto":
        return preference
    import torch

    if torch.cuda.is_available():
        return "cuda"
    if torch.backends.mps.is_available():
        return "mps"
    return "cpu"


def chunk_spans(n_tokens: int, max_len: int, overlap: int) -> list[tuple[int, int]]:
    """Token spans [start, end) covering 0..n_tokens, each at most max_len long.

    Consecutive spans share `overlap` tokens. The spans are sized evenly, so a text
    slightly over the limit gives two similar halves instead of a full chunk plus a
    few leftover tokens that would count as much as the full chunk in the mean.
    """
    if max_len <= 0:
        raise ValueError(f"max_len must be positive, got {max_len}")
    if not 0 <= overlap < max_len:
        raise ValueError(f"overlap must be in [0, max_len), got overlap={overlap}, max_len={max_len}")
    if n_tokens <= 0:
        return []
    if n_tokens <= max_len:
        return [(0, n_tokens)]
    k = math.ceil((n_tokens - overlap) / (max_len - overlap))
    size = math.ceil((n_tokens + (k - 1) * overlap) / k)
    step = size - overlap
    return [(i * step, min(i * step + size, n_tokens)) for i in range(k)]


def passage_spans(n_tokens: int, rule: str, max_len: int, overlap: int) -> list[tuple[int, int]]:
    if rule == "chunk":
        return chunk_spans(n_tokens, max_len, overlap)
    if rule == "truncate":
        if max_len <= 0:
            raise ValueError(f"max_len must be positive, got {max_len}")
        return [(0, min(n_tokens, max_len))] if n_tokens > 0 else []
    raise ValueError(f"unknown long-text rule {rule!r}; expected one of {RULES}")


def l2_normalize(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    norms = np.linalg.norm(x, axis=-1, keepdims=True)
    if np.any(norms == 0) or not np.all(np.isfinite(norms)):
        raise ValueError("cannot L2-normalize a zero or non-finite vector")
    return x / norms


def pool_passages(passage_embeddings: np.ndarray, owner: Sequence[int], n_docs: int) -> np.ndarray:
    """Mean of each document's passage embeddings, then L2-normalized."""
    passage_embeddings = np.asarray(passage_embeddings, dtype=np.float64)
    owner = np.asarray(owner, dtype=np.int64)
    if passage_embeddings.shape[0] != owner.shape[0]:
        raise ValueError(f"{passage_embeddings.shape[0]} passage embeddings but {owner.shape[0]} owners")
    counts = np.bincount(owner, minlength=n_docs)
    if len(counts) != n_docs or np.any(counts == 0):
        raise ValueError(f"documents without any passage: {np.flatnonzero(counts[:n_docs] == 0).tolist()[:5]}")
    sums = np.zeros((n_docs, passage_embeddings.shape[1]))
    np.add.at(sums, owner, passage_embeddings)
    return l2_normalize(sums / counts[:, None])


class Encoder(Protocol):
    max_content_tokens: int

    def tokenize(self, text: str) -> list[int]: ...

    def encode_ids(self, passages: Sequence[Sequence[int]], batch_size: int) -> np.ndarray: ...


class SentenceTransformerEncoder:
    def __init__(self, name: str, device: str = "auto", max_tokens: int | None = None):
        from sentence_transformers import SentenceTransformer

        self.name = name
        self.device = detect_device(device)
        self.model = SentenceTransformer(name, device=self.device)
        self.model.eval()
        self.tokenizer = self.model.tokenizer
        self.n_special = self.tokenizer.num_special_tokens_to_add()
        if self.n_special != 2 or self.tokenizer.cls_token_id is None or self.tokenizer.sep_token_id is None:
            raise ValueError(f"{name}: expected a tokenizer that wraps passages as [CLS] ... [SEP]")
        self.max_seq_length = int(self.model.max_seq_length)
        if max_tokens is not None:
            self.max_seq_length = min(self.max_seq_length, int(max_tokens))
        self.max_content_tokens = self.max_seq_length - self.n_special

    def tokenize(self, text: str) -> list[int]:
        return self.tokenizer(text, add_special_tokens=False, truncation=False, verbose=False)["input_ids"]

    def encode_ids(self, passages: Sequence[Sequence[int]], batch_size: int = 16) -> np.ndarray:
        import torch

        tok = self.tokenizer
        wrapped = [[tok.cls_token_id, *ids, tok.sep_token_id] for ids in passages]
        too_long = [len(w) for w in wrapped if len(w) > self.max_seq_length]
        if too_long:
            raise ValueError(f"{len(too_long)} passages exceed max_seq_length={self.max_seq_length} (longest {max(too_long)})")
        order = sorted(range(len(wrapped)), key=lambda i: -len(wrapped[i]))
        out: np.ndarray | None = None
        for start in range(0, len(order), batch_size):
            idx = order[start:start + batch_size]
            width = len(wrapped[idx[0]])
            input_ids = torch.full((len(idx), width), tok.pad_token_id, dtype=torch.long)
            mask = torch.zeros((len(idx), width), dtype=torch.long)
            for row, i in enumerate(idx):
                input_ids[row, :len(wrapped[i])] = torch.tensor(wrapped[i], dtype=torch.long)
                mask[row, :len(wrapped[i])] = 1
            with torch.no_grad():
                emb = self.model({"input_ids": input_ids.to(self.device), "attention_mask": mask.to(self.device)})["sentence_embedding"]
            emb = emb.detach().cpu().numpy().astype(np.float64)
            if out is None:
                out = np.zeros((len(wrapped), emb.shape[1]))
            out[idx] = emb
            if (start // batch_size) % 10 == 0:
                print(f"[embed]   {min(start + batch_size, len(order))}/{len(order)} passages")
        if out is None:
            raise ValueError("no passages to encode")
        return out


@dataclass(frozen=True)
class EmbeddingResult:
    embeddings: np.ndarray
    n_passages: list[int]
    n_tokens: list[int]
    n_tokens_used: list[int]


def embed_texts(texts: Sequence[str], encoder: Encoder, rule: str, overlap: int = 32, batch_size: int = 16) -> EmbeddingResult:
    passages: list[list[int]] = []
    owner: list[int] = []
    n_passages: list[int] = []
    n_total: list[int] = []
    n_used: list[int] = []
    for doc, text in enumerate(texts):
        ids = encoder.tokenize(text)
        spans = passage_spans(len(ids), rule, encoder.max_content_tokens, overlap)
        if not spans:
            raise ValueError(f"document {doc} has no tokens; empty texts cannot be embedded")
        for start, end in spans:
            passages.append(list(ids[start:end]))
            owner.append(doc)
        n_passages.append(len(spans))
        n_total.append(len(ids))
        n_used.append(spans[-1][1])
    pooled = pool_passages(encoder.encode_ids(passages, batch_size), owner, len(texts))
    return EmbeddingResult(pooled.astype(np.float32), n_passages, n_total, n_used)


@dataclass(frozen=True)
class EmbeddingPaths:
    array: Path
    species: Path
    meta: Path

    def all(self) -> list[Path]:
        return [self.array, self.species, self.meta]


def embedding_paths(directory: Path, model: str, rule: str, mask_level: str) -> EmbeddingPaths:
    stem = f"{model_slug(model)}__{rule}__{mask_level}"
    directory = Path(directory)
    return EmbeddingPaths(directory / f"{stem}.npy", directory / f"{stem}.species.txt", directory / f"{stem}.meta.json")


def write_species_order(path: Path, species: Sequence[str]) -> Path:
    return atomic_write_text(Path(path), "\n".join(species) + "\n")


def read_species_order(path: Path) -> list[str]:
    return [line for line in Path(path).read_text(encoding="utf-8").splitlines() if line]


def save_embeddings(paths: EmbeddingPaths, embeddings: np.ndarray, species: Sequence[str], **meta: Any) -> None:
    if embeddings.shape[0] != len(species):
        raise ValueError(f"{embeddings.shape[0]} embedding rows but {len(species)} species")
    if len(set(species)) != len(species):
        raise ValueError("duplicate species in embedding order")
    paths.array.parent.mkdir(parents=True, exist_ok=True)
    tmp = paths.array.with_name(f".{paths.array.name}.tmp.npy")
    np.save(tmp, embeddings)
    tmp.replace(paths.array)
    write_species_order(paths.species, species)
    write_json(paths.meta, {"file": paths.array.name, "sha256": sha256(paths.array), "shape": list(embeddings.shape), **meta})


def load_embeddings(paths: EmbeddingPaths) -> tuple[np.ndarray, list[str]]:
    embeddings = np.load(paths.array)
    species = read_species_order(paths.species)
    if embeddings.shape[0] != len(species):
        raise ValueError(f"{paths.array.name}: {embeddings.shape[0]} rows but {len(species)} species in {paths.species.name}")
    return embeddings, species


def cache_is_current(paths: EmbeddingPaths, params: dict[str, Any]) -> tuple[bool, str]:
    if not is_cached(paths.all()):
        return False, "not cached"
    meta = read_json(paths.meta)
    for key, value in params.items():
        if meta.get(key) != value:
            return False, f"{key} changed ({meta.get(key)!r} -> {value!r})"
    return True, "cached"
