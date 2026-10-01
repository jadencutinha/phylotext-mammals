import hashlib

import numpy as np
import pytest

from lineage_charisma.embed import (
    cache_is_current,
    chunk_spans,
    embed_texts,
    embedding_paths,
    l2_normalize,
    load_embeddings,
    model_slug,
    passage_spans,
    pool_passages,
    save_embeddings,
)


class FakeEncoder:
    """Whitespace tokens; a passage embeds as the sum of per-token hash vectors."""

    def __init__(self, max_content_tokens=8, dim=16):
        self.max_content_tokens = max_content_tokens
        self.dim = dim
        self.seen: list[list[int]] = []

    def tokenize(self, text):
        return [int.from_bytes(hashlib.sha256(w.encode()).digest()[:4], "little") for w in text.split()]

    def _vector(self, token):
        return np.random.default_rng(token).normal(size=self.dim)

    def encode_ids(self, passages, batch_size=16):
        self.seen.extend(list(p) for p in passages)
        return np.stack([np.sum([self._vector(t) for t in p], axis=0) for p in passages])


@pytest.mark.parametrize("max_len,overlap", [(8, 0), (8, 2), (10, 3), (510, 32), (382, 32), (5, 4)])
@pytest.mark.parametrize("n", [1, 5, 8, 9, 10, 17, 100, 511, 1021, 2609])
def test_chunk_spans_never_exceed_max_length_and_cover_all_tokens(n, max_len, overlap):
    spans = chunk_spans(n, max_len, overlap)
    assert all(0 < end - start <= max_len for start, end in spans)
    covered = set()
    for start, end in spans:
        covered.update(range(start, end))
    assert covered == set(range(n))
    assert spans[0][0] == 0 and spans[-1][1] == n
    assert all(a[0] < b[0] and a[1] < b[1] for a, b in zip(spans, spans[1:]))


def test_chunk_spans_short_text_is_one_chunk():
    assert chunk_spans(8, 8, 2) == [(0, 8)]
    assert chunk_spans(3, 8, 2) == [(0, 3)]
    assert chunk_spans(0, 8, 2) == []


def test_chunk_spans_overlap_and_balance():
    spans = chunk_spans(1000, 510, 32)
    assert len(spans) == 3
    assert all(prev[1] - nxt[0] == 32 for prev, nxt in zip(spans, spans[1:]))
    sizes = [end - start for start, end in spans]
    assert max(sizes) - min(sizes) <= len(spans)


def test_chunk_spans_uses_fewest_chunks():
    assert len(chunk_spans(510, 510, 32)) == 1
    assert len(chunk_spans(511, 510, 32)) == 2
    assert len(chunk_spans(988, 510, 32)) == 2
    assert len(chunk_spans(989, 510, 32)) == 3


@pytest.mark.parametrize("max_len,overlap", [(0, 0), (8, 8), (8, -1), (8, 9)])
def test_chunk_spans_rejects_bad_parameters(max_len, overlap):
    with pytest.raises(ValueError):
        chunk_spans(20, max_len, overlap)


def test_truncate_keeps_only_the_first_tokens():
    assert passage_spans(100, "truncate", 8, 2) == [(0, 8)]
    assert passage_spans(5, "truncate", 8, 2) == [(0, 5)]
    assert passage_spans(20, "chunk", 8, 2) == chunk_spans(20, 8, 2)
    with pytest.raises(ValueError):
        passage_spans(20, "mean", 8, 2)


def test_l2_normalize_gives_unit_norm():
    x = np.random.default_rng(0).normal(size=(7, 12)) * 50
    out = l2_normalize(x)
    assert np.allclose(np.linalg.norm(out, axis=1), 1.0)
    assert np.allclose(out[0] * np.linalg.norm(x[0]), x[0])
    with pytest.raises(ValueError):
        l2_normalize(np.zeros((2, 3)))


def test_pool_passages_is_mean_then_unit_norm():
    emb = np.array([[1.0, 0.0], [0.0, 1.0], [3.0, 0.0], [0.0, 2.0], [0.0, 4.0]])
    pooled = pool_passages(emb, [0, 0, 1, 2, 2], 3)
    assert np.allclose(np.linalg.norm(pooled, axis=1), 1.0)
    assert np.allclose(pooled[0], np.array([0.5, 0.5]) / np.linalg.norm([0.5, 0.5]))
    assert np.allclose(pooled[1], [1.0, 0.0])
    assert np.allclose(pooled[2], [0.0, 1.0])


def test_pool_passages_rejects_document_without_passages():
    with pytest.raises(ValueError):
        pool_passages(np.ones((2, 3)), [0, 2], 3)
    with pytest.raises(ValueError):
        pool_passages(np.ones((2, 3)), [0], 1)


def test_embed_texts_chunk_rule_is_unit_norm_and_covers_every_token():
    enc = FakeEncoder(max_content_tokens=8)
    texts = [" ".join(f"w{i}_{j}" for j in range(n)) for i, n in enumerate([3, 8, 9, 30, 61])]
    result = embed_texts(texts, enc, "chunk", overlap=2)
    assert result.embeddings.shape == (5, enc.dim)
    assert np.allclose(np.linalg.norm(result.embeddings, axis=1), 1.0, atol=1e-6)
    assert result.n_passages[:2] == [1, 1] and all(n > 1 for n in result.n_passages[2:])
    assert result.n_tokens == result.n_tokens_used == [3, 8, 9, 30, 61]
    assert all(len(p) <= 8 for p in enc.seen)
    assert {t for p in enc.seen for t in p} == {t for text in texts for t in enc.tokenize(text)}


def test_embed_texts_truncate_rule_ignores_text_past_the_limit():
    head = " ".join(f"h{j}" for j in range(8))
    enc = FakeEncoder(max_content_tokens=8)
    result = embed_texts([head, head + " tail one two three", head + " something else entirely"], enc, "truncate")
    assert result.n_passages == [1, 1, 1]
    assert result.n_tokens_used == [8, 8, 8] and result.n_tokens == [8, 12, 11]
    assert np.allclose(result.embeddings[0], result.embeddings[1])
    assert np.allclose(result.embeddings[0], result.embeddings[2])
    chunked = embed_texts([head, head + " tail one two three"], FakeEncoder(max_content_tokens=8), "chunk", overlap=2)
    assert not np.allclose(chunked.embeddings[0], chunked.embeddings[1])


def test_embed_texts_rejects_empty_text():
    with pytest.raises(ValueError):
        embed_texts(["some words", ""], FakeEncoder(), "chunk", overlap=2)


def test_model_slug_and_paths(tmp_path):
    assert model_slug("BAAI/bge-large-en-v1.5") == "baai_bge_large_en_v1_5"
    paths = embedding_paths(tmp_path, "sentence-transformers/all-mpnet-base-v2", "truncate", "taxonomy")
    assert paths.array.name == "sentence_transformers_all_mpnet_base_v2__truncate__taxonomy.npy"
    assert paths.species.name == "sentence_transformers_all_mpnet_base_v2__truncate__taxonomy.species.txt"


def test_save_and_load_embeddings_round_trip_and_cache_check(tmp_path):
    paths = embedding_paths(tmp_path, "m/x", "chunk", "none")
    emb = l2_normalize(np.random.default_rng(1).normal(size=(3, 4))).astype(np.float32)
    species = ["Canis lupus", "Panthera leo", "Ursus arctos"]
    params = {"model": "m/x", "rule": "chunk", "chunk_overlap": 32}
    assert cache_is_current(paths, params) == (False, "not cached")
    save_embeddings(paths, emb, species, **params)
    loaded, order = load_embeddings(paths)
    assert order == species and np.array_equal(loaded, emb)
    assert cache_is_current(paths, params)[0]
    assert not cache_is_current(paths, {**params, "chunk_overlap": 64})[0]


def test_save_embeddings_rejects_misaligned_species(tmp_path):
    paths = embedding_paths(tmp_path, "m/x", "chunk", "none")
    with pytest.raises(ValueError):
        save_embeddings(paths, np.ones((3, 4)), ["a", "b"])
    with pytest.raises(ValueError):
        save_embeddings(paths, np.ones((2, 4)), ["a", "a"])


def test_load_embeddings_rejects_row_count_mismatch(tmp_path):
    paths = embedding_paths(tmp_path, "m/x", "chunk", "none")
    save_embeddings(paths, np.ones((3, 4)), ["a", "b", "c"])
    paths.species.write_text("a\nb\n")
    with pytest.raises(ValueError):
        load_embeddings(paths)
