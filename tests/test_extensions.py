from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pytest

from lineage_charisma import extensions as ext
from lineage_charisma.extensions import (
    EXCLUDED,
    PROVIDERS,
    ApiEmbedder,
    ApiError,
    Registered,
    Tokenizer,
    assign_groups,
    classify_heading,
    cost_estimate,
    embed_with_provider,
    group_texts,
    normalize_heading,
    registered_tests,
    split_passages,
    tfidf_distance,
)
from lineage_charisma.stats import mantel, partial_mantel
from lineage_charisma.wiki import split_sections

ROOT = Path(__file__).resolve().parents[1]
AMENDMENT = ROOT / "docs" / "preregistration_amendment_extensions.md"


# ---------------------------------------------------------------- the code matches the amendment

def test_settings_match_the_amendment():
    doc = AMENDMENT.read_text()
    s = ext.TFIDF_SETTINGS
    assert s["min_df"] == 5 and s["max_df"] == 1.0 and s["ngram_range"] == (1, 1) and s["sublinear_tf"] and s["norm"] == "l2"
    assert f"`{s['token_pattern']}`" in doc
    assert "`min_df=5`" in doc and "`sublinear_tf=True`" in doc
    for p in PROVIDERS.values():
        assert f"`{p.model}`" in doc and f"`{p.key_env}`" in doc
    assert ext.MIN_TOKENS == 50 and ext.MIN_GROUP_SPECIES == 50 and "at least 50 tokens" in doc and "fewer than **50 species**" in doc


def test_heading_table_in_the_amendment_matches_the_code():
    doc = AMENDMENT.read_text()
    for group, headings in ext.EXACT_HEADINGS.items():
        line = next(l for l in doc.splitlines() if l.strip().startswith(f"| {group} |"))
        assert [h.strip() for h in line.split("|")[2].split(",")] == list(headings), group
    for group, words in ext.KEYWORDS.items():
        line = next(l for l in doc.splitlines() if l.strip().startswith(f"- {group}:"))
        assert [w.strip() for w in line.split(":", 1)[1].split(",")] == list(words), group


def test_english_stopword_list_has_the_registered_size():
    from sklearn.feature_extraction.text import ENGLISH_STOP_WORDS

    assert len(ENGLISH_STOP_WORDS) == 318


# ---------------------------------------------------------------- B. TF-IDF

def test_tfidf_ignores_mask_token_stopwords_digits_and_rare_terms():
    texts = ["The [TAXON] eats fish 12 times", "A [TAXON] eats fish and crabs", "The [TAXON] hunts fish", "x [TAXON] eats crabs"]
    from sklearn.feature_extraction.text import TfidfVectorizer

    vec = TfidfVectorizer(**{**ext.TFIDF_SETTINGS, "min_df": 2}).fit([t.replace("[TAXON]", " ") for t in texts])
    assert set(vec.vocabulary_) == {"eats", "fish", "crabs"}       # no taxon, the, a, and, 12, x; hunts is in 1 document
    res = tfidf_distance(texts, min_df=2)
    assert res.vocabulary_size == 3 and res.dropped == []
    d = res.distance
    assert np.allclose(d, d.T) and np.allclose(np.diag(d), 0) and d.min() >= 0
    assert d[0, 2] < d[0, 3]                                          # 0 and 2 share fish; 0 and 3 share eats only, with crabs unshared


def test_tfidf_reports_documents_with_no_vocabulary():
    res = tfidf_distance(["fish fish", "fish crabs", "the [TAXON] and"], min_df=1)
    assert res.dropped == [2] and res.kept == [0, 1] and res.distance.shape == (2, 2)


# ---------------------------------------------------------------- C. section groups

@pytest.mark.parametrize("heading, group, rule", [
    ("Description", "appearance", "exact"),
    ("Behaviour and ecology", "behavior", "exact"),
    ("Distribution & habitat", "distribution_habitat", "exact"),
    ("Hunting and diet", "diet", "exact"),
    ("External links", EXCLUDED, "exact"),
    ("Relationship with humans", "other", "exact"),
    ("Diet and behaviour", "diet", "keyword"),                       # first keyword wins
    ("Behaviour and diet", "behavior", "keyword"),
    ("Habitat use and diet", "distribution_habitat", "keyword"),
    ("Notes on ecology", EXCLUDED, "keyword"),
    ("Colouration of the pelage", "appearance", "keyword"),
    ("New Zealand", None, "none"),
])
def test_classify_heading(heading, group, rule):
    assert classify_heading(heading) == (group, rule)


def test_normalize_heading():
    assert normalize_heading("  Behaviour   &  Ecology ") == "behavior and ecology"
    assert normalize_heading("Colouration") == "coloration"


EXTRACT = """The animal is a small mammal of [[somewhere]].

== Description ==
It is brown and small.

=== Size ===
Head and body 10 cm.

== Ecology ==
It lives alone.

=== New Zealand ===
Introduced populations persist.

=== Diet ===
It eats insects.

== History ==
Described in 1900.

=== Early accounts ===
Old notes.

== References ==
Ref one.

=== Online ===
Link.
"""


def test_subsections_inherit_unless_their_own_heading_matches():
    groups = [g for g, _ in assign_groups(split_sections(EXTRACT))]
    assert groups == ["appearance", "appearance", "behavior", "behavior", "diet", "other", "other", EXCLUDED, EXCLUDED]


def test_group_texts_keep_the_lead_as_its_own_row_and_leave_out_excluded_sections():
    g = group_texts(EXTRACT)
    assert set(g) == {"appearance", "behavior", "diet", "other", "lead"}
    assert g["lead"] == "The animal is a small mammal of somewhere."
    assert g["appearance"] == "It is brown and small.\n\nHead and body 10 cm."
    assert g["behavior"] == "It lives alone.\n\nIntroduced populations persist."
    assert all("small mammal" not in v for k, v in g.items() if k != "lead") and "Ref one" not in " ".join(g.values())


# ---------------------------------------------------------------- A. commercial embeddings

WORDS = Tokenizer(lambda t: list(range(len(t.split()))), lambda ids: " ".join(f"w{i}" for i in ids), "words")


class FakeApi:
    def __init__(self, failures=0, dim=4):
        self.calls, self.failures, self.dim = [], failures, dim

    def __call__(self, url, headers, body):
        self.calls.append((url, headers, body))
        if self.failures:
            self.failures -= 1
            return 429, {"error": "rate limited"}
        data = [{"index": i, "embedding": [float(len(t)), 1.0] + [0.0] * (self.dim - 2)} for i, t in enumerate(body["input"])]
        return 200, {"data": data[::-1], "model": body["model"], "usage": {"total_tokens": 10 * len(data)}}


@pytest.fixture
def key(monkeypatch):
    monkeypatch.setenv("OPENAI_API_KEY", "sk-test-secret")


def test_embedder_batches_caches_and_orders_by_index(tmp_path, key):
    api = FakeApi()
    provider = PROVIDERS["openai"]
    emb = ApiEmbedder(provider, tmp_path, post=api, sleep=lambda s: None)
    passages = [f"text {'x' * i}" for i in range(200)]
    out = emb.embed(passages, [5] * 200)
    assert out.shape == (200, 4) and out[3, 0] == len(passages[3])
    assert [len(c[2]["input"]) for c in api.calls] == [96, 96, 8]
    assert api.calls[0][2] == {"model": "text-embedding-3-large", "input": passages[:96], "encoding_format": "float"}
    emb.embed(passages, [5] * 200)
    assert len(api.calls) == 3                                     # everything served from the cache
    cached = json.loads(emb.cache_path(passages[0]).read_text())
    assert cached["model"] == "text-embedding-3-large" and "sk-test" not in json.dumps(cached)


def test_embedder_splits_batches_by_tokens(tmp_path, key):
    api = FakeApi()
    ApiEmbedder(PROVIDERS["openai"], tmp_path, post=api, sleep=lambda s: None).embed(["a", "b", "c"], [60_000, 30_000, 30_000])
    assert [len(c[2]["input"]) for c in api.calls] == [2, 1]


def test_embedder_retries_then_stops(tmp_path, key):
    waits = []
    api = FakeApi(failures=2)
    ApiEmbedder(PROVIDERS["openai"], tmp_path, post=api, sleep=waits.append).embed(["a"], [1])
    assert waits == [2.0, 4.0] and len(api.calls) == 3
    with pytest.raises(ApiError, match="no valid response after 5 retries"):
        ApiEmbedder(PROVIDERS["openai"], tmp_path / "b", post=FakeApi(failures=99), sleep=lambda s: None).embed(["a"], [1])


def test_client_errors_are_not_retried_and_never_show_the_key(tmp_path, key):
    calls = []

    def post(url, headers, body):
        calls.append(1)
        return 401, {"error": {"message": "invalid key"}}

    with pytest.raises(ApiError) as e:
        ApiEmbedder(PROVIDERS["openai"], tmp_path, post=post, sleep=lambda s: None).embed(["a"], [1])
    assert len(calls) == 1 and "sk-test-secret" not in str(e.value)


def test_missing_key_sends_nothing(tmp_path, monkeypatch):
    monkeypatch.delenv("VOYAGE_API_KEY", raising=False)
    api = FakeApi()
    with pytest.raises(ApiError, match="VOYAGE_API_KEY is not set"):
        ApiEmbedder(PROVIDERS["voyage"], tmp_path, post=api).embed(["a"], [1])
    assert api.calls == []


def test_voyage_request_body():
    assert PROVIDERS["voyage"].body(["a"]) == {"model": "voyage-3-large", "input": ["a"], "input_type": None, "truncation": False,
                                               "output_dtype": "float"}


def test_long_texts_are_chunked_at_the_provider_limit_and_pooled(tmp_path, key):
    texts = ["one two three", " ".join(["w"] * 25)]
    passages, owner, counts = split_passages(texts, WORDS, limit=10, overlap=2)
    assert owner == [0, 1, 1, 1] and passages[0] == "one two three" and max(counts) <= 10
    res = embed_with_provider(texts, ApiEmbedder(PROVIDERS["openai"], tmp_path, post=FakeApi(), sleep=lambda s: None), WORDS, 10, 2)
    assert res.n_passages == [1, 3] and np.allclose(np.linalg.norm(res.embeddings, axis=1), 1, atol=1e-6)


def test_cost_estimate():
    rows = cost_estimate({"pilot": 1_000_000}, PROVIDERS["openai"]) + cost_estimate({"pilot": 1_000_000}, PROVIDERS["voyage"])
    assert [round(r["usd"], 4) for r in rows] == [0.13, 0.18]
    assert [round(r["usd_pessimistic"], 4) for r in rows] == [0.26, 0.45]


# ---------------------------------------------------------------- shared tests

def _dist(points):
    return np.abs(points[:, None] - points[None, :]) if points.ndim == 1 else np.linalg.norm(points[:, None] - points[None, :], axis=-1)


def test_registered_tests_match_the_pilot_implementation():
    rng = np.random.default_rng(1)
    n = 30
    order = [f"s{i:02d}" for i in range(n)]
    phylo = _dist(rng.normal(size=(n, 2)))
    text = phylo + _dist(rng.normal(size=(n, 3)))
    keep = order[2:]
    idx = np.arange(2, n)
    eco, length = _dist(rng.normal(size=(n - 2, 2))), _dist(rng.normal(size=n - 2))
    reg = Registered(order, phylo, keep, eco, length)
    h1, h2 = registered_tests(text, reg, 199, 7, 0.05, workers=1)
    ref1 = mantel(text, phylo, permutations=199, seed=7)
    sub = np.ix_(idx, idx)
    ref2 = partial_mantel(phylo[sub], text[sub], [eco, length], permutations=199, seed=7)
    assert h1["r"] == pytest.approx(ref1.r) and h1["p"] == pytest.approx(ref1.p) and h1["outcome"] == "positive"
    assert h2["r"] == pytest.approx(ref2.r) and h2["p"] == pytest.approx(ref2.p) and h2["raw_r"] == pytest.approx(ref2.r_marginal)
    assert h2["rule_a"] == "survives" and h2["n_species"] == n - 2


def test_restrict_keeps_matrices_aligned():
    order = ["a", "b", "c", "d"]
    m = _dist(np.arange(4.0))
    reg = Registered(order, m, ["b", "c", "d"], m[1:, 1:], m[1:, 1:]).restrict(["a", "c", "d"])
    assert reg.order == ["a", "c", "d"] and reg.h2_species == ["c", "d"]
    assert reg.phylo[0, 1] == 2.0 and reg.ecology[0, 1] == 1.0


def test_guard_refuses_an_uncommitted_file(tmp_path):
    (tmp_path / "docs").mkdir()
    (tmp_path / ext.AMENDMENT).write_text("draft")
    ok, why = ext.amendment_is_committed(tmp_path)
    assert not ok and "must be committed" in why
    assert ext.amendment_is_committed(tmp_path / "nowhere")[0] is False


def test_scripts_import_without_side_effects():
    import importlib.util

    for name in ("36_api_embeddings.py", "37_tfidf_baseline.py", "38_section_groups.py"):
        spec = importlib.util.spec_from_file_location(re.sub(r"\W", "_", name), ROOT / "scripts" / name)
        spec.loader.exec_module(importlib.util.module_from_spec(spec))
