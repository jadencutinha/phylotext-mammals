import numpy as np
import pytest
import requests

from lineage_charisma.attention import (
    PageviewsClient,
    abs_difference_matrix,
    attention_matrices,
    attention_variables,
    combined_distance_matrix,
    month_range,
    pairwise_mean_matrix,
    resolve_articles,
    standardize,
    summed_monthly_views,
)
from lineage_charisma.distances import validate_distance_matrix
from lineage_charisma.wiki import PageInfo


class FakeResponse:
    def __init__(self, status, items=None, headers=None):
        self.status_code, self._items, self.headers = status, items or [], headers or {}

    def json(self):
        return {"items": self._items}

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError(str(self.status_code))


class FakeSession:
    def __init__(self, responses):
        self.responses, self.headers, self.urls = list(responses), {}, []

    def get(self, url, timeout=None):
        self.urls.append(url)
        r = self.responses.pop(0)
        if isinstance(r, Exception):
            raise r
        return r


def items(*pairs):
    return [{"timestamp": f"{month.replace('-', '')}0100", "views": views} for month, views in pairs]


def make_client(tmp_path, responses, **kw):
    session = FakeSession(responses)
    client = PageviewsClient("https://pv/per-article", "test-agent/0.1", cache_dir=tmp_path, min_interval_seconds=0, backoff_seconds=0, session=session, **kw)
    return client, session


def test_month_range_is_inclusive():
    months = month_range("2023-01", "2025-12")
    assert len(months) == 36 and months[0] == "2023-01" and months[-1] == "2025-12"


def test_url_encodes_title_and_window(tmp_path):
    client, _ = make_client(tmp_path, [])
    assert client.url("Rüppell's fox", "2023-01", "2024-02") == (
        "https://pv/per-article/en.wikipedia/all-access/user/R%C3%BCppell%27s_fox/monthly/2023010100/2024022900")
    assert "AC%2FDC" in client.url("AC/DC", "2023-01", "2023-01")


def test_client_sets_user_agent_caches_and_treats_404_as_no_views(tmp_path):
    client, session = make_client(tmp_path, [FakeResponse(200, items(("2023-01", 10), ("2023-02", 30))), FakeResponse(404)])
    assert session.headers["User-Agent"] == "test-agent/0.1"
    assert client.monthly("Lion", "2023-01", "2023-02") == {"2023-01": 10, "2023-02": 30}
    assert client.monthly("Lion", "2023-01", "2023-02") == {"2023-01": 10, "2023-02": 30}
    assert client.monthly("Old lion title", "2023-01", "2023-02") == {}
    assert client.monthly("Old lion title", "2023-01", "2023-02") == {}
    assert client.network_calls == 2


def test_client_retries_then_fails(tmp_path):
    client, _ = make_client(tmp_path, [FakeResponse(429), requests.ConnectionError(), FakeResponse(200, items(("2023-01", 5)))])
    assert client.monthly("Lion", "2023-01", "2023-01") == {"2023-01": 5}
    failing, _ = make_client(tmp_path / "b", [FakeResponse(503)] * 3, max_retries=2)
    with pytest.raises(RuntimeError, match="failed after 3 attempts"):
        failing.monthly("Tiger", "2023-01", "2023-01")
    offline, _ = make_client(tmp_path / "c", [], offline=True)
    with pytest.raises(RuntimeError, match="offline"):
        offline.monthly("Tiger", "2023-01", "2023-01")


def test_summed_views_add_redirects_and_fill_empty_months(tmp_path):
    client, _ = make_client(tmp_path, [FakeResponse(200, items(("2023-02", 100), ("2023-03", 120))), FakeResponse(200, items(("2023-01", 90), ("2023-02", 5))), FakeResponse(404)])
    total = summed_monthly_views(client, ["New title", "Old title", "Misspelling"], "2023-01", "2023-04")
    assert total.to_dict() == {"2023-01": 90, "2023-02": 105, "2023-03": 120, "2023-04": 0}
    assert total.median() == 97.5


class FakeWiki:
    api_url = "https://api"

    def __init__(self):
        self.calls = []

    def query_titles(self, titles):
        return {"Old name": PageInfo("Old name", "New name", exists=True, redirected=True), "Lion": PageInfo("Lion", "Lion", exists=True),
                "Gone": PageInfo("Gone", "Gone", exists=False)}

    def get(self, url, params):
        self.calls.append(params)
        if params["titles"] == "New name" and "rdcontinue" not in params:
            return {"continue": {"rdcontinue": "x", "continue": "||"}, "query": {"pages": [{"title": "New name", "redirects": [{"title": "Old name"}]}]}}
        if params["titles"] == "New name":
            return {"query": {"pages": [{"title": "New name", "redirects": [{"title": "Older name"}]}]}}
        return {"query": {"pages": [{"title": "Lion"}]}}


def test_resolve_articles_follows_renames_and_collects_redirects():
    wiki = FakeWiki()
    out = resolve_articles(wiki, ["Old name", "Lion"])
    assert out == {"Old name": ("New name", ["Old name", "Older name"]), "Lion": ("Lion", [])}
    assert all(c["rdnamespace"] == 0 for c in wiki.calls)
    with pytest.raises(RuntimeError, match="no longer exists"):
        resolve_articles(wiki, ["Gone"])


def test_attention_variables():
    v = attention_variables([0, 99, 9999], [50, 500, 2000])
    np.testing.assert_allclose(v["log_pageviews"], np.log([1, 100, 10000]))
    np.testing.assert_allclose(v["log_length"], np.log([50, 500, 2000]))
    for bad in (([-1, 2], [5, 5]), ([1, 2], [0, 5]), ([np.nan, 2], [5, 5])):
        with pytest.raises(ValueError):
            attention_variables(*bad)


def test_distance_matrices_known_answer():
    np.testing.assert_allclose(abs_difference_matrix([1.0, 4.0, 2.0]), [[0, 3, 1], [3, 0, 2], [1, 2, 0]])
    z = standardize([1.0, 2.0, 3.0])
    assert z.mean() == pytest.approx(0) and z.std() == pytest.approx(1)
    with pytest.raises(ValueError):
        standardize([2.0, 2.0])
    # two variables that differ only in scale: each contributes the same standardized difference
    a = [0.0, 1.0, 2.0]
    combined = combined_distance_matrix(a, [10 * x for x in a])
    gap = abs(standardize(a)[1] - standardize(a)[0])
    assert combined[0, 1] == pytest.approx(np.sqrt(2) * gap) and combined[0, 2] == pytest.approx(2 * np.sqrt(2) * gap)
    mean = pairwise_mean_matrix([1.0, 3.0])
    np.testing.assert_allclose(mean, [[1, 2], [2, 3]])


def test_attention_matrices_are_valid_and_in_input_order():
    rng = np.random.default_rng(0)
    v = attention_variables(rng.integers(0, 50000, size=30), rng.integers(50, 2500, size=30))
    mats = attention_matrices(v)
    assert set(mats) == {"attention_pageviews", "attention_length", "attention_combined"}
    for m in mats.values():
        assert m.shape == (30, 30)
        validate_distance_matrix(m)
    assert mats["attention_pageviews"][3, 7] == pytest.approx(abs(v["log_pageviews"][3] - v["log_pageviews"][7]))
    # reordering the species reorders the matrices the same way
    perm = rng.permutation(30)
    shuffled = attention_matrices(v.iloc[perm].reset_index(drop=True))
    for name, m in mats.items():
        np.testing.assert_allclose(shuffled[name], m[np.ix_(perm, perm)])
