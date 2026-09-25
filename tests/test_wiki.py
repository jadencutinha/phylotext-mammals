import json
from pathlib import Path

import pytest
import requests

from lineage_charisma.wiki import (
    PageInfo,
    WikiClient,
    classify_article,
    clean_text,
    find_description,
    page_from_response,
    parse_entities,
    parse_title_query,
    preprocess,
    resolve_species_titles,
    split_sections,
)

FIXTURES = Path(__file__).parent / "fixtures"
HEADINGS = ["description", "physical description", "characteristics", "appearance"]


def load(name):
    return json.loads((FIXTURES / name).read_text())


def test_parse_title_query_follows_normalization_and_redirects():
    requested = ["Canis lupaster", "Canis anthus", "felis catus", "Leopardus colocola", "Zzyzx nonexistens", "Jaguar (disambiguation)"]
    info = parse_title_query(load("title_query.json"), requested)
    assert info["Canis lupaster"].title == "African wolf" and info["Canis lupaster"].redirected
    assert info["Canis anthus"].wikibase_item == info["Canis lupaster"].wikibase_item == "Q20744349"
    assert info["felis catus"].title == "Cat"
    assert info["Leopardus colocola"].title == "Pampas cat"
    assert not info["Zzyzx nonexistens"].exists
    assert info["Jaguar (disambiguation)"].disambiguation


def test_parse_entities_extracts_taxon_name_rank_and_sitelink():
    ents = parse_entities(load("entities.json"))
    assert ents["Q20744349"]["taxon_names"] == ["Canis anthus"]
    assert "Q7432" in ents["Q20744349"]["ranks"]
    assert ents["Q192233"]["taxon_names"] == ["Pardofelis temminckii"]
    assert ents["Q146"]["taxon_names"] == []
    assert ents["Q20980826"]["taxon_names"] == ["Felis catus"]
    assert ents["Q23772769"]["enwiki"] is None


@pytest.mark.parametrize(
    "species,other,taxon,ranks,source,redirected,expected",
    [
        ("Panthera leo", set(), ["Panthera leo"], ["Q7432"], "scientific_name", False, "exact"),
        ("Canis rufus", set(), ["Canis lupus rufus"], ["Q68947"], "scientific_name", True, "subspecies_article"),
        ("Urva javanica", {"Herpestes javanicus"}, ["Herpestes javanicus"], ["Q7432"], "common_name", True, "synonym"),
        ("Catopuma temminckii", set(), ["Pardofelis temminckii"], ["Q7432"], "scientific_name", True, "genus_transfer"),
        ("Canis lupaster", set(), ["Canis anthus"], ["Q7432"], "scientific_name", True, "name_redirect"),
        ("Canis lupaster", set(), ["Canis anthus"], ["Q7432"], "common_name", True, None),
        ("Canis lupaster", set(), ["Canis anthus"], ["Q7432"], "scientific_name", False, None),
        ("Canis lupaster", set(), ["Canis aureus lupaster", "Canis lupus"], ["Q68947"], "synonym", True, "subspecies_article"),
        ("Leopardus garleppi", set(), ["Leopardus colocola"], ["Q7432"], "common_name", True, None),
        ("Panthera leo", set(), [], [], "scientific_name", True, None),
    ],
)
def test_classify_article(species, other, taxon, ranks, source, redirected, expected):
    assert classify_article(species, other, taxon, ranks, source, redirected) == expected


class FakeClient:
    def __init__(self, pages, entities, search=None):
        self.pages, self.ents, self.search = pages, entities, search or {}
        self.title_queries = []

    def query_titles(self, titles):
        titles = list(titles)
        self.title_queries.append(titles)
        return {t: self.pages.get(t, PageInfo(t, t, exists=False)) for t in titles}

    def entities(self, qids):
        return {q: self.ents[q] for q in qids if q in self.ents}

    def search_taxon(self, name):
        return self.search.get(name, [])


def page(query, title, qid, redirected=True, disambiguation=False):
    return PageInfo(query=query, title=title, exists=True, redirected=redirected, wikibase_item=qid, disambiguation=disambiguation)


def ent(names, rank="Q7432", enwiki=None):
    return {"taxon_names": names, "ranks": [rank], "enwiki": enwiki}


def test_resolve_prefers_exact_and_split_species_loses_parent_article():
    pages = {
        "Leopardus colocola": page("Leopardus colocola", "Pampas cat", "Q1"),
        "Leopardus pajeros": page("Leopardus pajeros", "Pampas cat", "Q1"),
        "Pampas cat": page("Pampas cat", "Pampas cat", "Q1", redirected=False),
    }
    client = FakeClient(pages, {"Q1": ent(["Leopardus colocola"])})
    cands = {
        "Leopardus colocola": [("Leopardus colocola", "scientific_name")],
        "Leopardus pajeros": [("Leopardus pajeros", "scientific_name"), ("Pampas cat", "common_name")],
    }
    res, notes = resolve_species_titles(cands, {s: set() for s in cands}, client)
    assert res["Leopardus colocola"]["wiki_title"] == "Pampas cat"
    assert "Leopardus pajeros" not in res
    assert "Leopardus pajeros" in notes


def test_resolve_name_redirect_when_wikidata_lags_taxonomy():
    pages = {"Canis lupaster": page("Canis lupaster", "African wolf", "Q2")}
    client = FakeClient(pages, {"Q2": ent(["Canis anthus"])})
    res, _ = resolve_species_titles({"Canis lupaster": [("Canis lupaster", "scientific_name")]}, {"Canis lupaster": set()}, client)
    assert res["Canis lupaster"]["wiki_title"] == "African wolf"
    assert res["Canis lupaster"]["wiki_validation"] == "name_redirect"


def test_resolve_skips_disambiguation_and_uses_next_candidate():
    pages = {
        "Genus species": page("Genus species", "Thing", "Q9", disambiguation=True),
        "Common thing": page("Common thing", "Common thing", "Q3", redirected=False),
    }
    client = FakeClient(pages, {"Q3": ent(["Genus species"])})
    cands = {"Genus species": [("Genus species", "scientific_name"), ("Common thing", "common_name")]}
    res, _ = resolve_species_titles(cands, {"Genus species": set()}, client)
    assert res["Genus species"]["wiki_title"] == "Common thing"
    assert res["Genus species"]["wiki_source"] == "common_name"


def test_resolve_wikidata_fallback_canonicalizes_sitelink_redirect():
    pages = {
        "Felis catus": page("Felis catus", "Cat", "Q146"),
    }
    ents = {"Q146": ent([], rank=None), "Q20980826": ent(["Felis catus"], enwiki="Felis catus")}
    client = FakeClient(pages, ents, search={"Felis catus": ["Q20980826"]})
    res, _ = resolve_species_titles({"Felis catus": [("Felis catus", "scientific_name")]}, {"Felis catus": set()}, client)
    assert res["Felis catus"]["wiki_title"] == "Cat"
    assert res["Felis catus"]["wiki_source"] == "wikidata_search"
    assert res["Felis catus"]["wiki_qid"] == "Q20980826"


def test_resolve_reports_items_without_english_article():
    client = FakeClient({}, {"Q5": ent(["Catopuma temminckii"], enwiki=None)}, search={"Catopuma temminckii": ["Q5"]})
    res, notes = resolve_species_titles({"Catopuma temminckii": [("Catopuma temminckii", "scientific_name")]}, {"Catopuma temminckii": set()}, client)
    assert res == {}
    assert "no English Wikipedia article" in notes["Catopuma temminckii"]


def test_page_from_response_real_fixture():
    p = page_from_response(load("page_aardwolf.json"))
    assert p["title"] == "Aardwolf"
    assert p["revid"] and p["pageid"] and p["wikibase_item"]
    assert len(p["extract"]) > 1000


def test_page_from_response_missing_raises():
    with pytest.raises(ValueError):
        page_from_response({"query": {"pages": [{"title": "X", "missing": True}]}})


def test_preprocess_real_page_with_description():
    p = page_from_response(load("page_aardwolf.json"))
    out = preprocess(p["extract"], HEADINGS)
    assert out["has_description"] and out["description_heading"] == "Description"
    assert out["lead"].startswith("The aardwolf")
    assert "==" not in out["text"]
    assert "Distribution and habitat" not in out["text"]
    assert out["text"].startswith(out["lead"])
    assert out["n_words"] == len(out["text"].split())


def test_preprocess_falls_back_through_heading_list():
    p = page_from_response(load("page_sokoke.json"))
    out = preprocess(p["extract"], HEADINGS)
    assert out["description_heading"] == "Characteristics"
    assert out["lead"]


def test_preprocess_lead_only_when_no_description():
    out = preprocess("Lead text here.\n\n== Taxonomy ==\nStuff.\n", HEADINGS)
    assert out == {**out, "lead": "Lead text here.", "description": "", "has_description": False, "description_heading": None, "text": "Lead text here."}


def test_find_description_includes_subsections_and_stops_at_same_level():
    text = "Lead.\n== Description ==\nBody.\n=== Size ===\nBig.\n=== Colour ===\nRed.\n== Ecology ==\nEats.\n"
    heading, body = find_description(split_sections(text), HEADINGS)
    assert heading == "Description"
    assert "Body." in body and "Big." in body and "Red." in body and "Eats." not in body


def test_find_description_matches_nested_heading():
    text = "Lead.\n== Biology ==\n=== Description ===\nFur.\n=== Diet ===\nMeat.\n"
    heading, body = find_description(split_sections(text), HEADINGS)
    assert heading == "Description" and "Fur." in body and "Meat." not in body


def test_find_description_respects_priority_order():
    text = "Lead.\n== Characteristics ==\nA.\n== Description ==\nB.\n"
    assert find_description(split_sections(text), HEADINGS) == ("Description", "\nB.\n")


@pytest.mark.parametrize(
    "raw,expected",
    [
        ("The lion[1] is big.[citation needed]", "The lion is big."),
        ("Text[a] more[note 2] end", "Text more end"),
        ("A <ref>Smith 2000</ref>cat", "A cat"),
        ("A {{convert|3|m}} tail", "A tail"),
        ("See [[Panthera|panthers]] and [[Felidae]]", "See panthers and Felidae"),
        ("'''Bold''' and ''italic''", "Bold and italic"),
        ("Aardwolf ( ; ) is", "Aardwolf is"),
        ("Lion (listen) roars", "Lion roars"),
        ("a b   c , d", "a b c, d"),
        ("x\n\n\n\ny", "x\n\ny"),
    ],
)
def test_clean_text(raw, expected):
    assert clean_text(raw) == expected


class FakeResponse:
    def __init__(self, status, payload=None, headers=None):
        self.status_code, self._payload, self.headers = status, payload or {}, headers or {}

    def json(self):
        return self._payload

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError(str(self.status_code))


class FakeSession:
    def __init__(self, responses):
        self.responses = list(responses)
        self.headers = {}
        self.calls = []

    def get(self, url, params=None, timeout=None):
        self.calls.append(params)
        r = self.responses.pop(0)
        if isinstance(r, Exception):
            raise r
        return r


def make_client(tmp_path, responses, **kw):
    session = FakeSession(responses)
    client = WikiClient("https://api", "https://wd", "test-agent/0.1", cache_dir=tmp_path, min_interval_seconds=0, backoff_seconds=0, session=session, **kw)
    return client, session


def test_client_sets_user_agent_and_caches(tmp_path):
    client, session = make_client(tmp_path, [FakeResponse(200, {"ok": 1})])
    assert session.headers["User-Agent"] == "test-agent/0.1"
    assert client.get("https://api", {"a": 1}) == {"ok": 1}
    assert client.get("https://api", {"a": 1}) == {"ok": 1}
    assert len(session.calls) == 1 and client.network_calls == 1
    assert session.calls[0]["maxlag"] == 5


def test_client_retries_on_429_5xx_maxlag_and_connection_errors(tmp_path):
    responses = [
        FakeResponse(429, headers={"Retry-After": "0"}),
        FakeResponse(503),
        requests.ConnectionError("boom"),
        FakeResponse(200, {"error": {"code": "maxlag"}}),
        FakeResponse(200, {"ok": 2}),
    ]
    client, session = make_client(tmp_path, responses)
    assert client.get("https://api", {"b": 1}) == {"ok": 2}
    assert len(session.calls) == 5


def test_client_gives_up_after_max_retries(tmp_path):
    client, _ = make_client(tmp_path, [FakeResponse(500)] * 3, max_retries=2)
    with pytest.raises(RuntimeError):
        client.get("https://api", {"c": 1})


def test_client_raises_on_api_error(tmp_path):
    client, _ = make_client(tmp_path, [FakeResponse(200, {"error": {"code": "badvalue"}})])
    with pytest.raises(RuntimeError, match="badvalue"):
        client.get("https://api", {"d": 1})


def test_client_offline_requires_cache(tmp_path):
    client, session = make_client(tmp_path, [], offline=True)
    with pytest.raises(RuntimeError, match="offline"):
        client.get("https://api", {"e": 1})
    assert session.calls == []


def test_query_titles_batches_by_50(tmp_path):
    titles = [f"T{i}" for i in range(120)]
    payload = {"query": {"pages": [{"title": t, "pageid": 1} for t in titles]}}
    client, session = make_client(tmp_path, [FakeResponse(200, payload)] * 3)
    out = client.query_titles(titles)
    assert len(session.calls) == 3
    assert [len(c["titles"].split("|")) for c in session.calls] == [50, 50, 20]
    assert all(out[t].exists for t in titles)
