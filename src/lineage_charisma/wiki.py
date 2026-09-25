from __future__ import annotations

import hashlib
import json
import re
import time
import unicodedata
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import requests

from .io_utils import atomic_write_text
from .taxonomy import elevated_name, epithet_stem, name_tokens, normalize_name, split_binomial

SPECIES_RANK_QID = "Q7432"
VALIDATION_GRADE = {"exact": 5, "subspecies_article": 4, "synonym": 3, "genus_transfer": 2, "name_redirect": 1}
TRUSTED_QUERY_SOURCES = {"scientific_name", "tree_name", "legacy_name"}


@dataclass
class PageInfo:
    query: str
    title: str
    exists: bool
    redirected: bool = False
    fragment: str | None = None
    wikibase_item: str | None = None
    disambiguation: bool = False


class WikiClient:
    def __init__(
        self,
        api_url: str,
        wikidata_api_url: str,
        user_agent: str,
        cache_dir: Path | None = None,
        min_interval_seconds: float = 0.5,
        max_retries: int = 5,
        backoff_seconds: float = 2.0,
        timeout_seconds: float = 30,
        session: requests.Session | None = None,
        offline: bool = False,
    ):
        self.api_url = api_url
        self.wikidata_api_url = wikidata_api_url
        self.cache_dir = Path(cache_dir) if cache_dir else None
        self.min_interval = min_interval_seconds
        self.max_retries = max_retries
        self.backoff = backoff_seconds
        self.timeout = timeout_seconds
        self.offline = offline
        self.session = session or requests.Session()
        self.session.headers.update({"User-Agent": user_agent, "Accept-Encoding": "gzip"})
        self._last_request = 0.0
        self.network_calls = 0

    @classmethod
    def from_config(cls, cfg, **overrides: Any) -> "WikiClient":
        w = cfg["wikipedia"]
        kwargs = dict(
            api_url=w["api_url"],
            wikidata_api_url=w["wikidata_api_url"],
            user_agent=w["user_agent"],
            cache_dir=cfg.api_cache_dir(),
            min_interval_seconds=w.get("min_interval_seconds", 0.5),
            max_retries=w.get("max_retries", 5),
            backoff_seconds=w.get("backoff_seconds", 2.0),
            timeout_seconds=w.get("timeout_seconds", 30),
        )
        kwargs.update(overrides)
        return cls(**kwargs)

    def _cache_file(self, url: str, params: dict[str, Any]) -> Path | None:
        if not self.cache_dir:
            return None
        key = hashlib.sha1((url + "?" + json.dumps(params, sort_keys=True)).encode()).hexdigest()
        return self.cache_dir / key[:2] / f"{key}.json"

    def _throttle(self) -> None:
        wait = self.min_interval - (time.monotonic() - self._last_request)
        if wait > 0:
            time.sleep(wait)
        self._last_request = time.monotonic()

    def get(self, url: str, params: dict[str, Any]) -> dict[str, Any]:
        params = {"format": "json", "formatversion": 2, "maxlag": 5, **params}
        cache = self._cache_file(url, params)
        if cache and cache.exists():
            return json.loads(cache.read_text())
        if self.offline:
            raise RuntimeError(f"offline mode and no cached response for {params}")
        last_error: Exception | None = None
        for attempt in range(self.max_retries + 1):
            self._throttle()
            try:
                self.network_calls += 1
                resp = self.session.get(url, params=params, timeout=self.timeout)
                if resp.status_code == 429 or resp.status_code >= 500:
                    retry_after = resp.headers.get("Retry-After")
                    delay = float(retry_after) if retry_after and retry_after.isdigit() else self.backoff * 2**attempt
                    last_error = requests.HTTPError(f"HTTP {resp.status_code}")
                    time.sleep(delay)
                    continue
                resp.raise_for_status()
                data = resp.json()
                if "error" in data:
                    if data["error"].get("code") == "maxlag":
                        last_error = RuntimeError("maxlag")
                        time.sleep(self.backoff * 2**attempt)
                        continue
                    raise RuntimeError(f"API error: {data['error']}")
                if cache:
                    cache.parent.mkdir(parents=True, exist_ok=True)
                    atomic_write_text(cache, json.dumps(data, ensure_ascii=False))
                return data
            except (requests.ConnectionError, requests.Timeout) as exc:
                last_error = exc
                time.sleep(self.backoff * 2**attempt)
        raise RuntimeError(f"request failed after {self.max_retries + 1} attempts: {last_error}")

    def query_titles(self, titles: Iterable[str], batch_size: int = 50) -> dict[str, PageInfo]:
        titles = list(dict.fromkeys(t for t in titles if t))
        out: dict[str, PageInfo] = {}
        for i in range(0, len(titles), batch_size):
            batch = titles[i : i + batch_size]
            data = self.get(
                self.api_url,
                {"action": "query", "titles": "|".join(batch), "redirects": 1, "prop": "pageprops", "ppprop": "wikibase_item|disambiguation"},
            )
            out.update(parse_title_query(data, batch))
        return out

    def entities(self, qids: Iterable[str], batch_size: int = 50) -> dict[str, dict[str, Any]]:
        qids = sorted({q for q in qids if q})
        out: dict[str, dict[str, Any]] = {}
        for i in range(0, len(qids), batch_size):
            batch = qids[i : i + batch_size]
            data = self.get(
                self.wikidata_api_url,
                {"action": "wbgetentities", "ids": "|".join(batch), "props": "claims|sitelinks", "sitefilter": "enwiki"},
            )
            out.update(parse_entities(data))
        return out

    def search_taxon(self, name: str, limit: int = 5) -> list[str]:
        data = self.get(
            self.wikidata_api_url,
            {"action": "query", "list": "search", "srsearch": f'haswbstatement:"P225={name}"', "srlimit": limit, "srnamespace": 0},
        )
        return [hit["title"] for hit in data.get("query", {}).get("search", [])]

    def fetch_page(self, title: str) -> dict[str, Any]:
        return self.get(
            self.api_url,
            {
                "action": "query",
                "titles": title,
                "redirects": 1,
                "prop": "extracts|revisions|pageprops",
                "explaintext": 1,
                "exsectionformat": "wiki",
                "rvprop": "ids|timestamp",
                "ppprop": "wikibase_item",
            },
        )


def parse_title_query(data: dict[str, Any], requested: Iterable[str]) -> dict[str, PageInfo]:
    q = data.get("query", {})
    normalized = {n["from"]: n["to"] for n in q.get("normalized", [])}
    redirects = {r["from"]: (r["to"], r.get("tofragment")) for r in q.get("redirects", [])}
    pages = {p["title"]: p for p in q.get("pages", [])}
    out = {}
    for title in requested:
        cur, fragment, redirected, seen = normalized.get(title, title), None, False, set()
        while cur in redirects and cur not in seen:
            seen.add(cur)
            cur, frag = redirects[cur]
            fragment = frag or fragment
            redirected = True
        page = pages.get(cur, {})
        props = page.get("pageprops", {}) or {}
        out[title] = PageInfo(
            query=title,
            title=cur,
            exists=bool(page) and not page.get("missing") and not page.get("invalid"),
            redirected=redirected,
            fragment=fragment,
            wikibase_item=props.get("wikibase_item"),
            disambiguation="disambiguation" in props,
        )
    return out


def _claim_values(claims: dict[str, Any], prop: str) -> list[Any]:
    values = []
    for claim in claims.get(prop, []):
        snak = claim.get("mainsnak", {})
        if snak.get("snaktype") != "value":
            continue
        v = snak.get("datavalue", {}).get("value")
        values.append(v.get("id") if isinstance(v, dict) and "id" in v else v)
    return values


def parse_entities(data: dict[str, Any]) -> dict[str, dict[str, Any]]:
    out = {}
    for qid, ent in data.get("entities", {}).items():
        if "missing" in ent:
            continue
        claims = ent.get("claims", {})
        sitelinks = ent.get("sitelinks", {})
        out[qid] = {
            "taxon_names": [v for v in _claim_values(claims, "P225") if isinstance(v, str)],
            "ranks": _claim_values(claims, "P105"),
            "enwiki": sitelinks.get("enwiki", {}).get("title"),
        }
    return out


def classify_article(
    species: str,
    other_names: set[str],
    taxon_names: Iterable[str],
    ranks: Iterable[str] = (),
    query_source: str | None = None,
    redirected: bool = False,
) -> str | None:
    parsed = [(normalize_name(t), elevated_name(t), len(name_tokens(t)) >= 3) for t in taxon_names]
    if any(b == species and not tri for b, _, tri in parsed):
        return "exact"
    if any(e == species and tri for _, e, tri in parsed):
        return "subspecies_article"
    if any(b in other_names and not tri for b, _, tri in parsed):
        return "synonym"
    if query_source not in TRUSTED_QUERY_SOURCES:
        return None
    species_rank = SPECIES_RANK_QID in set(ranks)
    stem = epithet_stem(split_binomial(species)[1])
    if species_rank and any(b and not tri and epithet_stem(split_binomial(b)[1]) == stem for b, _, tri in parsed):
        return "genus_transfer"
    if species_rank and redirected and query_source == "scientific_name" and parsed and not any(tri for _, _, tri in parsed):
        return "name_redirect"
    return None


def _settle_conflicts(results: dict[str, dict[str, Any]], banned: dict[str, set[str]]) -> list[str]:
    by_title: dict[str, list[str]] = defaultdict(list)
    for sp, r in results.items():
        by_title[r["wiki_title"]].append(sp)
    losers = []
    for title, sps in by_title.items():
        if len(sps) < 2:
            continue
        ranked = sorted(sps, key=lambda s: (-VALIDATION_GRADE[results[s]["wiki_validation"]], results[s]["wiki_rank"]))
        top = ranked[0]
        tie = [s for s in ranked[1:] if VALIDATION_GRADE[results[s]["wiki_validation"]] == VALIDATION_GRADE[results[top]["wiki_validation"]]]
        losers_here = ranked if tie else ranked[1:]
        for sp in losers_here:
            banned[sp].add(title)
            losers.append(sp)
    for sp in losers:
        results.pop(sp, None)
    return losers


def resolve_species_titles(
    candidates: dict[str, list[tuple[str, str]]],
    other_names: dict[str, set[str]],
    client: WikiClient,
    search_names: dict[str, list[str]] | None = None,
) -> tuple[dict[str, dict[str, Any]], dict[str, str]]:
    queues = {sp: list(enumerate(c)) for sp, c in candidates.items()}
    results: dict[str, dict[str, Any]] = {}
    banned: dict[str, set[str]] = defaultdict(set)
    notes: dict[str, str] = {}

    def try_batch(batch: dict[str, tuple[int, tuple[str, str]]]) -> None:
        infos = client.query_titles(title for _, (title, _) in batch.values())
        ents = client.entities(i.wikibase_item for i in infos.values() if i.exists and not i.disambiguation)
        for sp, (rank, (title, source)) in batch.items():
            info = infos.get(title)
            if info is None or not info.exists:
                continue
            if info.disambiguation:
                notes.setdefault(sp, f"{title!r} is a disambiguation page")
                continue
            if info.title in banned[sp]:
                continue
            ent = ents.get(info.wikibase_item or "", {})
            verdict = classify_article(
                sp, other_names.get(sp, set()), ent.get("taxon_names", []), ent.get("ranks", []), source, info.redirected
            )
            if verdict is None:
                notes.setdefault(sp, f"{title!r} -> {info.title!r} describes {ent.get('taxon_names') or 'no taxon'}")
                continue
            results[sp] = {
                "wiki_title": info.title,
                "wiki_qid": info.wikibase_item,
                "wiki_query": title,
                "wiki_source": source,
                "wiki_validation": verdict,
                "wiki_rank": rank,
            }

    while True:
        batch = {}
        for sp, q in queues.items():
            if sp in results:
                continue
            while q and q[0][1][0] in banned[sp]:
                q.pop(0)
            if q:
                batch[sp] = q.pop(0)
        if not batch:
            break
        try_batch(batch)
        _settle_conflicts(results, banned)

    unresolved = [sp for sp in candidates if sp not in results]
    if unresolved:
        found: dict[str, list[str]] = {}
        for sp in unresolved:
            names = (search_names or {}).get(sp) or [sp]
            qids: list[str] = []
            for name in names:
                qids.extend(client.search_taxon(name))
            found[sp] = list(dict.fromkeys(qids))
        ents = client.entities(q for qs in found.values() for q in qs)
        canonical = client.query_titles(e["enwiki"] for e in ents.values() if e.get("enwiki"))
        for sp, qids in found.items():
            for qid in qids:
                ent = ents.get(qid, {})
                info = canonical.get(ent.get("enwiki") or "")
                if info is None or not info.exists or info.disambiguation:
                    continue
                title = info.title
                if title in banned[sp]:
                    continue
                verdict = classify_article(sp, other_names.get(sp, set()), ent.get("taxon_names", []))
                if verdict:
                    results[sp] = {
                        "wiki_title": title,
                        "wiki_qid": qid,
                        "wiki_query": f"wikidata:{qid}",
                        "wiki_source": "wikidata_search",
                        "wiki_validation": verdict,
                        "wiki_rank": 10_000,
                    }
                    break
            else:
                if not qids:
                    notes.setdefault(sp, "no Wikidata item with this taxon name")
                elif all(not ents.get(q, {}).get("enwiki") for q in qids):
                    notes[sp] = f"Wikidata item(s) {qids} exist but have no English Wikipedia article"
        _settle_conflicts(results, banned)
    for sp in results:
        notes.pop(sp, None)
    return results, notes


HEADING_RE = re.compile(r"^(={2,6})\s*(.*?)\s*\1\s*$", re.M)


@dataclass
class Section:
    level: int
    heading: str
    body: str


def split_sections(text: str) -> list[Section]:
    sections = []
    matches = list(HEADING_RE.finditer(text))
    lead_end = matches[0].start() if matches else len(text)
    sections.append(Section(0, "", text[:lead_end]))
    for i, m in enumerate(matches):
        end = matches[i + 1].start() if i + 1 < len(matches) else len(text)
        sections.append(Section(len(m.group(1)), m.group(2).strip(), text[m.end() : end]))
    return sections


def find_description(sections: list[Section], headings: Iterable[str]) -> tuple[str | None, str]:
    wanted = [h.lower() for h in headings]
    body_sections = sections[1:]
    for target in wanted:
        for idx, sec in enumerate(body_sections):
            if sec.heading.lower() != target:
                continue
            parts = [sec.body]
            for sub in body_sections[idx + 1 :]:
                if sub.level <= sec.level:
                    break
                parts.append(sub.body)
            return sec.heading, "\n".join(parts)
    return None, ""


CITATION_RE = re.compile(
    r"\[(?:\d+|[a-z]|[ivx]+|note \d+|nb \d+|citation needed|clarification needed|dubious[^\]]*|when\?|who\?|which\?|verification needed|failed verification|page needed)\]",
    re.I,
)


def clean_text(text: str) -> str:
    t = unicodedata.normalize("NFC", text)
    t = re.sub(r"<!--.*?-->", " ", t, flags=re.S)
    t = re.sub(r"<ref[^>]*/>|<ref[^>]*>.*?</ref>", " ", t, flags=re.S | re.I)
    t = re.sub(r"<[^>]+>", " ", t)
    prev = None
    while prev != t:
        prev = t
        t = re.sub(r"\{\{[^{}]*\}\}", " ", t)
    t = re.sub(r"\{\\displaystyle[^}]*\}", " ", t)
    t = CITATION_RE.sub("", t)
    t = re.sub(r"\[\[(?:[^|\]]*\|)?([^\]]*)\]\]", r"\1", t)
    t = re.sub(r"'{2,}", "", t)
    t = re.sub(r"\(\s*listen\s*\)", "", t, flags=re.I)
    t = re.sub(r"\(\s*[;,:/\s]*\)", "", t)
    t = t.replace(" ", " ")
    t = re.sub(r"[ \t]+", " ", t)
    t = re.sub(r" +([,.;:)])", r"\1", t)
    t = re.sub(r"\( +", "(", t)
    t = "\n".join(line.strip() for line in t.split("\n"))
    t = re.sub(r"\n{3,}", "\n\n", t)
    return t.strip()


def page_from_response(data: dict[str, Any]) -> dict[str, Any]:
    pages = data.get("query", {}).get("pages", [])
    if not pages or pages[0].get("missing"):
        raise ValueError("page missing from response")
    page = pages[0]
    rev = (page.get("revisions") or [{}])[0]
    return {
        "title": page.get("title"),
        "pageid": page.get("pageid"),
        "revid": rev.get("revid"),
        "rev_timestamp": rev.get("timestamp"),
        "wikibase_item": (page.get("pageprops") or {}).get("wikibase_item"),
        "extract": page.get("extract") or "",
    }


def preprocess(extract: str, headings: Iterable[str]) -> dict[str, Any]:
    sections = split_sections(extract)
    lead = clean_text(sections[0].body)
    heading, desc_raw = find_description(sections, headings)
    description = clean_text(desc_raw)
    text = "\n\n".join(p for p in (lead, description) if p)
    return {
        "lead": lead,
        "description": description,
        "description_heading": heading,
        "has_description": bool(description),
        "text": text,
        "n_chars": len(text),
        "n_words": len(text.split()),
    }
