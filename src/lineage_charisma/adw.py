from __future__ import annotations

import time
from dataclasses import dataclass, field
from html.parser import HTMLParser
from pathlib import Path
from typing import Iterable

import requests

from lineage_charisma.io_utils import read_json, write_json
from lineage_charisma.taxonomy import elevated_name, normalize_name
from lineage_charisma.wiki import clean_text

# docs/preregistration_amendment_secondary.md, section A and Change 1
BASE_URL = "https://animaldiversity.org/accounts/{}/"
CRAWL_DELAY_SECONDS = 8.0
SECTIONS = ("geographic_range", "habitat", "physical_description", "behavior", "food_habits")   # in this order
SKIP_TAGS = {"refs", "ref", "figure", "figcaption", "table", "aside", "script", "style", "button", "nav"}


@dataclass
class ParsedAccount:
    title: str = ""
    rank: str = ""                      # from the h1 class, e.g. "Species" or "Genus"
    not_found: bool = False
    sections: dict[str, str] = field(default_factory=dict)

    @property
    def is_species_account(self) -> bool:
        return self.rank == "Species" and not self.not_found


class _AccountParser(HTMLParser):
    """Collects the page's h1 (title and rank) and the paragraph text of each <section> keyed by its h2 id."""

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.out = ParsedAccount()
        self._in_h1 = self._in_title = False
        self._section: str | None = None
        self._in_p = 0
        self._skip = 0
        self._buf: list[str] = []
        self._title_buf: list[str] = []

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if tag in SKIP_TAGS:
            self._skip += 1
        elif tag == "title":
            self._in_title = True
        elif tag == "h1" and "rank--" in (a.get("class") or ""):
            self._in_h1 = True
            self.out.rank = a["class"].split("rank--", 1)[1].split()[0]
        elif tag == "section":
            self._section = None
        elif tag == "h2" and a.get("id"):
            self._section = a["id"]
            self.out.sections.setdefault(self._section, "")
        elif tag == "p" and self._section is not None:
            self._in_p += 1
            self._buf = []

    def handle_endtag(self, tag):
        if tag in SKIP_TAGS:
            self._skip = max(0, self._skip - 1)
        elif tag == "title":
            self._in_title = False
        elif tag == "h1":
            self._in_h1 = False
        elif tag == "p" and self._in_p:
            self._in_p -= 1
            text = " ".join("".join(self._buf).split())
            if text:
                prev = self.out.sections[self._section]
                self.out.sections[self._section] = f"{prev}\n{text}" if prev else text
        elif tag == "section":
            self._section = None

    def handle_data(self, data):
        if self._in_title and "Page not found" in data:
            self.out.not_found = True
        if self._in_h1:
            self.out.title += data
        elif self._in_p and not self._skip:
            self._buf.append(data)


def parse_account(html: str) -> ParsedAccount:
    parser = _AccountParser()
    parser.feed(html)
    parser.out.title = " ".join(parser.out.title.split())
    return parser.out


def extract_text(account: ParsedAccount) -> str:
    """The registered text rule: paragraph text of the five sections, in order, through the Wikipedia cleaning step."""
    parts = [account.sections[s] for s in SECTIONS if account.sections.get(s)]
    return clean_text("\n\n".join(parts))


def query_names(species: str, tree_binomial: str | None, synonym_combinations: Iterable[str]) -> list[str]:
    """Names to try, in the registered order: current MDD name, Upham tip binomial, MDD synonyms in MDD order."""
    names: list[str] = []
    for n in [species, normalize_name(tree_binomial) if tree_binomial else None, *(elevated_name(c) for c in synonym_combinations)]:
        if n and n not in names:
            names.append(n)
    return names


class ADWClient:
    """Fetches ADW pages at most once per CRAWL_DELAY_SECONDS, caching every response (hits and misses) on disk."""

    def __init__(self, cache_dir: Path, user_agent: str, delay: float = CRAWL_DELAY_SECONDS, offline: bool = False):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.delay, self.offline = delay, offline
        self.session = requests.Session()
        self.session.headers.update({"User-Agent": user_agent})
        self._last = 0.0
        self.requests_made = 0

    def _paths(self, name: str) -> tuple[Path, Path]:
        slug = name.replace(" ", "_")
        return self.cache_dir / f"{slug}.html", self.cache_dir / f"{slug}.json"

    def fetch(self, name: str) -> tuple[str, dict]:
        html_path, meta_path = self._paths(name)
        if html_path.exists() and meta_path.exists():
            return html_path.read_text(encoding="utf-8"), read_json(meta_path)
        if self.offline:
            raise RuntimeError(f"offline: no cached ADW page for {name}")
        wait = self.delay - (time.monotonic() - self._last)
        if wait > 0:
            time.sleep(wait)
        url = BASE_URL.format(name.replace(" ", "_"))
        for attempt in range(4):
            try:
                resp = self.session.get(url, timeout=60)
                self._last = time.monotonic()
                if resp.status_code in (429, 500, 502, 503, 504):
                    time.sleep(max(self.delay, float(resp.headers.get("Retry-After", 0) or 0)) * (attempt + 1))
                    continue
                break
            except requests.RequestException:
                self._last = time.monotonic()
                if attempt == 3:
                    raise
                time.sleep(self.delay * (attempt + 1))
        self.requests_made += 1
        meta = {"name": name, "url": url, "final_url": resp.url, "status": resp.status_code, "fetched_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
        html_path.write_text(resp.text, encoding="utf-8")
        write_json(meta_path, meta)
        return resp.text, meta
