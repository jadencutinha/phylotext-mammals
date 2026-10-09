from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path
from typing import Any, Iterable, Sequence
from urllib.parse import quote

import numpy as np
import pandas as pd
import requests

from lineage_charisma.io_utils import atomic_write_text


def month_range(start: str, end: str) -> list[str]:
    """Every month from `start` to `end` inclusive, as YYYY-MM."""
    return [str(p) for p in pd.period_range(start, end, freq="M")]


class PageviewsClient:
    """Monthly per-article pageviews from the Wikimedia REST API, cached on disk, rate-limited, retried with backoff."""

    def __init__(
        self,
        api_url: str,
        user_agent: str,
        project: str = "en.wikipedia",
        access: str = "all-access",
        agent: str = "user",
        cache_dir: Path | None = None,
        min_interval_seconds: float = 0.1,
        max_retries: int = 5,
        backoff_seconds: float = 2.0,
        timeout_seconds: float = 30,
        session: requests.Session | None = None,
        offline: bool = False,
    ):
        self.api_url, self.project, self.access, self.agent = api_url.rstrip("/"), project, access, agent
        self.cache_dir = Path(cache_dir) if cache_dir else None
        self.min_interval, self.max_retries, self.backoff, self.timeout = min_interval_seconds, max_retries, backoff_seconds, timeout_seconds
        self.offline = offline
        self.session = session or requests.Session()
        self.session.headers.update({"User-Agent": user_agent, "Accept-Encoding": "gzip"})
        self._last_request = 0.0
        self.network_calls = 0

    @classmethod
    def from_config(cls, cfg, **overrides: Any) -> "PageviewsClient":
        a, w = cfg["attention"], cfg["wikipedia"]
        kwargs = dict(
            api_url=a["pageviews_api_url"], user_agent=w["user_agent"], project=a["project"], access=a["access"], agent=a["agent"],
            cache_dir=cfg.api_cache_dir() / "pageviews", min_interval_seconds=a.get("min_interval_seconds", 0.1),
            max_retries=w.get("max_retries", 5), backoff_seconds=w.get("backoff_seconds", 2.0), timeout_seconds=w.get("timeout_seconds", 30),
        )
        kwargs.update(overrides)
        return cls(**kwargs)

    def url(self, title: str, start: str, end: str) -> str:
        first = pd.Period(start, freq="M").strftime("%Y%m0100")
        last = pd.Period(end, freq="M").end_time.strftime("%Y%m%d00")
        article = quote(title.replace(" ", "_"), safe="")
        return f"{self.api_url}/{self.project}/{self.access}/{self.agent}/{article}/monthly/{first}/{last}"

    def _cache_file(self, url: str) -> Path | None:
        if not self.cache_dir:
            return None
        key = hashlib.sha1(url.encode()).hexdigest()
        return self.cache_dir / key[:2] / f"{key}.json"

    def _throttle(self) -> None:
        wait = self.min_interval - (time.monotonic() - self._last_request)
        if wait > 0:
            time.sleep(wait)
        self._last_request = time.monotonic()

    def monthly(self, title: str, start: str, end: str) -> dict[str, int]:
        """Views per month (YYYY-MM) for one title. A title with no views in the window gives an empty dict."""
        url = self.url(title, start, end)
        cache = self._cache_file(url)
        if cache and cache.exists():
            return json.loads(cache.read_text())["views"]
        if self.offline:
            raise RuntimeError(f"offline mode and no cached pageviews for {title!r}")
        last_error: Exception | None = None
        for attempt in range(self.max_retries + 1):
            self._throttle()
            try:
                self.network_calls += 1
                resp = self.session.get(url, timeout=self.timeout)
                if resp.status_code == 404:      # the API answers 404 when a title has no recorded views in the window
                    views: dict[str, int] = {}
                elif resp.status_code == 429 or resp.status_code >= 500:
                    retry_after = resp.headers.get("Retry-After")
                    last_error = requests.HTTPError(f"HTTP {resp.status_code}")
                    time.sleep(float(retry_after) if retry_after and retry_after.isdigit() else self.backoff * 2**attempt)
                    continue
                else:
                    resp.raise_for_status()
                    views = {f"{item['timestamp'][:4]}-{item['timestamp'][4:6]}": int(item["views"]) for item in resp.json()["items"]}
                if cache:
                    cache.parent.mkdir(parents=True, exist_ok=True)
                    atomic_write_text(cache, json.dumps({"title": title, "url": url, "views": views}, ensure_ascii=False))
                return views
            except (requests.ConnectionError, requests.Timeout) as exc:
                last_error = exc
                time.sleep(self.backoff * 2**attempt)
        raise RuntimeError(f"pageviews request for {title!r} failed after {self.max_retries + 1} attempts: {last_error}")


def resolve_articles(wiki_client, titles: Iterable[str]) -> dict[str, tuple[str, list[str]]]:
    """For each title, the article it now resolves to and every main-namespace redirect pointing at that article.

    Pageviews are recorded under the title a reader asked for, so an article that was renamed has its
    earlier views under the old title, which is now a redirect.
    """
    titles = list(dict.fromkeys(titles))
    infos = wiki_client.query_titles(titles)
    out: dict[str, tuple[str, list[str]]] = {}
    for title in titles:
        info = infos.get(title)
        if info is None or not info.exists:
            raise RuntimeError(f"Wikipedia article {title!r} no longer exists")
        redirects, cont = [], {}
        while True:
            data = wiki_client.get(wiki_client.api_url, {"action": "query", "titles": info.title, "prop": "redirects", "rdnamespace": 0,
                                                         "rdlimit": "max", "rdprop": "title", **cont})
            for page in data.get("query", {}).get("pages", []):
                redirects += [r["title"] for r in page.get("redirects", [])]
            if "continue" not in data:
                break
            cont = data["continue"]
        out[title] = (info.title, sorted(set(redirects) - {info.title}))
    return out


def summed_monthly_views(client: PageviewsClient, titles: Sequence[str], start: str, end: str) -> pd.Series:
    """Views per month summed over `titles`, with 0 for months in the window that have none."""
    total = pd.Series(0, index=month_range(start, end), dtype=np.int64)
    for title in titles:
        for month, views in client.monthly(title, start, end).items():
            if month in total.index:
                total[month] += views
    return total


def attention_variables(median_pageviews: Sequence[float], n_tokens: Sequence[float]) -> pd.DataFrame:
    views, tokens = np.asarray(median_pageviews, dtype=np.float64), np.asarray(n_tokens, dtype=np.float64)
    if (views < 0).any() or (tokens <= 0).any() or not (np.isfinite(views).all() and np.isfinite(tokens).all()):
        raise ValueError("pageviews must be finite and non-negative, and token counts finite and positive")
    return pd.DataFrame({"log_pageviews": np.log1p(views), "log_length": np.log(tokens)})


def abs_difference_matrix(values: Sequence[float]) -> np.ndarray:
    v = np.asarray(values, dtype=np.float64)
    return np.abs(v[:, None] - v[None, :])


def standardize(values: Sequence[float]) -> np.ndarray:
    """Mean 0 and standard deviation 1 (population standard deviation)."""
    v = np.asarray(values, dtype=np.float64)
    sd = v.std()
    if sd == 0:
        raise ValueError("cannot standardize a constant variable")
    return (v - v.mean()) / sd


def combined_distance_matrix(*variables: Sequence[float]) -> np.ndarray:
    """Euclidean distance between species on the standardized variables."""
    z = np.column_stack([standardize(v) for v in variables])
    dist = np.sqrt(((z[:, None, :] - z[None, :, :]) ** 2).sum(axis=-1))
    np.fill_diagonal(dist, 0.0)
    return dist


def pairwise_mean_matrix(values: Sequence[float]) -> np.ndarray:
    """Mean of the two species' values for each pair. High when both are high; this is not a distance (its diagonal is the species' own value)."""
    v = np.asarray(values, dtype=np.float64)
    return (v[:, None] + v[None, :]) / 2.0


def attention_matrices(variables: pd.DataFrame) -> dict[str, np.ndarray]:
    return {
        "attention_pageviews": abs_difference_matrix(variables["log_pageviews"]),
        "attention_length": abs_difference_matrix(variables["log_length"]),
        "attention_combined": combined_distance_matrix(variables["log_pageviews"], variables["log_length"]),
    }
