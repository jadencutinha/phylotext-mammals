from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = PROJECT_ROOT / "config.yaml"


@dataclass(frozen=True)
class Config:
    data: dict[str, Any]
    root: Path = field(default=PROJECT_ROOT)

    def __getitem__(self, key: str) -> Any:
        return self.data[key]

    def get(self, key: str, default: Any = None) -> Any:
        return self.data.get(key, default)

    @property
    def clade_rank(self) -> str | None:
        return (self.data.get("target_clade") or {}).get("rank")

    @property
    def clade_name(self) -> str | None:
        return (self.data.get("target_clade") or {}).get("name")

    @property
    def clade_slug(self) -> str:
        name = self.clade_name or "all"
        return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")

    @property
    def output_subdir(self) -> str | None:
        """Set by a `clades` overlay (e.g. "mammalia"): interim, processed and table outputs go one level down."""
        return self.data.get("output_subdir")

    def path(self, key: str, shared: bool = False) -> Path:
        p = Path(self.data["paths"][key])
        p = p if p.is_absolute() else self.root / p
        if self.output_subdir and not shared and key in ("interim", "processed"):
            p = p / self.output_subdir
        p.mkdir(parents=True, exist_ok=True)
        return p

    @property
    def raw(self) -> Path:
        return self.path("raw")

    @property
    def interim(self) -> Path:
        return self.path("interim")

    @property
    def processed(self) -> Path:
        return self.path("processed")

    def mdd_species_path(self) -> Path:
        return self.raw / f"mdd_species_{self.clade_slug}.csv"

    def mdd_synonyms_path(self) -> Path:
        return self.raw / f"mdd_synonyms_{self.clade_slug}.csv"

    def tree_path(self) -> Path:
        return self.raw / "tree" / self["sources"]["tree"]["mcc_filename"]

    def tips_path(self) -> Path:
        return self.interim / f"tree_tips_{self.clade_slug}.csv"

    def matched_path(self) -> Path:
        return self.interim / "matched.csv"

    def residual_path(self) -> Path:
        return self.interim / "residual_review.csv"

    def excluded_species_path(self) -> Path:
        return self.interim / "excluded_species.csv"

    def unused_tips_path(self) -> Path:
        return self.interim / "unused_tree_tips.csv"

    def wiki_dir(self) -> Path:
        p = self.path("interim", shared=True) / "wiki"     # shared across clades
        p.mkdir(parents=True, exist_ok=True)
        return p

    def api_cache_dir(self) -> Path:
        p = self.path("interim", shared=True) / "api_cache"   # shared across clades
        p.mkdir(parents=True, exist_ok=True)
        return p

    def pruned_tree_path(self) -> Path:
        return self.interim / "tree_pruned_mdd_names.nwk"

    def stub_species_path(self) -> Path:
        return self.interim / "stub_species.csv"

    def corpus_path(self) -> Path:
        return self.processed / "corpus.parquet"

    def embeddings_dir(self) -> Path:
        p = self.processed / "embeddings"
        p.mkdir(parents=True, exist_ok=True)
        return p

    def matrices_dir(self, min_tokens: int | None = None) -> Path:
        p = self.processed / ("matrices" if min_tokens is None else f"matrices_min{int(min_tokens)}")
        p.mkdir(parents=True, exist_ok=True)
        return p

    def species_order_path(self, min_tokens: int | None = None) -> Path:
        return self.matrices_dir(min_tokens) / "species_order.txt"

    def mask_terms_path(self) -> Path:
        return self.interim / "mask_terms.csv"

    @property
    def mask_levels(self) -> list[str]:
        return list((self.data.get("masking") or {}).get("levels") or ["none"])

    def combos(self) -> list[tuple[str, str, str]]:
        """Every (model, rule, mask level); the first is the primary combination."""
        emb = self.data["embedding"]
        return [(m, r, k) for k in self.mask_levels for m in emb["models"] for r in emb["rules"]]

    @property
    def reports(self) -> Path:
        return self.path("reports")

    def figures_dir(self) -> Path:
        p = self.reports / "figures" / (self.output_subdir or "")
        p.mkdir(parents=True, exist_ok=True)
        return p

    def week4_figures_dir(self) -> Path:
        # a subdirectory, because scripts/08_sanity_report.py clears reports/figures/*.png on every run
        p = self.figures_dir() / "week4"
        p.mkdir(parents=True, exist_ok=True)
        return p

    def week5_figures_dir(self) -> Path:
        p = self.figures_dir() / "week5"
        p.mkdir(parents=True, exist_ok=True)
        return p

    def week6_figures_dir(self) -> Path:
        p = self.reports / "figures" / "week6"
        p.mkdir(parents=True, exist_ok=True)
        return p

    def tables_dir(self) -> Path:
        p = self.reports / "tables" / (self.output_subdir or "")
        p.mkdir(parents=True, exist_ok=True)
        return p

    @property
    def primary_combo(self) -> tuple[str, str, str]:
        return self.combos()[0]


def _deep_merge(base: dict[str, Any], over: dict[str, Any]) -> dict[str, Any]:
    out = dict(base)
    for k, v in over.items():
        out[k] = _deep_merge(out[k], v) if isinstance(v, dict) and isinstance(out.get(k), dict) else v
    return out


def load_config(path: str | os.PathLike | None = None, clade: str | None = None) -> Config:
    """Load config.yaml; `clade` (or LC_CLADE) names an entry of `clades` that is deep-merged over the defaults."""
    cfg_path = Path(path) if path else Path(os.environ.get("LC_CONFIG", DEFAULT_CONFIG))
    with open(cfg_path) as fh:
        data = yaml.safe_load(fh)
    clade = clade or os.environ.get("LC_CLADE")
    overlays = data.pop("clades", None) or {}
    if clade:
        if clade not in overlays:
            raise KeyError(f"clade {clade!r} not in config clades: {sorted(overlays)}")
        data = _deep_merge(data, overlays[clade])
    return Config(data=data, root=cfg_path.resolve().parent)


def base_parser(description: str):
    import argparse

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--config", default=None, help="path to config.yaml")
    parser.add_argument("--force", action="store_true", help="rebuild outputs even if cached")
    parser.add_argument("--clade", default=None, help="entry of `clades` in config.yaml, e.g. mammalia (default: LC_CLADE, else the pilot)")
    return parser
