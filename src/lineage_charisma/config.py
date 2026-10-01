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

    def path(self, key: str) -> Path:
        p = Path(self.data["paths"][key])
        p = p if p.is_absolute() else self.root / p
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
        p = self.interim / "wiki"
        p.mkdir(parents=True, exist_ok=True)
        return p

    def api_cache_dir(self) -> Path:
        p = self.interim / "api_cache"
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
        p = self.reports / "figures"
        p.mkdir(parents=True, exist_ok=True)
        return p


def load_config(path: str | os.PathLike | None = None) -> Config:
    cfg_path = Path(path) if path else Path(os.environ.get("LC_CONFIG", DEFAULT_CONFIG))
    with open(cfg_path) as fh:
        data = yaml.safe_load(fh)
    return Config(data=data, root=cfg_path.resolve().parent)


def base_parser(description: str):
    import argparse

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--config", default=None, help="path to config.yaml")
    parser.add_argument("--force", action="store_true", help="rebuild outputs even if cached")
    return parser
