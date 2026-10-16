from __future__ import annotations

from typing import Sequence

import numpy as np

from lineage_charisma.config import Config
from lineage_charisma.distances import load_matrix, phylo_matrix_path, text_matrix_path
from lineage_charisma.io_utils import read_csv
from lineage_charisma.stats import mantel, mrm, partial_mantel, submatrix
from lineage_charisma.traits import TRAITS

# every name and model below is fixed in docs/preregistration_amendment_h2.md
ECOLOGY_FILES = {"ecology": "ecology.npy", "ecology_unweighted": "ecology_unweighted.npy"}
ATTENTION_FILES = {"length": "attention_length.npy", "pageviews": "attention_pageviews.npy"}
PRIMARY_CONTROLS = ("ecology", "length")
MODELS = {
    "primary": PRIMARY_CONTROLS,
    "S1": ("ecology",),
    "S2": ("length",),
    "S3": ("pageviews",),
    "S4": ("ecology", "length", "pageviews"),
}
MRM_PREDICTORS = ("phylogeny", "ecology", "length", "pageviews")
CONFOUNDED_SHARE = 0.5


def load_h2_matrices(cfg: Config) -> tuple[dict[str, np.ndarray], list[str], list[str]]:
    """Phylogeny and every control matrix cut to the species with ecology data: (matrices, those species, the full order)."""
    out_dir, order_path = cfg.matrices_dir(), cfg.species_order_path()
    absent = [f for f in (*ECOLOGY_FILES.values(), *ATTENTION_FILES.values()) if not (out_dir / f).exists()]
    if absent:
        raise SystemExit(f"missing matrices {absent}; run make ecology and make attention first")
    matrices, species = {}, None
    for name, filename in ECOLOGY_FILES.items():
        matrices[name], eco_order = load_matrix(out_dir / filename, out_dir / "ecology_species.txt")
        if species is not None and eco_order != species:
            raise ValueError(f"{filename} is in a different species order from the other ecology matrix")
        species = eco_order
    phylo, order = load_matrix(phylo_matrix_path(out_dir), order_path)
    if species != [s for s in order if s in set(species)]:
        raise ValueError("ecology species are not in the order of species_order.txt")
    matrices["phylogeny"] = submatrix(phylo, order, species)
    for name, filename in ATTENTION_FILES.items():
        full, _ = load_matrix(out_dir / filename, order_path)
        matrices[name] = submatrix(full, order, species)
    return matrices, species, order


def load_text(cfg: Config, combo: tuple[str, str, str], species: Sequence[str] | None = None) -> np.ndarray:
    text, order = load_matrix(text_matrix_path(cfg.matrices_dir(), *combo), cfg.species_order_path())
    return text if species is None else submatrix(text, order, species)


def complete_trait_species(cfg: Config) -> list[str]:
    """Species with no missing value in any trait column, in the order of species_order.txt."""
    table = read_csv(cfg.processed / "traits.csv")
    return table.loc[table[[t.column for t in TRAITS]].notna().all(axis=1), "species_id"].tolist()


def cut(matrices: dict[str, np.ndarray], species: Sequence[str], keep: Sequence[str]) -> dict[str, np.ndarray]:
    return {name: submatrix(m, species, keep) for name, m in matrices.items()}


def share_of_raw(raw_r: float, partial_r: float) -> float:
    """Partial r as a share of raw r; not defined when the raw r is at or below zero."""
    return partial_r / raw_r if raw_r > 0 else float("nan")


def decision(raw_r: float, partial_r: float, partial_p: float, alpha: float) -> dict[str, bool | None]:
    """The two pre-registered rules: (a) the lineage signal survives, (b) the signal is largely confounded."""
    share = share_of_raw(raw_r, partial_r)
    return {"survives": bool(partial_r > 0 and partial_p < alpha), "confounded": None if np.isnan(share) else bool(share < CONFOUNDED_SHARE)}


def partial_result(text: np.ndarray, phylo: np.ndarray, controls: Sequence[np.ndarray], method: str, permutations: int, seed: int) -> dict[str, float]:
    """Raw and partial text-phylogeny correlation on the same species. Both tests reorder the phylogenetic matrix."""
    raw = mantel(text, phylo, method=method, permutations=permutations, seed=seed)
    partial = partial_mantel(phylo, text, list(controls), method=method, permutations=permutations, seed=seed)
    return {"raw_r": raw.r, "raw_p": raw.p, "partial_r": partial.r, "partial_p": partial.p, "share_of_raw": share_of_raw(raw.r, partial.r),
            "change": raw.r - partial.r, "n_species": raw.n, "n_pairs": raw.n_pairs}


def mrm_rows(text: np.ndarray, matrices: dict[str, np.ndarray], method: str, permutations: int, seed: int,
             predictors: Sequence[str] = MRM_PREDICTORS) -> list[dict]:
    fit = mrm(text, {name: matrices[name] for name in predictors}, method=method, permutations=permutations, seed=seed)
    return [{"predictor": name, "coef_std": fit.coef[j], "coef_raw": fit.coef_raw[j], "t": fit.t[j], "p": fit.p[j], "r2": fit.r2, "p_r2": fit.p_r2,
             "n_species": fit.n, "n_pairs": fit.n_pairs} for j, name in enumerate(fit.names)]
