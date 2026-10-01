from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any, Sequence

import dendropy
import numpy as np

from lineage_charisma.embed import l2_normalize, model_slug, read_species_order
from lineage_charisma.io_utils import read_meta, slugify_species, write_meta
from lineage_charisma.phylo import cophenetic_matrix


class SpeciesOrderError(AssertionError):
    pass


def order_hash(species: Sequence[str]) -> str:
    return hashlib.sha256("\n".join(species).encode("utf-8")).hexdigest()


def canonical_order(species: Sequence[str]) -> list[str]:
    species = list(species)
    if len(set(species)) != len(species):
        raise SpeciesOrderError("duplicate species ids; cannot build a canonical order")
    return sorted(species)


def assert_same_order(actual: Sequence[str], expected: Sequence[str], what: str = "matrix") -> None:
    actual, expected = list(actual), list(expected)
    if actual == expected:
        return
    if set(actual) != set(expected):
        only_actual, only_expected = sorted(set(actual) - set(expected)), sorted(set(expected) - set(actual))
        raise SpeciesOrderError(
            f"{what}: species set differs from the canonical order "
            f"({len(only_actual)} extra, e.g. {only_actual[:3]}; {len(only_expected)} missing, e.g. {only_expected[:3]})"
        )
    first = next(i for i, (a, e) in enumerate(zip(actual, expected)) if a != e)
    raise SpeciesOrderError(f"{what}: same species but different order; first mismatch at row {first}: {actual[first]!r} != {expected[first]!r}")


def reindex_rows(values: np.ndarray, current: Sequence[str], target: Sequence[str], what: str = "array") -> np.ndarray:
    position = {s: i for i, s in enumerate(current)}
    if len(position) != len(current):
        raise SpeciesOrderError(f"{what}: duplicate species in its order")
    missing = [s for s in target if s not in position]
    if missing:
        raise SpeciesOrderError(f"{what}: {len(missing)} species in the canonical order are absent, e.g. {missing[:3]}")
    return values[[position[s] for s in target]]


def reindex_matrix(matrix: np.ndarray, current: Sequence[str], target: Sequence[str], what: str = "matrix") -> np.ndarray:
    if matrix.shape != (len(current), len(current)):
        raise SpeciesOrderError(f"{what}: shape {matrix.shape} does not match {len(current)} species")
    idx = reindex_rows(np.arange(len(current)), current, target, what)
    return matrix[np.ix_(idx, idx)]


def validate_distance_matrix(matrix: np.ndarray, name: str = "matrix", atol: float = 1e-8) -> None:
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"{name}: not square, shape {matrix.shape}")
    if not np.all(np.isfinite(matrix)):
        raise ValueError(f"{name}: contains {int((~np.isfinite(matrix)).sum())} NaN or infinite values")
    asym = float(np.abs(matrix - matrix.T).max())
    if asym > atol:
        raise ValueError(f"{name}: not symmetric (max |D - D^T| = {asym:.3g})")
    diag = float(np.abs(np.diag(matrix)).max())
    if diag > atol:
        raise ValueError(f"{name}: diagonal is not zero (max |D_ii| = {diag:.3g})")
    if matrix.min() < 0:
        raise ValueError(f"{name}: contains negative distances (min = {matrix.min():.3g})")


def cosine_distance_matrix(embeddings: np.ndarray, norm_atol: float = 1e-3) -> np.ndarray:
    emb = np.asarray(embeddings, dtype=np.float64)
    norms = np.linalg.norm(emb, axis=1)
    if not np.allclose(norms, 1.0, atol=norm_atol):
        raise ValueError(f"embeddings are not L2-normalized (norms range {norms.min():.4f}..{norms.max():.4f})")
    emb = l2_normalize(emb)
    dist = 1.0 - emb @ emb.T
    dist = (dist + dist.T) / 2.0
    np.fill_diagonal(dist, 0.0)
    # float rounding can leave identical vectors at about -1e-16
    return np.clip(dist, 0.0, 2.0)


def phylo_distance_matrix(tree: dendropy.Tree, species_order: Sequence[str]) -> np.ndarray:
    dist = cophenetic_matrix(tree)
    labels = [slugify_species(s) for s in species_order]
    if len(set(labels)) != len(labels):
        raise SpeciesOrderError("phylo: two species ids map to the same tip label")
    missing = [s for s, lab in zip(species_order, labels) if lab not in dist.index]
    if missing:
        raise SpeciesOrderError(f"phylo: {len(missing)} species are not tips of the tree, e.g. {missing[:3]}")
    return dist.loc[labels, labels].to_numpy(dtype=np.float64)


def tree_tip_order(tree: dendropy.Tree, species_order: Sequence[str]) -> list[str]:
    by_label = {slugify_species(s): s for s in species_order}
    return [by_label[leaf.taxon.label] for leaf in tree.leaf_node_iter() if leaf.taxon.label in by_label]


def save_matrix(path: Path, matrix: np.ndarray, species_order: Sequence[str], **meta: Any) -> Path:
    path = Path(path)
    if matrix.shape != (len(species_order), len(species_order)):
        raise SpeciesOrderError(f"{path.name}: shape {matrix.shape} does not match {len(species_order)} species")
    validate_distance_matrix(matrix, path.name)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp.npy")
    np.save(tmp, matrix)
    tmp.replace(path)
    write_meta(path, shape=list(matrix.shape), species_order_sha256=order_hash(species_order), **meta)
    return path


def load_matrix(path: Path, order_path: Path) -> tuple[np.ndarray, list[str]]:
    path = Path(path)
    species = read_species_order(order_path)
    matrix = np.load(path)
    if matrix.shape != (len(species), len(species)):
        raise SpeciesOrderError(f"{path.name}: shape {matrix.shape} does not match {len(species)} species in {Path(order_path).name}")
    recorded = read_meta(path).get("species_order_sha256")
    if recorded != order_hash(species):
        raise SpeciesOrderError(f"{path.name} was built for a different species order than {Path(order_path).name}; rebuild the matrices")
    return matrix, species


def text_matrix_path(directory: Path, model: str, rule: str, mask_level: str) -> Path:
    return Path(directory) / f"text__{model_slug(model)}__{rule}__{mask_level}.npy"


def phylo_matrix_path(directory: Path) -> Path:
    return Path(directory) / "phylo.npy"
