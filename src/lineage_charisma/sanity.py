from __future__ import annotations

from typing import Sequence

import numpy as np
import pandas as pd


def off_diagonal(matrix: np.ndarray) -> np.ndarray:
    return matrix[np.triu_indices(matrix.shape[0], k=1)]


def nearest_neighbors(matrix: np.ndarray, order: Sequence[str], species: str, k: int = 5) -> list[tuple[str, float]]:
    i = list(order).index(species)
    row = matrix[i].copy()
    row[i] = np.inf
    # stable sort: ties (common on an ultrametric tree) keep the canonical species order
    nearest = np.argsort(row, kind="stable")[:k]
    return [(order[j], float(row[j])) for j in nearest]


def knn_indices(matrix: np.ndarray, k: int) -> np.ndarray:
    work = matrix.copy()
    np.fill_diagonal(work, np.inf)
    return np.argsort(work, axis=1, kind="stable")[:, :k]


def knn_label_agreement(matrix: np.ndarray, labels: Sequence[str], k: int) -> float:
    """Mean share of each species' k nearest neighbours that carry its own label."""
    labels = np.asarray(labels)
    return float((labels[knn_indices(matrix, k)] == labels[:, None]).mean())


def knn_overlap(a: np.ndarray, b: np.ndarray, k: int) -> float:
    """Mean number of shared species between the k nearest neighbours under two matrices."""
    na, nb = knn_indices(a, k), knn_indices(b, k)
    return float(np.mean([len(set(x) & set(y)) for x, y in zip(na, nb)]))


def family_table(matrix: np.ndarray, families: Sequence[str], min_size: int = 2) -> pd.DataFrame:
    families = np.asarray(families)
    rows = []
    for family in sorted(set(families)):
        inside = families == family
        n = int(inside.sum())
        if n < min_size or n == len(families):
            continue
        within = off_diagonal(matrix[np.ix_(inside, inside)])
        between = matrix[np.ix_(inside, ~inside)]
        rows.append({
            "family": family,
            "n": n,
            "within": float(within.mean()),
            "between": float(between.mean()),
            "ratio": float(within.mean() / between.mean()),
        })
    return pd.DataFrame(rows, columns=["family", "n", "within", "between", "ratio"]).sort_values("n", ascending=False).reset_index(drop=True)


def pooled_within_between(matrix: np.ndarray, labels: Sequence[str]) -> tuple[float, float]:
    labels = np.asarray(labels)
    same = labels[:, None] == labels[None, :]
    upper = np.triu(np.ones_like(same, dtype=bool), k=1)
    return float(matrix[same & upper].mean()), float(matrix[~same & upper].mean())


def near_neighbor_counts(matrix: np.ndarray, quantile: float) -> tuple[float, np.ndarray]:
    radius = float(np.quantile(off_diagonal(matrix), quantile))
    work = matrix.copy()
    np.fill_diagonal(work, np.inf)
    return radius, (work <= radius).sum(axis=1)


def mean_knn_distance(matrix: np.ndarray, k: int) -> np.ndarray:
    return np.take_along_axis(matrix, knn_indices(matrix, k), axis=1).mean(axis=1)


def binned_median(x: np.ndarray, y: np.ndarray, edges: np.ndarray, min_count: int = 1) -> pd.DataFrame:
    which = np.digitize(x, edges[1:-1])
    rows = []
    for b in range(len(edges) - 1):
        mask = which == b
        if mask.sum() >= min_count:
            rows.append({"lo": float(edges[b]), "hi": float(edges[b + 1]), "n": int(mask.sum()), "x": float(np.median(x[mask])), "y": float(np.median(y[mask]))})
    return pd.DataFrame(rows, columns=["lo", "hi", "n", "x", "y"])


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    rx, ry = pd.Series(x).rank().to_numpy(), pd.Series(y).rank().to_numpy()
    return float(np.corrcoef(rx, ry)[0, 1])
