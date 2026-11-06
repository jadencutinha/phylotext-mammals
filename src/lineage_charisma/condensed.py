from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from typing import Iterator, Mapping, Sequence

import numpy as np
from scipy.stats import rankdata

from lineage_charisma.distances import validate_distance_matrix
from lineage_charisma.stats import METHODS, MantelResult, MRMResult, PartialMantelResult, _p_value

# docs/preregistration_amendment_fullscale.md, section 5: at full scale the permutation tests work on condensed
# upper-triangle vectors, rank each matrix once, and permute through index arrays in chunks. The square-matrix code in
# stats.py is the reference; tests/test_condensed.py checks that both give the same statistics and null distributions.
DEFAULT_CHUNK = 256


def default_workers() -> int:
    return max(1, min(8, (os.cpu_count() or 2) - 2))


@dataclass(frozen=True)
class Prepared:
    """One distance matrix as a centred, unit-norm upper-triangle vector (of ranks, for Spearman)."""

    vector: np.ndarray        # float64, length n(n-1)/2, row-major upper triangle
    n: int
    method: str
    sd: float                 # standard deviation of the values before standardizing (ranks for Spearman)
    mean: float

    def square32(self) -> np.ndarray:
        """Symmetric float32 square holding the vector off the diagonal, for gathering permuted values."""
        square = np.zeros((self.n, self.n), dtype=np.float32)
        rows, cols = np.triu_indices(self.n, k=1)
        square[rows, cols] = self.vector
        square[cols, rows] = self.vector
        return square


def condensed(matrix: np.ndarray) -> np.ndarray:
    return np.asarray(matrix)[np.triu_indices(matrix.shape[0], k=1)]


def prepare(matrix: np.ndarray, method: str = "spearman", name: str = "matrix", check: bool = True) -> Prepared:
    if method not in METHODS:
        raise ValueError(f"unknown method {method!r}; expected one of {METHODS}")
    matrix = np.asarray(matrix, dtype=np.float64)
    if check:
        validate_distance_matrix(matrix, name)
    return prepare_vector(condensed(matrix), matrix.shape[0], method, name)


def prepare_vector(values: np.ndarray, n: int, method: str = "spearman", name: str = "matrix") -> Prepared:
    """As `prepare`, from a condensed upper triangle already extracted (e.g. a posterior tree cached on disk)."""
    if len(values) != n * (n - 1) // 2:
        raise ValueError(f"{name}: {len(values)} values is not a condensed matrix of {n} species")
    values = np.asarray(values, dtype=np.float64)
    if method == "spearman":
        values = rankdata(values)
    mean, sd = float(values.mean()), float(values.std())
    centred = values - mean
    norm = np.linalg.norm(centred)
    if norm == 0:
        raise ValueError(f"{name}: all off-diagonal distances are equal, so the correlation is undefined")
    return Prepared(vector=centred / norm, n=n, method=method, sd=sd, mean=mean)


def iter_permutations(n: int, permutations: int, seed: int | None, chunk_size: int = DEFAULT_CHUNK) -> Iterator[np.ndarray]:
    """The permutations of stats.draw_permutations, drawn chunk by chunk from one generator (same values for any chunk size)."""
    rng = np.random.default_rng(seed)
    for start in range(0, permutations, chunk_size):
        size = min(chunk_size, permutations - start)
        yield rng.permuted(np.tile(np.arange(n), (size, 1)), axis=1).astype(np.int32)


def null_correlations(moving: Prepared, fixed: Sequence[Prepared], permutations: int, seed: int | None,
                      chunk_size: int = DEFAULT_CHUNK, workers: int | None = None) -> np.ndarray:
    """Correlation of each `fixed` vector with `moving` after reordering moving's species by each permutation: permutations x len(fixed).

    Each permutation's statistic is computed on its own (gather the permuted vector, then one float32
    matrix-vector product), so the result does not depend on the chunk size or the number of threads.
    """
    return correlations_for_permutations(moving, fixed, iter_permutations(moving.n, permutations, seed, chunk_size), permutations, workers)


def correlations_for_permutations(moving: Prepared, fixed: Sequence[Prepared], chunks, permutations: int, workers: int | None = None) -> np.ndarray:
    """As `null_correlations`, for permutations supplied by the caller as an iterable of (chunk x n) int arrays."""
    n = moving.n
    if any(f.n != n for f in fixed):
        raise ValueError("all matrices must have the same number of species")
    flat = moving.square32().ravel()
    stack = np.stack([f.vector.astype(np.float32) for f in fixed])
    rows, cols = (a.astype(np.int32) for a in np.triu_indices(n, k=1))

    def one(p: np.ndarray) -> np.ndarray:
        idx = p[rows].astype(np.int32)
        idx *= n
        idx += p[cols]
        return stack @ flat.take(idx)

    out = np.empty((permutations, len(fixed)))
    done = 0
    with ThreadPoolExecutor(workers or default_workers()) as pool:
        for chunk in chunks:
            out[done:done + len(chunk)] = np.array(list(pool.map(one, chunk)), dtype=np.float64)
            done += len(chunk)
    if done != permutations:
        raise ValueError(f"expected {permutations} permutations, got {done}")
    return out


def correlation(a: Prepared, b: Prepared) -> float:
    return float(a.vector @ b.vector)


def _inverse(corr: np.ndarray, names: Sequence[str]) -> np.ndarray:
    if np.linalg.cond(corr) > 1e10:
        raise ValueError(f"the matrices {list(names)} are collinear, so their separate contributions are undefined")
    return np.linalg.inv(corr)


def _corr_matrix(items: Sequence[Prepared]) -> np.ndarray:
    k = len(items)
    return np.array([[1.0 if i == j else correlation(items[i], items[j]) for j in range(k)] for i in range(k)])


def partial_from_correlations(r_xy: np.ndarray, r_xz: np.ndarray, r_yz: np.ndarray, inv: np.ndarray) -> np.ndarray:
    """Partial correlation of x and y given controls z, from the correlations alone (as in stats.partial_mantel)."""
    y_unexplained = 1.0 - r_yz @ inv @ r_yz
    if y_unexplained <= 1e-12:
        raise ValueError("y is fully explained by the controls, so the partial correlation is undefined")
    x_unexplained = 1.0 - np.einsum("bi,ij,bj->b", r_xz, inv, r_xz)
    return (r_xy - r_xz @ inv @ r_yz) / np.sqrt(x_unexplained * y_unexplained)


def mantel(x: Prepared, y: Prepared, permutations: int = 9999, seed: int | None = None, alternative: str = "greater",
           null: np.ndarray | None = None, **kw) -> MantelResult:
    """Mantel test reordering y, as stats.mantel. `null` may hold correlations already computed by `null_correlations` with y moving."""
    r = correlation(x, y)
    if null is None:
        null = null_correlations(y, [x], permutations, seed, **kw)[:, 0]
    return MantelResult(r=r, p=_p_value(null, r, alternative), null=null, n_pairs=len(x.vector), n=x.n, method=x.method,
                        alternative=alternative, permutations=len(null), seed=seed)


def partial_mantel(x: Prepared, y: Prepared, controls: Sequence[Prepared], permutations: int = 9999, seed: int | None = None,
                   alternative: str = "greater", null_corr: np.ndarray | None = None, **kw) -> PartialMantelResult:
    """Partial Mantel test reordering x, as stats.partial_mantel.

    `null_corr`, if given, holds the permuted correlations of [y, *controls] with x moving (columns in that order).
    """
    if not controls:
        raise ValueError("partial_mantel needs at least one control matrix; use mantel for none")
    inv = _inverse(_corr_matrix(controls), [f"control {i}" for i in range(len(controls))])
    r_yz = np.array([correlation(y, z) for z in controls])
    r_xy = correlation(x, y)
    r = float(partial_from_correlations(np.array([r_xy]), np.array([[correlation(x, z) for z in controls]]), r_yz, inv)[0])
    if null_corr is None:
        null_corr = null_correlations(x, [y, *controls], permutations, seed, **kw)
    null = partial_from_correlations(null_corr[:, 0], null_corr[:, 1:], r_yz, inv)
    return PartialMantelResult(r=r, p=_p_value(null, r, alternative), null=null, r_marginal=r_xy, n_pairs=len(x.vector), n=x.n,
                               n_controls=len(controls), method=x.method, alternative=alternative, permutations=len(null), seed=seed)


def mrm(response: Prepared, predictors: Mapping[str, Prepared], permutations: int = 9999, seed: int | None = None,
        null_corr: np.ndarray | None = None, **kw) -> MRMResult:
    """MRM reordering the response, as stats.mrm. `null_corr` may hold the permuted correlations of the predictors with the response moving."""
    names = tuple(predictors)
    xs = list(predictors.values())
    k, n_pairs = len(xs), len(response.vector)
    inv = _inverse(_corr_matrix(xs), names)
    dof = n_pairs - k - 1

    def fit(r_xy: np.ndarray):
        coef = r_xy @ inv
        r2 = np.einsum("bi,bi->b", coef, r_xy)
        se = np.sqrt(np.clip(1.0 - r2, 0.0, None)[:, None] / dof * np.diag(inv)[None, :])
        with np.errstate(divide="ignore", invalid="ignore"):
            t = coef / se
        return coef, t, r2

    coef, t, r2 = fit(np.array([[correlation(x, response) for x in xs]]))
    if null_corr is None:
        null_corr = null_correlations(response, xs, permutations, seed, **kw)
    _, null_t, null_r2 = fit(null_corr)
    p = np.array([_p_value(null_t[:, j], t[0, j], "two-sided") for j in range(k)])
    sds, means = np.array([x.sd for x in xs]), np.array([x.mean for x in xs])
    coef_raw = coef[0] * response.sd / sds
    intercept = float(response.mean - coef_raw @ means)
    return MRMResult(names=names, coef=coef[0], coef_raw=coef_raw, intercept=intercept, t=t[0], p=p, r2=float(r2[0]),
                     p_r2=_p_value(null_r2, float(r2[0]), "greater"), null_t=null_t, null_r2=null_r2, n_pairs=n_pairs, n=response.n,
                     method=response.method, permutations=len(null_corr), seed=seed)


def subset(matrix: np.ndarray, order: Sequence[str], keep: Sequence[str]) -> np.ndarray:
    """Rows and columns for the species in `keep`, in the order they appear in `order` (as stats.submatrix)."""
    wanted = set(keep)
    missing = wanted - set(order)
    if missing:
        raise ValueError(f"{len(missing)} species are not in the matrix order, e.g. {sorted(missing)[:3]}")
    idx = np.array([i for i, s in enumerate(order) if s in wanted])
    return np.asarray(matrix)[np.ix_(idx, idx)]


def mantel_on_pairs(x: np.ndarray, y: np.ndarray, keep: np.ndarray, method: str = "spearman", permutations: int = 999, seed: int | None = None,
                    alternative: str = "greater", chunk_size: int = DEFAULT_CHUNK, workers: int | None = None) -> MantelResult:
    """As stats.mantel_on_pairs, with the pair subset given as a boolean mask over the condensed upper triangle.

    y's rows and columns are reordered together while the subset of pairs stays where it is; for
    Spearman the permuted values in the subset are re-ranked every time. y is gathered in float64 so
    that ties in it are kept exactly.
    """
    x, y = np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)
    n = x.shape[0]
    rows, cols = np.triu_indices(n, k=1)
    rows, cols = rows[keep].astype(np.int32), cols[keep].astype(np.int32)
    if len(rows) < 3:
        raise ValueError(f"need at least 3 pairs, got {len(rows)}")

    def standardized(values: np.ndarray) -> np.ndarray:
        if method == "spearman":
            values = rankdata(values)
        centred = values - values.mean()
        norm = np.linalg.norm(centred)
        if norm == 0:
            raise ValueError("all distances in the subset are equal, so the correlation is undefined")
        return centred / norm

    xs = standardized(x[rows, cols])
    r = float(xs @ standardized(y[rows, cols]))
    flat = y.ravel()

    def one(p: np.ndarray) -> float:
        idx = p[rows].astype(np.int64)
        idx *= n
        idx += p[cols]
        return float(standardized(flat.take(idx)) @ xs)

    null = np.empty(permutations)
    done = 0
    with ThreadPoolExecutor(workers or default_workers()) as pool:
        for chunk in iter_permutations(n, permutations, seed, chunk_size):
            null[done:done + len(chunk)] = list(pool.map(one, chunk))
            done += len(chunk)
    return MantelResult(r=r, p=_p_value(null, r, alternative), null=null, n_pairs=len(rows), n=n, method=method,
                        alternative=alternative, permutations=permutations, seed=seed)
