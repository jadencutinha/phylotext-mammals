from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from scipy.stats import rankdata

from lineage_charisma.distances import validate_distance_matrix

METHODS = ("spearman", "pearson")
ALTERNATIVES = ("greater", "less", "two-sided")


@dataclass(frozen=True)
class MantelResult:
    r: float
    p: float
    null: np.ndarray          # the statistic under each permutation, in the order drawn
    n_pairs: int
    n: int
    method: str
    alternative: str
    permutations: int
    seed: int | None


def upper_triangle(matrix: np.ndarray) -> np.ndarray:
    return matrix[np.triu_indices(matrix.shape[0], k=1)]


def _standardized_square(values: np.ndarray, n: int, name: str) -> np.ndarray:
    """Symmetric matrix whose upper triangle is `values` centred and scaled to unit norm; zero diagonal."""
    centred = values - values.mean()
    norm = np.linalg.norm(centred)
    if norm == 0:
        raise ValueError(f"{name}: all off-diagonal distances are equal, so the correlation is undefined")
    square = np.zeros((n, n))
    square[np.triu_indices(n, k=1)] = centred / norm
    return square + square.T


def draw_permutations(n: int, permutations: int, seed: int | None) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.permuted(np.tile(np.arange(n), (permutations, 1)), axis=1)


def permuted_statistics(x_square: np.ndarray, y_square: np.ndarray, perms: np.ndarray, batch_size: int = 128) -> np.ndarray:
    """Correlation of x with y after reordering y's rows and columns by each permutation.

    Both inputs come from `_standardized_square`, so the correlation over the upper triangle
    is half the sum of the elementwise product over the whole matrix. Permuting does not
    change the mean or norm of y's upper triangle, so nothing is re-standardized.
    """
    out = np.empty(len(perms))
    for start in range(0, len(perms), batch_size):
        p = perms[start:start + batch_size]
        permuted = y_square[p[:, :, None], p[:, None, :]]
        out[start:start + len(p)] = np.einsum("ij,bij->b", x_square, permuted) / 2.0
    return out


def mantel(
    x: np.ndarray,
    y: np.ndarray,
    method: str = "spearman",
    permutations: int = 9999,
    seed: int | None = None,
    alternative: str = "greater",
    batch_size: int = 128,
) -> MantelResult:
    """Mantel test between two distance matrices in the same species order.

    The statistic is the correlation between the two upper triangles. The null distribution
    reorders the rows and columns of `y` together. p = (1 + permuted statistics at least as
    extreme as the observed one) / (1 + permutations). For Spearman the distances are ranked
    once (average ranks for ties); permuting species only moves ranks around.
    """
    if method not in METHODS:
        raise ValueError(f"unknown method {method!r}; expected one of {METHODS}")
    if alternative not in ALTERNATIVES:
        raise ValueError(f"unknown alternative {alternative!r}; expected one of {ALTERNATIVES}")
    x, y = np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)
    validate_distance_matrix(x, "x")
    validate_distance_matrix(y, "y")
    if x.shape != y.shape:
        raise ValueError(f"matrices differ in shape: {x.shape} and {y.shape}")
    n = x.shape[0]
    if n < 3:
        raise ValueError(f"a Mantel test needs at least 3 species, got {n}")

    xv, yv = upper_triangle(x), upper_triangle(y)
    if method == "spearman":
        xv, yv = rankdata(xv), rankdata(yv)
    x_square, y_square = _standardized_square(xv, n, "x"), _standardized_square(yv, n, "y")

    r = float((x_square * y_square).sum() / 2.0)
    null = permuted_statistics(x_square, y_square, draw_permutations(n, permutations, seed), batch_size)

    return MantelResult(r=r, p=_p_value(null, r, alternative), null=null, n_pairs=len(xv), n=n, method=method,
                        alternative=alternative, permutations=permutations, seed=seed)


def submatrix(matrix: np.ndarray, order: Sequence[str], keep: Sequence[str]) -> np.ndarray:
    """Rows and columns of `matrix` for the species in `keep`, in the order they appear in `order`."""
    wanted = set(keep)
    missing = wanted - set(order)
    if missing:
        raise ValueError(f"{len(missing)} species are not in the matrix order, e.g. {sorted(missing)[:3]}")
    idx = [i for i, s in enumerate(order) if s in wanted]
    return matrix[np.ix_(idx, idx)]


@dataclass(frozen=True)
class PartialMantelResult:
    r: float                  # partial correlation of x and y given the controls
    p: float
    null: np.ndarray
    r_marginal: float         # plain Mantel correlation of x and y on the same species
    n_pairs: int
    n: int
    n_controls: int
    method: str
    alternative: str
    permutations: int
    seed: int | None


@dataclass(frozen=True)
class MRMResult:
    names: tuple[str, ...]
    coef: np.ndarray          # standardized regression coefficients, one per predictor
    t: np.ndarray             # pseudo-t of each coefficient
    p: np.ndarray             # two-sided permutation p-value of each pseudo-t
    r2: float
    p_r2: float
    null_t: np.ndarray        # permutations x predictors
    null_r2: np.ndarray
    n_pairs: int
    n: int
    method: str
    permutations: int
    seed: int | None


def _prepared_squares(matrices: Sequence[np.ndarray], names: Sequence[str], method: str) -> tuple[list[np.ndarray], int, int]:
    if method not in METHODS:
        raise ValueError(f"unknown method {method!r}; expected one of {METHODS}")
    arrays = [np.asarray(m, dtype=np.float64) for m in matrices]
    for a, name in zip(arrays, names):
        validate_distance_matrix(a, name)
        if a.shape != arrays[0].shape:
            raise ValueError(f"matrices differ in shape: {name} is {a.shape}, {names[0]} is {arrays[0].shape}")
    n = arrays[0].shape[0]
    if n < 3:
        raise ValueError(f"a Mantel test needs at least 3 species, got {n}")
    vectors = [upper_triangle(a) for a in arrays]
    if method == "spearman":
        vectors = [rankdata(v) for v in vectors]
    return [_standardized_square(v, n, name) for v, name in zip(vectors, names)], n, len(vectors[0])


def permuted_correlations(fixed: Sequence[np.ndarray], moving: np.ndarray, perms: np.ndarray, batch_size: int = 128) -> np.ndarray:
    """Correlation of each matrix in `fixed` with `moving` after reordering `moving` by each permutation: permutations x len(fixed)."""
    stack = np.stack(fixed)
    out = np.empty((len(perms), len(stack)))
    for start in range(0, len(perms), batch_size):
        p = perms[start:start + batch_size]
        permuted = moving[p[:, :, None], p[:, None, :]]
        out[start:start + len(p)] = np.einsum("kij,bij->bk", stack, permuted) / 2.0
    return out


def _correlation(a: np.ndarray, b: np.ndarray) -> float:
    return float((a * b).sum() / 2.0)


def _control_inverse(controls: Sequence[np.ndarray], names: Sequence[str]) -> np.ndarray:
    k = len(controls)
    corr = np.array([[1.0 if i == j else _correlation(controls[i], controls[j]) for j in range(k)] for i in range(k)])
    if np.linalg.cond(corr) > 1e10:
        raise ValueError(f"the matrices {list(names)} are collinear, so their separate contributions are undefined")
    return np.linalg.inv(corr)


def _p_value(null: np.ndarray, observed: float, alternative: str) -> float:
    # the tolerance keeps a permutation that reproduces the observed value from being lost to rounding
    tol = 1e-12
    if alternative == "greater":
        extreme = null >= observed - tol
    elif alternative == "less":
        extreme = null <= observed + tol
    else:
        extreme = np.abs(null) >= abs(observed) - tol
    return (1 + int(extreme.sum())) / (1 + len(null))


def partial_mantel(
    x: np.ndarray,
    y: np.ndarray,
    controls: Sequence[np.ndarray],
    method: str = "spearman",
    permutations: int = 9999,
    seed: int | None = None,
    alternative: str = "greater",
    batch_size: int = 128,
) -> PartialMantelResult:
    """Partial Mantel test: the correlation of x and y after the linear effect of the control matrices is removed from both.

    The statistic is the partial correlation of the upper triangles (of their ranks, for Spearman).
    The null distribution reorders the rows and columns of `x` together and recomputes the partial
    correlation each time, leaving `y` and the controls as they are (the scheme of vegan's
    mantel.partial, extended to several controls).
    """
    if alternative not in ALTERNATIVES:
        raise ValueError(f"unknown alternative {alternative!r}; expected one of {ALTERNATIVES}")
    controls = list(controls)
    if not controls:
        raise ValueError("partial_mantel needs at least one control matrix; use mantel for none")
    names = ["x", "y", *[f"control {i}" for i in range(len(controls))]]
    squares, n, n_pairs = _prepared_squares([x, y, *controls], names, method)
    xs, ys, zs = squares[0], squares[1], squares[2:]
    inv = _control_inverse(zs, names[2:])
    r_yz = np.array([_correlation(ys, z) for z in zs])
    y_unexplained = 1.0 - r_yz @ inv @ r_yz
    if y_unexplained <= 1e-12:
        raise ValueError("y is fully explained by the controls, so the partial correlation is undefined")

    def partial(r_xy: np.ndarray, r_xz: np.ndarray) -> np.ndarray:
        x_unexplained = 1.0 - np.einsum("bi,ij,bj->b", r_xz, inv, r_xz)
        return (r_xy - r_xz @ inv @ r_yz) / np.sqrt(x_unexplained * y_unexplained)

    r_xy = _correlation(xs, ys)
    r = float(partial(np.array([r_xy]), np.array([[_correlation(xs, z) for z in zs]]))[0])
    corr = permuted_correlations([ys, *zs], xs, draw_permutations(n, permutations, seed), batch_size)
    null = partial(corr[:, 0], corr[:, 1:])
    return PartialMantelResult(r=r, p=_p_value(null, r, alternative), null=null, r_marginal=r_xy, n_pairs=n_pairs, n=n, n_controls=len(zs),
                               method=method, alternative=alternative, permutations=permutations, seed=seed)


def mrm(
    response: np.ndarray,
    predictors: dict[str, np.ndarray],
    method: str = "spearman",
    permutations: int = 9999,
    seed: int | None = None,
    batch_size: int = 128,
) -> MRMResult:
    """Multiple regression on distance matrices (Lichstein 2007).

    The response's upper triangle is regressed on the predictors' upper triangles, all standardized
    (and ranked first, for Spearman), so the coefficients are standardized. Significance comes from
    reordering the rows and columns of the response together: each coefficient is tested two-sided
    on its pseudo-t, and the whole model on R squared.
    """
    if not predictors:
        raise ValueError("mrm needs at least one predictor matrix")
    names = tuple(predictors)
    squares, n, n_pairs = _prepared_squares([response, *predictors.values()], ["response", *names], method)
    ys, xs = squares[0], squares[1:]
    k = len(xs)
    inv = _control_inverse(xs, names)
    dof = n_pairs - k - 1

    def fit(r_xy: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        coef = r_xy @ inv
        r2 = np.einsum("bi,bi->b", coef, r_xy)
        se = np.sqrt(np.clip(1.0 - r2, 0.0, None)[:, None] / dof * np.diag(inv)[None, :])
        with np.errstate(divide="ignore", invalid="ignore"):
            t = coef / se
        return coef, t, r2

    coef, t, r2 = fit(np.array([[_correlation(x, ys) for x in xs]]))
    _, null_t, null_r2 = fit(permuted_correlations(xs, ys, draw_permutations(n, permutations, seed), batch_size))
    p = np.array([_p_value(null_t[:, j], t[0, j], "two-sided") for j in range(k)])
    return MRMResult(names=names, coef=coef[0], t=t[0], p=p, r2=float(r2[0]), p_r2=_p_value(null_r2, float(r2[0]), "greater"),
                     null_t=null_t, null_r2=null_r2, n_pairs=n_pairs, n=n, method=method, permutations=permutations, seed=seed)
