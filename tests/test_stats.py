import time

import numpy as np
import pytest
from scipy.stats import kstest, pearsonr, spearmanr

from lineage_charisma.stats import draw_permutations, mantel, submatrix, upper_triangle


def random_distance_matrix(n, rng, dim=3):
    points = rng.normal(size=(n, dim))
    return np.linalg.norm(points[:, None, :] - points[None, :, :], axis=-1)


@pytest.fixture
def rng():
    return np.random.default_rng(0)


def test_identical_matrices_give_r_one_and_smallest_p(rng):
    d = random_distance_matrix(30, rng)
    for method in ("spearman", "pearson"):
        res = mantel(d, d, method=method, permutations=999, seed=1)
        assert res.r == pytest.approx(1.0)
        assert res.p == pytest.approx(1 / 1000)
        assert res.n_pairs == 30 * 29 // 2
        assert res.null.shape == (999,)


def test_monotone_transform_keeps_spearman_at_one_but_not_pearson(rng):
    d = random_distance_matrix(25, rng)
    assert mantel(d, d ** 3, method="spearman", permutations=99, seed=1).r == pytest.approx(1.0)
    assert mantel(d, d ** 3, method="pearson", permutations=99, seed=1).r < 0.99


def test_statistic_matches_scipy_on_upper_triangles(rng):
    a = random_distance_matrix(20, rng)
    b = np.round(random_distance_matrix(20, rng), 1)   # rounding makes ties
    np.fill_diagonal(b, 0.0)
    assert mantel(a, b, method="spearman", permutations=9, seed=1).r == pytest.approx(spearmanr(upper_triangle(a), upper_triangle(b))[0])
    assert mantel(a, b, method="pearson", permutations=9, seed=1).r == pytest.approx(pearsonr(upper_triangle(a), upper_triangle(b))[0])


def test_null_matches_explicit_permutation_loop(rng):
    a, b = random_distance_matrix(12, rng), random_distance_matrix(12, rng)
    perms = draw_permutations(12, 50, seed=7)
    for method, corr in (("spearman", spearmanr), ("pearson", pearsonr)):
        res = mantel(a, b, method=method, permutations=50, seed=7, batch_size=8)
        slow = [corr(upper_triangle(a), upper_triangle(b[np.ix_(p, p)]))[0] for p in perms]
        np.testing.assert_allclose(res.null, slow, atol=1e-10)


def test_permutations_are_permutations():
    perms = draw_permutations(15, 40, seed=3)
    assert perms.shape == (40, 15)
    assert all(sorted(p) == list(range(15)) for p in perms)
    assert len({tuple(p) for p in perms}) > 35


def test_seed_reproducibility(rng):
    a, b = random_distance_matrix(20, rng), random_distance_matrix(20, rng)
    first, again, other = (mantel(a, b, permutations=199, seed=s) for s in (5, 5, 6))
    np.testing.assert_array_equal(first.null, again.null)
    assert first.p == again.p
    assert not np.array_equal(first.null, other.null)
    assert first.r == other.r


def test_batch_size_does_not_change_the_result(rng):
    a, b = random_distance_matrix(20, rng), random_distance_matrix(20, rng)
    np.testing.assert_allclose(mantel(a, b, permutations=100, seed=5, batch_size=7).null,
                               mantel(a, b, permutations=100, seed=5, batch_size=100).null)


def test_independent_matrices_give_roughly_uniform_p(rng):
    ps = np.array([mantel(random_distance_matrix(15, rng), random_distance_matrix(15, rng), permutations=199, seed=i).p
                   for i in range(300)])
    assert kstest(ps, "uniform").pvalue > 0.01
    assert 0.02 < (ps <= 0.05).mean() < 0.09


def test_related_matrices_are_detected_and_alternatives_are_consistent(rng):
    points = rng.normal(size=(40, 3))
    a = np.linalg.norm(points[:, None] - points[None], axis=-1)
    noisy = points + rng.normal(scale=0.5, size=points.shape)
    b = np.linalg.norm(noisy[:, None] - noisy[None], axis=-1)
    greater = mantel(a, b, permutations=999, seed=1)
    less = mantel(a, b, permutations=999, seed=1, alternative="less")
    two = mantel(a, b, permutations=999, seed=1, alternative="two-sided")
    assert greater.r > 0.5 and greater.p == pytest.approx(0.001)
    assert less.p == pytest.approx(1.0)
    assert two.p == pytest.approx(0.001)


def test_agrees_with_scikit_bio(rng):
    skbio_mantel = pytest.importorskip("skbio.stats.distance").mantel
    points = rng.normal(size=(40, 2))
    a = np.linalg.norm(points[:, None] - points[None], axis=-1)
    for scale in (1.5, 6.0):   # a clear effect and a marginal one
        noisy = points + rng.normal(scale=scale, size=points.shape)
        b = np.linalg.norm(noisy[:, None] - noisy[None], axis=-1)
        for method in ("spearman", "pearson"):
            ours = mantel(a, b, method=method, permutations=999, seed=1)
            r, p, n = skbio_mantel(a, b, method=method, permutations=999, alternative="greater")
            assert ours.r == pytest.approx(r, abs=1e-9)
            assert n == ours.n
            # two independent 999-permutation estimates of the same p
            assert abs(ours.p - p) < 5 * np.sqrt(max(p * (1 - p), 0.001) * 2 / 999) + 2 / 1000


def test_rejects_bad_input(rng):
    d = random_distance_matrix(10, rng)
    with pytest.raises(ValueError, match="shape"):
        mantel(d, d[:9, :9])
    with pytest.raises(ValueError, match="unknown method"):
        mantel(d, d, method="kendall")
    with pytest.raises(ValueError, match="NaN"):
        bad = d.copy(); bad[0, 1] = bad[1, 0] = np.nan
        mantel(d, bad)
    with pytest.raises(ValueError, match="undefined"):
        mantel(d, np.ones((10, 10)) - np.eye(10))


def test_submatrix_keeps_matrix_order():
    order = ["a", "b", "c", "d"]
    m = np.arange(16.0).reshape(4, 4)
    np.testing.assert_array_equal(submatrix(m, order, ["d", "b"]), m[np.ix_([1, 3], [1, 3])])
    with pytest.raises(ValueError):
        submatrix(m, order, ["z"])


def test_full_size_run_is_fast(rng):
    a, b = random_distance_matrix(284, rng), random_distance_matrix(284, rng)
    start = time.perf_counter()
    mantel(a, b, permutations=9999, seed=1)
    assert time.perf_counter() - start < 60


# --- partial Mantel and MRM (Week 5 machinery; tested on synthetic matrices only) ---

from scipy.stats import rankdata  # noqa: E402

from lineage_charisma.stats import mrm, partial_mantel  # noqa: E402


def distance(points):
    return np.linalg.norm(points[:, None] - points[None], axis=-1)


def residual_correlation(x, y, controls, ranked):
    vecs = [upper_triangle(m) for m in (x, y, *controls)]
    if ranked:
        vecs = [rankdata(v) for v in vecs]
    design = np.column_stack([np.ones(len(vecs[0])), *vecs[2:]])
    resid = [v - design @ np.linalg.lstsq(design, v, rcond=None)[0] for v in vecs[:2]]
    return pearsonr(*resid)[0]


@pytest.fixture
def confounded(rng):
    """x and y each follow z and are otherwise unrelated; w is unrelated to everything."""
    z = rng.normal(size=(40, 2))
    return (distance(z + rng.normal(scale=0.7, size=z.shape)), distance(z + rng.normal(scale=0.7, size=z.shape)),
            distance(z), distance(rng.normal(size=(40, 2))))


def test_partial_mantel_matches_residual_correlation(confounded):
    x, y, z, w = confounded
    for method in ("spearman", "pearson"):
        for controls in ([z], [z, w]):
            res = partial_mantel(x, y, controls, method=method, permutations=49, seed=1)
            assert res.r == pytest.approx(residual_correlation(x, y, controls, ranked=method == "spearman"), abs=1e-10)
            assert res.r_marginal == pytest.approx(mantel(x, y, method=method, permutations=9, seed=1).r)
            assert res.n_controls == len(controls) and res.n_pairs == 40 * 39 // 2


def test_partial_mantel_one_control_matches_textbook_formula(confounded):
    x, y, z, _ = confounded
    rxy, rxz, ryz = (mantel(a, b, method="pearson", permutations=9, seed=1).r for a, b in ((x, y), (x, z), (y, z)))
    expected = (rxy - rxz * ryz) / np.sqrt((1 - rxz ** 2) * (1 - ryz ** 2))
    assert partial_mantel(x, y, [z], method="pearson", permutations=9, seed=1).r == pytest.approx(expected)


def test_partial_mantel_null_matches_explicit_loop(confounded):
    x, y, z, w = confounded
    x, y, z, w = (m[:12, :12] for m in (x, y, z, w))
    perms = draw_permutations(12, 30, seed=4)
    res = partial_mantel(x, y, [z, w], method="spearman", permutations=30, seed=4, batch_size=7)
    slow = [residual_correlation(x[np.ix_(p, p)], y, [z, w], ranked=True) for p in perms]
    np.testing.assert_allclose(res.null, slow, atol=1e-10)
    again = partial_mantel(x, y, [z, w], method="spearman", permutations=30, seed=4)
    np.testing.assert_allclose(res.null, again.null)
    assert res.p == again.p


def test_partial_mantel_removes_a_shared_driver(confounded):
    x, y, z, _ = confounded
    marginal = mantel(x, y, permutations=499, seed=1)
    partial = partial_mantel(x, y, [z], permutations=499, seed=1)
    assert marginal.r > 0.2 and marginal.p <= 0.01
    # the partial can overshoot below zero: distances are not linear in the shared driver, and the residuals are
    assert partial.r < marginal.r / 2 and partial.p > 0.05


def test_partial_mantel_keeps_a_direct_link(rng):
    z, own = rng.normal(size=(40, 2)), rng.normal(size=(40, 2))
    x = distance(np.hstack([z, own]) + rng.normal(scale=0.3, size=(40, 4)))
    y = distance(np.hstack([z, own]) + rng.normal(scale=0.3, size=(40, 4)))
    res = partial_mantel(x, y, [distance(z)], permutations=499, seed=1)
    assert res.r > 0.3 and res.p == pytest.approx(1 / 500)


def test_partial_mantel_p_is_roughly_uniform_when_only_the_control_links_them(rng):
    ps = []
    for i in range(200):
        z = rng.normal(size=(15, 2))
        x, y = distance(z + rng.normal(scale=0.7, size=z.shape)), distance(z + rng.normal(scale=0.7, size=z.shape))
        ps.append(partial_mantel(x, y, [distance(z)], permutations=199, seed=i).p)
    ps = np.array(ps)
    # partial Mantel tests are known to be somewhat liberal; this guards against gross miscalibration only
    assert (ps <= 0.05).mean() < 0.15 and 0.3 < np.median(ps) < 0.7


def test_partial_mantel_rejects_bad_input(confounded):
    x, y, z, _ = confounded
    with pytest.raises(ValueError, match="at least one control"):
        partial_mantel(x, y, [])
    with pytest.raises(ValueError, match="collinear"):
        partial_mantel(x, y, [z, z * 2], method="pearson")
    with pytest.raises(ValueError, match="fully explained"):
        partial_mantel(x, y, [y], method="pearson")
    with pytest.raises(ValueError, match="shape"):
        partial_mantel(x, y, [z[:10, :10]])


def test_mrm_matches_least_squares(confounded):
    x, y, z, w = confounded
    for method in ("spearman", "pearson"):
        res = mrm(x, {"y": y, "z": z, "w": w}, method=method, permutations=49, seed=1)
        vecs = [upper_triangle(m) for m in (x, y, z, w)]
        if method == "spearman":
            vecs = [rankdata(v) for v in vecs]
        vecs = [(v - v.mean()) / v.std() for v in vecs]
        design = np.column_stack(vecs[1:])
        beta, rss = np.linalg.lstsq(design, vecs[0], rcond=None)[:2]
        np.testing.assert_allclose(res.coef, beta, atol=1e-10)
        r2 = 1 - rss[0] / (vecs[0] ** 2).sum()
        assert res.r2 == pytest.approx(r2)
        dof = len(vecs[0]) - 3 - 1
        se = np.sqrt(rss[0] / dof * np.diag(np.linalg.inv(design.T @ design)))
        np.testing.assert_allclose(res.t, beta / se, rtol=1e-8)
        assert res.names == ("y", "z", "w") and res.null_t.shape == (49, 3) and res.null_r2.shape == (49,)


def test_mrm_with_one_predictor_is_the_mantel_test(confounded):
    x, y, _, _ = confounded
    res = mrm(x, {"y": y}, method="spearman", permutations=299, seed=3)
    plain = mantel(y, x, method="spearman", permutations=299, seed=3, alternative="two-sided")
    assert res.coef[0] == pytest.approx(plain.r) and res.r2 == pytest.approx(plain.r ** 2)
    assert res.p[0] == plain.p == res.p_r2
    skbio_mantel = pytest.importorskip("skbio.stats.distance").mantel
    assert res.coef[0] == pytest.approx(skbio_mantel(x, y, method="spearman", permutations=0)[0], abs=1e-9)


def test_mrm_null_matches_explicit_loop(confounded):
    x, y, z, _ = (m[:12, :12] for m in confounded)
    perms = draw_permutations(12, 25, seed=2)
    res = mrm(x, {"y": y, "z": z}, method="pearson", permutations=25, seed=2, batch_size=6)
    design = np.column_stack([np.ones(66), upper_triangle(y), upper_triangle(z)])
    slow = []
    for p in perms:
        v = upper_triangle(x[np.ix_(p, p)])
        fitted = design @ np.linalg.lstsq(design, v, rcond=None)[0]
        slow.append(1 - ((v - fitted) ** 2).sum() / ((v - v.mean()) ** 2).sum())
    np.testing.assert_allclose(res.null_r2, slow, atol=1e-10)


def test_mrm_separates_a_real_predictor_from_an_irrelevant_one(confounded):
    x, _, z, w = confounded
    res = mrm(x, {"z": z, "w": w}, permutations=499, seed=1)
    assert res.coef[0] > 0.3 and res.p[0] == pytest.approx(1 / 500) and res.p_r2 == pytest.approx(1 / 500)
    assert abs(res.coef[1]) < 0.15 and res.p[1] > 0.05
    with pytest.raises(ValueError, match="collinear"):
        mrm(x, {"z": z, "z2": z * 3}, method="pearson")
    with pytest.raises(ValueError, match="at least one predictor"):
        mrm(x, {})
