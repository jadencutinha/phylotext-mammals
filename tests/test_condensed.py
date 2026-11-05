import numpy as np
import pytest

from lineage_charisma import condensed as cd
from lineage_charisma.stats import draw_permutations, mantel, mrm, partial_mantel


def random_distance_matrix(n, rng, dim=3, decimals=None):
    points = rng.normal(size=(n, dim))
    d = np.linalg.norm(points[:, None, :] - points[None, :, :], axis=-1)
    if decimals is not None:   # rounding makes ties, as in patristic distances
        d = np.round(d, decimals)
        np.fill_diagonal(d, 0.0)
    return d


@pytest.fixture
def mats():
    rng = np.random.default_rng(3)
    base = random_distance_matrix(40, rng)
    noisy = base + random_distance_matrix(40, rng)
    return {"x": base, "y": np.round(noisy, 1) * (1 - np.eye(40)), "z1": random_distance_matrix(40, rng) + 0.3 * base,
            "z2": random_distance_matrix(40, rng, decimals=1)}


@pytest.mark.parametrize("method", ["spearman", "pearson"])
def test_condensed_mantel_matches_square_mantel(mats, method):
    square = mantel(mats["x"], mats["y"], method=method, permutations=999, seed=20261004)
    fast = cd.mantel(cd.prepare(mats["x"], method), cd.prepare(mats["y"], method), permutations=999, seed=20261004)
    assert fast.r == pytest.approx(square.r, abs=1e-12)
    np.testing.assert_allclose(fast.null, square.null, atol=1e-6)
    assert fast.p == square.p
    assert fast.n_pairs == square.n_pairs


@pytest.mark.parametrize("method", ["spearman", "pearson"])
def test_condensed_partial_mantel_matches_square_version(mats, method):
    controls = [mats["z1"], mats["z2"]]
    square = partial_mantel(mats["x"], mats["y"], controls, method=method, permutations=999, seed=7)
    fast = cd.partial_mantel(cd.prepare(mats["x"], method), cd.prepare(mats["y"], method), [cd.prepare(c, method) for c in controls],
                             permutations=999, seed=7)
    assert fast.r == pytest.approx(square.r, abs=1e-12)
    assert fast.r_marginal == pytest.approx(square.r_marginal, abs=1e-12)
    np.testing.assert_allclose(fast.null, square.null, atol=1e-5)
    assert fast.p == square.p


def test_condensed_mrm_matches_square_version(mats):
    preds = {"x": mats["x"], "z1": mats["z1"], "z2": mats["z2"]}
    square = mrm(mats["y"], preds, permutations=499, seed=11)
    fast = cd.mrm(cd.prepare(mats["y"]), {k: cd.prepare(v) for k, v in preds.items()}, permutations=499, seed=11)
    for field in ("coef", "coef_raw", "t"):
        np.testing.assert_allclose(getattr(fast, field), getattr(square, field), rtol=1e-10)
    assert fast.intercept == pytest.approx(square.intercept)
    assert fast.r2 == pytest.approx(square.r2)
    np.testing.assert_array_equal(fast.p, square.p)
    assert fast.p_r2 == square.p_r2


@pytest.mark.parametrize("chunk", [1, 7, 64, 1000])
def test_chunked_permutations_are_the_unchunked_draw(chunk):
    expected = draw_permutations(25, 300, seed=20261004)
    got = np.concatenate(list(cd.iter_permutations(25, 300, 20261004, chunk_size=chunk)))
    np.testing.assert_array_equal(got, expected)


def test_fewer_permutations_are_a_prefix_of_more_with_the_same_seed():
    # lets a 999-permutation robustness test reuse the first 999 permuted statistics of a 9,999-permutation pass
    np.testing.assert_array_equal(draw_permutations(30, 99, seed=5), draw_permutations(30, 999, seed=5)[:99])


@pytest.mark.parametrize("chunk,workers", [(1, 1), (13, 3), (256, 8)])
def test_null_distribution_is_identical_for_any_chunk_size_and_thread_count(mats, chunk, workers):
    x, y, z = (cd.prepare(mats[k]) for k in ("x", "y", "z1"))
    reference = cd.null_correlations(y, [x, z], 300, seed=1, chunk_size=300, workers=1)
    np.testing.assert_array_equal(cd.null_correlations(y, [x, z], 300, seed=1, chunk_size=chunk, workers=workers), reference)


def test_shared_null_pass_equals_separate_tests(mats):
    x, y, z1, z2 = (cd.prepare(mats[k]) for k in ("x", "y", "z1", "z2"))
    corr = cd.null_correlations(x, [y, z1, z2], 199, seed=4)
    separate = cd.partial_mantel(x, y, [z1, z2], permutations=199, seed=4)
    shared = cd.partial_mantel(x, y, [z1, z2], seed=4, null_corr=corr)
    np.testing.assert_array_equal(shared.null, separate.null)
    raw = cd.mantel(y, x, permutations=199, seed=4)
    assert cd.mantel(y, x, seed=4, null=corr[:, 0]).p == raw.p


def test_subset_matches_index_order():
    rng = np.random.default_rng(0)
    d = random_distance_matrix(6, rng)
    order = list("abcdef")
    np.testing.assert_array_equal(cd.subset(d, order, ["e", "b"]), d[np.ix_([1, 4], [1, 4])])


def test_condensed_mantel_on_pairs_matches_square_version(mats):
    from lineage_charisma.stats import mantel_on_pairs

    n = mats["x"].shape[0]
    rng = np.random.default_rng(2)
    keep = rng.random(n * (n - 1) // 2) < 0.3
    square = np.zeros((n, n), dtype=bool)
    square[np.triu_indices(n, k=1)] = keep
    square |= square.T
    ref = mantel_on_pairs(mats["x"], mats["y"], square, permutations=299, seed=9)
    got = cd.mantel_on_pairs(mats["x"], mats["y"], keep, permutations=299, seed=9, workers=3, chunk_size=50)
    assert got.r == pytest.approx(ref.r, abs=1e-12)
    np.testing.assert_allclose(got.null, ref.null, atol=1e-12)
    assert got.p == ref.p
