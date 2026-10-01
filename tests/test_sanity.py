import numpy as np
import pytest

from lineage_charisma.sanity import (
    binned_median,
    family_table,
    knn_label_agreement,
    knn_overlap,
    mean_knn_distance,
    near_neighbor_counts,
    nearest_neighbors,
    off_diagonal,
    pooled_within_between,
    spearman,
)

ORDER = ["a1", "a2", "a3", "b1", "b2"]
FAMILIES = ["A", "A", "A", "B", "B"]


@pytest.fixture
def matrix():
    m = np.full((5, 5), 10.0)
    m[:3, :3] = [[0, 1, 2], [1, 0, 3], [2, 3, 0]]
    m[3:, 3:] = [[0, 4], [4, 0]]
    return m


def test_nearest_neighbors_excludes_self_and_breaks_ties_by_order(matrix):
    assert nearest_neighbors(matrix, ORDER, "a1", 2) == [("a2", 1.0), ("a3", 2.0)]
    assert nearest_neighbors(matrix, ORDER, "b1", 3) == [("b2", 4.0), ("a1", 10.0), ("a2", 10.0)]


def test_family_table_within_and_between(matrix):
    table = family_table(matrix, FAMILIES).set_index("family")
    assert table.loc["A", "n"] == 3 and table.loc["A", "within"] == pytest.approx(2.0)
    assert table.loc["B", "within"] == pytest.approx(4.0)
    assert table.loc["A", "between"] == table.loc["B", "between"] == pytest.approx(10.0)
    assert table.loc["A", "ratio"] == pytest.approx(0.2)
    assert family_table(matrix, FAMILIES, min_size=3)["family"].tolist() == ["A"]


def test_pooled_within_between(matrix):
    within, between = pooled_within_between(matrix, FAMILIES)
    assert within == pytest.approx((1 + 2 + 3 + 4) / 4)
    assert between == pytest.approx(10.0)


def test_knn_label_agreement_and_overlap(matrix):
    assert knn_label_agreement(matrix, FAMILIES, 1) == 1.0
    assert knn_label_agreement(matrix, FAMILIES, 2) == pytest.approx(0.8)
    assert knn_overlap(matrix, matrix, 2) == 2.0
    other = matrix[::-1, ::-1].copy()
    assert 0.0 <= knn_overlap(matrix, other, 2) <= 2.0


def test_near_neighbor_counts_and_mean_knn_distance(matrix):
    radius, counts = near_neighbor_counts(matrix, 0.4)
    assert radius == pytest.approx(np.quantile(off_diagonal(matrix), 0.4))
    assert counts.tolist() == [2, 2, 2, 1, 1]
    assert mean_knn_distance(matrix, 2).tolist() == [1.5, 2.0, 2.5, 7.0, 7.0]


def test_binned_median_and_spearman():
    x, y = np.array([1.0, 2.0, 6.0, 7.0, 8.0]), np.array([10.0, 20.0, 1.0, 2.0, 3.0])
    out = binned_median(x, y, np.array([0.0, 5.0, 10.0]))
    assert out["n"].tolist() == [2, 3] and out["y"].tolist() == [15.0, 2.0]
    assert binned_median(x, y, np.array([0.0, 5.0, 10.0]), min_count=3)["n"].tolist() == [3]
    assert spearman(np.array([1.0, 2.0, 3.0]), np.array([10.0, 100.0, 1000.0])) == pytest.approx(1.0)
