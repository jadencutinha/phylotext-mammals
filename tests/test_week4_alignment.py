import numpy as np
import pytest

from lineage_charisma.config import load_config
from lineage_charisma.distances import load_matrix, validate_distance_matrix
from lineage_charisma.embed import read_species_order
from lineage_charisma.io_utils import read_csv

FULL_ORDER = ["attention_pageviews.npy", "attention_length.npy", "attention_combined.npy", "attention_mean_log_pageviews.npy"]
ECOLOGY = ["ecology.npy", "ecology_unweighted.npy"]


@pytest.fixture(scope="module")
def cfg():
    return load_config()


def built(cfg, name):
    path = cfg.matrices_dir() / name
    if not path.exists():
        pytest.skip(f"{name} not built")
    return path


@pytest.mark.parametrize("name", FULL_ORDER)
def test_attention_matrices_follow_species_order(cfg, name):
    matrix, order = load_matrix(built(cfg, name), cfg.species_order_path())   # checks shape and the recorded order hash
    assert order == read_species_order(cfg.species_order_path())
    table = read_csv(cfg.processed / "attention.csv")
    assert table["species_id"].tolist() == order
    i, j = 0, len(order) - 1
    if name == "attention_mean_log_pageviews.npy":
        assert matrix[i, j] == pytest.approx((table["log_pageviews"][i] + table["log_pageviews"][j]) / 2)
        np.testing.assert_allclose(np.diag(matrix), table["log_pageviews"])
    else:
        validate_distance_matrix(matrix, name)
    if name == "attention_pageviews.npy":
        assert matrix[i, j] == pytest.approx(abs(table["log_pageviews"][i] - table["log_pageviews"][j]))
    if name == "attention_length.npy":
        assert matrix[i, j] == pytest.approx(abs(table["log_length"][i] - table["log_length"][j]))


@pytest.mark.parametrize("name", ECOLOGY)
def test_ecology_matrices_follow_their_species_list(cfg, name):
    path = built(cfg, name)
    eco_order_path = cfg.matrices_dir() / "ecology_species.txt"
    matrix, eco_order = load_matrix(path, eco_order_path)
    validate_distance_matrix(matrix, name)
    assert matrix.max() <= 1.0
    full = read_species_order(cfg.species_order_path())
    # the ecology species are the full order with any species lacking an ecological distance removed, in the same sequence
    assert eco_order == [s for s in full if s in set(eco_order)]
    traits = read_csv(cfg.processed / "traits.csv")
    assert traits["species_id"].tolist() == full
