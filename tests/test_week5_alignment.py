import numpy as np
import pytest

from lineage_charisma.config import load_config
from lineage_charisma.distances import load_matrix, phylo_matrix_path, text_matrix_path, validate_distance_matrix
from lineage_charisma.embed import read_species_order
from lineage_charisma.h2 import ATTENTION_FILES, ECOLOGY_FILES, complete_trait_species, cut, decision, load_h2_matrices, load_text, share_of_raw
from lineage_charisma.io_utils import read_csv


@pytest.fixture(scope="module")
def cfg():
    config = load_config()
    needed = [config.matrices_dir() / f for f in (*ECOLOGY_FILES.values(), *ATTENTION_FILES.values())]
    if not all(p.exists() for p in needed):
        pytest.skip("ecology and attention matrices not built")
    return config


@pytest.fixture(scope="module")
def loaded(cfg):
    return load_h2_matrices(cfg)


def test_ecology_subset_is_the_full_order_less_the_species_without_ecology(cfg, loaded):
    matrices, species, order = loaded
    assert order == read_species_order(cfg.species_order_path())
    assert species == read_species_order(cfg.matrices_dir() / "ecology_species.txt")
    assert species == [s for s in order if s in set(species)]
    assert len(set(species)) == len(species)


def test_every_h2_matrix_is_a_distance_matrix_on_the_subset(loaded):
    matrices, species, _ = loaded
    assert set(matrices) == {"phylogeny", "ecology", "ecology_unweighted", "length", "pageviews"}
    for name, matrix in matrices.items():
        assert matrix.shape == (len(species), len(species)), name
        validate_distance_matrix(matrix, name)


def test_subset_entries_are_the_full_matrix_entries_for_the_same_species(cfg, loaded):
    matrices, species, order = loaded
    position = {s: i for i, s in enumerate(order)}
    full = {"phylogeny": load_matrix(phylo_matrix_path(cfg.matrices_dir()), cfg.species_order_path())[0]}
    for name, filename in ATTENTION_FILES.items():
        full[name] = load_matrix(cfg.matrices_dir() / filename, cfg.species_order_path())[0]
    rng = np.random.default_rng(0)
    for i, j in rng.integers(0, len(species), size=(200, 2)):
        a, b = position[species[i]], position[species[j]]
        for name, matrix in full.items():
            assert matrices[name][i, j] == matrix[a, b], name


@pytest.mark.parametrize("which", ["first", "all"])
def test_text_matrices_are_cut_to_the_same_species(cfg, loaded, which):
    _, species, order = loaded
    position = {s: i for i, s in enumerate(order)}
    combos = cfg.combos() if which == "all" else [cfg.primary_combo]
    for combo in combos:
        path = text_matrix_path(cfg.matrices_dir(), *combo)
        if not path.exists():
            pytest.skip(f"{path.name} not built")
        full, cut_down = load_text(cfg, combo), load_text(cfg, combo, species)
        assert cut_down.shape == (len(species), len(species))
        idx = [position[s] for s in species]
        np.testing.assert_array_equal(cut_down, full[np.ix_(idx, idx)])


def test_complete_trait_species_are_a_subset_with_no_missing_trait(cfg, loaded):
    matrices, species, _ = loaded
    complete = complete_trait_species(cfg)
    assert set(complete) <= set(species)
    assert complete == [s for s in species if s in set(complete)]
    traits = read_csv(cfg.processed / "traits.csv").set_index("species_id")
    assert not traits.loc[complete].isna().any().any()
    assert traits.drop(index=complete).isna().any(axis=1).all()
    smaller = cut(matrices, species, complete)
    assert all(m.shape == (len(complete), len(complete)) for m in smaller.values())
    i, j = species.index(complete[3]), species.index(complete[-2])
    assert smaller["ecology"][3, len(complete) - 2] == matrices["ecology"][i, j]


def test_decision_rules():
    assert decision(0.10, 0.08, 0.01, 0.05) == {"survives": True, "confounded": False}
    assert decision(0.10, 0.04, 0.01, 0.05) == {"survives": True, "confounded": True}
    assert decision(0.10, 0.06, 0.20, 0.05) == {"survives": False, "confounded": False}
    assert decision(0.10, -0.02, 0.70, 0.05) == {"survives": False, "confounded": True}
    assert decision(-0.01, 0.02, 0.30, 0.05) == {"survives": False, "confounded": None}
    assert np.isnan(share_of_raw(0.0, 0.02))
    assert share_of_raw(0.2, 0.05) == pytest.approx(0.25)
