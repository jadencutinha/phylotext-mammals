import numpy as np
import pytest

from lineage_charisma.config import load_config
from lineage_charisma.distances import load_matrix, order_hash, phylo_matrix_path, text_matrix_path
from lineage_charisma.embed import read_species_order
from lineage_charisma.io_utils import read_json, read_meta

# docs/preregistration_amendment_fullscale.md, section 3: every full-scale matrix uses the frozen species_order.txt
FROZEN_SHA256 = "4730f40785d21f82ed7a93d06e6d8f71a3e75e290fddcbac63aa0d9872763f0e"


@pytest.fixture(scope="module")
def cfg():
    config = load_config(clade="mammalia")
    if not phylo_matrix_path(config.matrices_dir()).exists():
        pytest.skip("full-scale matrices not built")
    return config


def test_frozen_order_is_the_committed_list(cfg):
    import hashlib

    assert hashlib.sha256(cfg.species_order_path().read_bytes()).hexdigest() == FROZEN_SHA256


def test_every_frozen_set_matrix_is_aligned_to_species_order(cfg):
    order = read_species_order(cfg.species_order_path())
    paths = [phylo_matrix_path(cfg.matrices_dir()), *(text_matrix_path(cfg.matrices_dir(), *c) for c in cfg.combos()),
             *(cfg.matrices_dir() / f for f in ("attention_length.npy", "attention_pageviews.npy", "attention_mean_log_pageviews.npy"))]
    for path in paths:
        assert read_meta(path)["species_order_sha256"] == order_hash(order), path.name
        matrix, loaded = load_matrix(path, cfg.species_order_path())
        assert loaded == order and matrix.shape == (len(order), len(order)), path.name


def test_ecology_matrices_are_aligned_to_an_ordered_subset(cfg):
    order = read_species_order(cfg.species_order_path())
    eco_path = cfg.matrices_dir() / "ecology_species.txt"
    species = read_species_order(eco_path)
    assert species == [s for s in order if s in set(species)]
    for name in ("ecology.npy", "ecology_unweighted.npy", "ecology_no_volant.npy"):
        matrix, loaded = load_matrix(cfg.matrices_dir() / name, eco_path)
        assert loaded == species and matrix.shape == (len(species), len(species)), name


def test_text_matrix_rows_follow_the_embedding_species(cfg):
    from lineage_charisma.distances import cosine_distance_matrix, reindex_rows
    from lineage_charisma.embed import embedding_paths, load_embeddings

    order = read_species_order(cfg.species_order_path())
    combo = cfg.primary_combo
    emb, emb_species = load_embeddings(embedding_paths(cfg.embeddings_dir(), *combo))
    rng = np.random.default_rng(0)
    pick = sorted(rng.choice(len(order), size=40, replace=False))
    sub = reindex_rows(emb, emb_species, [order[i] for i in pick], "embeddings")
    text, _ = load_matrix(text_matrix_path(cfg.matrices_dir(), *combo), cfg.species_order_path())
    np.testing.assert_allclose(text[np.ix_(pick, pick)], cosine_distance_matrix(sub), atol=1e-6)


def test_posterior_tree_caches_match_the_frozen_order(cfg):
    caches = sorted((cfg.interim / "posterior").glob("*.json"))
    if not caches:
        pytest.skip("posterior trees not cached")
    order = read_species_order(cfg.species_order_path())
    n = len(order)
    for key in caches:
        assert read_json(key)["species_order_sha256"] == order_hash(order), key.name
        assert len(np.load(key.with_suffix(".npy"), mmap_mode="r")) == n * (n - 1) // 2
