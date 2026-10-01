import dendropy
import numpy as np
import pandas as pd
import pytest

from lineage_charisma.corpus import build_corpus, read_corpus, stub_table, write_corpus
from lineage_charisma.distances import (
    SpeciesOrderError,
    assert_same_order,
    canonical_order,
    cosine_distance_matrix,
    load_matrix,
    phylo_distance_matrix,
    reindex_matrix,
    reindex_rows,
    save_matrix,
    tree_tip_order,
    validate_distance_matrix,
)
from lineage_charisma.embed import embed_texts, embedding_paths, load_embeddings, save_embeddings, write_species_order
from test_embed import FakeEncoder

NEWICK = "(((Panthera_leo:2,Panthera_tigris:2):3,Felis_catus:5):5,(Canis_lupus:4,Vulpes_vulpes:4):6);"

DESCRIPTIONS = {
    "Vulpes vulpes": "the red fox is a small canid with a bushy tail and a red coat",
    "Panthera leo": "the lion is a large cat with a tawny coat and the male has a mane",
    "Canis lupus": "the wolf is a large canid with a long muzzle and a grey coat " * 3,
    "Felis catus": "the cat is a small cat",
    "Panthera tigris": "the tiger is a large cat with an orange coat and black stripes " * 2,
}


@pytest.fixture
def tree():
    return dendropy.Tree.get(data=NEWICK, schema="newick", preserve_underscores=True, rooting="default-rooted")


@pytest.fixture
def corpus():
    matched = pd.DataFrame({
        "species": list(DESCRIPTIONS),
        "family": ["Canidae", "Felidae", "Canidae", "Felidae", "Felidae"],
    })
    return build_corpus(matched, DESCRIPTIONS, lambda texts: [len(t.split()) for t in texts], min_tokens=7)


def test_corpus_is_keyed_by_species_and_flags_stubs(corpus, tmp_path):
    assert list(corpus.columns) == ["species_id", "family", "text", "n_tokens", "is_stub"]
    assert corpus["species_id"].tolist() == sorted(DESCRIPTIONS)
    assert len(corpus) == 5
    assert corpus.set_index("species_id")["is_stub"].to_dict() == {
        "Canis lupus": False, "Felis catus": True, "Panthera leo": False, "Panthera tigris": False, "Vulpes vulpes": False,
    }
    assert stub_table(corpus)["species_id"].tolist() == ["Felis catus"]
    path = write_corpus(corpus, tmp_path / "corpus.parquet")
    assert read_corpus(path).equals(corpus)


def test_build_corpus_rejects_missing_text():
    matched = pd.DataFrame({"species": ["Panthera leo", "Lynx lynx"], "family": ["Felidae", "Felidae"]})
    with pytest.raises(KeyError):
        build_corpus(matched, DESCRIPTIONS, lambda texts: [1] * len(texts), min_tokens=1)


def test_validate_accepts_a_proper_distance_matrix():
    validate_distance_matrix(np.array([[0.0, 1.0, 2.0], [1.0, 0.0, 3.0], [2.0, 3.0, 0.0]]))


@pytest.mark.parametrize("bad,message", [
    (np.array([[0.0, 1.0], [2.0, 0.0]]), "symmetric"),
    (np.array([[0.5, 1.0], [1.0, 0.0]]), "diagonal"),
    (np.array([[0.0, np.nan], [np.nan, 0.0]]), "NaN"),
    (np.array([[0.0, -1.0], [-1.0, 0.0]]), "negative"),
    (np.zeros((2, 3)), "square"),
])
def test_validate_rejects_bad_matrices(bad, message):
    with pytest.raises(ValueError, match=message):
        validate_distance_matrix(bad)


def test_cosine_distance_matrix_is_valid_and_correct():
    emb = np.array([[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [1.0, 0.0]])
    dist = cosine_distance_matrix(emb)
    validate_distance_matrix(dist)
    assert np.array_equal(dist, dist.T) and np.all(np.diag(dist) == 0)
    assert np.allclose(dist[0], [0.0, 1.0, 2.0, 0.0])


def test_cosine_distance_matrix_requires_normalized_embeddings():
    with pytest.raises(ValueError, match="normalized"):
        cosine_distance_matrix(np.array([[2.0, 0.0], [0.0, 1.0]]))


def test_phylo_distance_matrix_follows_the_requested_order(tree):
    order = ["Vulpes vulpes", "Panthera leo", "Felis catus"]
    dist = phylo_distance_matrix(tree, order)
    validate_distance_matrix(dist)
    assert dist.tolist() == [[0.0, 20.0, 20.0], [20.0, 0.0, 10.0], [20.0, 10.0, 0.0]]
    with pytest.raises(SpeciesOrderError):
        phylo_distance_matrix(tree, ["Panthera leo", "Lynx lynx"])


def test_tree_tip_order_uses_species_ids(tree):
    assert tree_tip_order(tree, sorted(DESCRIPTIONS)) == ["Panthera leo", "Panthera tigris", "Felis catus", "Canis lupus", "Vulpes vulpes"]


def test_assert_same_order_reports_set_and_order_mismatches():
    assert_same_order(["a", "b"], ["a", "b"])
    with pytest.raises(SpeciesOrderError, match="different order"):
        assert_same_order(["b", "a"], ["a", "b"])
    with pytest.raises(SpeciesOrderError, match="species set differs"):
        assert_same_order(["a", "c"], ["a", "b"])
    with pytest.raises(SpeciesOrderError):
        assert_same_order(["a"], ["a", "b"])


def test_reindex_matrix_moves_rows_and_columns_together():
    matrix = np.array([[0.0, 1.0, 2.0], [1.0, 0.0, 3.0], [2.0, 3.0, 0.0]])
    out = reindex_matrix(matrix, ["a", "b", "c"], ["c", "a", "b"])
    assert out.tolist() == [[0.0, 2.0, 3.0], [2.0, 0.0, 1.0], [3.0, 1.0, 0.0]]
    with pytest.raises(SpeciesOrderError):
        reindex_matrix(matrix, ["a", "b", "c"], ["a", "b", "d"])
    with pytest.raises(SpeciesOrderError):
        reindex_matrix(matrix, ["a", "b"], ["a", "b"])


def test_text_and_phylo_matrices_align_on_one_species_order(tree, corpus, tmp_path):
    # embeddings are cached in a deliberately different order from the canonical one
    shuffled = corpus.sample(frac=1.0, random_state=3)
    assert shuffled["species_id"].tolist() != sorted(DESCRIPTIONS)
    result = embed_texts(shuffled["text"].tolist(), FakeEncoder(max_content_tokens=8), "chunk", overlap=2)
    paths = embedding_paths(tmp_path / "embeddings", "fake/model", "chunk", "none")
    save_embeddings(paths, result.embeddings, shuffled["species_id"].tolist())

    order = canonical_order(corpus["species_id"])
    order_path = write_species_order(tmp_path / "matrices" / "species_order.txt", order)
    embeddings, emb_species = load_embeddings(paths)
    text = cosine_distance_matrix(reindex_rows(embeddings, emb_species, order))
    phylo = phylo_distance_matrix(tree, order)
    text_path = save_matrix(tmp_path / "matrices" / "text.npy", text, order)
    phylo_path = save_matrix(tmp_path / "matrices" / "phylo.npy", phylo, order)

    text_loaded, text_order = load_matrix(text_path, order_path)
    phylo_loaded, phylo_order = load_matrix(phylo_path, order_path)
    assert text_order == phylo_order == order == sorted(DESCRIPTIONS)
    assert text_loaded.shape == phylo_loaded.shape == (5, 5)

    i, j, k = order.index("Panthera leo"), order.index("Panthera tigris"), order.index("Canis lupus")
    assert phylo_loaded[i, j] == 4.0 and phylo_loaded[i, k] == 20.0
    lion = embed_texts([DESCRIPTIONS["Panthera leo"], DESCRIPTIONS["Panthera tigris"]], FakeEncoder(max_content_tokens=8), "chunk", overlap=2).embeddings
    assert text_loaded[i, j] == pytest.approx(1.0 - float(lion[0].astype(float) @ lion[1].astype(float)), abs=1e-5)


def test_load_matrix_fails_loudly_when_the_species_order_changes(tree, tmp_path):
    order = canonical_order(DESCRIPTIONS)
    order_path = write_species_order(tmp_path / "species_order.txt", order)
    path = save_matrix(tmp_path / "phylo.npy", phylo_distance_matrix(tree, order), order)
    load_matrix(path, order_path)
    write_species_order(order_path, list(reversed(order)))
    with pytest.raises(SpeciesOrderError, match="different species order"):
        load_matrix(path, order_path)
    write_species_order(order_path, order[:-1])
    with pytest.raises(SpeciesOrderError, match="does not match"):
        load_matrix(path, order_path)


def test_save_matrix_rejects_invalid_or_misaligned_input(tmp_path):
    with pytest.raises(SpeciesOrderError):
        save_matrix(tmp_path / "m.npy", np.zeros((3, 3)), ["a", "b"])
    with pytest.raises(ValueError):
        save_matrix(tmp_path / "m.npy", np.array([[0.0, 1.0], [2.0, 0.0]]), ["a", "b"])
