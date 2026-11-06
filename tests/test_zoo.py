from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from lineage_charisma import zoo

FIXTURE = Path(__file__).parent / "fixtures" / "fake_signage.csv"
MATCHING = {"fuzzy_review_threshold": 85, "fuzzy_accept_threshold": 92, "fuzzy_min_margin": 3}


@pytest.fixture
def mdd():
    species = pd.DataFrame({
        "species": ["Panthera uncia", "Gorilla gorilla", "Ailurus fulgens", "Ailurus styani", "Prionailurus bengalensis"],
        "genus": ["Panthera", "Gorilla", "Ailurus", "Ailurus", "Prionailurus"],
        "family": ["Felidae", "Hominidae", "Ailuridae", "Ailuridae", "Felidae"],
        "mdd_id": [1, 2, 3, 4, 5],
        "msw3_name": ["Uncia uncia", "Gorilla gorilla", "Ailurus fulgens", None, "Prionailurus bengalensis"],
        "cmw_name": [None, None, None, None, None],
        "original_combination": [None, None, None, None, None],
    })
    synonyms = pd.DataFrame({"mdd_id": [5], "normalized_original_combination": ["Felis bengalensis"], "original_combination": ["Felis bengalensis"],
                             "root_name": ["bengalensis"]})
    return species, synonyms


def test_fake_rows_are_read_in_the_template_format():
    signs = zoo.read_signage(FIXTURE)
    assert list(signs.columns) == ["row", *zoo.TEMPLATE_COLUMNS]
    assert len(signs) == 3
    assert signs["transcription"].str.startswith("FAKE TEST ROW").all()


def test_template_columns_are_enforced(tmp_path):
    bad = tmp_path / "bad.csv"
    bad.write_text("species_common,transcription\nX,Y\n")
    with pytest.raises(ValueError, match="missing template columns"):
        zoo.read_signage(bad)


def test_template_file_has_exactly_the_template_header():
    template = Path(__file__).resolve().parents[1] / "data" / "raw" / "zoo_signage" / "template.csv"
    if not template.exists():
        pytest.skip("template not created")
    assert pd.read_csv(template).columns.tolist() == zoo.TEMPLATE_COLUMNS


def test_fake_rows_resolve_to_mdd_species(mdd):
    resolved = zoo.resolve_signage(zoo.read_signage(FIXTURE), *mdd, MATCHING)
    assert resolved["species_id"].tolist() == ["Panthera uncia", "Gorilla gorilla", "Ailurus fulgens"]
    assert resolved["resolved"].all()
    assert set(resolved["method"]) == {"exact"}           # "Gorilla gorilla gorilla" is read as its binomial


def test_legacy_synonym_fuzzy_and_failures(mdd):
    lookup, sources = zoo.name_index(*mdd)
    assert zoo.resolve_name("Uncia uncia", lookup, sources, MATCHING).method == "legacy_name"
    syn = zoo.resolve_name("Felis bengalensis", lookup, sources, MATCHING)
    assert (syn.species, syn.method) == ("Prionailurus bengalensis", "synonym")
    fuzzy = zoo.resolve_name("Panthera unica", lookup, sources, MATCHING)
    assert (fuzzy.species, fuzzy.method) == ("Panthera uncia", "fuzzy")
    assert zoo.resolve_name("Leo", lookup, sources, MATCHING).method == "unparsed"
    assert zoo.resolve_name("Canis lupus", lookup, sources, MATCHING).species is None


def test_similarity_table_on_fake_vectors():
    rng = np.random.default_rng(0)
    species = ["A a", "B b", "C c"]
    wiki = {s: rng.normal(size=8) for s in species}
    sign = {s: wiki[s] + 0.1 * rng.normal(size=8) for s in species}     # each sign close to its own article
    table = zoo.similarity_table(species, sign, wiki)
    assert table["sign_wiki_above_baseline"].all()
    assert table["sign_adw"].isna().all()
    summary = zoo.similarity_summary(table)
    assert summary["n_signed_species"] == 3 and summary["share_sign_wiki_above_baseline"] == 1.0
