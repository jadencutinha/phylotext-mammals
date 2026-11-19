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


# ---------------------------------------------------------------- signs with no printed scientific name (amendment 5, section D)

COMMON = pd.DataFrame({
    "species": ["Panthera uncia", "Ailurus fulgens", "Ailurus styani", "Elaphurus davidianus", "Gorilla gorilla"],
    "common_name": ["Snow Leopard", "Himalayan Red Panda", "Chinese Red Panda", "Père David's Deer", "Western Gorilla"],
    "other_common_names": ["Ounce", "Red Panda|Lesser Panda", "Red Panda", None, "Western Lowland Gorilla|Lowland Gorilla"],
})


@pytest.mark.parametrize("printed, expected", [
    ("Snow Leopard", "snow leopard"),
    ("  SNOW-LEOPARD ", "snow leopard"),
    ("Snow—Leopard", "snow leopard"),
    ("Père David’s Deer", "père davids deer"),          # accents kept, curly apostrophe deleted
    ("Pere David's Deer (Milu)", "pere davids deer milu"),
])
def test_common_names_are_normalized(printed, expected):
    assert zoo.normalize_common_name(printed) == expected


def test_common_name_match_is_exact_and_ambiguity_is_excluded():
    index = zoo.common_name_index(COMMON)
    hit = zoo.resolve_common_name("western lowland-gorilla", index)
    assert (hit.species, hit.method) == ("Gorilla gorilla", "common_name")
    assert zoo.resolve_common_name("Ounce", index).species == "Panthera uncia"              # other common names count
    amb = zoo.resolve_common_name("Red Panda", index)
    assert (amb.species, amb.method) == (None, "ambiguous_common_name") and "Ailurus fulgens" in amb.note
    assert zoo.resolve_common_name("Snow Leopards", index).method == "no_common_name_match"   # no fuzzy or plural matching
    assert zoo.resolve_common_name("Pere David's Deer", index).method == "no_common_name_match"  # accents are not folded


def test_rows_route_by_scientific_name_when_printed_else_common_name(mdd):
    species = mdd[0].merge(COMMON, on="species", how="left")
    signs = pd.DataFrame({"row": [1, 2, 3, 4, 5],
                          "species_common": ["Snow Leopard", "Western Lowland Gorilla", "Red Panda", "Snow Leopard", "Tiger"],
                          "species_scientific": ["Panthera uncia", "", "", "", "Panthera nonexistens"],
                          "exhibit": [""] * 5, "photo_filename": [""] * 5, "transcription": ["t"] * 5, "notes": [""] * 5})
    out = zoo.resolve_signage(signs, species, mdd[1], MATCHING)
    assert out["route"].tolist() == [zoo.SCIENTIFIC, zoo.COMMON, zoo.COMMON, zoo.COMMON, zoo.SCIENTIFIC]
    assert out["resolved"].tolist() == [True, True, False, True, False]   # row 5: an unresolved scientific name gets no common-name fallback
    assert out.loc[out["resolved"], "species_id"].tolist() == ["Panthera uncia", "Gorilla gorilla", "Panthera uncia"]
    assert zoo.route_counts(out) == {"n_species_by_scientific_name": 1, "n_species_by_common_name": 1,
                                     "n_common_name_rows_ambiguous": 1, "n_common_name_rows_no_match": 0}
