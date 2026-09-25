import pandas as pd
import pytest

from lineage_charisma.phylo import parse_tip_label
from lineage_charisma.taxonomy import (
    apply_bridges,
    elevated_name,
    epithet_stem,
    filter_clade,
    mdd_species_table,
    normalize_name,
    reconcile_tree,
    unused_tip_table,
)

MOCK_MDD = [
    (1, "Panthera_leo", "Felidae", "Panthera", "leo", "Lion", "Panthera_leo", None),
    (2, "Panthera_tigris", "Felidae", "Panthera", "tigris", "Tiger", "Panthera_tigris", None),
    (3, "Urva_javanica", "Herpestidae", "Urva", "javanica", "Javan Mongoose", "Herpestes_javanicus", None),
    (4, "Lycalopex_grisea", "Canidae", "Lycalopex", "grisea", "South American Gray Fox", "Lycalopex_griseus", None),
    (5, "Otaria_flavescens", "Otariidae", "Otaria", "flavescens", "South American Sea Lion", "Otaria_flavescens", None),
    (6, "Canis_lupus", "Canidae", "Canis", "lupus", "Gray Wolf", "Canis_lupus", None),
    (7, "Canis_lupaster", "Canidae", "Canis", "lupaster", "African Golden Wolf", None, None),
    (8, "Genetta_fieldiana", "Viverridae", "Genetta", "fieldiana", "Rusty-spotted Genet", None, None),
    (9, "Leopardus_colocola", "Felidae", "Leopardus", "colocola", "Pampas Cat", "Leopardus_colocolo", None),
    (10, "Leopardus_garleppi", "Felidae", "Leopardus", "garleppi", "Garlepp's Pampas Cat", None, None),
    (11, "Nasua_olivacea", "Procyonidae", "Nasua", "olivacea", "Western Mountain Coati", "Nasuella_olivacea", None),
    (12, "Mustela_nivalis", "Mustelidae", "Mustela", "nivalis", "Least Weasel", "Mustela_nivalis", None),
    (13, "Neogale_vison", "Mustelidae", "Neogale", "vison", "American Mink", "Neovison_vison", None),
    (14, "Proteles_cristatus", "Hyaenidae", "Proteles", "cristatus", "Aardwolf", None, None),
    (15, "Felis_catus", "Felidae", "Felis", "catus", "Domestic Cat", None, None),
    (16, "Mungos_mungo", "Herpestidae", "Mungos", "mungo", "Banded Mongoose", None, None),
]

MOCK_SYNONYMS = [
    (5, "Otaria flavescens", "byronia", "Phoca byronia", "Phoca byronia", "species", "synonym"),
    (8, "Genetta fieldiana", "maculata", "Viverra maculata", "Viverra maculata", "species", "synonym"),
    (11, "Nasua olivacea", "meridensis", "Nasuella meridensis", "Nasuella meridensis", "species", "synonym"),
    (12, "Mustela nivalis", "subpalmata", "Mustela subpalmata", "Mustela subpalmata", "species", "synonym"),
    (7, "Canis lupaster", "algirensis", "Canis aureus algirensis", "Canis aureus algirensis", "subspecies", "synonym"),
]

MOCK_TIPS = [
    "Panthera_leo_FELIDAE_CARNIVORA",
    "Panthera_tigris_FELIDAE_CARNIVORA",
    "Herpestes_javanicus_HERPESTIDAE_CARNIVORA",
    "Pseudalopex_griseus_CANIDAE_CARNIVORA",
    "Otaria_bryonia_OTARIIDAE_CARNIVORA",
    "Canis_lupus_CANIDAE_CARNIVORA",
    "Canis_anthus_CANIDAE_CARNIVORA",
    "Canis_aureus_CANIDAE_CARNIVORA",
    "Genetta_maculata_VIVERRIDAE_CARNIVORA",
    "Leopardus_colocolo_FELIDAE_CARNIVORA",
    "Nasuella_olivacea_PROCYONIDAE_CARNIVORA",
    "Nasuella_meridensis_PROCYONIDAE_CARNIVORA",
    "Mustela_nivalis_MUSTELIDAE_CARNIVORA",
    "Mustela_subpalmata_MUSTELIDAE_CARNIVORA",
    "Neovison_vison_MUSTELIDAE_CARNIVORA",
    "Proteles_cristata_HYAENIDAE_CARNIVORA",
    "Felis_catus_FELIDAE_CARNIVORA",
    "Mungos_mungu_HERPESTIDAE_CARNIVORA",
    "Smilodon_populator_FELIDAE_CARNIVORA",
]


@pytest.fixture
def species():
    raw = pd.DataFrame(
        MOCK_MDD,
        columns=["id", "sciName", "family", "genus", "specificEpithet", "mainCommonName", "MSW3_sciName", "CMW_sciName"],
    )
    raw["order"] = "Carnivora"
    raw["extinct"] = 0
    raw["domestic"] = (raw["sciName"] == "Felis_catus").astype(int)
    return mdd_species_table(raw)


@pytest.fixture
def synonyms():
    return pd.DataFrame(
        MOCK_SYNONYMS,
        columns=["mdd_id", "species", "root_name", "original_combination", "normalized_original_combination", "original_rank", "validity"],
    )


@pytest.fixture
def tips():
    return pd.DataFrame([parse_tip_label(t) for t in MOCK_TIPS])


@pytest.fixture
def result(species, synonyms, tips):
    return reconcile_tree(species, tips, synonyms, review_threshold=85, accept_threshold=92, min_margin=3)


@pytest.mark.parametrize(
    "raw,expected",
    [
        ("Panthera_leo", "Panthera leo"),
        ("  panthera   LEO ", "Panthera leo"),
        ("Canis lupus familiaris", "Canis lupus"),
        ("Canis cf. lupus", "Canis lupus"),
        ("Canis aff. lupus", "Canis lupus"),
        ("Canis sp.", None),
        ("Canis sp. 1", None),
        ("Mustela (Lutreola) lutreola", "Mustela lutreola"),
        ("incertae_sedis incertae_sedis", None),
        ("Canis", None),
        (None, None),
        (float("nan"), None),
    ],
)
def test_normalize_name(raw, expected):
    assert normalize_name(raw) == expected


def test_elevated_name_uses_last_epithet():
    assert elevated_name("Canis aureus algirensis") == "Canis algirensis"
    assert elevated_name("Phoca byronia") == "Phoca byronia"


@pytest.mark.parametrize(
    "a,b",
    [("griseus", "grisea"), ("cristatus", "cristata"), ("edwardsi", "edwardsii"), ("colocolo", "colocola"), ("brachyurus", "brachyura"), ("vitticollis", "vitticolla")],
)
def test_epithet_stem_handles_gender_and_spelling(a, b):
    assert epithet_stem(a) == epithet_stem(b)


def test_epithet_stem_keeps_distinct_names_apart():
    assert epithet_stem("lupus") != epithet_stem("lupaster")
    assert epithet_stem("leo") == "leo"


def test_filter_clade_is_case_insensitive():
    df = pd.DataFrame({"order": ["CARNIVORA", "Primates", "carnivora"], "x": [1, 2, 3]})
    assert filter_clade(df, "order", "Carnivora")["x"].tolist() == [1, 3]
    assert len(filter_clade(df, None, None)) == 3


def test_mdd_table_normalizes_names(species):
    assert species.loc[0, "species"] == "Panthera leo"
    assert species.set_index("species").loc["Urva javanica", "msw3_name"] == "Herpestes javanicus"


EXPECTED = {
    "Panthera leo": ("exact", "Panthera_leo_FELIDAE_CARNIVORA"),
    "Panthera tigris": ("exact", "Panthera_tigris_FELIDAE_CARNIVORA"),
    "Canis lupus": ("exact", "Canis_lupus_CANIDAE_CARNIVORA"),
    "Mustela nivalis": ("exact", "Mustela_nivalis_MUSTELIDAE_CARNIVORA"),
    "Felis catus": ("exact", "Felis_catus_FELIDAE_CARNIVORA"),
    "Urva javanica": ("legacy_name", "Herpestes_javanicus_HERPESTIDAE_CARNIVORA"),
    "Leopardus colocola": ("legacy_name", "Leopardus_colocolo_FELIDAE_CARNIVORA"),
    "Nasua olivacea": ("legacy_name", "Nasuella_olivacea_PROCYONIDAE_CARNIVORA"),
    "Neogale vison": ("legacy_name", "Neovison_vison_MUSTELIDAE_CARNIVORA"),
    "Genetta fieldiana": ("synonym", "Genetta_maculata_VIVERRIDAE_CARNIVORA"),
    "Lycalopex grisea": ("epithet_in_family", "Pseudalopex_griseus_CANIDAE_CARNIVORA"),
    "Proteles cristatus": ("epithet_in_family", "Proteles_cristata_HYAENIDAE_CARNIVORA"),
    "Otaria flavescens": ("fuzzy", "Otaria_bryonia_OTARIIDAE_CARNIVORA"),
}


@pytest.mark.parametrize("sp", sorted(EXPECTED))
def test_reconcile_methods(result, sp):
    a = result.assignments[sp]
    method, tip = EXPECTED[sp]
    assert a.status == "matched"
    assert a.match_method == method
    assert a.tree_tip == tip


def test_fuzzy_confidence_reflects_score(result):
    a = result.assignments["Otaria flavescens"]
    assert 0.92 <= a.confidence < 1.0
    assert a.matched_on == "Otaria byronia"


def test_trinomial_synonym_does_not_grab_parent_species_tip(result):
    a = result.assignments["Canis lupaster"]
    assert a.status == "no_tip"
    assert a.tree_tip is None


def test_split_species_without_tip_is_no_tip(result):
    assert result.assignments["Leopardus garleppi"].status == "no_tip"


def test_near_miss_goes_to_review(result):
    a = result.assignments["Mungos mungo"]
    assert a.status == "review"
    assert a.candidate_tip == "Mungos_mungu_HERPESTIDAE_CARNIVORA"
    assert 85 <= a.candidate_score < 92


def test_each_tip_claimed_at_most_once(result):
    tips = [a.tree_tip for a in result.assignments.values() if a.status == "matched"]
    assert len(tips) == len(set(tips))


def test_bridge_links_leftover_species_to_leftover_tip(result):
    apply_bridges(result, {"Canis lupaster": "Canis_anthus_CANIDAE_CARNIVORA"})
    a = result.assignments["Canis lupaster"]
    assert (a.status, a.match_method, a.tree_tip) == ("matched", "wiki_bridge", "Canis_anthus_CANIDAE_CARNIVORA")


def test_bridge_never_steals_a_claimed_tip(result):
    apply_bridges(result, {"Leopardus garleppi": "Leopardus_colocolo_FELIDAE_CARNIVORA"})
    assert result.assignments["Leopardus garleppi"].status == "no_tip"


def test_unused_tips_explained(result):
    unused = unused_tip_table(result).set_index("tree_binomial")
    assert unused.loc["Nasuella meridensis", "reason"] == "junior_synonym_of_mdd_species"
    assert unused.loc["Mustela subpalmata", "reason"] == "junior_synonym_of_mdd_species"
    assert unused.loc["Smilodon populator", "reason"] == "no_mdd_counterpart"
    assert "Canis anthus" in unused.index
