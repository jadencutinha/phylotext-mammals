import numpy as np
import pandas as pd
import pytest

from lineage_charisma.distances import validate_distance_matrix
from lineage_charisma.traits import BINARY, CATEGORICAL, CONTINUOUS, TRAITS, Trait, as_tips, build_trait_table, complete_species, gower_distance, missing_rates, reconcile_traits, trait_weights

MIXED = (
    Trait("mass", "mass", CONTINUOUS, "x", "mass", ""),
    Trait("stratum", "stratum", CATEGORICAL, "x", "stratum", ""),
    Trait("d1", "diet", CONTINUOUS, "x", "d1", ""),
    Trait("d2", "diet", CONTINUOUS, "x", "d2", ""),
    Trait("marine", "water", BINARY, "x", "marine", ""),
)


@pytest.fixture
def table():
    return pd.DataFrame({
        "mass": [1.0, 3.0, 5.0, np.nan],
        "stratum": ["G", "G", "Ar", "M"],
        "d1": [0.0, 100.0, 50.0, 0.0],
        "d2": [100.0, 0.0, 50.0, np.nan],
        "marine": [0.0, 0.0, 0.0, 1.0],
    }, index=["a", "b", "c", "d"])


def test_gower_mixed_types_known_answer(table):
    dist, shared = gower_distance(table, MIXED, weighting="column")
    # a-b: mass 2/4, stratum 0, d1 1, d2 1, marine 0
    assert dist[0, 1] == pytest.approx((0.5 + 0 + 1 + 1 + 0) / 5)
    # a-c: mass 1, stratum 1, d1 0.5, d2 0.5, marine 0
    assert dist[0, 2] == pytest.approx((1 + 1 + 0.5 + 0.5 + 0) / 5)
    assert shared[0, 1] == pytest.approx(1.0)


def test_gower_group_weighting_shares_weight_inside_a_group(table):
    np.testing.assert_allclose(trait_weights(MIXED, "group"), [1, 1, 0.5, 0.5, 1])
    dist, _ = gower_distance(table, MIXED, weighting="group")
    assert dist[0, 1] == pytest.approx((0.5 + 0 + 0.5 * 1 + 0.5 * 1 + 0) / 4)
    with pytest.raises(ValueError):
        trait_weights(MIXED, "equal")


def test_gower_missing_values_use_only_shared_traits(table):
    dist, shared = gower_distance(table, MIXED, weighting="column")
    # a-d: mass and d2 are missing for d, leaving stratum 1, d1 0, marine 1
    assert dist[0, 3] == pytest.approx((1 + 0 + 1) / 3)
    assert shared[0, 3] == pytest.approx(3 / 5)
    grouped, grouped_shared = gower_distance(table, MIXED, weighting="group")
    assert grouped[0, 3] == pytest.approx((1 + 0.5 * 0 + 1) / 2.5)
    assert grouped_shared[0, 3] == pytest.approx(2.5 / 4)


def test_gower_is_a_valid_distance_matrix(table):
    for weighting in ("group", "column"):
        dist, _ = gower_distance(table, MIXED, weighting)
        validate_distance_matrix(dist)
        assert dist.max() <= 1.0
        np.testing.assert_array_equal(dist, dist.T)
        np.testing.assert_array_equal(np.diag(dist), 0.0)


def test_gower_does_not_impute(table):
    # changing a value in a column where d is missing must not change any of d's distances
    before, _ = gower_distance(table, MIXED, "column")
    table.loc["a", "mass"] = 2.0
    after, _ = gower_distance(table, MIXED, "column")
    np.testing.assert_array_equal(before[3, 1:], after[3, 1:])
    assert before[0, 3] == after[0, 3]


def test_pairs_with_no_shared_trait_are_nan_and_dropped():
    t = pd.DataFrame({"mass": [1.0, 2.0, np.nan, np.nan], "stratum": ["G", "G", "M", None], "d1": [np.nan] * 4, "d2": [np.nan] * 4,
                      "marine": [np.nan] * 4}, index=list("abcd"))
    dist, shared = gower_distance(t, MIXED, "column")
    assert np.isnan(dist[0, 3]) and np.isnan(dist[3, 3]) and shared[0, 3] == 0
    assert dist[0, 2] == 1.0
    assert complete_species(dist, list("abcd")) == ["a", "b", "c"]


def test_constant_continuous_trait_contributes_zero():
    t = pd.DataFrame({"mass": [2.0, 2.0, 2.0], "stratum": ["G", "G", "M"], "d1": [0.0, 0, 0], "d2": [0.0, 0, 0], "marine": [0.0, 0, 1]}, index=list("abc"))
    dist, _ = gower_distance(t, MIXED, "column")
    assert dist[0, 1] == 0.0 and dist[0, 2] == pytest.approx(2 / 5)


def test_registered_traits_match_the_preregistration():
    groups = pd.Series([t.group for t in TRAITS]).value_counts().to_dict()
    assert groups == {"diet": 10, "activity time": 3, "terrestrial or aquatic": 3, "foraging stratum": 1, "habitat breadth": 1, "body mass": 1}
    assert trait_weights(TRAITS, "group").sum() == pytest.approx(6.0)
    assert len({t.column for t in TRAITS}) == len(TRAITS)


def test_reconcile_and_build_trait_table():
    species = pd.DataFrame({"mdd_id": [1, 2, 3], "species": ["Urva edwardsii", "Panthera leo", "Canis lupaster"], "genus": ["Urva", "Panthera", "Canis"],
                            "family": ["Herpestidae", "Felidae", "Canidae"], "msw3_name": ["Herpestes edwardsii", "Panthera leo", None], "cmw_name": [None, None, None]})
    elton = pd.DataFrame({"source_name": ["Herpestes edwardsii", "Panthera leo", "Mus musculus"], "family": ["Herpestidae", "Felidae", "Muridae"],
                          **{t.source_column: [10.0, 0.0, 5.0] for t in TRAITS if t.source == "elton" and t.kind != CATEGORICAL}, "ForStrat-Value": ["G", "G", "G"]})
    combine = pd.DataFrame({"source_name": ["Urva edwardsii", "Panthera leo", "Canis lupaster"], "family": ["Herpestidae", "Felidae", "Canidae"],
                            "habitat_breadth_n": [3.0, 5.0, np.nan], "terrestrial_non-volant": [1.0, 1.0, 1.0], "marine": [0.0, 0.0, 0.0], "freshwater": [0.0, 0.0, 0.0],
                            "adult_mass_g": [1000.0, 100000.0, 0.0]})
    matching = dict(fuzzy_review_threshold=85, fuzzy_accept_threshold=92, fuzzy_min_margin=3)
    assert list(as_tips(elton, species["family"].unique())["binomial"]) == ["Herpestes edwardsii", "Panthera leo"]
    emap, cmap = reconcile_traits(species, None, elton, matching), reconcile_traits(species, None, combine, matching)
    methods = emap.set_index("species")["match_method"]
    assert methods["Urva edwardsii"] == "legacy_name" and methods["Panthera leo"] == "exact" and pd.isna(methods["Canis lupaster"])
    assert (cmap["status"] == "matched").all()

    order = ["Canis lupaster", "Panthera leo", "Urva edwardsii"]
    table = build_trait_table(order, elton, combine, emap, cmap)
    assert list(table.index) == order and list(table.columns) == [t.column for t in TRAITS]
    assert table.loc["Urva edwardsii", "diet_inv"] == 10.0 and table.loc["Panthera leo", "log_body_mass"] == pytest.approx(5.0)
    assert table.loc["Canis lupaster", [t.column for t in TRAITS if t.source == "elton"]].isna().all()
    assert np.isnan(table.loc["Canis lupaster", "log_body_mass"])   # a mass of 0 is not a value
    rates = missing_rates(table).set_index("trait")
    assert rates.loc["diet_inv", "n_missing"] == 1 and rates.loc["marine", "n_missing"] == 0
    dist, _ = gower_distance(table)
    validate_distance_matrix(dist)
