import re
from pathlib import Path

import numpy as np
import pytest

from lineage_charisma.charisma import (
    DIMENSIONS,
    OUTPUT_SCHEMA,
    SYSTEM_PROMPT,
    InvalidResponse,
    charisma_distance,
    user_message,
    validate_scores,
    validation_sample,
)

AMENDMENT = Path(__file__).resolve().parents[1] / "docs" / "preregistration_amendment_secondary.md"


def test_valid_response_is_parsed():
    assert validate_scores('{"warmth": 0, "anthropomorphism": 3, "aesthetic": 1, "human_relevance": 2}') == {
        "warmth": 0, "anthropomorphism": 3, "aesthetic": 1, "human_relevance": 2}


@pytest.mark.parametrize("raw", [
    "not json",
    '{"warmth": 1, "anthropomorphism": 1, "aesthetic": 1',                                  # truncated
    "[1, 2, 3, 4]",                                                                          # not an object
    '{"warmth": 1, "anthropomorphism": 1, "aesthetic": 1}',                                  # missing key
    '{"warmth": 1, "anthropomorphism": 1, "aesthetic": 1, "human_relevance": 1, "x": 0}',   # extra key
    '{"warmth": 4, "anthropomorphism": 1, "aesthetic": 1, "human_relevance": 1}',           # above range
    '{"warmth": -1, "anthropomorphism": 1, "aesthetic": 1, "human_relevance": 1}',          # below range
    '{"warmth": 1.5, "anthropomorphism": 1, "aesthetic": 1, "human_relevance": 1}',         # not an integer
    '{"warmth": "2", "anthropomorphism": 1, "aesthetic": 1, "human_relevance": 1}',         # string
    '{"warmth": true, "anthropomorphism": 1, "aesthetic": 1, "human_relevance": 1}',        # boolean
])
def test_malformed_or_out_of_range_responses_are_rejected(raw):
    with pytest.raises(InvalidResponse):
        validate_scores(raw)


def test_schema_matches_the_validator():
    assert OUTPUT_SCHEMA["required"] == list(DIMENSIONS)
    assert OUTPUT_SCHEMA["additionalProperties"] is False
    assert all(OUTPUT_SCHEMA["properties"][d]["enum"] == [0, 1, 2, 3] for d in DIMENSIONS)


def test_prompt_text_is_the_registered_text():
    blocks = re.findall(r"```\n(.*?)\n```", AMENDMENT.read_text(), flags=re.S)
    system, user = blocks[0], blocks[1]
    assert SYSTEM_PROMPT == system
    assert user_message("{TEXT}") == user


def test_validation_sample_is_seeded_and_fixed():
    order = [f"S{i:03d}" for i in range(284)]
    a, b = validation_sample(order, 20261004), validation_sample(order, 20261004)
    assert a == b and len(set(a)) == 50
    idx = np.random.default_rng(20261004).choice(284, size=50, replace=False)
    assert a == [order[i] for i in idx]


def test_charisma_distance_is_symmetric_with_zero_diagonal():
    rng = np.random.default_rng(0)
    scores = rng.integers(0, 4, size=(30, 3)).astype(float)
    d = charisma_distance(scores)
    assert d.shape == (30, 30)
    np.testing.assert_array_equal(d, d.T)
    np.testing.assert_array_equal(np.diag(d), 0.0)
    assert (d >= 0).all()


def test_charisma_distance_is_euclidean_on_standardized_scores():
    scores = np.array([[0.0, 1.0], [2.0, 1.5], [3.0, 3.0], [1.0, 0.0]])
    z = (scores - scores.mean(0)) / scores.std(0)
    d = charisma_distance(scores)
    assert d[0, 2] == pytest.approx(np.linalg.norm(z[0] - z[2]))
    one = charisma_distance(scores[:, 0])
    assert one[1, 3] == pytest.approx(abs(z[1, 0] - z[3, 0]))


def test_charisma_distance_follows_species_order():
    from lineage_charisma.config import load_config
    from lineage_charisma.embed import read_species_order

    cfg = load_config()
    if not cfg.species_order_path().exists():
        pytest.skip("pilot species order not built")
    order = read_species_order(cfg.species_order_path())
    rng = np.random.default_rng(1)
    by_species = {s: rng.random(2) * 3 for s in reversed(order)}   # built in a different order on purpose
    d = charisma_distance(np.stack([by_species[s] for s in order]))
    z = np.stack([by_species[s] for s in order])
    z = (z - z.mean(0)) / z.std(0)
    i, j = order.index("Panthera leo"), order.index("Canis lupus")
    assert d.shape == (len(order), len(order))
    assert d[i, j] == pytest.approx(np.linalg.norm(z[i] - z[j]))


def test_no_variance_or_no_dimension_is_an_error():
    with pytest.raises(ValueError):
        charisma_distance(np.ones((5, 2)))
    with pytest.raises(ValueError):
        charisma_distance(np.zeros((5, 0)))
