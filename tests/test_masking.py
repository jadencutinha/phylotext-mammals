from collections import Counter

import pandas as pd
import pytest

from lineage_charisma.config import load_config
from lineage_charisma.corpus import build_corpus
from lineage_charisma.masking import (
    MASK_LEVELS,
    Masker,
    build_vocabulary,
    derived_terms,
    head_nouns,
    plural_forms,
    term_leaks,
    text_column,
)

T = "[TAXON]"

SPECIES = pd.DataFrame({
    "species": ["Panthera leo", "Panthera pardus", "Lynx lynx", "Canis lupus", "Vulpes vulpes", "Poecilogale albinucha", "Urva javanica"],
    "order": ["Carnivora"] * 7,
    "suborder": ["Feliformia", "Feliformia", "Feliformia", "Caniformia", "Caniformia", "Caniformia", "Feliformia"],
    "family": ["Felidae", "Felidae", "Felidae", "Canidae", "Canidae", "Mustelidae", "Herpestidae"],
    "subfamily": ["Pantherinae", "Pantherinae", "Felinae", "Caninae", "Caninae", "Ictonychinae", "Herpestinae"],
    "genus": ["Panthera", "Panthera", "Lynx", "Canis", "Vulpes", "Poecilogale", "Urva"],
    "specific_epithet": ["leo", "pardus", "lynx", "lupus", "vulpes", "albinucha", "javanica"],
    "common_name": ["Lion", "Leopard", "Eurasian Lynx", "Gray Wolf", "Red Fox", "African Striped Weasel", "Javan Mongoose"],
    "other_common_names": ["African Lion", None, None, "Grey Wolf|Timber Wolf", None, "White-naped Weasel", "Small Asian Mongoose"],
})
SYNONYMS = pd.DataFrame({
    "root_name": ["leo", "persica", "major", "javanicus", "vulgaris"],
    "normalized_original_combination": ["Felis leo", "Felis leo persicus", "Canis major", "Herpestes javanicus", "Wolf vulgaris"],
})

TEXTS = {
    "Panthera leo": "The lion (Panthera leo) is a large cat of the genus Panthera and the family Felidae. P. leo and P. l. persica are felids; "
                    "lions were once placed in Felis. Its canine teeth are long. The major prey of lions is zebra.",
    "Panthera pardus": "The leopard (Panthera pardus) is a pantherine. Leopards are smaller than Panthera leo, and P. pardus climbs well.",
    "Lynx lynx": "The Eurasian lynx (Lynx lynx) is a medium-sized cat. Lynxes hunt hares; the lynx has tufted ears.",
    "Canis lupus": "The wolf (Canis lupus), also known as the gray wolf, is a canine native to Eurasia. Wolves have large canines and are the "
                   "largest of the Caninae. C. lupus is a caniform carnivoran.",
    "Vulpes vulpes": "The red fox (Vulpes vulpes) is the largest of the true foxes. Like other canids, V. vulpes has a bushy tail.",
    "Poecilogale albinucha": "The African striped weasel (Poecilogale albinucha) is a mustelid, the lone member of Poecilogale. "
                             "It may be related to the extinct Propoecilogale bolti.",
    "Urva javanica": "The Javan mongoose (Urva javanica), formerly Herpestes javanicus, has a skull with a slight curvature.",
}


@pytest.fixture(scope="module")
def masker():
    vocabulary = build_vocabulary(SPECIES, SYNONYMS, extra_terms=["carnivoran"], keep_terms=["canines"],
                                  keep_when_followed_by={"canine": ["tooth", "teeth"]})
    return Masker(vocabulary, T)


def test_helpers():
    assert derived_terms("Felidae") == {"felid", "feline"}
    assert derived_terms("Lutrinae") == {"lutrine"}
    assert derived_terms("Caniformia") == {"caniform"}
    assert plural_forms("wolf") >= {"wolf", "wolves"} and plural_forms("fox") >= {"fox", "foxes"}
    assert head_nouns(["Short-eared Dog", "Gray Wolf|Timber Wolf", "Vietnam Ferret-badger", None]) == {"dog", "wolf", "ferret-badger", "badger"}
    assert text_column("none") == "text" and text_column("strict") == "text_strict"
    with pytest.raises(ValueError):
        text_column("partial")


def test_none_level_returns_the_text_unchanged(masker):
    assert masker.mask(TEXTS["Panthera leo"], "none") == TEXTS["Panthera leo"]
    with pytest.raises(ValueError):
        masker.mask("x", "partial")


def test_taxonomy_masks_scientific_names_at_every_rank(masker):
    out = masker.mask(TEXTS["Panthera leo"], "taxonomy")
    assert out.startswith(f"The lion ({T}) is a large cat of the genus {T} and the family {T}.")
    for gone in ["Panthera", "leo", "Felidae", "felids", "Felis", "persica", "P. l."]:
        assert gone not in out
    assert "lion" in out and "cat" in out


def test_abbreviated_binomials_and_trinomials_are_masked(masker):
    assert masker.mask("P. leo and P. l. persica differ.", "taxonomy") == f"{T} and {T} differ."
    assert masker.mask("C. lupus is widespread.", "taxonomy") == f"{T} is widespread."
    assert masker.mask("In fig. B. the skull is shown.", "taxonomy") == "In fig. B. the skull is shown."


def test_matching_is_case_insensitive_and_covers_plurals_and_derived_terms(masker):
    assert masker.mask("PANTHERA and panthera; felids, a feline, a canid, Caninae, carnivorans.", "taxonomy") == f"{T} and {T}; {T}, a {T}, a {T}, {T}, {T}."
    assert masker.mask("Lynxes hunt; the lynx waits.", "taxonomy") == f"{T} hunt; the {T} waits."


def test_synonym_epithets_are_masked_only_inside_a_scientific_name(masker):
    assert masker.mask("The major prey of Canis major is deer.", "taxonomy") == f"The major prey of {T} is deer."


def test_a_name_hidden_behind_another_name_is_still_masked(masker):
    assert masker.mask("Felidae includes felids and Panthera.", "taxonomy") == f"{T} includes {T} and {T}."


def test_keep_rules_spare_teeth(masker):
    out = masker.mask(TEXTS["Canis lupus"], "taxonomy")
    assert "large canines" in out
    assert f"is a {T} native to Eurasia" in out
    assert masker.mask("Its canine teeth are long.", "taxonomy") == "Its canine teeth are long."


def test_compound_scientific_names_are_masked_whole(masker):
    out = masker.mask(TEXTS["Poecilogale albinucha"], "taxonomy")
    assert "poecilogale" not in out.lower()
    assert masker.mask(TEXTS["Urva javanica"], "taxonomy").endswith("a slight curvature.")


def test_strict_also_masks_common_name_head_nouns(masker):
    taxonomy = masker.mask(TEXTS["Canis lupus"], "taxonomy")
    strict = masker.mask(TEXTS["Canis lupus"], "strict")
    assert "wolf" in taxonomy and "Wolves" in taxonomy
    assert "wolf" not in strict.lower() and "wolves" not in strict.lower()
    assert strict.startswith(f"The {T} ({T}), also known as the gray {T}, is a {T} native to Eurasia.")
    # a historical genus spelled like a current common name is left to the strict level
    assert "felis" in masker.vocabulary.names and "wolf" not in masker.vocabulary.names
    assert masker.mask("The old genus Wolf and the wolf.", "taxonomy") == "The old genus Wolf and the wolf."
    assert masker.mask("The old genus Wolf and the wolf.", "strict") == f"The old genus {T} and the {T}."


def test_adjacent_masks_collapse_and_are_counted(masker):
    counts = Counter()
    out = masker.mask("Panthera leo persica (Felis leo)", "taxonomy", counts)
    assert out == f"{T} ({T})"
    assert masker.count(out) == 2
    assert counts[("binomial", "panthera leo persica")] == 1


@pytest.mark.parametrize("level", ["taxonomy", "strict"])
def test_no_masked_text_contains_its_own_genus_name(masker, level):
    for species, text in TEXTS.items():
        genus, epithet = species.split()
        masked = masker.mask(text, level)
        assert term_leaks(text, genus), species
        assert not term_leaks(masked, genus), (species, masked)
        assert not term_leaks(masked, epithet), (species, masked)


def test_term_leaks_definition():
    assert term_leaks("the pumas ran", "Puma") and term_leaks("Propoecilogale bolti", "Poecilogale")
    assert not term_leaks("a slight curvature", "Urva")
    assert not term_leaks("Canisius College", "Canis")


def test_build_corpus_adds_masked_columns_and_keeps_stub_flags_on_unmasked_text(masker):
    matched = SPECIES[["species", "family"]]
    corpus = build_corpus(matched, TEXTS, lambda texts: [len(t.split()) for t in texts], 12, masker, ["taxonomy", "none", "strict"])
    assert {"text", "text_taxonomy", "text_strict", "n_masked_taxonomy", "n_masked_strict"} <= set(corpus.columns)
    assert corpus["text"].tolist() == [TEXTS[s] for s in sorted(TEXTS)]
    assert (corpus["n_masked_strict"] >= corpus["n_masked_taxonomy"]).all() and (corpus["n_masked_taxonomy"] > 0).all()
    assert corpus.set_index("species_id")["n_tokens"].to_dict() == {s: len(t.split()) for s, t in TEXTS.items()}
    with pytest.raises(ValueError):
        build_corpus(matched, TEXTS, lambda texts: [1] * len(texts), 1, None, ["taxonomy"])


@pytest.mark.parametrize("level", [k for k in MASK_LEVELS if k != "none"])
def test_real_corpus_has_no_own_genus_in_masked_text(level):
    cfg = load_config()
    path = cfg.root / cfg["paths"]["processed"] / "corpus.parquet"
    if not path.exists():
        pytest.skip("corpus.parquet not built")
    corpus = pd.read_parquet(path)
    if text_column(level) not in corpus.columns:
        pytest.skip(f"corpus has no {level} column")
    leaking = [s for s, t in zip(corpus["species_id"], corpus[text_column(level)]) if term_leaks(t, s.split()[0])]
    assert leaking == []
    assert sum(term_leaks(t, s.split()[0]) for s, t in zip(corpus["species_id"], corpus["text"])) > 0.9 * len(corpus)
