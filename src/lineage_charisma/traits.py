from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

from lineage_charisma.taxonomy import normalize_name, reconcile_tree

CONTINUOUS, CATEGORICAL, BINARY = "continuous", "categorical", "binary"
WEIGHTINGS = ("group", "column")


@dataclass(frozen=True)
class Trait:
    column: str        # column in the trait table
    group: str
    kind: str
    source: str        # "elton" or "combine"
    source_column: str
    description: str


ELTON_DIET = {
    "Diet-Inv": "invertebrates",
    "Diet-Vend": "mammals and birds",
    "Diet-Vect": "reptiles and amphibians",
    "Diet-Vfish": "fish",
    "Diet-Vunk": "vertebrates, type unknown",
    "Diet-Scav": "carrion",
    "Diet-Fruit": "fruit",
    "Diet-Nect": "nectar and pollen",
    "Diet-Seed": "seeds",
    "Diet-PlantO": "other plant material",
}

# the pre-registered trait list (docs/preregistration.md, section 7); nothing else enters the ecological distance
TRAITS: tuple[Trait, ...] = (
    *(Trait(f"diet_{c.split('-')[1].lower()}", "diet", CONTINUOUS, "elton", c, f"percentage of diet: {what}") for c, what in ELTON_DIET.items()),
    Trait("foraging_stratum", "foraging stratum", CATEGORICAL, "elton", "ForStrat-Value", "M marine, G ground, S scansorial, Ar arboreal, A aerial"),
    Trait("activity_nocturnal", "activity time", BINARY, "elton", "Activity-Nocturnal", "active at night"),
    Trait("activity_crepuscular", "activity time", BINARY, "elton", "Activity-Crepuscular", "active at dawn and dusk"),
    Trait("activity_diurnal", "activity time", BINARY, "elton", "Activity-Diurnal", "active by day"),
    Trait("habitat_breadth", "habitat breadth", CONTINUOUS, "combine", "habitat_breadth_n", "number of IUCN habitat types used"),
    Trait("terrestrial", "terrestrial or aquatic", BINARY, "combine", "terrestrial_non-volant", "lives on land"),
    Trait("marine", "terrestrial or aquatic", BINARY, "combine", "marine", "lives in the sea"),
    Trait("freshwater", "terrestrial or aquatic", BINARY, "combine", "freshwater", "lives in fresh water"),
    Trait("log_body_mass", "body mass", CONTINUOUS, "combine", "adult_mass_g", "log10 of adult body mass in grams"),
)


def md5(path: Path) -> str:
    return hashlib.md5(Path(path).read_bytes()).hexdigest()


def load_elton(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", encoding="latin-1")
    df = df[df["Scientific"].notna()].copy()   # the file ends with blank rows
    df["source_name"] = df["Scientific"].str.strip()
    df["family"] = df["MSWFamilyLatin"].str.strip()
    return df.reset_index(drop=True)


def load_combine(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, low_memory=False)
    # COMBINE has one row per PHYLACINE name, so an IUCN species that lumps several of them appears on several
    # rows carrying the same trait values; taxa IUCN does not list (mostly extinct) keep their PHYLACINE name
    iucn = df["iucn2020_binomial"].str.strip()
    df["source_name"] = iucn.where(iucn != "Not recognised", df["phylacine_binomial"].str.strip())
    return df.drop_duplicates("source_name").reset_index(drop=True)


def as_tips(table: pd.DataFrame, families: Sequence[str]) -> pd.DataFrame:
    """A trait table's names in the shape `reconcile_tree` expects for tree tips, limited to the clade's families."""
    wanted = {f.upper() for f in families}
    rows = []
    for name, family in zip(table["source_name"], table["family"]):
        binomial = normalize_name(name)
        if binomial is None or str(family).upper() not in wanted:
            continue
        genus, epithet = binomial.split(" ", 1)
        rows.append({"tip_label": name, "binomial": binomial, "genus": genus, "epithet": epithet, "family": str(family).upper()})
    tips = pd.DataFrame(rows, columns=["tip_label", "binomial", "genus", "epithet", "family"])
    return tips.drop_duplicates("binomial").reset_index(drop=True)


def reconcile_traits(species: pd.DataFrame, synonyms: pd.DataFrame | None, table: pd.DataFrame, matching: dict) -> pd.DataFrame:
    """Match MDD species to one trait database's names with the Week 2 engine (exact, legacy name, synonym, epithet in family, fuzzy)."""
    tips = as_tips(table, species["family"].unique())
    result = reconcile_tree(
        species, tips, synonyms,
        review_threshold=matching["fuzzy_review_threshold"],
        accept_threshold=matching["fuzzy_accept_threshold"],
        min_margin=matching["fuzzy_min_margin"],
    )
    out = result.table().rename(columns={"tree_tip": "source_name", "candidate_tip": "candidate_name"})
    out["note"] = out["note"].str.replace("no tree tip", "no name in the trait database", regex=False)
    return out[["species", "status", "source_name", "match_method", "confidence", "matched_on", "candidate_name", "candidate_score", "note"]]


def build_trait_table(order: Sequence[str], elton: pd.DataFrame, combine: pd.DataFrame, elton_map: pd.DataFrame, combine_map: pd.DataFrame) -> pd.DataFrame:
    """One row per species in `order`, one column per pre-registered trait; NaN where a database has no record or no value."""
    sources = {"elton": (elton.set_index("source_name"), elton_map), "combine": (combine.set_index("source_name"), combine_map)}
    out = pd.DataFrame(index=pd.Index(order, name="species_id"))
    for trait in TRAITS:
        table, mapping = sources[trait.source]
        matched = mapping[mapping["status"] == "matched"].set_index("species")["source_name"]
        names = matched.reindex(order)
        values = pd.Series(table[trait.source_column].reindex(names.to_numpy()).to_numpy(), index=order)
        if trait.column == "log_body_mass":
            values = np.log10(pd.to_numeric(values, errors="coerce").where(lambda v: v > 0))
        elif trait.kind != CATEGORICAL:
            values = pd.to_numeric(values, errors="coerce")
        out[trait.column] = values.to_numpy()
    return out


def trait_weights(traits: Sequence[Trait], weighting: str) -> np.ndarray:
    """`group`: every trait group counts once and its columns share that weight. `column`: every column counts once."""
    if weighting not in WEIGHTINGS:
        raise ValueError(f"unknown weighting {weighting!r}; expected one of {WEIGHTINGS}")
    if weighting == "column":
        return np.ones(len(traits))
    sizes = pd.Series([t.group for t in traits]).value_counts()
    return np.array([1.0 / sizes[t.group] for t in traits])


def gower_distance(table: pd.DataFrame, traits: Sequence[Trait] = TRAITS, weighting: str = "group") -> tuple[np.ndarray, np.ndarray]:
    """Gower distance with the pairwise-available rule.

    Continuous traits contribute |difference| / range, categorical and binary traits 0 for a match
    and 1 for a mismatch. A pair's distance is the weighted mean over the traits both species have.
    Returns the distance matrix (NaN where a pair shares no trait) and the matrix of shared weight
    as a fraction of the total.
    """
    weights = trait_weights(traits, weighting)
    n = len(table)
    total, shared = np.zeros((n, n)), np.zeros((n, n))
    for trait, w in zip(traits, weights):
        col = table[trait.column]
        have = col.notna().to_numpy()
        both = have[:, None] & have[None, :]
        if trait.kind == CONTINUOUS:
            v = col.to_numpy(dtype=np.float64)
            span = np.nanmax(v) - np.nanmin(v) if have.any() else 0.0
            diff = np.abs(v[:, None] - v[None, :]) / span if span > 0 else np.zeros((n, n))
        else:
            v = col.to_numpy(dtype=object)
            diff = (v[:, None] != v[None, :]).astype(np.float64)
        total += np.where(both, w * diff, 0.0)
        shared += np.where(both, w, 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        dist = np.where(shared > 0, total / shared, np.nan)
    np.fill_diagonal(dist, np.where(np.diag(shared) > 0, 0.0, np.nan))
    return dist, shared / weights.sum()


def missing_rates(table: pd.DataFrame, traits: Sequence[Trait] = TRAITS) -> pd.DataFrame:
    return pd.DataFrame([{"trait": t.column, "group": t.group, "type": t.kind, "source": t.source, "source_column": t.source_column,
                          "n_missing": int(table[t.column].isna().sum()), "missing_rate": float(table[t.column].isna().mean())} for t in traits])


def complete_species(dist: np.ndarray, order: Sequence[str]) -> list[str]:
    """The largest set of species, found by dropping the worst first, whose pairwise distances are all defined."""
    keep = np.ones(len(order), dtype=bool)
    bad = np.isnan(dist)
    while bad[np.ix_(keep, keep)].any():
        counts = np.where(keep, (bad & keep[None, :]).sum(axis=1), -1)
        keep[int(counts.argmax())] = False
    return [s for s, k in zip(order, keep) if k]
