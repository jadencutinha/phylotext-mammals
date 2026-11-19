from __future__ import annotations

import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd
from rapidfuzz import fuzz

from lineage_charisma.taxonomy import build_candidate_index, normalize_name, synonym_lookup

# docs/preregistration_amendment_secondary.md, section C (descriptive only); the common-name route for signs with no printed
# scientific name is docs/preregistration_amendment_extensions.md, section D
TEMPLATE_COLUMNS = ["species_common", "species_scientific", "exhibit", "photo_filename", "transcription", "notes"]
SOURCE_ORDER = {"exact": 0, "legacy_name": 1, "synonym": 2}
MIN_H1_SPECIES = 10
COMMON_NAME_COLUMNS = ("common_name", "other_common_names")
SCIENTIFIC, COMMON = "scientific name", "common name"
_DASHES = re.compile("[\u002d\u2010-\u2015\u2212]")


def read_signage(path: Path) -> pd.DataFrame:
    """Read a transcription sheet in the template format, with every field stripped and whitespace in transcriptions collapsed."""
    df = pd.read_csv(path, dtype=str, keep_default_na=False)
    missing = [c for c in TEMPLATE_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"{Path(path).name} is missing template columns {missing}; expected {TEMPLATE_COLUMNS}")
    df = df[TEMPLATE_COLUMNS].copy()
    for c in TEMPLATE_COLUMNS:
        df[c] = df[c].str.strip()
    df["transcription"] = df["transcription"].str.split().str.join(" ")
    df.insert(0, "row", range(1, len(df) + 1))
    return df


@dataclass(frozen=True)
class Resolution:
    species: str | None
    method: str
    matched_on: str | None
    score: float | None
    note: str


def name_index(species: pd.DataFrame, synonyms: pd.DataFrame | None) -> tuple[dict[str, set[str]], dict[str, dict[str, str]]]:
    """MDD names with the Week 2 engine's candidate rules: name -> species, and species -> {name: source}."""
    candidates = build_candidate_index(species, synonyms)
    return synonym_lookup(candidates), {sp: dict(c) for sp, c in candidates.items()}


def resolve_name(printed: str, lookup: Mapping[str, set[str]], sources: Mapping[str, Mapping[str, str]], matching: Mapping[str, float]) -> Resolution:
    """Resolve one printed scientific name to one MDD species: current name, then legacy and synonym names, then fuzzy at the engine's thresholds."""
    name = normalize_name(printed)
    if name is None:
        return Resolution(None, "unparsed", None, None, f"{printed!r} is not a binomial")
    hits = lookup.get(name, set())
    if hits:
        ranked = sorted(hits, key=lambda sp: SOURCE_ORDER.get(sources[sp].get(name, "synonym"), 3))
        best = SOURCE_ORDER.get(sources[ranked[0]].get(name, "synonym"), 3)
        tied = [sp for sp in ranked if SOURCE_ORDER.get(sources[sp].get(name, "synonym"), 3) == best]
        if len(tied) > 1:
            return Resolution(None, "ambiguous", name, None, f"{name} names several MDD species: {tied}")
        return Resolution(tied[0], sources[tied[0]].get(name, "synonym"), name, 100.0, "")
    scores = sorted(((fuzz.ratio(name.lower(), cand.lower()), cand) for cand in lookup), reverse=True)
    best, cand = scores[0]
    second = next((s for s, c in scores[1:] if lookup[c] != lookup[cand]), 0.0)
    if best >= matching["fuzzy_accept_threshold"] and best - second >= matching["fuzzy_min_margin"] and len(lookup[cand]) == 1:
        return Resolution(next(iter(lookup[cand])), "fuzzy", cand, round(best, 1), "")
    status = "review" if best >= matching["fuzzy_review_threshold"] else "no_match"
    return Resolution(None, status, cand, round(best, 1), f"best fuzzy candidate {cand!r} at {best:.1f}, runner-up {second:.1f}")


def normalize_common_name(name: str) -> str:
    """Amendment 5, section D: lowercase, hyphens and dashes to spaces, other non-letter/digit/space characters deleted, spaces collapsed.
    Accents are letters and are kept."""
    t = _DASHES.sub(" ", name.lower())
    t = "".join(ch for ch in t if ch.isalnum() or ch.isspace())
    return " ".join(t.split())


def common_name_index(species: pd.DataFrame) -> dict[str, set[str]]:
    """Normalized MDD main and other common names -> the species carrying them (other names are `|`-separated)."""
    index: dict[str, set[str]] = defaultdict(set)
    for row in species[["species", *COMMON_NAME_COLUMNS]].itertuples(index=False):
        names = [row.common_name] if isinstance(row.common_name, str) else []
        names += row.other_common_names.split("|") if isinstance(row.other_common_names, str) else []
        for name in names:
            key = normalize_common_name(name)
            if key:
                index[key].add(row.species)
    return dict(index)


def resolve_common_name(printed: str, index: Mapping[str, set[str]]) -> Resolution:
    """Exact match after normalization; a name carried by several species is ambiguous, with no precedence between main and other names."""
    key = normalize_common_name(printed)
    hits = sorted(index.get(key, set()))
    if len(hits) == 1:
        return Resolution(hits[0], "common_name", key, 100.0, "")
    if hits:
        return Resolution(None, "ambiguous_common_name", key, None, f"{printed!r} is a common name of several MDD species: {hits}")
    return Resolution(None, "no_common_name_match", key, None, f"{printed!r} equals no MDD common name after normalization")


def resolve_signage(signs: pd.DataFrame, species: pd.DataFrame, synonyms: pd.DataFrame | None, matching: Mapping[str, float]) -> pd.DataFrame:
    """Each row by its printed scientific name; a row with none printed by its printed common name (amendment 5, section D)."""
    lookup, sources = name_index(species, synonyms)
    common = common_name_index(species) if (signs["species_scientific"] == "").any() else {}
    rows = []
    for sci, com in zip(signs["species_scientific"], signs["species_common"]):
        if sci:
            rows.append({**vars(resolve_name(sci, lookup, sources, matching)), "route": SCIENTIFIC})
        else:
            rows.append({**vars(resolve_common_name(com, common)), "route": COMMON})
    out = pd.concat([signs.reset_index(drop=True), pd.DataFrame(rows).rename(columns={"species": "species_id"})], axis=1)
    out["resolved"] = out["species_id"].notna()
    return out


def route_counts(resolved: pd.DataFrame) -> dict[str, int]:
    """Amendment 5, section D: signed species matched by each route, and the common-name rows excluded."""
    ok = resolved[resolved["resolved"]]
    by_species = ok.groupby("species_id")["route"].agg(lambda r: SCIENTIFIC if SCIENTIFIC in set(r) else COMMON)
    com = resolved[resolved["route"] == COMMON]
    return {"n_species_by_scientific_name": int((by_species == SCIENTIFIC).sum()),
            "n_species_by_common_name": int((by_species == COMMON).sum()),
            "n_common_name_rows_ambiguous": int((com["method"] == "ambiguous_common_name").sum()),
            "n_common_name_rows_no_match": int((com["method"] == "no_common_name_match").sum())}


def cosine(a: np.ndarray, b: np.ndarray) -> float:
    return float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b)))


def similarity_table(species: Sequence[str], sign: Mapping[str, np.ndarray], wiki: Mapping[str, np.ndarray],
                     adw: Mapping[str, np.ndarray] | None = None) -> pd.DataFrame:
    """Within-species cross-register cosine similarity and each sign's between-species baseline (section C, item 1)."""
    adw = adw or {}
    with_wiki = [s for s in species if s in wiki]
    rows = []
    for s in species:
        others = [o for o in with_wiki if o != s]
        rows.append({
            "species_id": s,
            "sign_wiki": cosine(sign[s], wiki[s]) if s in wiki else np.nan,
            "sign_adw": cosine(sign[s], adw[s]) if s in adw else np.nan,
            "wiki_adw": cosine(wiki[s], adw[s]) if s in wiki and s in adw else np.nan,
            "baseline_sign_other_wiki": float(np.mean([cosine(sign[s], wiki[o]) for o in others])) if others else np.nan,
        })
    table = pd.DataFrame(rows)
    table["sign_wiki_above_baseline"] = table["sign_wiki"] > table["baseline_sign_other_wiki"]
    return table


def similarity_summary(table: pd.DataFrame) -> dict[str, float | int]:
    defined = table.dropna(subset=["sign_wiki", "baseline_sign_other_wiki"])
    return {"n_signed_species": len(table), "n_with_wikipedia": int(table["sign_wiki"].notna().sum()),
            "n_with_adw": int(table["sign_adw"].notna().sum()),
            "median_sign_wiki": float(table["sign_wiki"].median()), "median_baseline": float(table["baseline_sign_other_wiki"].median()),
            "median_sign_adw": float(table["sign_adw"].median()), "median_wiki_adw": float(table["wiki_adw"].median()),
            "share_sign_wiki_above_baseline": float(defined["sign_wiki_above_baseline"].mean()) if len(defined) else float("nan")}
