from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd
from rapidfuzz import fuzz

from lineage_charisma.taxonomy import build_candidate_index, normalize_name, synonym_lookup

# docs/preregistration_amendment_secondary.md, section C (descriptive only)
TEMPLATE_COLUMNS = ["species_common", "species_scientific", "exhibit", "photo_filename", "transcription", "notes"]
SOURCE_ORDER = {"exact": 0, "legacy_name": 1, "synonym": 2}
MIN_H1_SPECIES = 10


def read_signage(path: Path) -> pd.DataFrame:
    """Read a transcription sheet in the template format; rows without a transcription or a scientific name are dropped and counted."""
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


def resolve_signage(signs: pd.DataFrame, species: pd.DataFrame, synonyms: pd.DataFrame | None, matching: Mapping[str, float]) -> pd.DataFrame:
    lookup, sources = name_index(species, synonyms)
    rows = [vars(resolve_name(p, lookup, sources, matching)) for p in signs["species_scientific"]]
    out = pd.concat([signs.reset_index(drop=True), pd.DataFrame(rows).rename(columns={"species": "species_id"})], axis=1)
    out["resolved"] = out["species_id"].notna()
    return out


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
