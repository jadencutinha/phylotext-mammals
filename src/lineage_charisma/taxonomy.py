from __future__ import annotations

import math
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import pandas as pd
from rapidfuzz import fuzz

MDD_COLUMNS = {
    "id": "mdd_id",
    "sciName": "species",
    "order": "order",
    "suborder": "suborder",
    "family": "family",
    "subfamily": "subfamily",
    "genus": "genus",
    "specificEpithet": "specific_epithet",
    "mainCommonName": "common_name",
    "otherCommonNames": "other_common_names",
    "originalNameCombination": "original_combination",
    "nominalNames": "nominal_names",
    "MSW3_sciName": "msw3_name",
    "CMW_sciName": "cmw_name",
    "iucnStatus": "iucn_status",
    "extinct": "extinct",
    "domestic": "domestic",
}

SYNONYM_COLUMNS = {
    "MDD_species_id": "mdd_id",
    "MDD_species": "species",
    "MDD_root_name": "root_name",
    "MDD_original_combination": "original_combination",
    "MDD_normalized_original_combination": "normalized_original_combination",
    "MDD_original_rank": "original_rank",
    "MDD_validity": "validity",
    "MDD_nomenclature_status": "nomenclature_status",
}

QUALIFIERS = {"cf", "aff", "sp", "spp", "ssp", "subsp", "var", "nr", "n", "gen", "indet"}
LATIN_ENDINGS = ("ae", "us", "um", "is", "er", "a", "e", "i", "o")

METHOD_CONFIDENCE = {
    "exact": 1.0,
    "legacy_name": 0.98,
    "synonym": 0.95,
    "epithet_in_family": 0.9,
    "wiki_bridge": 0.85,
}


def _clean_token(token: str) -> str:
    return token.strip(".,;:?!\"'()[]")


def name_tokens(name: object) -> list[str]:
    if name is None or (isinstance(name, float) and math.isnan(name)):
        return []
    text = str(name).replace("_", " ")
    text = re.sub(r"\([^)]*\)", " ", text)
    out = []
    for raw in text.split():
        tok = _clean_token(raw)
        if not tok or tok.lower().rstrip(".") in QUALIFIERS or tok == "×":
            continue
        out.append(tok)
    return out


def _valid_parts(genus: str, epithet: str) -> bool:
    if genus.lower() in {"incertae", "unknown", "nan", "na"}:
        return False
    return genus.isalpha() and epithet.replace("-", "").isalpha()


def normalize_name(name: object) -> str | None:
    toks = name_tokens(name)
    if len(toks) < 2 or not _valid_parts(toks[0], toks[1]):
        return None
    return f"{toks[0].capitalize()} {toks[1].lower()}"


def elevated_name(name: object) -> str | None:
    toks = name_tokens(name)
    if len(toks) < 2 or not _valid_parts(toks[0], toks[-1]):
        return None
    return f"{toks[0].capitalize()} {toks[-1].lower()}"


def epithet_stem(epithet: str) -> str:
    e = epithet.lower().replace("-", "")
    if e.endswith("ii"):
        e = e[:-1]
    for ending in LATIN_ENDINGS:
        if e.endswith(ending) and len(e) - len(ending) >= 3:
            return e[: -len(ending)]
    return e


def split_binomial(name: str) -> tuple[str, str]:
    genus, epithet = name.split(" ", 1)
    return genus, epithet


def load_mdd(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, low_memory=False)


def filter_clade(df: pd.DataFrame, rank: str | None, name: str | None) -> pd.DataFrame:
    if not rank or not name:
        return df.copy()
    if rank not in df.columns:
        raise KeyError(f"rank column {rank!r} not in table")
    return df[df[rank].astype(str).str.lower() == str(name).lower()].copy()


def mdd_species_table(df: pd.DataFrame) -> pd.DataFrame:
    keep = [c for c in MDD_COLUMNS if c in df.columns]
    out = df[keep].rename(columns=MDD_COLUMNS).copy()
    out["species"] = out["species"].map(normalize_name)
    for col in ("msw3_name", "cmw_name"):
        if col in out:
            out[col] = out[col].map(normalize_name)
    for col in ("extinct", "domestic"):
        if col in out:
            out[col] = out[col].fillna(0).astype(int)
    return out.reset_index(drop=True)


def load_mdd_synonyms(path: Path, mdd_ids: Iterable[int]) -> pd.DataFrame:
    ids = set(int(i) for i in mdd_ids)
    df = pd.read_csv(path, usecols=list(SYNONYM_COLUMNS), low_memory=False)
    df = df[df["MDD_species_id"].isin(ids)]
    return df.rename(columns=SYNONYM_COLUMNS).reset_index(drop=True)


def name_candidates(row: pd.Series, synonyms: pd.DataFrame | None = None) -> list[tuple[str, str]]:
    seen: dict[str, str] = {}

    def add(name: str | None, source: str) -> None:
        if isinstance(name, str) and name and name not in seen:
            seen[name] = source

    add(row["species"], "exact")
    add(row.get("msw3_name"), "legacy_name")
    add(row.get("cmw_name"), "legacy_name")
    genus = row["genus"]
    if synonyms is not None and len(synonyms):
        for syn in synonyms.itertuples(index=False):
            for combo in (syn.normalized_original_combination, syn.original_combination):
                add(elevated_name(combo), "synonym")
            if isinstance(syn.root_name, str):
                add(normalize_name(f"{genus} {syn.root_name}"), "synonym")
    add(elevated_name(row.get("original_combination")), "synonym")
    return list(seen.items())


def build_candidate_index(species: pd.DataFrame, synonyms: pd.DataFrame | None) -> dict[str, list[tuple[str, str]]]:
    by_id: dict[int, pd.DataFrame] = {}
    if synonyms is not None and len(synonyms):
        by_id = {int(k): g for k, g in synonyms.groupby("mdd_id")}
    return {
        row["species"]: name_candidates(row, by_id.get(int(row["mdd_id"])) if "mdd_id" in row else None)
        for _, row in species.iterrows()
    }


def synonym_lookup(candidates: dict[str, list[tuple[str, str]]]) -> dict[str, set[str]]:
    lookup: dict[str, set[str]] = defaultdict(set)
    for sp, cands in candidates.items():
        for name, _ in cands:
            lookup[name].add(sp)
    return lookup


@dataclass
class Assignment:
    species: str
    status: str = "unmatched"
    tree_tip: str | None = None
    tree_binomial: str | None = None
    match_method: str | None = None
    confidence: float | None = None
    matched_on: str | None = None
    candidate_tip: str | None = None
    candidate_score: float | None = None
    note: str | None = None


@dataclass
class ReconcileResult:
    assignments: dict[str, Assignment]
    tips: pd.DataFrame
    candidates: dict[str, list[tuple[str, str]]] = field(default_factory=dict)

    def claimed_tips(self) -> set[str]:
        return {a.tree_tip for a in self.assignments.values() if a.status == "matched" and a.tree_tip}

    def unclaimed_tips(self) -> pd.DataFrame:
        return self.tips[~self.tips["tip_label"].isin(self.claimed_tips())]

    def table(self) -> pd.DataFrame:
        return pd.DataFrame([vars(a) for a in self.assignments.values()])

    def assign(self, species: str, tip_row: pd.Series, method: str, matched_on: str, confidence: float | None = None) -> None:
        a = self.assignments[species]
        a.status = "matched"
        a.tree_tip = tip_row["tip_label"]
        a.tree_binomial = tip_row["binomial"]
        a.match_method = method
        a.matched_on = matched_on
        a.confidence = METHOD_CONFIDENCE.get(method, 0.0) if confidence is None else confidence
        a.note = None


def _pending(result: ReconcileResult) -> list[str]:
    return [s for s, a in result.assignments.items() if a.status != "matched"]


def _resolve_claims(result: ReconcileResult, claims: dict[str, list[tuple[str, str]]], method: str, tips_by_label: dict[str, pd.Series]) -> None:
    by_species: dict[str, set[str]] = defaultdict(set)
    for tip, claimants in claims.items():
        for sp, _ in claimants:
            by_species[sp].add(tip)
    for tip, claimants in claims.items():
        species = {sp for sp, _ in claimants}
        if len(species) == 1:
            sp = next(iter(species))
            if len(by_species[sp]) == 1:
                matched_on = claimants[0][1]
                result.assign(sp, tips_by_label[tip], method, matched_on)
            else:
                a = result.assignments[sp]
                a.status = "review"
                a.note = f"{method}: multiple tips {sorted(by_species[sp])}"
        else:
            for sp in species:
                a = result.assignments[sp]
                if a.status != "matched":
                    a.status = "review"
                    a.candidate_tip = tip
                    a.note = f"{method}: tip {tip} also claimed by {sorted(species - {sp})}"


def _name_pass(result: ReconcileResult, sources: set[str], method: str) -> None:
    free = result.unclaimed_tips()
    by_binomial = {r["binomial"]: r for _, r in free.iterrows()}
    tips_by_label = {r["tip_label"]: r for _, r in free.iterrows()}
    claims: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for sp in _pending(result):
        for name, src in result.candidates.get(sp, []):
            if src in sources and name in by_binomial:
                claims[by_binomial[name]["tip_label"]].append((sp, name))
    _resolve_claims(result, claims, method, tips_by_label)


def _epithet_pass(result: ReconcileResult, species: pd.DataFrame) -> None:
    free = result.unclaimed_tips()
    tips_by_label = {r["tip_label"]: r for _, r in free.iterrows()}
    info = species.set_index("species")
    claims: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for sp in _pending(result):
        row = info.loc[sp]
        names = [n for n, src in result.candidates.get(sp, []) if src in {"exact", "legacy_name"}]
        stems = {epithet_stem(split_binomial(n)[1]) for n in names}
        genera = {split_binomial(n)[0] for n in names}
        family = str(row.get("family", "")).upper()
        hits = free[
            free["epithet"].map(epithet_stem).isin(stems)
            & ((free["family"].str.upper() == family) | free["genus"].isin(genera))
        ]
        for _, tip in hits.iterrows():
            claims[tip["tip_label"]].append((sp, tip["binomial"]))
    _resolve_claims(result, claims, "epithet_in_family", tips_by_label)


def _fuzzy_pass(result: ReconcileResult, species: pd.DataFrame, review_threshold: float, accept_threshold: float, min_margin: float) -> None:
    info = species.set_index("species")
    proposals: dict[str, tuple[str, float, str]] = {}
    for sp in _pending(result):
        free = result.unclaimed_tips()
        family = str(info.loc[sp].get("family", "")).upper()
        pool = free[free["family"].str.upper() == family]
        if pool.empty:
            pool = free
        if pool.empty:
            continue
        scores: dict[str, tuple[float, str]] = {}
        for name, _ in result.candidates.get(sp, []):
            for _, tip in pool.iterrows():
                s = fuzz.ratio(name.lower(), tip["binomial"].lower())
                if s > scores.get(tip["tip_label"], (-1, ""))[0]:
                    scores[tip["tip_label"]] = (s, name)
        ranked = sorted(scores.items(), key=lambda kv: kv[1][0], reverse=True)
        best_tip, (best, best_name) = ranked[0]
        second = ranked[1][1][0] if len(ranked) > 1 else 0.0
        a = result.assignments[sp]
        a.candidate_tip, a.candidate_score = best_tip, round(best, 1)
        if best >= accept_threshold and best - second >= min_margin:
            proposals[sp] = (best_tip, best, best_name)
        elif best >= review_threshold:
            a.status = "review"
            a.note = f"fuzzy: best {best:.1f} via {best_name!r}, runner-up {second:.1f}"
    by_tip: dict[str, list[str]] = defaultdict(list)
    for sp, (tip, _, _) in proposals.items():
        by_tip[tip].append(sp)
    tips_by_label = {r["tip_label"]: r for _, r in result.tips.iterrows()}
    for tip, sps in by_tip.items():
        if len(sps) == 1:
            sp = sps[0]
            _, score, name = proposals[sp]
            result.assign(sp, tips_by_label[tip], "fuzzy", name, confidence=round(score / 100, 3))
        else:
            for sp in sps:
                result.assignments[sp].status = "review"
                result.assignments[sp].note = f"fuzzy: tip {tip} proposed for {sps}"


def reconcile_tree(
    species: pd.DataFrame,
    tips: pd.DataFrame,
    synonyms: pd.DataFrame | None = None,
    review_threshold: float = 85,
    accept_threshold: float = 92,
    min_margin: float = 3,
) -> ReconcileResult:
    candidates = build_candidate_index(species, synonyms)
    result = ReconcileResult(
        assignments={sp: Assignment(sp) for sp in species["species"]},
        tips=tips.reset_index(drop=True),
        candidates=candidates,
    )
    _name_pass(result, {"exact"}, "exact")
    _name_pass(result, {"legacy_name"}, "legacy_name")
    _name_pass(result, {"synonym"}, "synonym")
    _epithet_pass(result, species)
    _fuzzy_pass(result, species, review_threshold, accept_threshold, min_margin)
    for a in result.assignments.values():
        if a.status == "unmatched":
            a.status = "no_tip"
            a.note = a.note or "no tree tip within fuzzy review threshold"
    return result


def apply_bridges(result: ReconcileResult, bridges: dict[str, str], method: str = "wiki_bridge") -> None:
    tips_by_label = {r["tip_label"]: r for _, r in result.tips.iterrows()}
    claimed = result.claimed_tips()
    for sp, tip in bridges.items():
        a = result.assignments[sp]
        if a.status == "matched" or tip in claimed:
            continue
        result.assign(sp, tips_by_label[tip], method, tips_by_label[tip]["binomial"])
        claimed.add(tip)


def unused_tip_table(result: ReconcileResult) -> pd.DataFrame:
    lookup = synonym_lookup(result.candidates)
    matched_by_species = {a.species: a.tree_tip for a in result.assignments.values() if a.status == "matched"}
    rows = []
    for _, tip in result.unclaimed_tips().iterrows():
        owners = sorted(lookup.get(tip["binomial"], set()))
        if owners:
            reason = "junior_synonym_of_mdd_species"
            detail = "; ".join(f"{o} -> {matched_by_species.get(o) or 'unmatched'}" for o in owners)
        else:
            reason = "no_mdd_counterpart"
            detail = "not a current MDD species or listed synonym in the target clade (often a fossil or tree-only taxon)"
        rows.append({"tree_tip": tip["tip_label"], "tree_binomial": tip["binomial"], "family": tip["family"], "reason": reason, "detail": detail})
    return pd.DataFrame(rows, columns=["tree_tip", "tree_binomial", "family", "reason", "detail"])
