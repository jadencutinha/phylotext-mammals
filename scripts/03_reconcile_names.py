from __future__ import annotations

import pandas as pd

from lineage_charisma.config import base_parser, load_config
from lineage_charisma.io_utils import read_csv, skip_if_cached, write_csv
from lineage_charisma.taxonomy import apply_bridges, reconcile_tree, unused_tip_table
from lineage_charisma.wiki import WikiClient, resolve_species_titles

MATCHED_COLUMNS = [
    "species", "tree_tip", "wiki_title", "match_method", "confidence",
    "wiki_validation", "wiki_source", "wiki_query", "wiki_qid", "matched_on", "tree_binomial",
    "mdd_id", "family", "genus", "common_name", "extinct", "domestic",
]
RESIDUAL_COLUMNS = [
    "species", "issue", "candidate_tip", "candidate_score", "wiki_title", "note", "mdd_id", "family", "common_name",
]
EXCLUDED_COLUMNS = [
    "species", "reason", "nearest_tip", "nearest_score", "wiki_title", "note", "mdd_id", "family", "common_name", "extinct", "domestic",
]


def title_candidates(species: pd.DataFrame, result) -> tuple[dict, dict, dict]:
    candidates, other_names, search_names = {}, {}, {}
    for _, row in species.iterrows():
        sp = row["species"]
        a = result.assignments[sp]
        names = result.candidates.get(sp, [])
        items: list[tuple[str, str]] = [(sp, "scientific_name")]
        if a.tree_binomial:
            items.append((a.tree_binomial, "tree_name"))
        items += [(n, "legacy_name") for n, s in names if s == "legacy_name"]
        if isinstance(row.get("common_name"), str):
            items.append((row["common_name"], "common_name"))
        if isinstance(row.get("other_common_names"), str):
            items += [(c.strip(), "other_common_name") for c in row["other_common_names"].split("|") if c.strip()]
        items += [(n, "synonym") for n, s in names if s == "synonym"]
        candidates[sp] = list(dict.fromkeys(items))
        other_names[sp] = {n for n, _ in names} | ({a.tree_binomial} if a.tree_binomial else set())
        search_names[sp] = list(dict.fromkeys([sp] + ([a.tree_binomial] if a.tree_binomial else [])))
    return candidates, other_names, search_names


def find_bridges(result, species: pd.DataFrame, wiki: dict, client: WikiClient) -> tuple[dict[str, str], dict[str, str]]:
    family = species.set_index("species")["family"].str.upper().to_dict()
    leftovers = [s for s, a in result.assignments.items() if a.status != "matched" and s in wiki]
    free = result.unclaimed_tips()
    if not leftovers or free.empty:
        return {}, {}
    infos = client.query_titles(free["binomial"])
    tip_article = {
        r["tip_label"]: infos[r["binomial"]].title
        for _, r in free.iterrows()
        if infos.get(r["binomial"]) and infos[r["binomial"]].exists and not infos[r["binomial"]].disambiguation
    }
    tip_family = free.set_index("tip_label")["family"].to_dict()
    bridges, rejections = {}, {}
    for sp in leftovers:
        title = wiki[sp]["wiki_title"]
        hits = [t for t, art in tip_article.items() if art == title and tip_family.get(t) == family.get(sp)]
        if len(hits) == 1:
            bridges[sp] = hits[0]
        cand = result.assignments[sp].candidate_tip
        if not hits and cand in tip_article and tip_article[cand] != title:
            rejections[sp] = f"candidate tip {cand} resolves to a different article ({tip_article[cand]!r})"
    return bridges, rejections


def main() -> None:
    parser = base_parser("Reconcile MDD names with tree tips and English Wikipedia titles.")
    parser.add_argument("--offline", action="store_true", help="use only cached API responses")
    args = parser.parse_args()
    cfg = load_config(args.config)
    outputs = [cfg.matched_path(), cfg.residual_path(), cfg.excluded_species_path(), cfg.unused_tips_path()]
    if skip_if_cached(outputs, args.force):
        return

    species = read_csv(cfg.mdd_species_path())
    synonyms = read_csv(cfg.mdd_synonyms_path())
    tips = read_csv(cfg.tips_path())
    m = cfg["matching"]
    result = reconcile_tree(
        species, tips, synonyms,
        review_threshold=m["fuzzy_review_threshold"],
        accept_threshold=m["fuzzy_accept_threshold"],
        min_margin=m["fuzzy_min_margin"],
    )
    counts = result.table()["status"].value_counts().to_dict()
    print(f"[tree-match] {counts}")

    client = WikiClient.from_config(cfg, offline=args.offline)
    candidates, other_names, search_names = title_candidates(species, result)
    wiki, wiki_notes = resolve_species_titles(candidates, other_names, client, search_names)
    print(f"[wiki] resolved {len(wiki)}/{len(species)} titles ({client.network_calls} network calls)")

    bridges, rejections = find_bridges(result, species, wiki, client)
    apply_bridges(result, bridges)
    for sp, tip in bridges.items():
        print(f"[bridge] {sp} -> {tip} (shared article {wiki[sp]['wiki_title']!r})")
    for sp, why in rejections.items():
        a = result.assignments[sp]
        if a.status == "review":
            a.status, a.note = "no_tip", why
            print(f"[auto-reject] {sp}: {why}")

    info = species.set_index("species")
    matched, residual, excluded = [], [], []
    for sp, a in result.assignments.items():
        row = info.loc[sp]
        w = wiki.get(sp, {})
        base = dict(species=sp, mdd_id=row["mdd_id"], family=row["family"], common_name=row.get("common_name"))
        if a.status == "matched" and w:
            matched.append({**base, **{k: v for k, v in w.items() if k != "wiki_rank"}, "tree_tip": a.tree_tip, "tree_binomial": a.tree_binomial,
                            "match_method": a.match_method, "confidence": a.confidence, "matched_on": a.matched_on,
                            "genus": row["genus"], "extinct": row["extinct"], "domestic": row["domestic"]})
        elif a.status == "matched":
            residual.append({**base, "issue": "no_validated_wikipedia_article", "candidate_tip": a.tree_tip, "note": wiki_notes.get(sp)})
        elif a.status == "review":
            residual.append({**base, "issue": "ambiguous_tree_match", "candidate_tip": a.candidate_tip, "candidate_score": a.candidate_score,
                             "wiki_title": w.get("wiki_title"), "note": a.note})
        else:
            reason = "extinct_no_tree_tip" if row["extinct"] else ("domestic_no_tree_tip" if row["domestic"] else "no_tree_tip")
            excluded.append({**base, "reason": reason, "nearest_tip": a.candidate_tip, "nearest_score": a.candidate_score,
                             "wiki_title": w.get("wiki_title"), "note": a.note, "extinct": row["extinct"], "domestic": row["domestic"]})

    matched_df = pd.DataFrame(matched, columns=MATCHED_COLUMNS).sort_values(["family", "species"])
    residual_df = pd.DataFrame(residual, columns=RESIDUAL_COLUMNS).sort_values("species")
    excluded_df = pd.DataFrame(excluded, columns=EXCLUDED_COLUMNS).sort_values(["reason", "species"])
    unused_df = unused_tip_table(result).sort_values(["reason", "tree_tip"])

    meta = dict(clade=f"{cfg.clade_rank}={cfg.clade_name}", matching=m)
    write_csv(matched_df, cfg.matched_path(), **meta)
    write_csv(residual_df, cfg.residual_path(), **meta)
    write_csv(excluded_df, cfg.excluded_species_path(), **meta)
    write_csv(unused_df, cfg.unused_tips_path(), **meta)

    print(f"[done] matched={len(matched_df)} residual={len(residual_df)} excluded_species={len(excluded_df)} unused_tips={len(unused_df)}")
    print(f"[done] match methods: {matched_df['match_method'].value_counts().to_dict()}")
    print(f"[done] wiki validation: {matched_df['wiki_validation'].value_counts().to_dict()}")
    if len(residual_df):
        print("[review] residual items:")
        print(residual_df[["species", "issue", "candidate_tip", "note"]].to_string(index=False))


if __name__ == "__main__":
    main()
