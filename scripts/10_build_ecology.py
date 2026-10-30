from __future__ import annotations

import numpy as np
import pandas as pd

from lineage_charisma.config import base_parser, load_config
from lineage_charisma.distances import save_matrix
from lineage_charisma.embed import read_species_order, write_species_order
from lineage_charisma.io_utils import atomic_write_text, download, is_cached, read_csv, write_csv
from lineage_charisma.stats import submatrix
from lineage_charisma.traits import (
    TRAITS,
    Trait,
    build_trait_table,
    complete_species,
    gower_distance,
    load_combine,
    load_elton,
    md5,
    missing_rates,
    reconcile_traits,
    trait_list,
    trait_weights,
)

DATABASES = {"elton": ("elton_traits", load_elton, "EltonTraits 1.0"), "combine": ("combine", load_combine, "COMBINE")}
MATRICES = {"group": "ecology.npy", "column": "ecology_unweighted.npy"}
VOLANT_COLUMN = "volant"
NO_VOLANT_MATRIX = "ecology_no_volant.npy"   # group-weighted without `volant`; written only when traits.include_volant is set


def fetch(cfg, key: str, force: bool):
    src = cfg["sources"][key]
    path = cfg.raw / "traits" / src["filename"]
    if force or not is_cached(path):
        print(f"[download] {src['url']}")
        download(src["url"], path, user_agent=cfg["wikipedia"]["user_agent"])
    else:
        print(f"[cache] {path} exists")
    found = md5(path)
    if found != src["md5"]:
        raise SystemExit(f"{path.name}: md5 {found} does not match the recorded {src['md5']}. The file at {src['url']} is not the one "
                         "this project was built against; look at it before going on.")
    return path


def traits_doc(table: pd.DataFrame, rates: pd.DataFrame, shared: dict[str, np.ndarray], n_complete: int, traits: tuple[Trait, ...] = TRAITS) -> str:
    weights = {w: trait_weights(traits, w) / trait_weights(traits, w).sum() for w in MATRICES}
    desc = {t.column: t.description for t in traits}
    out = ["# Ecological traits", "",
           "Written by `scripts/10_build_ecology.py`; do not edit by hand. The trait list is fixed in `docs/preregistration.md` (section 7 and amendment A1.1) "
           "and in `TRAITS` in `src/lineage_charisma/traits.py`.", "",
           "## Sources", "",
           "- **EltonTraits 1.0** (Wilman et al. 2014, *Ecology* 95:2027), mammal file `MamFuncDat.txt`: diet, foraging stratum, activity time.",
           "- **COMBINE** (Soria et al. 2021, *Ecology* 102:e03344), `trait_data_reported.csv`: habitat breadth, terrestrial or aquatic, body mass. "
           "This is the file of reported values. COMBINE's imputed file is not used.", "",
           "## Traits", "",
           f"{len(traits)} columns in {len(set(t.group for t in traits))} groups, for {len(table)} species. "
           "\"Weight\" is each column's share of the distance when a pair has every trait: equal by group in the primary matrix, equal by column in the robustness variant.", "",
           "| column | group | type | source | source column | meaning | group weight | column weight | missing | missing rate |",
           "|---|---|---|---|---|---|---:|---:|---:|---:|"]
    for i, r in rates.iterrows():
        out.append(f"| `{r['trait']}` | {r['group']} | {r['type']} | {r['source']} | `{r['source_column']}` | {desc[r['trait']]} | "
                   f"{weights['group'][i]:.3f} | {weights['column'][i]:.3f} | {r['n_missing']} | {r['missing_rate']:.1%} |")
    out += ["", "## How the distance is computed", "",
            "- Continuous traits contribute |difference| / range, where the range is taken over the species with a value. "
            "Categorical and binary traits contribute 0 for a match and 1 for a mismatch.",
            "- A pair's distance is the weighted mean of those contributions over the traits both species have (Gower's pairwise-available rule). Nothing is imputed.",
            "- `ecology.npy` weights the six groups equally (primary). `ecology_unweighted.npy` weights every column equally (robustness variant).",
            "- Body mass is log10 of `adult_mass_g`. The three activity flags are not exclusive: a species can be both nocturnal and crepuscular.",
            ("- `volant` is COMBINE's `terrestrial_volant` (full scale only, docs/preregistration_amendment_fullscale.md section 4). "
             f"`{NO_VOLANT_MATRIX}` is the group-weighted matrix without it, a declared robustness variant."
             if VOLANT_COLUMN in desc else "- `terrestrial_volant` in COMBINE is 0 for every Carnivora species that has a value and is not used."), "",
            "## Coverage", "",
            f"{n_complete} of {len(table)} species have an ecological distance to each other and are in the matrices (`ecology_species.txt`). "
            f"{int(table.notna().all(axis=1).sum())} species have every trait."]
    for w, name in MATRICES.items():
        frac = shared[w][np.triu_indices(len(table), k=1)]
        out.append(f"- `{name}`: over all {len(frac):,} pairs of the {len(table)} species, a pair shares on average {frac.mean():.1%} of the total trait weight; "
                   f"{(frac < 1 - 1e-9).mean():.1%} of pairs are missing at least one trait and {int((frac == 0).sum())} pairs share none.")
    return "\n".join(out) + "\n"


def main() -> None:
    parser = base_parser("Download EltonTraits and COMBINE, reconcile names, and build the ecological (Gower) distance matrices.")
    parser.add_argument("--accept-residual", action="store_true",
                        help="go on even if more than traits.max_unmatched species lack a record in a trait database (after reviewing traits_residual_review.csv)")
    args = parser.parse_args()
    cfg = load_config(args.config, args.clade)
    order = read_species_order(cfg.species_order_path())
    traits = trait_list(bool(cfg["traits"].get("include_volant", False)))
    species = read_csv(cfg.mdd_species_path())
    species = species[species["species"].isin(order)].reset_index(drop=True)
    synonyms = read_csv(cfg.mdd_synonyms_path())

    tables, maps = {}, {}
    for name, (key, loader, _) in DATABASES.items():
        tables[name] = loader(fetch(cfg, key, args.force))
        maps[name] = reconcile_traits(species, synonyms, tables[name], cfg["matching"])
        counts = maps[name]["match_method"].value_counts().to_dict()
        print(f"[traits] {name}: {int((maps[name]['status'] == 'matched').sum())}/{len(order)} species matched {counts}")

    name_map = pd.concat([m.assign(database=name) for name, m in maps.items()], ignore_index=True)
    write_csv(name_map[["species", "database", *[c for c in name_map.columns if c not in ("species", "database")]]], cfg.interim / "traits_name_map.csv")
    info = species.set_index("species")
    residual = name_map[name_map["status"] != "matched"].copy()
    residual["family"] = info.loc[residual["species"], "family"].to_numpy()
    residual["common_name"] = info.loc[residual["species"], "common_name"].to_numpy()
    residual["issue"] = residual["status"].map({"review": "ambiguous_match", "no_tip": "no_record"})
    residual = residual[["species", "database", "issue", "candidate_name", "candidate_score", "note", "family", "common_name"]].sort_values(["database", "species"])
    residual_path = cfg.interim / "traits_residual_review.csv"
    review = name_map[(name_map["status"] == "matched") & ~name_map["match_method"].isin(["exact", "legacy_name"])]
    if len(review):
        print("[traits] matches made by synonym, epithet or fuzzy rule; check them in traits_name_map.csv:")
        print(review[["species", "database", "source_name", "match_method", "confidence"]].to_string(index=False))

    reviewed, limit = set(cfg["traits"].get("reviewed_missing") or []), int(cfg["traits"]["max_unmatched"])
    residual["reviewed"] = residual["species"].isin(reviewed)
    write_csv(residual, residual_path)
    unreviewed = sorted(set(residual["species"]) - reviewed)
    if len(residual):
        print(f"[traits] {residual['species'].nunique()} species lack a record in at least one database, {len(unreviewed)} of them not yet reviewed "
              f"({residual_path.relative_to(cfg.root)}):")
        print(residual.to_string(index=False))
    if len(unreviewed) > limit and not args.accept_residual:
        raise SystemExit(f"[traits] stopping: {len(unreviewed)} unreviewed unmatched species is more than traits.max_unmatched = {limit}: {unreviewed}. "
                         "Review the list, then add them to traits.reviewed_missing in config.yaml or rerun with --accept-residual (make ecology ACCEPT_RESIDUAL=1).")

    table = build_trait_table(order, tables["elton"], tables["combine"], maps["elton"], maps["combine"], traits)
    rates = missing_rates(table, traits)
    write_csv(table.reset_index(), cfg.processed / "traits.csv")
    write_csv(rates, cfg.tables_dir() / "ecology_missing.csv", n_species=len(order))

    dist, shared = {}, {}
    for weighting in MATRICES:
        dist[weighting], shared[weighting] = gower_distance(table, traits, weighting)
    # which pairs share a trait does not depend on the weights, so one species set serves both matrices
    complete = complete_species(dist["group"], order)
    dropped = sorted(set(order) - set(complete))
    order_path = cfg.matrices_dir() / "ecology_species.txt"
    write_species_order(order_path, complete)
    if dropped:
        print(f"[ecology] {len(dropped)} species have no ecological distance to some other species and are left out: {dropped}")
    for weighting, filename in MATRICES.items():
        matrix = submatrix(dist[weighting], order, complete)
        save_matrix(cfg.matrices_dir() / filename, matrix, complete, kind="gower", weighting=weighting, traits=[t.column for t in traits], dropped_species=dropped)
        off = matrix[np.triu_indices(len(complete), k=1)]
        print(f"[ecology] {filename}: {len(complete)} species, mean distance {off.mean():.3f}, max {off.max():.3f}")
    if traits != TRAITS:
        variant, _ = gower_distance(table, TRAITS, "group")
        matrix = submatrix(variant, order, complete)
        off = matrix[np.triu_indices(len(complete), k=1)]
        if np.isnan(off).any():
            raise SystemExit(f"[ecology] {NO_VOLANT_MATRIX}: {int(np.isnan(off).sum())} pairs of ecology_species.txt share no trait once `volant` is left out")
        save_matrix(cfg.matrices_dir() / NO_VOLANT_MATRIX, matrix, complete, kind="gower", weighting="group", traits=[t.column for t in TRAITS],
                    dropped_species=dropped, role="robustness variant: ecology without volant")
        print(f"[ecology] {NO_VOLANT_MATRIX}: {len(complete)} species, mean distance {off.mean():.3f}, max {off.max():.3f}")

    doc_name = f"ecology_traits_{cfg.output_subdir}.md" if cfg.output_subdir else "ecology_traits.md"
    atomic_write_text(cfg.root / "docs" / doc_name, traits_doc(table, rates, shared, len(complete), traits))
    print("[ecology] missing data per trait:")
    print(rates[["trait", "n_missing", "missing_rate"]].to_string(index=False, float_format=lambda v: f"{v:.1%}"))


if __name__ == "__main__":
    main()
