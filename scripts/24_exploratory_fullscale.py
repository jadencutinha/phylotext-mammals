from __future__ import annotations

import numpy as np
import pandas as pd

from lineage_charisma import condensed as cd
from lineage_charisma.config import base_parser, load_config
from lineage_charisma.distances import load_matrix, phylo_matrix_path, text_matrix_path
from lineage_charisma.embed import read_species_order
from lineage_charisma.h2 import PRIMARY_CONTROLS, share_of_raw
from lineage_charisma.io_utils import read_csv, write_csv
from lineage_charisma.plotstyle import combo_label

# everything here is exploratory: docs/preregistration_amendment_fullscale.md, section 9
PERMS = 999
MIN_ORDER_SIZE = 50
CAP_MY = 30.0
TRANSFORMS = {"none": lambda d: d, "log": np.log1p, "capped": lambda d: np.minimum(d, CAP_MY)}
# as in the pilot: a log transform keeps the order of the distances, so it is run with Pearson only
CHECKS = [("none", "spearman"), ("capped", "spearman"), ("none", "pearson"), ("log", "pearson"), ("capped", "pearson")]
QUARTILE_NAMES = ["Q1 (least viewed)", "Q2", "Q3", "Q4 (most viewed)"]


def quartile_rows(group, text, phylo, views_matrix, seed):
    views, tv, pv = cd.condensed(views_matrix), cd.condensed(text), cd.condensed(phylo)
    quartile = pd.qcut(views, 4, labels=False)
    rows = []
    for q, name in enumerate(QUARTILE_NAMES):
        inside = quartile == q
        row = {"status": "exploratory", "group": group, "n_species": text.shape[0], "quartile": name, "n_pairs": int(inside.sum()),
               "min_mean_log_pageviews": views[inside].min(), "max_mean_log_pageviews": views[inside].max(),
               "mean_text_distance": tv[inside].mean(), "median_text_distance": np.median(tv[inside]), "mean_phylo_distance": pv[inside].mean(),
               "r_text_phylo": np.nan, "p": np.nan}
        if q in (0, len(QUARTILE_NAMES) - 1):
            res = cd.mantel_on_pairs(text, phylo, inside, "spearman", PERMS, seed)
            row.update(r_text_phylo=res.r, p=res.p)
        rows.append(row)
        print(f"[exploratory] {group:14s} {name:18s} {row['n_pairs']:9d} pairs  mean text {row['mean_text_distance']:.4f}  mean phylo {row['mean_phylo_distance']:.1f} My"
              + ("" if np.isnan(row["r_text_phylo"]) else f"  r = {row['r_text_phylo']:+.4f}  p = {row['p']:.4f}"), flush=True)
    return rows


def main() -> None:
    args = base_parser("Full scale, exploratory: pageview quartiles (all species and within orders) and the 30 My cap.").parse_args()
    cfg = load_config(args.config, args.clade or "mammalia")
    seed = int(cfg["stats"]["seed"])
    out_dir, order_path = cfg.matrices_dir(), cfg.species_order_path()
    primary = cfg.primary_combo
    order = read_species_order(order_path)
    text = load_matrix(text_matrix_path(out_dir, *primary), order_path)[0]
    phylo = load_matrix(phylo_matrix_path(out_dir), order_path)[0]
    views = load_matrix(out_dir / "attention_mean_log_pageviews.npy", order_path)[0]
    mdd_order = read_csv(cfg.mdd_species_path()).set_index("species")["order"]
    tables = cfg.tables_dir()

    # 1. pageview quartiles, across all species and within each order with at least 50 species
    rows = quartile_rows("all mammals", text, phylo, views, seed)
    counts = mdd_order.loc[order].value_counts()
    for name in counts[counts >= MIN_ORDER_SIZE].index:
        members = [s for s in order if mdd_order[s] == name]
        rows += quartile_rows(name, *(cd.subset(m, order, members) for m in (text, phylo, views)), seed)
    write_csv(pd.DataFrame(rows), tables / "exploratory_pageview_quartiles.csv", status="exploratory", permutations=PERMS, seed=seed, alternative="greater",
              permuted="phylogeny, with the set of pairs held fixed", text_matrix=combo_label(*primary))
    del views

    # 2. transformed phylogenetic distance, including the 30 My cap
    species = read_species_order(out_dir / "ecology_species.txt")
    controls = {"ecology": load_matrix(out_dir / "ecology.npy", out_dir / "ecology_species.txt")[0],
                "length": cd.subset(load_matrix(out_dir / "attention_length.npy", order_path)[0], order, species)}
    text_h2 = cd.subset(text, order, species)
    rows = []
    for transform, method in CHECKS:
        moved = TRANSFORMS[transform](phylo)
        h1 = cd.mantel(cd.prepare(text, method), cd.prepare(moved, method), PERMS, seed)
        x, y = cd.prepare(cd.subset(moved, order, species), method), cd.prepare(text_h2, method)
        z = [cd.prepare(controls[c], method) for c in PRIMARY_CONTROLS]
        null = cd.null_correlations(x, [y, *z], PERMS, seed)
        raw = cd.mantel(y, x, seed=seed, null=null[:, 0])
        part = cd.partial_mantel(x, y, z, seed=seed, null_corr=null)
        rows.append({"status": "exploratory", "transform": transform, "method": method, "h1_r": h1.r, "h1_p": h1.p, "h1_n_species": h1.n,
                     "raw_r": raw.r, "raw_p": raw.p, "partial_r": part.r, "partial_p": part.p, "share_of_raw": share_of_raw(raw.r, part.r),
                     "change": raw.r - part.r, "h2_n_species": part.n})
        print(f"[exploratory] {transform:7s} {method:8s} H1 r = {h1.r:+.4f} (p = {h1.p:.4f})  partial r = {part.r:+.4f} (p = {part.p:.4f})", flush=True)
    pv = cd.condensed(phylo)
    write_csv(pd.DataFrame(rows), tables / "exploratory_nonlinear.csv", status="exploratory", permutations=PERMS, seed=seed, alternative="greater",
              cap_my=CAP_MY, share_of_pairs_at_or_below_cap=float((pv <= CAP_MY).mean()), controls=list(PRIMARY_CONTROLS), text_matrix=combo_label(*primary))
    print(f"[exploratory] {float((pv <= CAP_MY).mean()):.2%} of pairs are at or below {CAP_MY:.0f} My")


if __name__ == "__main__":
    main()
