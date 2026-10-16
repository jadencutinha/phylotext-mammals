from __future__ import annotations

import numpy as np
import pandas as pd

from lineage_charisma.config import base_parser, load_config
from lineage_charisma.distances import load_matrix, phylo_matrix_path, text_matrix_path
from lineage_charisma.h2 import PRIMARY_CONTROLS, load_h2_matrices, partial_result
from lineage_charisma.io_utils import write_csv
from lineage_charisma.plotstyle import BLUE, INK, MUTED, ORANGE, SURFACE, apply_style, combo_label, plt, save
from lineage_charisma.stats import mantel, mantel_on_pairs, submatrix, upper_triangle

# everything here is exploratory: docs/preregistration_amendment_h2.md, section 9
CAP_MY = 30.0
TRANSFORMS = {
    "none": lambda d: d,
    "log": np.log1p,
    "capped": lambda d: np.minimum(d, CAP_MY),
}
TRANSFORM_NAMES = {"none": "untransformed", "log": "log(1 + distance)", "capped": f"capped at {CAP_MY:.0f} My"}
# a log transform keeps the order of the distances, so its Spearman result is the untransformed one; it is run with Pearson only
CHECKS = [("none", "spearman"), ("capped", "spearman"), ("none", "pearson"), ("log", "pearson"), ("capped", "pearson")]
QUARTILE_NAMES = ["Q1 (least viewed)", "Q2", "Q3", "Q4 (most viewed)"]


def pair_matrix(values: np.ndarray, n: int) -> np.ndarray:
    """Symmetric n x n matrix whose upper triangle is `values`."""
    out = np.zeros((n, n), dtype=values.dtype)
    out[np.triu_indices(n, k=1)] = values
    return out + out.T


def plot_quartiles(table: pd.DataFrame, label: str, path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.6), constrained_layout=True, sharey=True)
    ys = np.arange(len(table))[::-1]
    for ax, column, colour, xlabel in ((axes[0], "mean_text_distance", BLUE, "mean text distance"), (axes[1], "mean_phylo_distance", MUTED, "mean phylogenetic distance (My)")):
        values = table[column].to_numpy()
        ax.barh(ys, values, height=0.56, color=colour, edgecolor=SURFACE, linewidth=1)
        for y, v in zip(ys, values):
            ax.annotate(f"{v:.3f}" if column == "mean_text_distance" else f"{v:.1f}", xy=(v, y), xytext=(4, 0), textcoords="offset points", va="center", color=INK)
        ax.set_xlim(0, values.max() * 1.18)
        ax.set_xlabel(xlabel)
        ax.grid(True, axis="x")
        ax.set_axisbelow(True)
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)
    axes[0].set_yticks(ys, table["quartile"])
    fig.suptitle(f"Exploratory: species pairs by quartile of mean log pageviews ({label})", x=0.01, ha="left", fontweight="bold", fontsize=10)
    save(fig, path)


def plot_nonlinear(table: pd.DataFrame, path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 3.8), constrained_layout=True)
    ys = np.arange(len(table))[::-1]
    ax.hlines(ys, table[["h1_r", "partial_r"]].min(axis=1), table[["h1_r", "partial_r"]].max(axis=1), color="#c3c2b7", linewidth=2, zorder=1)
    ax.scatter(table["h1_r"], ys, s=64, color=BLUE, edgecolor=SURFACE, linewidth=1.5, zorder=3, label="H1: raw r, 284 species")
    ax.scatter(table["partial_r"], ys, s=64, color=ORANGE, edgecolor=SURFACE, linewidth=1.5, zorder=3, label="primary H2 model: partial r, 283 species")
    ax.axvline(0, color=MUTED, linewidth=0.8)
    ax.set_yticks(ys, [f"{TRANSFORM_NAMES[t]}, {m.capitalize()}" for t, m in zip(table["transform"], table["method"])])
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(-0.6, len(table) - 0.4)
    ax.set_xlabel("Mantel r, text distance against (transformed) phylogenetic distance")
    ax.set_title("Exploratory: transforming phylogenetic distance")
    ax.grid(True, axis="x")
    ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False)
    fig.legend(loc="outside lower center", ncol=2, handletextpad=0.3, columnspacing=1.6)
    save(fig, path)


def main() -> None:
    args = base_parser("Exploratory analyses: pairs of well-known species, and transformed phylogenetic distance.").parse_args()
    cfg = load_config(args.config)
    apply_style()
    scfg = cfg["stats"]
    perms, seed = int(scfg["permutations"]), int(scfg["seed"])
    primary, out_dir, order_path = cfg.primary_combo, cfg.matrices_dir(), cfg.species_order_path()
    label = combo_label(*primary)
    text, order = load_matrix(text_matrix_path(out_dir, *primary), order_path)
    phylo, _ = load_matrix(phylo_matrix_path(out_dir), order_path)
    n = len(order)

    # 1. "both well-known": pairs split by the mean of the two species' log pageviews
    mean_views, _ = load_matrix(out_dir / "attention_mean_log_pageviews.npy", order_path)
    views, tv, pv = upper_triangle(mean_views), upper_triangle(text), upper_triangle(phylo)
    quartile = pd.qcut(views, 4, labels=False)
    rows = []
    for q, name in enumerate(QUARTILE_NAMES):
        inside = quartile == q
        row = {"status": "exploratory", "quartile": name, "n_pairs": int(inside.sum()), "min_mean_log_pageviews": views[inside].min(),
               "max_mean_log_pageviews": views[inside].max(), "mean_text_distance": tv[inside].mean(), "median_text_distance": np.median(tv[inside]),
               "mean_phylo_distance": pv[inside].mean(), "r_text_phylo": np.nan, "p": np.nan}
        if q in (0, len(QUARTILE_NAMES) - 1):
            res = mantel_on_pairs(text, phylo, pair_matrix(inside, n), method="spearman", permutations=perms, seed=seed)
            row.update(r_text_phylo=res.r, p=res.p)
        rows.append(row)
        print(f"[exploratory] {name:18s} {row['n_pairs']:6d} pairs  mean text distance {row['mean_text_distance']:.4f}  "
              f"mean phylo {row['mean_phylo_distance']:.1f} My" + ("" if np.isnan(row["r_text_phylo"]) else f"  r = {row['r_text_phylo']:+.4f}  p = {row['p']:.4f}"))
    well_known = pd.DataFrame(rows)
    write_csv(well_known, cfg.tables_dir() / "exploratory_well_known.csv", status="exploratory", permutations=perms, seed=seed, alternative="greater",
              permuted="phylogeny, with the set of pairs held fixed", text_matrix=label, n_species=n)

    # 2. nonlinear phylogeny: H1 on all species and the primary H2 model on the ecology species, with phylogenetic distance transformed
    matrices, species, _ = load_h2_matrices(cfg)
    text_h2, controls = submatrix(text, order, species), [matrices[c] for c in PRIMARY_CONTROLS]
    rows = []
    for transform, method in CHECKS:
        moved = TRANSFORMS[transform](phylo)
        h1 = mantel(text, moved, method=method, permutations=perms, seed=seed)
        res = partial_result(text_h2, submatrix(moved, order, species), controls, method, perms, seed)
        rows.append({"status": "exploratory", "transform": transform, "method": method, "h1_r": h1.r, "h1_p": h1.p, "h1_n_species": h1.n,
                     "raw_r": res["raw_r"], "raw_p": res["raw_p"], "partial_r": res["partial_r"], "partial_p": res["partial_p"],
                     "share_of_raw": res["share_of_raw"], "change": res["change"], "h2_n_species": res["n_species"]})
        print(f"[exploratory] {TRANSFORM_NAMES[transform]:20s} {method:8s} H1 r = {h1.r:+.4f} (p = {h1.p:.4f})  "
              f"partial r = {res['partial_r']:+.4f} (p = {res['partial_p']:.4f})")
    nonlinear = pd.DataFrame(rows)
    write_csv(nonlinear, cfg.tables_dir() / "exploratory_nonlinear.csv", status="exploratory", permutations=perms, seed=seed, alternative="greater",
              cap_my=CAP_MY, share_of_pairs_at_or_below_cap=float((pv <= CAP_MY).mean()), controls=list(PRIMARY_CONTROLS), text_matrix=label)

    fig_dir = cfg.week5_figures_dir()
    plot_quartiles(well_known, label, fig_dir / "exploratory_well_known.png")
    plot_nonlinear(nonlinear, fig_dir / "exploratory_nonlinear.png")
    print(f"[exploratory] wrote exploratory_well_known.csv, exploratory_nonlinear.csv and 2 figures in {fig_dir.relative_to(cfg.root)}")


if __name__ == "__main__":
    main()
