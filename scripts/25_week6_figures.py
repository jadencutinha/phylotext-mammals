from __future__ import annotations

import numpy as np
import pandas as pd

from lineage_charisma.config import base_parser, load_config
from lineage_charisma.io_utils import read_csv
from lineage_charisma.plotstyle import AXIS, BLUE, INK, MUTED, ORANGE, SURFACE, apply_style, plt, save

PREDICTOR_NAMES = {"phylogeny": "phylogenetic distance", "ecology": "ecological distance", "length": "article-length difference", "pageviews": "pageview difference"}


def p_label(p: float, permutations: int) -> str:
    return f"p ≤ {1 / (permutations + 1):.4f}" if p <= 1.0001 / (permutations + 1) else f"p = {p:.4f}"


def hbar(ax, labels, values, colours, notes, xlabel):
    ys = np.arange(len(values))[::-1]
    values = np.asarray(values, dtype=float)
    ax.barh(ys, values, height=0.56, color=colours, edgecolor=SURFACE, linewidth=1)
    span = max(values.max(), 0) - min(values.min(), 0)
    for y, v, note in zip(ys, values, notes):
        ax.annotate(note, xy=(max(v, 0), y), xytext=(5, 0), textcoords="offset points", va="center", color=INK)
    ax.axvline(0, color=MUTED, linewidth=0.8)
    ax.set_xlim(min(values.min(), 0) - 0.05 * span, max(values.max(), 0) + 0.6 * span)
    ax.set_yticks(ys, labels)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel(xlabel)
    ax.grid(True, axis="x")
    ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False)


def plot_controls(partial: pd.DataFrame, path) -> None:
    by = partial.set_index("analysis")
    head = by.loc["primary"]
    items = [("no control (raw)", head["raw_r"], head["raw_p"], head["permutations"], MUTED)]
    for key, name in (("secondary S1", "ecology"), ("secondary S2", "length"), ("secondary S3", "pageviews"),
                      ("primary", "ecology + length\n(primary H2 model)"), ("secondary S4", "ecology + length\n+ pageviews"),
                      ("ecology_without_flight", "ecology without flight\n+ length (robustness)")):
        r = by.loc[key]
        items.append((name, r["partial_r"], r["partial_p"], r["permutations"], ORANGE if key == "ecology_without_flight" else BLUE))
    fig, ax = plt.subplots(figsize=(7.2, 4.4), constrained_layout=True)
    hbar(ax, [i[0] for i in items], [i[1] for i in items], [i[4] for i in items],
         [f"r = {i[1]:.3f}, {p_label(i[2], int(i[3]))}" for i in items], "Spearman Mantel r, text distance against phylogenetic distance")
    ax.set_title(f"All mammals: text–phylogeny correlation, raw and after each control\n{int(head['n_species']):,} species, one-sided p "
                 f"({int(head['permutations']):,} permutations; robustness variant 999)")
    save(fig, path)


def plot_mrm(mrm: pd.DataFrame, path) -> None:
    rows = mrm[mrm["is_primary"]]
    perms = int(rows["permutations"].iloc[0])
    fig, ax = plt.subplots(figsize=(7.2, 3.2), constrained_layout=True)
    hbar(ax, [PREDICTOR_NAMES[p] for p in rows["predictor"]], rows["coef_std"], BLUE,
         [f"β = {c:.3f}, {p_label(p, perms)}" for c, p in zip(rows["coef_std"], rows["p"])], "standardized coefficient on ranked distances")
    ax.set_title(f"All mammals, MRM: text distance on four predictors, two-sided p\n{int(rows['n_species'].iloc[0]):,} species, "
                 f"R² = {rows['r2'].iloc[0]:.4f}, {p_label(rows['p_r2'].iloc[0], perms)}")
    save(fig, path)


def plot_posterior(post: pd.DataFrame, h1: pd.Series, h2: pd.Series, path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.6), constrained_layout=True)
    panels = ((axes[0], [("h1_r", "H1 r", BLUE, h1["r"])], "H1 (frozen set)"),
              (axes[1], [("raw_r", "raw r", BLUE, h2["raw_r"]), ("partial_r", "partial r", ORANGE, h2["partial_r"])], "Primary H2 model (H2 species)"))
    for ax, series, title in panels:
        for k, (col, name, colour, mcc) in enumerate(series):
            ys = np.full(len(post), k) + np.linspace(-0.12, 0.12, len(post))
            ax.scatter(post[col], ys, s=18, color=colour, alpha=0.8, edgecolor="none")
            ax.scatter([mcc], [k], marker="D", s=46, color=INK, zorder=3)
            ax.annotate(f"MCC {mcc:.3f}; trees {post[col].min():.3f}–{post[col].max():.3f}", xy=(mcc, k), xytext=(0, 11), textcoords="offset points",
                        ha="center", fontsize=7.5, color=INK)
        ax.set_yticks(range(len(series)), [s[1] for s in series])
        ax.set_ylim(-0.6, len(series) - 0.3)
        ax.axvline(0, color=MUTED, linewidth=0.8)
        ax.set_xlabel("Spearman Mantel r")
        ax.set_title(title)
        ax.grid(True, axis="x")
        ax.set_axisbelow(True)
        ax.tick_params(axis="y", length=0)
        ax.spines["left"].set_visible(False)
    fig.suptitle(f"All mammals: {len(post)} posterior trees (dots) and the MCC tree (diamond), 999 permutations per tree", x=0.01, ha="left",
                 fontweight="bold", fontsize=10)
    save(fig, path)


def plot_orders(orders: pd.DataFrame, h1: pd.Series, path) -> None:
    orders = orders.sort_values("h1_r")
    fig, ax = plt.subplots(figsize=(7.2, 0.42 * len(orders) + 1.6), constrained_layout=True)
    ys = np.arange(len(orders))
    ax.hlines(ys, orders[["h1_r", "partial_r"]].min(axis=1), orders[["h1_r", "partial_r"]].max(axis=1), color=AXIS, linewidth=2, zorder=1)
    ax.scatter(orders["h1_r"], ys, s=60, color=BLUE, edgecolor=SURFACE, linewidth=1.4, zorder=3, label="H1 r")
    ax.scatter(orders["partial_r"], ys, s=60, color=ORANGE, edgecolor=SURFACE, linewidth=1.4, zorder=3, label="primary H2 partial r (ecology + length)")
    ax.axvline(0, color=MUTED, linewidth=0.8)
    ax.axvline(h1["r"], color=BLUE, linewidth=0.8, linestyle="--", label=f"all mammals, H1 r = {h1['r']:.3f}")
    ax.set_yticks(ys, [f"{o} (n = {n:,})" for o, n in zip(orders["order"], orders["n_species"])])
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Spearman Mantel r, text distance against phylogenetic distance, within the order")
    ax.set_title("Within-order correlations, orders with at least 50 species (descriptive)")
    ax.grid(True, axis="x")
    ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False)
    fig.legend(loc="outside lower center", ncol=3, handletextpad=0.3, columnspacing=1.2)
    save(fig, path)


def plot_quartiles(q: pd.DataFrame, path) -> None:
    ends = q[q["r_text_phylo"].notna()]
    groups = list(dict.fromkeys(ends["group"]))
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 0.4 * len(groups) + 1.9), constrained_layout=True, sharey=True)
    ys = np.arange(len(groups))[::-1]
    allq = q[q["group"] == "all mammals"]
    for ax, column, xlabel in ((axes[0], "mean_text_distance", "mean text distance"), (axes[1], "r_text_phylo", "Spearman r, text vs phylogeny")):
        for qname, colour, offset in (("Q1 (least viewed)", MUTED, 0.14), ("Q4 (most viewed)", BLUE, -0.14)):
            sub = q[q["quartile"] == qname].set_index("group").loc[groups]
            ax.scatter(sub[column], ys + offset, s=36, color=colour, edgecolor=SURFACE, linewidth=1, zorder=3, label=qname)
        ax.set_xlabel(xlabel)
        ax.grid(True, axis="x")
        ax.set_axisbelow(True)
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)
    axes[1].axvline(0, color=MUTED, linewidth=0.8)
    axes[0].set_yticks(ys, groups)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=2)
    fig.suptitle(f"Exploratory: least- and most-viewed quartiles of species pairs (by mean log pageviews)\n"
                 f"all mammals: {int(allq['n_pairs'].sum()):,} pairs; quartiles recomputed within each order", x=0.01, ha="left", fontweight="bold", fontsize=10)
    save(fig, path)


def main() -> None:
    args = base_parser("Week 6 figures from the full-scale result tables.").parse_args()
    cfg = load_config(args.config, args.clade or "mammalia")
    apply_style()
    t, out = cfg.tables_dir(), cfg.week6_figures_dir()
    h1 = read_csv(t / "h1_mantel.csv")
    partial = read_csv(t / "h2_partial_mantel.csv")
    head_h1, head_h2 = h1[h1["is_primary"]].iloc[0], partial[partial["is_primary"]].iloc[0]
    plot_controls(partial, out / "h2_controls.png")
    plot_mrm(read_csv(t / "h2_mrm.csv"), out / "h2_mrm.png")
    plot_orders(read_csv(t / "per_order.csv"), head_h1, out / "per_order.png")
    if (t / "posterior.csv").exists():
        plot_posterior(read_csv(t / "posterior.csv"), head_h1, head_h2, out / "posterior.png")
    if (t / "exploratory_pageview_quartiles.csv").exists():
        plot_quartiles(read_csv(t / "exploratory_pageview_quartiles.csv"), out / "pageview_quartiles.png")
    print(f"[figures] wrote {sorted(p.name for p in out.glob('*.png'))} in {out.relative_to(cfg.root)}")


if __name__ == "__main__":
    main()
