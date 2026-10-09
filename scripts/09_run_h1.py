from __future__ import annotations

import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm

from lineage_charisma.config import base_parser, load_config
from lineage_charisma.corpus import read_corpus
from lineage_charisma.distances import assert_same_order, load_matrix, phylo_matrix_path, text_matrix_path
from lineage_charisma.io_utils import write_csv
from lineage_charisma.plotstyle import AXIS, BLUE, BLUE_RAMP, HALO, INK, INK_2, MASK_COLORS, MASK_NAMES, MEDIAN_LINE, MUTED, SURFACE, apply_style, combo_label, plt, save
from lineage_charisma.sanity import binned_median
from lineage_charisma.stats import METHODS, mantel, submatrix, upper_triangle

N_FAMILIES = 3


def plot_scatter(text, phylo, r, p, label, path) -> None:
    x, y = upper_triangle(phylo), upper_triangle(text)
    medians = binned_median(x, y, np.arange(0, x.max() + 5, 5), min_count=30)
    fig, ax = plt.subplots(figsize=(7.2, 4.8), constrained_layout=True)
    hexes = ax.hexbin(x, y, gridsize=55, cmap=BLUE_RAMP, norm=LogNorm(), mincnt=1, linewidths=0)
    ax.plot(medians["x"], medians["y"], label="median text distance per 5 My bin", **MEDIAN_LINE)
    fig.legend(loc="outside lower center", handlelength=2.2)
    ax.set_xlabel("Phylogenetic distance (patristic, million years)")
    ax.set_ylabel("Text distance (1 − cosine similarity)")
    ax.set_title(f"Text against phylogenetic distance, {len(x):,} pairs ({label})\nSpearman Mantel r = {r:.3f}, one-sided p = {p:.4f}")
    ax.grid(True, axis="y")
    ax.set_axisbelow(True)
    bar = fig.colorbar(hexes, ax=ax, pad=0.02)
    bar.set_label("species pairs per hexagon (log scale)", color=INK_2)
    bar.outline.set_visible(False)
    save(fig, path)


def plot_null(result, label, path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 3.8), constrained_layout=True)
    lo, hi = min(result.null.min(), result.r), max(result.null.max(), result.r)
    pad = 0.04 * (hi - lo)
    ax.hist(result.null, bins=60, color=BLUE, edgecolor=SURFACE, linewidth=0.4)
    ax.axvline(result.r, color=INK, linewidth=1.6)
    side = "right" if result.r < (lo + hi) / 2 else "left"
    ax.annotate(f"observed r = {result.r:.3f}\np = {result.p:.4f}", xy=(result.r, 0.9), xycoords=("data", "axes fraction"),
                xytext=(6 if side == "right" else -6, 0), textcoords="offset points", ha="left" if side == "right" else "right",
                va="center", color=INK, path_effects=HALO)
    ax.set_xlim(lo - pad, hi + pad)
    ax.set_xlabel(f"Mantel r ({result.method}) under permutation of species")
    ax.set_ylabel("permutations")
    ax.set_title(f"Permutation null and observed r, {result.permutations:,} permutations ({label})")
    ax.grid(True, axis="y")
    ax.set_axisbelow(True)
    save(fig, path)


def plot_variants(table: pd.DataFrame, model_rules, masks, primary, path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 0.62 * len(model_rules) + 1.7), constrained_layout=True)
    rows = {mr: i for i, mr in enumerate(reversed(model_rules))}
    by_combo = table.set_index(["model", "rule", "mask"])
    for mr, row in rows.items():
        values = [by_combo.loc[(*mr, mask), "r"] for mask in masks]
        ax.hlines(row, min(values), max(values), color=AXIS, linewidth=2, zorder=1)
    for mask in masks:
        ax.scatter([by_combo.loc[(*mr, mask), "r"] for mr in rows], list(rows.values()), s=64, color=MASK_COLORS[mask],
                   edgecolor=SURFACE, linewidth=1.5, zorder=3, label=MASK_NAMES[mask])
    pr = by_combo.loc[primary, "r"]
    ax.scatter([pr], [rows[primary[:2]]], s=210, facecolor="none", edgecolor=INK, linewidth=1.4, zorder=4, label="primary specification")
    ax.annotate(f"primary, r = {pr:.3f}", xy=(pr, rows[primary[:2]]), xytext=(0, 13), textcoords="offset points",
                ha="center", color=INK, path_effects=HALO)
    ax.axvline(0, color=MUTED, linewidth=0.8)
    ax.set_ylim(-0.6, len(rows) - 0.3)
    ax.set_yticks(list(rows.values()), [combo_label(*mr) for mr in rows])
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Spearman Mantel r, text distance against phylogenetic distance")
    ax.set_title("H1 across all 12 variants")
    ax.grid(True, axis="x")
    ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False)
    fig.legend(loc="outside lower center", ncol=2, handletextpad=0.3, columnspacing=1.6)
    save(fig, path)


def main() -> None:
    args = base_parser("H1: Mantel tests of text distance against phylogenetic distance.").parse_args()
    cfg = load_config(args.config)
    apply_style()
    scfg, ecfg = cfg["stats"], cfg["embedding"]
    perms, seed, primary_method = int(scfg["permutations"]), int(scfg["seed"]), scfg["method"]
    combos, primary = cfg.combos(), cfg.primary_combo
    methods = [primary_method] + [m for m in METHODS if m != primary_method]

    order_path = cfg.species_order_path()
    phylo, order = load_matrix(phylo_matrix_path(cfg.matrices_dir()), order_path)
    text = {}
    for combo in combos:
        text[combo], combo_order = load_matrix(text_matrix_path(cfg.matrices_dir(), *combo), order_path)
        assert_same_order(combo_order, order, combo_label(*combo))

    rows, results = [], {}
    for combo in combos:
        for method in methods:
            res = mantel(text[combo], phylo, method=method, permutations=perms, seed=seed)
            results[(combo, method)] = res
            is_primary = combo == primary and method == primary_method
            rows.append({"model": combo[0], "rule": combo[1], "mask": combo[2], "method": method, "r": res.r, "p": res.p,
                         "n_pairs": res.n_pairs, "is_primary": is_primary})
            print(f"[h1] {combo_label(*combo):45s} {method:8s} r = {res.r:+.4f}  p = {res.p:.4f}{'   <- primary' if is_primary else ''}")
    table = pd.DataFrame(rows)
    write_csv(table, cfg.tables_dir() / "h1_mantel.csv", permutations=perms, seed=seed, alternative="greater", n_species=len(order))

    corpus = read_corpus(cfg.corpus_path()).set_index("species_id").loc[order]
    family_rows = []
    for family, n in corpus["family"].value_counts().head(N_FAMILIES).items():
        members = corpus.index[corpus["family"] == family].tolist()
        for method in methods:
            res = mantel(submatrix(text[primary], order, members), submatrix(phylo, order, members), method=method, permutations=perms, seed=seed)
            family_rows.append({"family": family, "n_species": int(n), "model": primary[0], "rule": primary[1], "mask": primary[2],
                                "method": method, "r": res.r, "p": res.p, "n_pairs": res.n_pairs})
            print(f"[h1] within {family} (n = {n}) {method:8s} r = {res.r:+.4f}  p = {res.p:.4f}")
    write_csv(pd.DataFrame(family_rows), cfg.tables_dir() / "h1_mantel_families.csv", permutations=perms, seed=seed, alternative="greater")

    fig_dir, head = cfg.week4_figures_dir(), results[(primary, primary_method)]
    model_rules = [(m, r) for m in ecfg["models"] for r in ecfg["rules"]]
    masks = [m for m in MASK_COLORS if m in cfg.mask_levels]
    plot_scatter(text[primary], phylo, head.r, head.p, combo_label(*primary), fig_dir / "h1_scatter_primary.png")
    plot_null(head, combo_label(*primary), fig_dir / "h1_null_primary.png")
    plot_variants(table[table["method"] == primary_method], model_rules, masks, primary, fig_dir / "h1_variants.png")
    print(f"[h1] wrote {cfg.tables_dir().relative_to(cfg.root)}/h1_mantel.csv, h1_mantel_families.csv and 3 figures in {fig_dir.relative_to(cfg.root)}")


if __name__ == "__main__":
    main()
