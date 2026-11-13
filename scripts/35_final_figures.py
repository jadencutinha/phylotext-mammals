from __future__ import annotations

import json
import re

import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from lineage_charisma.config import base_parser, load_config
from lineage_charisma.distances import load_matrix, text_matrix_path
from lineage_charisma.embed import read_species_order
from lineage_charisma.io_utils import read_csv, write_csv
from lineage_charisma.plotstyle import AXIS, BLUE, GRID, HALO, INK, INK_2, MUTED, ORANGE, SURFACE, apply_style, plt, save_final
from lineage_charisma.sanity import knn_label_agreement

# Week 8 publication figures, read from the saved result tables (and, for figure 2 only, the pilot distance matrices).
# One style for every figure: 7.2 in wide (a full text column), 9 pt text, blue = pilot (Carnivora) and orange = full scale
# (all mammals), grey = reference or no-control values. Filled markers mean one-sided p < 0.05, hollow markers p >= 0.05, so no
# reading depends on colour alone. Aqua is not used (it fails the deuteranopia check against grey). Exploratory and post-hoc
# figures say so in the title. Output: reports/figures/final/figNN_*.png (300 dpi) and .pdf.

WIDTH = 7.2
PILOT, FULL = BLUE, ORANGE
MASKS = [("none", "no masking"), ("taxonomy", "taxonomy\nmasked"), ("strict", "strict\n(primary)")]
PRIMARY_MODEL, PRIMARY_RULE = "BAAI/bge-large-en-v1.5", "chunk"
ALPHA = 0.05


def p_text(p: float, perms: int) -> str:
    return f"p ≤ {1 / (perms + 1):.4f}" if p <= 1.0001 / (perms + 1) else f"p = {p:.4f}"


def ann(ax, text, xy, **kw):
    """Annotation with a true minus sign and a surface halo, so labels stay legible over grid and zero lines."""
    kw.setdefault("path_effects", HALO)
    return ax.annotate(re.sub(r"-(?=\d)", "−", text), xy, **kw)


def dot(ax, x, y, colour, p=None, size=48, zorder=3, **kw):
    """Filled when p < 0.05 (or no p), hollow otherwise."""
    filled = p is None or p < ALPHA
    ax.scatter([x], [y], s=size, facecolor=colour if filled else SURFACE, edgecolor=colour, linewidth=1.6, zorder=zorder, **kw)


def tidy(ax, axis="x"):
    ax.grid(True, axis=axis)
    ax.set_axisbelow(True)
    if axis == "x":
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)


def zero(ax, axis="x"):
    (ax.axvline if axis == "x" else ax.axhline)(0, color=MUTED, linewidth=0.8, zorder=1)


def sig_legend(fig_or_ax, extra=(), **kw):
    handles = list(extra) + [
        Line2D([], [], marker="o", linestyle="", markersize=6.5, markerfacecolor=INK_2, markeredgecolor=INK_2, label="one-sided p < 0.05"),
        Line2D([], [], marker="o", linestyle="", markersize=6.5, markerfacecolor=SURFACE, markeredgecolor=INK_2, markeredgewidth=1.4,
               label="p ≥ 0.05")]
    fig_or_ax.legend(handles=handles, **kw)


def series_handle(colour, label):
    return Line2D([], [], marker="o", linestyle="", markersize=6.5, markerfacecolor=colour, markeredgecolor=colour, label=label)


# ---------------------------------------------------------------- 1. pipeline overview
def fig_pipeline(out):
    fig, ax = plt.subplots(figsize=(WIDTH, 4.6))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 62)
    ax.axis("off")
    cols = [
        (1, "Data sources", ["MDD v2 taxonomy\n+ synonyms", "Upham et al. 2019 tree\n(MCC + posterior)", "English Wikipedia\narticles",
                             "EltonTraits 1.0\n+ COMBINE traits", "Wikipedia pageviews\n(2023-01 to 2025-12)"]),
        (21, "Matching", ["names reconciled to\nMDD and tree tips\n(exact, synonym,\nreviewed fuzzy)", "article lead +\ndescription sections;\nstubs (< 50 tokens)\ndropped",
                          "pilot: 284 Carnivora\nfull: 4,611 mammals"]),
        (41, "Masking", ["strict (primary):\nscientific names +\ncommon-name\nhead nouns", "taxonomy:\nscientific names only", "none"]),
        (61, "Embeddings", ["bge-large-en-v1.5\n(primary);\nall-mpnet-base-v2", "chunk (primary):\nmean of overlapping\npassages; truncate:\nfirst passage only"]),
        (81, "Distance matrices", ["text: cosine\ndistance", "phylogeny: patristic\ndistance (My)", "ecology: Gower\ndistance, 19–20 traits",
                                   "length: |Δ log tokens|", "pageviews: |Δ log views|"]),
    ]
    box_w = 18
    centres = {}
    for x0, title, items in cols:
        ax.text(x0 + box_w / 2, 60, title, ha="center", va="top", fontsize=9.5, fontweight="bold", color=INK)
        n = len(items)
        h = min(8.6, (46 - 1.4 * (n - 1)) / n)
        top = 54
        for i, item in enumerate(items):
            y = top - i * (h + 1.4) - h
            ax.add_patch(FancyBboxPatch((x0, y), box_w, h, boxstyle="round,pad=0,rounding_size=1.2", facecolor="#eef4fc" if x0 < 81 else "#fdf0ea",
                                        edgecolor=AXIS, linewidth=0.8))
            ax.text(x0 + box_w / 2, y + h / 2, item, ha="center", va="center", fontsize=6.4, color=INK, linespacing=1.15)
        centres[x0] = (top - (n * h + (n - 1) * 1.4) / 2)
    for (a, _, _), (b, _, _) in zip(cols, cols[1:]):
        ax.add_patch(FancyArrowPatch((a + box_w + 0.3, 30), (b - 0.3, 30), arrowstyle="-|>", mutation_scale=10, color=INK_2, linewidth=1))
    tests = ("Tests (Spearman, 9,999 permutations, seed 20261004)\nH1: Mantel r(text, phylogeny)   ·   "
             "H2: partial Mantel r(text, phylogeny | ecology + length)   ·   MRM on all four predictors")
    ax.add_patch(FancyBboxPatch((1, 0.6), 98, 6.6, boxstyle="round,pad=0,rounding_size=1.2", facecolor="#f3f2ee", edgecolor=AXIS, linewidth=0.8))
    ax.text(50, 3.9, tests, ha="center", va="center", fontsize=7.0, color=INK, linespacing=1.4)
    ax.add_patch(FancyArrowPatch((90, 9.0), (90, 7.4), arrowstyle="-|>", mutation_scale=10, color=INK_2, linewidth=1))
    ax.set_title("Study pipeline: from species lists and articles to distance matrices and permutation tests", loc="left")
    save_final(fig, out / "fig01_pipeline")


# ---------------------------------------------------------------- 2. masking effect (pilot)
def masking_table(cfg) -> pd.DataFrame:
    """Share of species whose nearest text neighbour is in the same family, per mask level (week 3 sanity check, recomputed from
    the saved pilot matrices), with H1 r from h1_mantel.csv."""
    md, op = cfg.matrices_dir(), cfg.species_order_path()
    order = read_species_order(op)
    families = read_csv(cfg.matched_path()).set_index("species").loc[order, "family"].to_numpy()
    sizes = pd.Series(families).map(pd.Series(families).value_counts()).to_numpy()
    chance = float(((sizes - 1) / (len(order) - 1)).mean())
    h1 = read_csv(cfg.tables_dir() / "h1_mantel.csv").query("method == 'spearman'").set_index(["model", "rule", "mask"])
    rows = []
    for model in cfg["embedding"]["models"]:
        for rule in cfg["embedding"]["rules"]:
            for mask, _ in MASKS:
                matrix = load_matrix(text_matrix_path(md, model, rule, mask), op)[0]
                rows.append({"status": "descriptive (week 3 sanity check, recomputed)", "model": model, "rule": rule, "mask": mask,
                             "nn_same_family": knn_label_agreement(matrix, families, 1), "chance_same_family": chance,
                             "h1_r": h1.loc[(model, rule, mask), "r"], "h1_p": h1.loc[(model, rule, mask), "p"], "n_species": len(order)})
    return pd.DataFrame(rows)


def fig_masking(table, out):
    fig, axes = plt.subplots(1, 2, figsize=(WIDTH, 3.3), constrained_layout=True)
    x = np.arange(len(MASKS))
    prim = table[(table["model"] == PRIMARY_MODEL) & (table["rule"] == PRIMARY_RULE)].set_index("mask")
    others = table[~((table["model"] == PRIMARY_MODEL) & (table["rule"] == PRIMARY_RULE))]
    for ax, col, ylabel, fmt in ((axes[0], "nn_same_family", "species whose nearest text neighbour\nis in the same family (%)", "{:.0%}"),
                                 (axes[1], "h1_r", "H1 Mantel r (Spearman),\ntext vs phylogenetic distance", "{:.3f}")):
        scale = 100 if col == "nn_same_family" else 1
        vals = [prim.loc[m, col] for m, _ in MASKS]
        ax.bar(x, np.array(vals) * scale, width=0.56, color=PILOT, edgecolor=SURFACE, linewidth=1, zorder=2)
        for i, (m, _) in enumerate(MASKS):
            ann(ax, fmt.format(vals[i]), (x[i], vals[i] * scale), xytext=(0, 3), textcoords="offset points", ha="center", va="bottom", color=INK)
            o = others[others["mask"] == m][col] * scale
            ax.scatter(np.full(len(o), x[i] + 0.38), o, s=16, color=MUTED, zorder=3, label="other model / rule" if i == 0 else None)
        ax.set_xticks(x, [label for _, label in MASKS])
        ax.set_xlim(-0.5, len(MASKS) - 0.3)
        ax.set_ylabel(ylabel)
        ax.grid(True, axis="y")
        ax.set_axisbelow(True)
    axes[0].axhline(table["chance_same_family"].iloc[0] * 100, color=INK_2, linewidth=0.9, linestyle="--", zorder=1)
    ann(axes[0], f"chance ({table['chance_same_family'].iloc[0]:.0%})", (-0.48, table["chance_same_family"].iloc[0] * 100), xytext=(0, 3),
                     textcoords="offset points", ha="left", va="bottom", color=INK_2, fontsize=8)
    axes[0].set_ylim(0, 100)
    axes[0].set_title("Family structure in text space")
    axes[1].set_title("Text–phylogeny correlation (H1)")
    axes[1].legend(loc="upper left", handletextpad=0.2)
    n = int(table["n_species"].iloc[0])
    fig.suptitle(f"Pilot ({n} Carnivora): masking common names removes most of the family signal", x=0.01, ha="left", fontweight="bold", fontsize=10)
    save_final(fig, out / "fig02_masking_effect")


# ---------------------------------------------------------------- 3. pilot vs full-scale summary
def fig_summary(t, out):
    m = t / "mammalia"
    ph1, fh1 = read_csv(t / "h1_mantel.csv").query("is_primary").iloc[0], read_csv(m / "h1_mantel.csv").query("is_primary").iloc[0]
    ph2, fh2 = read_csv(t / "h2_partial_mantel.csv").query("is_primary").iloc[0], read_csv(m / "h2_partial_mantel.csv").query("is_primary").iloc[0]
    rows = [("H1: Mantel r\n(text vs phylogeny)", (ph1["r"], ph1["p"], 9999, 284), (fh1["r"], fh1["p"], int(fh1["permutations"]), int(fh1["n_species"]))),
            ("H2 baseline: raw r\n(species with ecology)", (ph2["raw_r"], ph2["raw_p"], 9999, int(ph2["n_species"])),
             (fh2["raw_r"], fh2["raw_p"], 9999, int(fh2["n_species"]))),
            ("H2: partial r\n(| ecology + length)", (ph2["partial_r"], ph2["partial_p"], 9999, int(ph2["n_species"])),
             (fh2["partial_r"], fh2["partial_p"], int(fh2["permutations"]), int(fh2["n_species"])))]
    decisions = {"pilot": ["H1 supported (r > 0, p < 0.05)", "", f"rule (a) not met; rule (b) not met (partial = {ph2['share_of_raw']:.1%} of raw)"],
                 "full": ["H1 not supported", "", "rule (a) not met; rule (b) undefined (raw r ≤ 0)"]}
    fig, ax = plt.subplots(figsize=(WIDTH, 3.5), constrained_layout=True)
    ys = np.arange(len(rows))[::-1]
    for y, (label, pilot, full) in zip(ys, rows):
        for (r, p, perms, n), colour, off, key in ((pilot, PILOT, 0.15, "pilot"), (full, FULL, -0.15, "full")):
            dot(ax, r, y + off, colour, p)
            dec = decisions[key][len(rows) - 1 - y]
            ann(ax, f"r = {r:.3f}, {p_text(p, perms)}, n = {n:,}" + (f"  —  {dec}" if dec else ""), (r, y + off), xytext=(9, 0),
                        textcoords="offset points", va="center", fontsize=7.6, color=INK)
    zero(ax)
    ax.set_yticks(ys, [r[0] for r in rows])
    ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.set_xlim(-0.06, 0.16)
    ax.set_xlabel("Spearman correlation (one-sided permutation p)")
    tidy(ax)
    sig_legend(fig, [series_handle(PILOT, "pilot: Carnivora"), series_handle(FULL, "full scale: all mammals")], loc="outside lower center", ncol=4,
               handletextpad=0.2, columnspacing=1.2)
    ax.set_title("Pre-registered decisions: pilot vs full scale (primary specification, MCC tree)")
    save_final(fig, out / "fig03_pilot_vs_full")


# ---------------------------------------------------------------- 4. H2 controls, pilot and full scale
def fig_controls(t, out):
    p, f = read_csv(t / "h2_partial_mantel.csv"), read_csv(t / "mammalia" / "h2_partial_mantel.csv")
    p_by = p[p["analysis"].isin(["primary", "secondary"]) & (p["mask"] == "strict") & (p["rule"] == "chunk") & (p["model"] == PRIMARY_MODEL)].set_index("controls")
    f_by = f[f["analysis"].isin(["primary", "secondary S1", "secondary S2", "secondary S3", "secondary S4"])].set_index("controls")
    items = [(None, "no control (raw r)"), ("ecology", "ecology"), ("length", "article length"), ("pageviews", "pageviews"),
             ("ecology + length", "ecology + length\n(primary H2 model)"), ("ecology + length + pageviews", "ecology + length\n+ pageviews")]
    fig, ax = plt.subplots(figsize=(WIDTH, 3.9), constrained_layout=True)
    ys = np.arange(len(items))[::-1]
    for y, (key, _) in zip(ys, items):
        for by, colour, off in ((p_by, PILOT, 0.14), (f_by, FULL, -0.14)):
            head = by.loc["ecology + length"]
            r, pv, perms = (head["raw_r"], head["raw_p"], 9999) if key is None else (by.loc[key, "partial_r"], by.loc[key, "partial_p"], 9999)
            dot(ax, r, y + off, colour, pv)
            ann(ax, f"{r:.3f} ({p_text(pv, perms)})", (r, y + off), xytext=(8, 0), textcoords="offset points", va="center", fontsize=7.4, color=INK_2)
    ax.axhspan(ys[4] - 0.42, ys[4] + 0.42, color=GRID, alpha=0.5, zorder=0, linewidth=0)
    zero(ax)
    ax.set_yticks(ys, [i[1] for i in items])
    ax.set_ylim(-0.6, len(items) - 0.4)
    ax.set_xlim(-0.05, 0.075)
    ax.set_xlabel("Spearman partial Mantel r, text vs phylogenetic distance given the controls")
    tidy(ax)
    n_p, n_f = int(p_by.loc["ecology + length", "n_species"]), int(f_by.loc["ecology + length", "n_species"])
    sig_legend(fig, [series_handle(PILOT, f"pilot ({n_p} Carnivora)"), series_handle(FULL, f"full scale ({n_f:,} mammals)")], loc="outside lower center",
               ncol=4, handletextpad=0.2, columnspacing=1.2)
    ax.set_title("H2: the text–phylogeny correlation before and after each control")
    save_final(fig, out / "fig04_h2_controls")


# ---------------------------------------------------------------- 5. per order, with stub share
def fig_orders(t, out):
    m = t / "mammalia"
    orders = read_csv(m / "per_order.csv").sort_values("h1_r")
    stubs = read_csv(m / "stubs_by_order.csv").set_index("order")
    allr = read_csv(m / "h1_mantel.csv").query("is_primary").iloc[0]["r"]
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(WIDTH, 4.0), constrained_layout=True, sharey=True, gridspec_kw={"width_ratios": [2.3, 1]})
    ys = np.arange(len(orders))
    for y, (_, r) in zip(ys, orders.iterrows()):
        ax.hlines(y, min(r["h1_r"], r["partial_r"]), max(r["h1_r"], r["partial_r"]), color=AXIS, linewidth=2, zorder=1)
        ax.scatter([r["h1_r"]], [y + 0.12], marker="o", s=46, facecolor=PILOT if r["h1_p"] < ALPHA else SURFACE, edgecolor=PILOT, linewidth=1.6, zorder=3)
        ax.scatter([r["partial_r"]], [y - 0.12], marker="s", s=40, facecolor=FULL if r["partial_p"] < ALPHA else SURFACE, edgecolor=FULL,
                   linewidth=1.6, zorder=3)
    zero(ax)
    ax.axvline(allr, color=INK_2, linewidth=0.9, linestyle="--")
    ann(ax, f"all mammals, H1 r = {allr:.3f}", (allr, -0.75), xytext=(4, 0), textcoords="offset points", ha="left", va="center",
                fontsize=7.4, color=INK_2)
    ax.set_yticks(ys, [f"{o} (n = {n:,})" for o, n in zip(orders["order"], orders["n_species"])])
    ax.set_ylim(-1.2, len(orders) - 0.4)
    ax.set_xlabel("Spearman Mantel r within the order")
    tidy(ax)
    share = stubs.loc[orders["order"], "stub_rate"].to_numpy() * 100
    bx.barh(ys, share, height=0.56, color=MUTED, edgecolor=SURFACE)
    for y, s_ in zip(ys, share):
        ann(bx, f"{s_:.0f}%" if s_ >= 0.5 or s_ == 0 else "<1%", (s_, y), xytext=(3, 0), textcoords="offset points", va="center", fontsize=7.4, color=INK_2)
    bx.set_xlim(0, max(share) * 1.3)
    bx.set_xlabel("articles dropped as stubs\n(< 50 tokens, %)")
    tidy(bx)
    handles = [Line2D([], [], marker="o", linestyle="", markersize=6.5, markerfacecolor=PILOT, markeredgecolor=PILOT, label="H1 r"),
               Line2D([], [], marker="s", linestyle="", markersize=6, markerfacecolor=FULL, markeredgecolor=FULL, label="H2 partial r (ecology + length)"),
               Line2D([], [], marker="o", linestyle="", markersize=6.5, markerfacecolor=SURFACE, markeredgecolor=INK_2, label="hollow: p ≥ 0.05")]
    fig.legend(handles=handles, loc="outside lower center", ncol=3, handletextpad=0.2, columnspacing=1.4)
    fig.suptitle("Full scale, within each order with at least 50 species (descriptive)", x=0.01, ha="left", fontweight="bold", fontsize=10)
    save_final(fig, out / "fig05_per_order")


# ---------------------------------------------------------------- 6. posterior trees
def fig_posterior(t, out):
    m = t / "mammalia"
    ph1, ph2 = read_csv(t / "h1_posterior.csv"), read_csv(t / "h2_posterior.csv")
    fp = read_csv(m / "posterior.csv")
    mcc = {"pilot": (read_csv(t / "h1_mantel.csv").query("is_primary").iloc[0]["r"], read_csv(t / "h2_partial_mantel.csv").query("is_primary").iloc[0]["partial_r"]),
           "full": (read_csv(m / "h1_mantel.csv").query("is_primary").iloc[0]["r"], read_csv(m / "h2_partial_mantel.csv").query("is_primary").iloc[0]["partial_r"])}
    rows = [("pilot H1 r", ph1["r"], ph1["p"], mcc["pilot"][0], PILOT), ("pilot H2 partial r", ph2["partial_r"], ph2["partial_p"], mcc["pilot"][1], PILOT),
            ("full-scale H1 r", fp["h1_r"], fp["h1_p"], mcc["full"][0], FULL), ("full-scale H2 partial r", fp["partial_r"], fp["partial_p"], mcc["full"][1], FULL)]
    rng = np.random.default_rng(20261004)   # jitter only
    fig, ax = plt.subplots(figsize=(WIDTH, 3.6), constrained_layout=True)
    ys = np.arange(len(rows))[::-1]
    for y, (label, r, p, m_, colour) in zip(ys, rows):
        jitter = rng.uniform(-0.17, 0.17, len(r))
        sig = p.to_numpy() < ALPHA
        ax.scatter(r[sig], y + jitter[sig], s=13, color=colour, alpha=0.75, edgecolor="none", zorder=2)
        ax.scatter(r[~sig], y + jitter[~sig], s=13, facecolor="none", edgecolor=colour, alpha=0.8, linewidth=0.8, zorder=2)
        ax.vlines(np.median(r), y - 0.3, y + 0.3, color=INK, linewidth=1.6, zorder=3)
        ax.scatter([m_], [y], marker="D", s=40, color=INK, edgecolor=SURFACE, linewidth=1, zorder=4)
        ann(ax, f"{len(r)} trees: median {np.median(r):.3f}, {int(sig.sum())}/{len(r)} with p < 0.05", (r.max(), y), xytext=(8, 0),
                    textcoords="offset points", va="center", fontsize=7.4, color=INK_2)
    zero(ax)
    ax.set_yticks(ys, [r[0] for r in rows])
    ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.set_xlim(-0.08, 0.27)
    ax.set_xlabel("Spearman r on each posterior tree (pilot 9,999 permutations per tree, full scale 999)")
    tidy(ax)
    handles = [Line2D([], [], marker="o", linestyle="", markersize=5, markerfacecolor=INK_2, markeredgecolor=INK_2, label="one tree, p < 0.05"),
               Line2D([], [], marker="o", linestyle="", markersize=5, markerfacecolor=SURFACE, markeredgecolor=INK_2, label="one tree, p ≥ 0.05"),
               Line2D([], [], color=INK, linewidth=1.6, label="median over trees"),
               Line2D([], [], marker="D", linestyle="", markersize=6, color=INK, label="MCC tree (the registered test)")]
    fig.legend(handles=handles, loc="outside lower center", ncol=4, handletextpad=0.3, columnspacing=1.2)
    ax.set_title("Robustness to tree uncertainty: the primary tests on trees drawn from the Upham et al. posterior")
    save_final(fig, out / "fig06_posterior_trees")


# ---------------------------------------------------------------- 7. pageview quartiles (exploratory)
def fig_quartiles(t, out):
    pil = read_csv(t / "exploratory_well_known.csv")
    ful = read_csv(t / "mammalia" / "exploratory_pageview_quartiles.csv").query("group == 'all mammals'")
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(WIDTH, 3.3), constrained_layout=True)
    x = np.arange(4)
    labels = ["Q1\nleast viewed", "Q2", "Q3", "Q4\nmost viewed"]
    for df, colour, name in ((pil, PILOT, "pilot (Carnivora)"), (ful, FULL, "full scale (all mammals)")):
        ax.plot(x, df["mean_text_distance"], color=colour, linewidth=2, marker="o", markersize=6, markeredgecolor=SURFACE, label=name)
    ax.set_xticks(x, labels)
    ax.set_ylabel("mean text distance (cosine)")
    ax.set_title("Popular pairs are described more alike")
    ax.grid(True, axis="y")
    ax.set_axisbelow(True)
    ax.legend(loc="lower left")
    ends = ["Q1 (least viewed)", "Q4 (most viewed)"]
    for df, colour, off, perms in ((pil, PILOT, -0.1, 9999), (ful, FULL, 0.1, 999)):
        sub = df.set_index("quartile").loc[ends]
        for i, (q, r) in enumerate(sub.iterrows()):
            bx.scatter([i + off], [r["r_text_phylo"]], s=48, facecolor=colour if r["p"] < ALPHA else SURFACE, edgecolor=colour, linewidth=1.6, zorder=3)
            ann(bx, f"{r['r_text_phylo']:.3f}", (i + off, r["r_text_phylo"]), xytext=(0, 6), textcoords="offset points", ha="center",
                        fontsize=7.4, color=INK_2)
        bx.plot([off, 1 + off], sub["r_text_phylo"], color=colour, linewidth=1, zorder=2)
    bx.axhline(0, color=MUTED, linewidth=0.8)
    bx.set_xticks([0, 1], ["Q1\nleast viewed", "Q4\nmost viewed"])
    bx.set_xlim(-0.5, 1.5)
    bx.set_ylim(-0.01, 0.155)
    bx.set_ylabel("Spearman r, text vs phylogeny\nwithin the quartile")
    bx.set_title("Obscure pairs track the tree more closely")
    bx.grid(True, axis="y")
    bx.set_axisbelow(True)
    fig.suptitle("Exploratory: species pairs split into quartiles by mean log pageviews", x=0.01, ha="left", fontweight="bold", fontsize=10)
    save_final(fig, out / "fig07_pageview_quartiles")


# ---------------------------------------------------------------- 8. length-matched deciles (post-hoc)
def fig_length_matched(t, out):
    m = t / "mammalia"
    bins = read_csv(m / "posthoc_length_matched_bins.csv")
    head = read_csv(m / "posthoc_length_matched.csv").iloc[0]
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(WIDTH, 3.4), constrained_layout=True, sharex=True)
    x = bins["bin"].to_numpy() + 1
    ax.bar(x, bins["r_spearman"], width=0.62, color=FULL, edgecolor=SURFACE, zorder=2)
    ax.axhline(head["r_lm"], color=INK, linewidth=1.2, linestyle="--", zorder=3)
    ax.axhline(head["unmatched_h1_r"], color=MUTED, linewidth=1.2, linestyle=":", zorder=3)
    ann(ax, f"r_LM = {head['r_lm']:+.4f}\ntwo-sided p = {head['p_two_sided']:.4f}", (5.5, 0.075), ha="center", fontsize=7.6, color=INK)
    ax.annotate("", (6.0, head["r_lm"]), xytext=(5.6, 0.07), arrowprops=dict(arrowstyle="-", color=INK, linewidth=0.8))
    ann(ax, f"unmatched H1 r = {head['unmatched_h1_r']:.3f}", (0.6, head["unmatched_h1_r"]), xytext=(0, -10), textcoords="offset points",
                fontsize=7.6, color=INK_2)
    ax.axhline(0, color=MUTED, linewidth=0.8)
    ax.set_ylabel("Spearman r, text vs phylogeny\nwithin the decile")
    ax.set_title("Correlation within each length decile")
    bx.bar(x, bins["different_minus_same"], width=0.62, color=FULL, edgecolor=SURFACE, zorder=2)
    bx.axhline(0, color=MUTED, linewidth=0.8)
    bx.set_ylabel("mean text distance: different-order\nminus same-order pairs")
    bx.set_title("Same-order pairs are closer in every decile")
    for a in (ax, bx):
        a.set_xticks(x, [str(int(v)) for v in x])
        a.set_xlabel("text-length decile")
        a.grid(True, axis="y")
        a.set_axisbelow(True)
    fig.supxlabel("text-length decile: 1 = shortest (50–64 tokens) to 10 = longest (555 tokens or more); equal species counts",
                  fontsize=8.5, color=INK_2)
    fig.suptitle(f"Post-hoc (amendment 4, section D): length-matched comparison, {int(head['n_species']):,} mammals in 10 deciles, "
                 f"{int(head['permutations']):,} within-decile permutations", x=0.01, ha="left", fontweight="bold", fontsize=10)
    save_final(fig, out / "fig08_length_matched_posthoc")


# ---------------------------------------------------------------- 9. ADW vs Wikipedia
def fig_register(t, out):
    s = t / "secondary"
    rc = read_csv(s / "register_comparison.csv").set_index(["test", "source"])
    bs = read_csv(s / "register_bootstrap.csv").set_index("test")
    tests = [("H1", "H1: Mantel r"), ("H2 partial", "H2: partial r\n(| ecology + length)")]
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(WIDTH, 2.9), constrained_layout=True, sharey=True, gridspec_kw={"width_ratios": [1.3, 1]})
    ys = np.arange(len(tests))[::-1]
    for y, (key, _) in zip(ys, tests):
        w, a = rc.loc[(key, "wikipedia")], rc.loc[(key, "adw")]
        ax.plot([w["r"], a["r"]], [y, y], color=AXIS, linewidth=2, zorder=1)
        dot(ax, w["r"], y, MUTED, w["p"])
        dot(ax, a["r"], y, PILOT, a["p"])
        for row, ha, dx in ((w, "right", -8), (a, "left", 8)):
            ann(ax, f"{row['r']:.3f}\n{p_text(row['p'], 9999)}", (row["r"], y), xytext=(dx, 0), textcoords="offset points", ha=ha, va="center",
                        fontsize=7.4, color=INK_2)
        b = bs.loc[key]
        bx.hlines(y, b["ci_low"], b["ci_high"], color=INK, linewidth=1.6)
        bx.scatter([b["delta"]], [y], s=44, color=INK, zorder=3)
        ann(bx, f"{b['delta']:+.3f} [{b['ci_low']:+.3f}, {b['ci_high']:+.3f}]", (b["delta"], y), xytext=(0, 8), textcoords="offset points",
                    ha="center", fontsize=7.4, color=INK)
    zero(ax)
    zero(bx)
    ax.set_yticks(ys, [f"{lab}\nn = {int(bs.loc[k, 'n_species'])}" for k, lab in tests])
    ax.set_ylim(-0.6, len(tests) - 0.4)
    ax.set_xlim(-0.04, 0.12)
    ax.set_xlabel("Spearman r, one-sided permutation p")
    ax.set_title("Each source on the same species")
    bx.set_xlim(-0.02, 0.12)
    bx.set_xlabel("Δ r (ADW − Wikipedia), 95% interval")
    bx.set_title(f"Paired species bootstrap ({int(bs['n_boot'].iloc[0]):,} draws)")
    tidy(ax)
    tidy(bx)
    sig_legend(fig, [series_handle(MUTED, "Wikipedia"), series_handle(PILOT, "Animal Diversity Web")], loc="outside lower center", ncol=4,
               handletextpad=0.2, columnspacing=1.2)
    fig.suptitle("Secondary analysis A: register comparison, Wikipedia vs Animal Diversity Web (Carnivora with usable ADW text)", x=0.01,
                 ha="left", fontweight="bold", fontsize=10)
    save_final(fig, out / "fig09_adw_vs_wikipedia")


# ---------------------------------------------------------------- 10. charisma (secondary B), once it has run
def fig_charisma(t, out) -> bool:
    s = t / "secondary"
    if not (s / "charisma_agreement.csv").exists():
        return False
    agree = read_csv(s / "charisma_agreement.csv")
    floor = json.loads((s / "charisma_agreement.csv.meta.json").read_text()).get("floor", 0.5)
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(WIDTH, 3.0), constrained_layout=True)
    ys = np.arange(len(agree))[::-1]
    for y, (_, r) in zip(ys, agree.iterrows()):
        ax.scatter([r["author_vs_model"]], [y + 0.12], s=46, marker="o", color=INK, zorder=3)
        ax.scatter([r["run1_vs_run2"]], [y - 0.12], s=40, marker="s", color=MUTED, zorder=3)
        ann(ax, "kept" if r["retained"] else "dropped", (1.0, y), xytext=(4, 0), textcoords="offset points", va="center", fontsize=7.4, color=INK_2)
    ax.axvline(floor, color=INK_2, linewidth=0.9, linestyle="--")
    ax.set_yticks(ys, [d.replace("_", " ") for d in agree["dimension"]])
    ax.set_xlim(min(0, agree[["author_vs_model", "run1_vs_run2"]].min().min() - 0.05), 1.0)
    ax.set_xlabel(f"Spearman agreement (drop rule: below {floor})")
    ax.set_title("Agreement per dimension")
    ax.set_ylim(-0.6, len(agree) - 0.4)
    fig.legend(handles=[series_handle(INK, f"author vs model (n = {int(agree['n_validation'].iloc[0])} texts)"),
                        Line2D([], [], marker="s", linestyle="", markersize=6, color=MUTED, label=f"run 1 vs run 2 (n = {int(agree['n_species'].iloc[0])} species)"),
                        Line2D([], [], marker="o", linestyle="", markersize=6.5, markerfacecolor=SURFACE, markeredgecolor=PILOT, markeredgewidth=1.4,
                               label="B tests: hollow = p ≥ 0.05")],
               loc="outside lower center", ncol=3, handletextpad=0.2, columnspacing=1.2)
    tidy(ax)
    if (s / "charisma_tests.csv").exists():
        tests = read_csv(s / "charisma_tests.csv")
        b1, b2 = tests.iloc[0], tests.iloc[1]
        items = [("B1: text ~ charisma\n| length + pageviews", b1["r"], b1["p"]), ("primary H2 partial r\n(same run)", b2["primary_h2_partial_r_same_run"], None),
                 ("B2: H2 partial r\n+ charisma control", b2["r"], b2["p"])]
        yb = np.arange(len(items))[::-1]
        for y, (label, r, p) in zip(yb, items):
            dot(bx, r, y, PILOT if p is not None else MUTED, p)
            ann(bx, f"{r:.3f}" + (f" ({p_text(p, 9999)})" if p is not None else ""), (r, y), xytext=(8, 0), textcoords="offset points",
                        va="center", fontsize=7.4, color=INK_2)
        zero(bx)
        bx.set_yticks(yb, [i[0] for i in items])
        bx.set_ylim(-0.6, len(items) - 0.4)
        bx.set_xlabel("Spearman partial Mantel r")
        bx.set_title("Tests B1 and B2")
        tidy(bx)
    fig.suptitle("Secondary analysis B: content charisma scored by claude-sonnet-5-5 (pilot)", x=0.01, ha="left", fontweight="bold", fontsize=10)
    save_final(fig, out / "fig10_charisma")
    return True


def main() -> None:
    args = base_parser("Week 8 publication figures (300 dpi PNG + PDF) in reports/figures/final/.").parse_args()
    cfg = load_config(args.config, None)
    apply_style()
    t, out = cfg.root / "reports" / "tables", cfg.reports / "figures" / "final"
    out.mkdir(parents=True, exist_ok=True)
    fig_pipeline(out)
    masking = masking_table(cfg)
    write_csv(masking, t / "final_masking_effect.csv", source="pilot text matrices and reports/tables/h1_mantel.csv", k=1)
    fig_masking(masking, out)
    fig_summary(t, out)
    fig_controls(t, out)
    fig_orders(t, out)
    fig_posterior(t, out)
    fig_quartiles(t, out)
    fig_length_matched(t, out)
    fig_register(t, out)
    made_b = fig_charisma(t, out)
    print(f"[final figures] wrote {len(list(out.glob('*.png')))} figures to {out.relative_to(cfg.root)}"
          + ("" if made_b else "; figure 10 (charisma) skipped: analysis B has not run"))


if __name__ == "__main__":
    main()
