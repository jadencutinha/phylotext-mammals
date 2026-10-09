from __future__ import annotations

from datetime import datetime, timezone

import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm

from lineage_charisma.config import base_parser, load_config
from lineage_charisma.corpus import length_summary, read_corpus
from lineage_charisma.distances import assert_same_order, load_matrix, phylo_matrix_path, reindex_matrix, reindex_rows, text_matrix_path, tree_tip_order
from lineage_charisma.embed import embedding_paths, load_embeddings, model_slug
from lineage_charisma.io_utils import atomic_write_text, read_csv, read_json
from lineage_charisma.masking import MASK_LEVELS, leak_counts, text_column
from lineage_charisma.phylo import load_tree
from lineage_charisma.plotstyle import AXIS, BLUE, BLUE_RAMP, HALO, INK, INK_2, MASK_COLORS, MASK_NAMES, MEDIAN_LINE, ORANGE, ORANGE_RAMP, SURFACE, apply_style, combo_label, plt, save
from lineage_charisma.sanity import (
    binned_median,
    family_table,
    knn_label_agreement,
    knn_overlap,
    mean_knn_distance,
    near_neighbor_counts,
    nearest_neighbors,
    off_diagonal,
    pooled_within_between,
    spearman,
)


# every figure this script writes starts with one of these; only they are cleared before a rebuild
FIGURE_PREFIXES = ("heatmaps", "text_vs_phylo", "family_within_between", "length_vs_neighbors", "family_ratio_by_mask")


def combo_slug(model: str, rule: str, mask: str | None = None) -> str:
    return f"{model_slug(model)}__{rule}" + (f"__{mask}" if mask else "")


def family_blocks(families: list[str]) -> list[tuple[str, int, int]]:
    blocks, start = [], 0
    for i in range(1, len(families) + 1):
        if i == len(families) or families[i] != families[start]:
            blocks.append((families[start], start, i))
            start = i
    return blocks


def plot_heatmaps(text, phylo, tip_families, label, path, min_label=8) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.9), constrained_layout=True)
    blocks = family_blocks(tip_families)
    off = off_diagonal(text)
    panels = [
        (axes[0], text, BLUE_RAMP, f"Text distance ({label})", "1 − cosine similarity", dict(vmin=np.quantile(off, 0.01), vmax=np.quantile(off, 0.99)), "both"),
        (axes[1], phylo, ORANGE_RAMP, "Phylogenetic distance", "patristic distance, million years", dict(vmin=0, vmax=phylo.max()), "neither"),
    ]
    for ax, matrix, cmap, title, unit, limits, extend in panels:
        image = ax.imshow(matrix, cmap=cmap, interpolation="nearest", **limits)
        for _, start, _ in blocks[1:]:
            ax.axhline(start - 0.5, color=SURFACE, linewidth=0.7)
            ax.axvline(start - 0.5, color=SURFACE, linewidth=0.7)
        big = [(name, (start + end - 1) / 2) for name, start, end in blocks if end - start >= min_label]
        ax.set_yticks([mid for _, mid in big], [name for name, _ in big])
        ax.set_xticks([mid for _, mid in big], [name for name, _ in big], rotation=45, ha="right")
        ax.tick_params(length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_title(title)
        bar = fig.colorbar(image, ax=ax, shrink=0.72, pad=0.02, extend=extend)
        bar.set_label(unit, color=INK_2)
        bar.outline.set_visible(False)
    fig.suptitle(f"Both matrices in tree tip order ({len(tip_families)} species; families with at least {min_label} species labelled)",
                 x=0.01, ha="left", fontsize=9, color=INK_2)
    save(fig, path)


def plot_text_vs_phylo(text, phylo, label, path) -> pd.DataFrame:
    x, y = off_diagonal(phylo), off_diagonal(text)
    medians = binned_median(x, y, np.arange(0, x.max() + 5, 5), min_count=30)
    fig, ax = plt.subplots(figsize=(7.2, 4.8), constrained_layout=True)
    hexes = ax.hexbin(x, y, gridsize=55, cmap=BLUE_RAMP, norm=LogNorm(), mincnt=1, linewidths=0)
    ax.plot(medians["x"], medians["y"], label="median per 5 My bin", **MEDIAN_LINE)
    fig.legend(loc="outside lower center", handlelength=2.2)
    ax.set_xlabel("Phylogenetic distance (patristic, million years)")
    ax.set_ylabel("Text distance (1 − cosine similarity)")
    ax.set_title(f"Text distance against phylogenetic distance, {len(x):,} species pairs ({label})")
    ax.grid(True, axis="y")
    ax.set_axisbelow(True)
    bar = fig.colorbar(hexes, ax=ax, pad=0.02)
    bar.set_label("species pairs per hexagon (log scale)", color=INK_2)
    bar.outline.set_visible(False)
    save(fig, path)
    return medians


def plot_family_gaps(table: pd.DataFrame, label, path) -> None:
    table = table.iloc[::-1].reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(7.2, 0.36 * len(table) + 1.5), constrained_layout=True)
    rows = np.arange(len(table))
    ax.hlines(rows, table["within"], table["between"], color=AXIS, linewidth=2, zorder=1)
    ax.scatter(table["within"], rows, s=64, color=BLUE, edgecolor=SURFACE, linewidth=1.5, zorder=3, label="within family")
    ax.scatter(table["between"], rows, s=64, color=ORANGE, edgecolor=SURFACE, linewidth=1.5, zorder=3, label="to other families")
    ax.set_yticks(rows, [f"{f}  (n={n})" for f, n in zip(table["family"], table["n"])])
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Mean text distance (1 − cosine similarity)")
    ax.set_title(f"Mean text distance within each family and to other families ({label})")
    ax.grid(True, axis="x")
    ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False)
    fig.legend(loc="outside lower center", ncol=2, handletextpad=0.3, columnspacing=1.6)
    save(fig, path)


def plot_mask_levels(tables: dict[str, pd.DataFrame], label, path) -> None:
    levels = list(tables)
    base = tables[levels[0]].iloc[::-1].reset_index(drop=True)
    ratios = {level: tables[level].set_index("family").loc[base["family"], "ratio"].to_numpy() for level in levels}
    stacked = np.vstack(list(ratios.values()))
    fig, ax = plt.subplots(figsize=(7.2, 0.36 * len(base) + 1.6), constrained_layout=True)
    rows = np.arange(len(base))
    ax.hlines(rows, stacked.min(axis=0), stacked.max(axis=0), color=AXIS, linewidth=2, zorder=1)
    for level in levels:
        ax.scatter(ratios[level], rows, s=64, color=MASK_COLORS[level], edgecolor=SURFACE, linewidth=1.5, zorder=3, label=MASK_NAMES[level])
    if stacked.max() > 0.97:
        ax.axvline(1.0, color=AXIS, linewidth=0.8, zorder=0)
        ax.annotate("1.0: family no tighter than its surroundings", (1.0, 1.0), xycoords=("data", "axes fraction"), xytext=(-4, -3),
                    textcoords="offset points", va="top", ha="right", color=INK_2, fontsize=8, path_effects=HALO)
    ax.set_yticks(rows, [f"{f}  (n={n})" for f, n in zip(base["family"], base["n"])])
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Within-family mean text distance ÷ mean distance to other families")
    ax.set_title(f"Family tightness in text space at each mask level ({label})")
    ax.grid(True, axis="x")
    ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False)
    fig.legend(loc="outside lower center", ncol=len(levels), handletextpad=0.3, columnspacing=1.6)
    save(fig, path)


def plot_length_effect(n_tokens, counts, knn_mean, radius, k, limit, note, names, highlight, label, path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.5), constrained_layout=True)
    edges = np.quantile(n_tokens, np.linspace(0, 1, 6))
    edges[-1] += 1
    panels = [
        (axes[0], counts, f"Species closer than {radius:.3f} in text distance", f"Number of near neighbours ({label})"),
        (axes[1], knn_mean, f"Mean text distance to the {k} nearest species", f"Distance to the {k} nearest neighbours ({label})"),
    ]
    for ax, values, ylabel, title in panels:
        ax.scatter(n_tokens, values, s=26, color=BLUE, alpha=0.6, edgecolor=SURFACE, linewidth=0.6, zorder=2)
        medians = binned_median(n_tokens.astype(float), values.astype(float), edges)
        ax.plot(medians["x"], medians["y"], zorder=3, label="median per length quintile" if ax is axes[0] else None, **MEDIAN_LINE)
        ax.axvline(limit, color=AXIS, linewidth=0.8, zorder=1)
        ax.annotate(f"{limit} tokens: {note}", (limit, 1.0), xycoords=("data", "axes fraction"), xytext=(4, -3), textcoords="offset points",
                    va="top", color=INK_2, fontsize=8, path_effects=HALO)
        for i in highlight:
            ax.annotate(names[i], (n_tokens[i], values[i]), xytext=(5, 4), textcoords="offset points", color=INK, fontsize=8, path_effects=HALO, zorder=4)
        ax.set_xscale("log")
        ax.set_xlabel("Description length (tokens, log scale)")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(True, axis="y")
        ax.set_axisbelow(True)
    fig.legend(loc="outside lower center", handlelength=2.2)
    save(fig, path)


def md_table(frame: pd.DataFrame, decimals: int = 3) -> str:
    cols = list(frame.columns)
    numeric = [pd.api.types.is_numeric_dtype(frame[c]) for c in cols]
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join("---:" if n else "---" for n in numeric) + "|"]
    for row in frame.itertuples(index=False):
        lines.append("| " + " | ".join(f"{v:.{decimals}f}" if isinstance(v, float) else str(v) for v in row) + " |")
    return "\n".join(lines)


def main() -> None:
    args = base_parser("Week 3 sanity report: nearest neighbours, family structure, mask levels, and figures.").parse_args()
    cfg = load_config(args.config)
    apply_style()
    ecfg, scfg = cfg["embedding"], cfg["week3_sanity"]
    k, quantile = scfg["k_neighbors"], scfg["near_neighbor_quantile"]
    combos = cfg.combos()
    primary = combos[0]
    primary_mask = primary[2]
    model_rules = [(m, r) for m in ecfg["models"] for r in ecfg["rules"]]
    masks = [level for level in MASK_LEVELS if level in cfg.mask_levels]
    fig_dir, report_path = cfg.figures_dir(), cfg.reports / "week3_sanity.md"
    for prefix in FIGURE_PREFIXES:
        for stale in fig_dir.glob(f"{prefix}__*.png"):
            stale.unlink()

    order_path = cfg.species_order_path()
    phylo, order = load_matrix(phylo_matrix_path(cfg.matrices_dir()), order_path)
    text = {}
    for combo in combos:
        text[combo], combo_order = load_matrix(text_matrix_path(cfg.matrices_dir(), *combo), order_path)
        assert_same_order(combo_order, order, combo_label(*combo))

    full_corpus = read_corpus(cfg.corpus_path())
    corpus = full_corpus.set_index("species_id")
    matched = read_csv(cfg.matched_path()).set_index("species")
    wiki = read_csv(cfg.interim / "wiki_texts.csv").set_index("species")
    in_matrices = corpus.loc[order]
    families = matched.loc[order, "family"].tolist()
    genera = matched.loc[order, "genus"].tolist()
    common = matched.loc[order, "wiki_title"].tolist()
    family_of, common_of = dict(zip(order, families)), dict(zip(order, common))

    tree = load_tree(cfg.pruned_tree_path())
    tip_order = tree_tip_order(tree, order)
    assert_same_order(sorted(tip_order), order, "tree tip order")
    tip_families = [family_of[s] for s in tip_order]

    def neighbor_cell(focal: str, other: str, value: float, digits: int) -> str:
        tag = "" if family_of[other] == family_of[focal] else f" [{family_of[other]}]"
        return f"{common_of[other]}{tag} ({value:.{digits}f})"

    showcase = [s for s in scfg["showcase_species"] if s in order]
    summary_rows, family_tables, length_rows, passage_rows, medians_by_combo, figures = [], {}, [], [], {}, {}
    for combo in combos:
        model, rule, mask = combo
        matrix, label = text[combo], combo_label(*combo)
        table = family_table(matrix, families, scfg["min_family_size"])
        family_tables[combo] = table
        within, between = pooled_within_between(matrix, families)
        summary_rows.append({
            "model / rule": combo_label(model, rule),
            "mask level": mask,
            "mean distance": float(off_diagonal(matrix).mean()),
            "within family": within,
            "between families": between,
            "ratio": within / between,
            "families with within < between": f"{int((table['within'] < table['between']).sum())} of {len(table)}",
            "nearest neighbour in same family": f"{knn_label_agreement(matrix, families, 1):.0%}",
            f"top {k} in same family": f"{knn_label_agreement(matrix, families, k):.0%}",
            "nearest neighbour in same genus": f"{knn_label_agreement(matrix, genera, 1):.0%}",
        })
        if mask != primary_mask:
            continue

        meta = read_json(embedding_paths(cfg.embeddings_dir(), *combo).meta)
        limit = meta["max_seq_length"] - 2
        n_tokens = np.array([meta["n_tokens"][s] for s in order])
        highlight = [order.index(s) for s in showcase[:1] + showcase[3:4] + showcase[8:9]] + [int(np.argmin(n_tokens))]
        tercile_of = np.digitize(n_tokens, np.quantile(n_tokens, [1 / 3, 2 / 3]))
        radius, counts = near_neighbor_counts(matrix, quantile)
        knn_mean = mean_knn_distance(matrix, k)
        for t, name in enumerate(["shortest third", "middle third", "longest third"]):
            group = tercile_of == t
            length_rows.append({
                "model / rule": combo_label(model, rule),
                "length group": name,
                "tokens": f"{int(n_tokens[group].min())}–{int(n_tokens[group].max())}",
                "species": int(group.sum()),
                "median near neighbours": f"{np.median(counts[group]):.1f}",
                f"median distance to {k} nearest": float(np.median(knn_mean[group])),
            })
        embeddings, emb_order = load_embeddings(embedding_paths(cfg.embeddings_dir(), *combo))
        embeddings = reindex_rows(embeddings, emb_order, order, label).astype(np.float64)
        centroid = embeddings.mean(axis=0)
        to_centroid = 1.0 - embeddings @ (centroid / np.linalg.norm(centroid))
        long = n_tokens > limit
        both_long, both_short = np.outer(long, long), np.outer(~long, ~long)
        upper = np.triu(np.ones_like(both_long), k=1)
        passage_rows.append({
            "model / rule": combo_label(model, rule),
            "passage limit (tokens)": int(limit),
            "species over the limit": int(long.sum()),
            "both over": float(matrix[both_long & upper].mean()),
            "one over": float(matrix[~both_long & ~both_short & upper].mean()),
            "neither over": float(matrix[both_short & upper].mean()),
            "to centroid, over": float(to_centroid[long].mean()),
            "to centroid, not over": float(to_centroid[~long].mean()),
        })
        slug = combo_slug(*combo)
        figures[combo] = {
            "heatmaps": fig_dir / f"heatmaps__{slug}.png",
            "scatter": fig_dir / f"text_vs_phylo__{slug}.png",
            "family": fig_dir / f"family_within_between__{slug}.png",
            "length": fig_dir / f"length_vs_neighbors__{slug}.png",
        }
        plot_heatmaps(reindex_matrix(matrix, order, tip_order), reindex_matrix(phylo, order, tip_order), tip_families, label, figures[combo]["heatmaps"])
        medians_by_combo[combo] = plot_text_vs_phylo(matrix, phylo, label, figures[combo]["scatter"])
        plot_family_gaps(table, label, figures[combo]["family"])
        note = "longer texts are split into chunks" if rule == "chunk" else "text beyond this is ignored"
        plot_length_effect(n_tokens, counts, knn_mean, radius, k, limit, note, common, highlight, label, figures[combo]["length"])
        print(f"[sanity] {label}: figures written")

    mask_figures = {}
    for model, rule in model_rules:
        mask_figures[(model, rule)] = fig_dir / f"family_ratio_by_mask__{combo_slug(model, rule)}.png"
        plot_mask_levels({mask: family_tables[(model, rule, mask)] for mask in masks}, combo_label(model, rule), mask_figures[(model, rule)])

    rel = lambda p: p.relative_to(report_path.parent).as_posix()
    stats = length_summary(full_corpus)
    stubs = corpus[corpus["is_stub"]].reset_index()
    excluded = [s for s in corpus.index if s not in set(order)]
    lead_only = int((~wiki.loc[corpus.index, "has_description"].astype(bool)).sum())
    out: list[str] = []
    out.append("# Week 3 sanity report: text and phylogenetic distance matrices")
    out.append(f"Generated {datetime.now(timezone.utc):%Y-%m-%d %H:%M} UTC by `scripts/08_sanity_report.py`. "
               f"Clade: {cfg.clade_rank} {cfg.clade_name}. Primary combination: {combo_label(*primary)}. "
               f"Everything here is descriptive; no Mantel or other test statistics are computed.")

    out.append("## 1. Corpus")
    out.append(md_table(pd.DataFrame([{
        "species": stats["n_species"], "min tokens": stats["min"], "lower quartile": stats["q1"], "median": stats["median"], "mean": stats["mean"],
        "upper quartile": stats["q3"], "max": stats["max"], "over 512 tokens": stats["over_512"], "over 384 tokens": stats["over_384"],
    }]), decimals=1))
    out.append(f"Token counts use the `{cfg['corpus']['tokenizer']}` tokenizer on the unmasked text, without special tokens. \"Over 512\" and \"over 384\" "
               f"count texts that do not fit one passage once the two special tokens are added; 384 is the limit of `all-mpnet-base-v2`. "
               f"{lead_only} species have no description section and contribute the article lead only.")
    out.append(f"**Stub species** (fewer than {cfg['corpus']['min_tokens']} tokens): {len(stubs)}. "
               + (f"Stubs are excluded from every matrix ({len(excluded)} dropped), leaving {len(order)} species. " if excluded else "They are flagged and are in every matrix. ")
               + "The full list is in `data/interim/stub_species.csv`.")
    if len(stubs):
        out.append(md_table(pd.DataFrame({
            "species": stubs["species_id"], "article": [matched.loc[s, "wiki_title"] for s in stubs["species_id"]],
            "family": stubs["family"], "tokens": stubs["n_tokens"], "text": stubs["text"],
        })))
    sensitivity = cfg.get("matrices", {}).get("sensitivity_min_tokens")
    short = corpus[(corpus["n_tokens"] < 100) & ~corpus["is_stub"]].sort_values("n_tokens").reset_index()
    if len(short):
        out.append(f"{len(short)} more species are under 100 tokens and stay in the primary matrices: "
                   + ", ".join(f"*{s}* ({n})" for s, n in zip(short["species_id"], short["n_tokens"])) + ". "
                   + (f"`scripts/07_build_matrices.py --sensitivity` drops everything under {sensitivity} tokens into a separate `matrices_min{sensitivity}/` directory." if sensitivity else ""))

    out.append("### Mask levels")
    out.append("- **none**: the cleaned text as fetched.\n"
               f"- **taxonomy**: scientific names replaced by `{cfg['masking']['token']}`: every MDD genus, species epithet, subfamily, family, suborder and order for the clade, "
               "MDD synonym names, abbreviated forms such as \"P. leo\", and derived terms such as \"felid\" and \"mustelid\".\n"
               "- **strict**: taxonomy masking plus the head noun of every MDD common name (\"cat\", \"bear\", \"mongoose\", \"seal\") and its plural.")
    genus_of = [s.split()[0] for s in in_matrices.index]
    epithet_of = [s.split()[1] for s in in_matrices.index]
    mask_rows = []
    for mask in masks:
        col = text_column(mask)
        n_masked = in_matrices[f"n_masked_{mask}"] if mask != "none" else pd.Series(0, index=in_matrices.index)
        mask_rows.append({
            "mask level": mask,
            "texts naming their own genus": leak_counts(in_matrices[col], genus_of),
            "own species epithet": leak_counts(in_matrices[col], epithet_of),
            "own family": leak_counts(in_matrices[col], in_matrices["family"]),
            "median replacements per text": float(n_masked.median()),
            "max replacements": int(n_masked.max()),
        })
    out.append(md_table(pd.DataFrame(mask_rows), decimals=1))
    terms = read_csv(cfg.mask_terms_path())
    top = lambda level, kinds, n: ", ".join(f"{s} ({c})" for s, c in terms[(terms["mask_level"] == level) & terms["kind"].isin(kinds)].groupby("surface")["count"].sum().sort_values(ascending=False).head(n).items())
    if "taxonomy" in masks:
        out.append(f"Most often replaced at the taxonomy level, apart from full binomials: {top('taxonomy', ['name', 'epithet', 'compound'], 15)}.")
    if "strict" in masks:
        out.append(f"Most often replaced common-name nouns at the strict level: {top('strict', ['head_noun'], 15)}.")
    out.append("Adjacent replacements collapse into one token, so a binomial becomes a single mask. Every replaced term and its count is in `data/interim/mask_terms.csv`.")

    out.append("## 2. Matrices")
    rows = [{"matrix": "phylo.npy", "what": "patristic distance, million years", "min off-diagonal": float(off_diagonal(phylo).min()),
             "mean": float(off_diagonal(phylo).mean()), "max": float(phylo.max())}]
    for combo in combos:
        off = off_diagonal(text[combo])
        rows.append({"matrix": text_matrix_path(cfg.matrices_dir(), *combo).name, "what": f"1 − cosine, {combo_label(*combo)}",
                     "min off-diagonal": float(off.min()), "mean": float(off.mean()), "max": float(off.max())})
    out.append(md_table(pd.DataFrame(rows)))
    out.append(f"All {len(rows)} matrices are {len(order)}×{len(order)}, share `data/processed/matrices/species_order.txt`, and passed the "
               f"validation in `save_matrix` (symmetric, zero diagonal, no NaN, non-negative). Patristic distance is the branch length between "
               f"two tips, which on this ultrametric tree is twice their divergence time.")

    out.append(f"## 3. Nearest neighbours, {combo_label(*primary)}")
    out.append(f"The {k} nearest species in text space and in phylogenetic space. A family in square brackets marks a neighbour from a different "
               f"family than the focal species. Phylogenetic ties are common (every species in a sister clade is equally far) and are broken alphabetically, "
               f"so each heading counts how many text neighbours are at least as close on the tree as the {k}th phylogenetic neighbour.")
    for species in showcase:
        tn = nearest_neighbors(text[primary], order, species, k)
        pn = nearest_neighbors(phylo, order, species, k)
        focal, kth = order.index(species), pn[-1][1]
        shared = sum(phylo[focal, order.index(s)] <= kth for s, _ in tn)
        out.append(f"**{common_of[species]}** (*{species}*, {family_of[species]}, {int(corpus.loc[species, 'n_tokens'])} tokens; "
                   f"{shared} of {k} text neighbours are within {kth:.1f} My)")
        out.append(md_table(pd.DataFrame({
            "rank": [str(i + 1) for i in range(k)],
            "text neighbour (distance)": [neighbor_cell(species, s, d, 3) for s, d in tn],
            "phylogenetic neighbour (My)": [neighbor_cell(species, s, d, 1) for s, d in pn],
        })))

    out.append("## 4. Family-level structure")
    out.append(f"\"Within family\" is the mean over all same-family pairs and \"between families\" the mean over all other pairs. The first "
               f"{len(model_rules)} rows are the primary mask level ({primary_mask}).")
    out.append(md_table(pd.DataFrame(summary_rows)))
    out.append(f"Per family for {combo_label(*primary)}; \"between\" is the mean distance from the family's members to every species outside it.")
    out.append(md_table(family_tables[primary]))
    out.append(f"![Mean text distance within each family and to other families]({rel(figures[primary]['family'])})")

    out.append("## 5. Mask levels compared on the family check")
    summary = pd.DataFrame(summary_rows)
    wide_rows = []
    for model, rule in model_rules:
        row = {"model / rule": combo_label(model, rule)}
        for mask in masks:
            rec = summary[(summary["model / rule"] == combo_label(model, rule)) & (summary["mask level"] == mask)].iloc[0]
            row[f"ratio, {mask}"] = float(rec["ratio"])
        for mask in masks:
            rec = summary[(summary["model / rule"] == combo_label(model, rule)) & (summary["mask level"] == mask)].iloc[0]
            row[f"nearest in family, {mask}"] = rec["nearest neighbour in same family"]
        wide_rows.append(row)
    out.append("Pooled within-to-between ratio (lower means families are tighter) and the share of species whose nearest text neighbour is in their own family:")
    out.append(md_table(pd.DataFrame(wide_rows)))
    out.append(f"Per-family within-to-between ratio at each mask level, {combo_label(*primary[:2])}:")
    per_family = family_tables[primary][["family", "n"]].copy()
    for mask in masks:
        per_family[mask] = family_tables[(*primary[:2], mask)].set_index("family").loc[per_family["family"], "ratio"].to_numpy()
    if "none" in masks and "strict" in masks:
        per_family["strict − none"] = per_family["strict"] - per_family["none"]
    out.append(md_table(per_family))
    out.append(f"![Family tightness at each mask level]({rel(mask_figures[primary[:2]])})")
    out.append("The same figure for the other models and rules: "
               + ", ".join(f"[{combo_label(*mr)}]({rel(path)})" for mr, path in mask_figures.items() if mr != primary[:2]) + ".")
    agree = []
    for model, rule in model_rules:
        for i, a in enumerate(masks):
            for b in masks[i + 1:]:
                agree.append({
                    "model / rule": combo_label(model, rule),
                    "mask levels": f"{a} vs {b}",
                    "rank correlation of pairwise distances": spearman(off_diagonal(text[(model, rule, a)]), off_diagonal(text[(model, rule, b)])),
                    f"shared species in top {k} (mean of {k})": knn_overlap(text[(model, rule, a)], text[(model, rule, b)], k),
                })
    out.append("How much the text matrices change between mask levels (text matrices compared with each other, not with the tree):")
    out.append(md_table(pd.DataFrame(agree), decimals=2))
    out.append(f"Top {k} text neighbours at each mask level, {combo_label(*primary[:2])}:")
    for species in showcase:
        out.append(f"**{common_of[species]}** (*{species}*)")
        out.append(md_table(pd.DataFrame({
            "rank": [str(i + 1) for i in range(k)],
            **{mask: [neighbor_cell(species, s, d, 3) for s, d in nearest_neighbors(text[(*primary[:2], mask)], order, species, k)] for mask in masks},
        })))

    out.append(f"## 6. Figures, {combo_label(*primary)}")
    out.append(f"![Text and phylogenetic distance heatmaps in tree tip order]({rel(figures[primary]['heatmaps'])})")
    out.append(f"![Text distance against phylogenetic distance]({rel(figures[primary]['scatter'])})")
    med = medians_by_combo[primary]
    out.append("Median text distance by phylogenetic-distance bin (the line in the figure above):")
    out.append(md_table(pd.DataFrame({
        "phylogenetic distance (My)": [f"{lo:.0f}–{hi:.0f}" for lo, hi in zip(med["lo"], med["hi"])],
        "species pairs": med["n"], "median text distance": med["y"],
    })))
    out.append(f"![Description length against near-neighbour count]({rel(figures[primary]['length'])})")
    out.append(f"A \"near neighbour\" is a species closer than the {quantile:.0%} quantile of all pairwise text distances for that model and rule. "
               f"Lengths are counted on the {primary_mask}-masked text with each model's own tokenizer, and length groups are thirds of the species by that count:")
    out.append(md_table(pd.DataFrame(length_rows)))
    out.append("Mean text distance by whether each text exceeds one passage (and is therefore chunked or truncated), and mean distance from "
               "each group's embeddings to the corpus centroid:")
    out.append(md_table(pd.DataFrame(passage_rows)))
    out.append(f"The same four figures exist for every model and rule at the {primary_mask} mask level:")
    out.append(md_table(pd.DataFrame([{
        "model / rule": combo_label(*combo[:2]), **{name: f"[{name}]({rel(path)})" for name, path in figures[combo].items()},
    } for combo in figures])))

    out.append(f"## 7. Second model and truncate rule, {primary_mask} mask level")
    at_primary = [(m, r, primary_mask) for m, r in model_rules]
    agree = []
    for i, a in enumerate(at_primary):
        for b in at_primary[i + 1:]:
            agree.append({
                "pair": f"{combo_label(*a[:2])} vs {combo_label(*b[:2])}",
                "rank correlation of pairwise distances": spearman(off_diagonal(text[a]), off_diagonal(text[b])),
                f"shared species in top {k} (mean of {k})": knn_overlap(text[a], text[b], k),
            })
    out.append("Agreement between text matrices (these compare text matrices with each other, not with the tree):")
    out.append(md_table(pd.DataFrame(agree), decimals=2))
    out.append(f"Top {k} text neighbours for each showcase species under every model and rule:")
    for species in showcase:
        out.append(f"**{common_of[species]}** (*{species}*)")
        out.append(md_table(pd.DataFrame({
            "rank": [str(i + 1) for i in range(k)],
            **{combo_label(*combo[:2]): [neighbor_cell(species, s, d, 3) for s, d in nearest_neighbors(text[combo], order, species, k)] for combo in at_primary},
        })))
    out.append("Per-family within-to-between ratio under every model and rule (below 1 means the family is tighter than its surroundings):")
    ratios = family_tables[primary][["family", "n"]].copy()
    for combo in at_primary:
        ratios[combo_label(*combo[:2])] = family_tables[combo].set_index("family").loc[ratios["family"], "ratio"].to_numpy()
    out.append(md_table(ratios, decimals=2))

    notes = cfg.reports / "week3_notes.md"
    if notes.exists():
        out.append(notes.read_text(encoding="utf-8").strip())
    atomic_write_text(report_path, "\n\n".join(out) + "\n")
    print(f"[sanity] report -> {report_path.relative_to(cfg.root)}")


if __name__ == "__main__":
    main()
