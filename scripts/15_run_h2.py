from __future__ import annotations

import numpy as np
import pandas as pd

from lineage_charisma.config import base_parser, load_config
from lineage_charisma.h2 import MODELS, MRM_PREDICTORS, PRIMARY_CONTROLS, complete_trait_species, cut, decision, load_h2_matrices, load_text, mrm_rows, partial_result
from lineage_charisma.io_utils import write_csv
from lineage_charisma.plotstyle import BLUE, INK, MASK_COLORS, MUTED, ORANGE, SURFACE, apply_style, combo_label, plot_variants, plt, save
from lineage_charisma.posterior import N_TREES, posterior_distance_matrices
from lineage_charisma.stats import submatrix

CONTROL_NAMES = {"ecology": "ecology", "ecology_unweighted": "ecology (columns equal)", "length": "length", "pageviews": "pageviews"}
PREDICTOR_NAMES = {"phylogeny": "phylogenetic distance", "ecology": "ecological distance", "length": "article-length difference", "pageviews": "pageview difference"}


def controls_label(controls) -> str:
    return " + ".join(controls)


def p_label(p: float, permutations: int) -> str:
    return f"p ≤ {1 / (permutations + 1):.4f}" if p <= 1.0001 / (permutations + 1) else f"p = {p:.4f}"


def plot_controls(primary_rows: pd.DataFrame, permutations: int, label: str, path) -> None:
    by_model = primary_rows.set_index("model_id")
    head = by_model.loc["primary"]
    bars = [("no control (raw)", head["raw_r"], head["raw_p"], MUTED)]
    for model_id, name in (("S1", "ecology"), ("S2", "length"), ("S3", "pageviews"), ("primary", "ecology + length\n(primary H2 model)"),
                           ("S4", "ecology + length\n+ pageviews")):
        bars.append((name, by_model.loc[model_id, "partial_r"], by_model.loc[model_id, "partial_p"], BLUE))
    fig, ax = plt.subplots(figsize=(7.2, 4.0), constrained_layout=True)
    ys = np.arange(len(bars))[::-1]
    values = [b[1] for b in bars]
    ax.barh(ys, values, height=0.56, color=[b[3] for b in bars], edgecolor=SURFACE, linewidth=1)
    span = max(max(values), 0) - min(min(values), 0)
    for y, (_, r, p, _) in zip(ys, bars):
        ax.annotate(f"r = {r:.3f}, {p_label(p, permutations)}", xy=(max(r, 0), y), xytext=(5, 0), textcoords="offset points", va="center", color=INK)
    ax.axvline(0, color=MUTED, linewidth=0.8)
    ax.set_xlim(min(min(values), 0) - 0.05 * span, max(max(values), 0) + 0.55 * span)
    ax.set_yticks(ys, [b[0] for b in bars])
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Spearman Mantel r, text distance against phylogenetic distance")
    ax.set_title(f"Text–phylogeny correlation, raw and after each control\n{int(head['n_species'])} species, one-sided p ({label})")
    ax.grid(True, axis="x")
    ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False)
    save(fig, path)


def plot_mrm(rows: pd.DataFrame, permutations: int, label: str, path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 3.2), constrained_layout=True)
    ys = np.arange(len(rows))[::-1]
    coefs = rows["coef_std"].to_numpy()
    ax.barh(ys, coefs, height=0.56, color=BLUE, edgecolor=SURFACE, linewidth=1)
    span = max(coefs.max(), 0) - min(coefs.min(), 0)
    for y, coef, p in zip(ys, coefs, rows["p"]):
        ax.annotate(f"β = {coef:.3f}, {p_label(p, permutations)}", xy=(max(coef, 0), y), xytext=(5, 0), textcoords="offset points", va="center", color=INK)
    ax.axvline(0, color=MUTED, linewidth=0.8)
    ax.set_xlim(min(coefs.min(), 0) - 0.05 * span, max(coefs.max(), 0) + 0.5 * span)
    ax.set_yticks(ys, [PREDICTOR_NAMES[p] for p in rows["predictor"]])
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("standardized coefficient on ranked distances")
    r2, p_r2 = rows["r2"].iloc[0], rows["p_r2"].iloc[0]
    ax.set_title(f"MRM: text distance on four predictors, two-sided p ({label})\nR² = {r2:.4f}, {p_label(p_r2, permutations)}")
    ax.grid(True, axis="x")
    ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False)
    save(fig, path)


def plot_posterior(table: pd.DataFrame, head: pd.Series, label: str, path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.0), constrained_layout=True)
    lo, hi = min(table["raw_r"].min(), table["partial_r"].min()), max(table["raw_r"].max(), table["partial_r"].max())
    bins = np.linspace(lo, hi, 26)
    ax.hist(table["partial_r"], bins=bins, color=ORANGE, edgecolor=SURFACE, linewidth=0.4, label="partial r (controlling for ecology and length)")
    ax.hist(table["raw_r"], bins=bins, histtype="step", color=BLUE, linewidth=2, label="raw r")
    ax.axvline(head["partial_r"], color=ORANGE, linewidth=1.4, linestyle="--", label=f"MCC tree, partial r = {head['partial_r']:.3f}")
    ax.axvline(head["raw_r"], color=BLUE, linewidth=1.4, linestyle="--", label=f"MCC tree, raw r = {head['raw_r']:.3f}")
    ax.axvline(0, color=INK, linewidth=0.6, linestyle=":")
    ax.set_xlabel("Spearman Mantel r, text distance against phylogenetic distance")
    ax.set_ylabel("posterior trees")
    ax.set_title(f"Raw and primary partial r on {len(table)} posterior trees, {int(head['n_species'])} species ({label})")
    ax.grid(True, axis="y")
    ax.set_axisbelow(True)
    fig.legend(loc="outside lower center", ncol=2, handlelength=2.2, columnspacing=1.6)
    save(fig, path)


def main() -> None:
    parser = base_parser("H2: partial Mantel tests and MRM of text distance on phylogeny, ecology and attention.")
    parser.add_argument("--skip-posterior", action="store_true", help="skip the 100 posterior trees (they need the 1.2 GB posterior archive)")
    args = parser.parse_args()
    cfg = load_config(args.config)
    apply_style()
    scfg, ecfg = cfg["stats"], cfg["embedding"]
    perms, seed, method, alpha = int(scfg["permutations"]), int(scfg["seed"]), scfg["method"], float(scfg["alpha"])
    combos, primary = cfg.combos(), cfg.primary_combo
    matrices, species, order = load_h2_matrices(cfg)
    print(f"[h2] {len(species)} species with ecology data; {perms:,} permutations, seed {seed}, {method}")

    partial, regression = [], []

    def add(analysis, model_id, combo, text, mats, controls, ecology="ecology", with_mrm=True):
        swap = {**mats, "ecology": mats[ecology]}    # the unweighted variant stands in wherever the model says ecology
        res = partial_result(text, swap["phylogeny"], [swap[c] for c in controls], method, perms, seed)
        is_primary = analysis == "primary" and model_id == "primary"
        base = {"analysis": analysis, "model_id": model_id, "model": combo[0], "rule": combo[1], "mask": combo[2], "ecology_matrix": ecology,
                "tree": "MCC", "method": method}
        partial.append({**base, "controls": controls_label(controls), **res, **decision(res["raw_r"], res["partial_r"], res["partial_p"], alpha),
                        "is_primary": is_primary})
        print(f"[h2] {analysis:16s} {model_id:8s} {combo_label(*combo):45s} given {controls_label(controls):30s} raw r = {res['raw_r']:+.4f}  "
              f"partial r = {res['partial_r']:+.4f}  p = {res['partial_p']:.4f}{'   <- primary' if is_primary else ''}")
        if with_mrm:
            regression.extend({**base, **row, "is_primary": is_primary} for row in mrm_rows(text, swap, method, perms, seed))
        return res

    text = load_text(cfg, primary, species)
    head = add("primary", "primary", primary, text, matrices, MODELS["primary"])
    for model_id in ("S1", "S2", "S3", "S4"):
        add("secondary", model_id, primary, text, matrices, MODELS[model_id], with_mrm=False)
    for combo in combos:
        if combo != primary:
            add("variant", "primary", combo, load_text(cfg, combo, species), matrices, PRIMARY_CONTROLS)
    add("unweighted_gower", "primary", primary, text, matrices, PRIMARY_CONTROLS, ecology="ecology_unweighted")
    complete = complete_trait_species(cfg)
    add("complete_traits", "primary", primary, submatrix(text, species, complete), cut(matrices, species, complete), PRIMARY_CONTROLS)

    partial_table, mrm_table = pd.DataFrame(partial), pd.DataFrame(regression)
    meta = dict(permutations=perms, seed=seed, alpha=alpha, text_matrix=combo_label(*primary), plan="docs/preregistration_amendment_h2.md")
    write_csv(partial_table, cfg.tables_dir() / "h2_partial_mantel.csv", alternative="greater", permuted="phylogeny", **meta)
    write_csv(mrm_table, cfg.tables_dir() / "h2_mrm.csv", alternative="two-sided on pseudo-t", permuted="text", predictors=list(MRM_PREDICTORS), **meta)

    fig_dir, label = cfg.week5_figures_dir(), combo_label(*primary)
    on_primary = partial_table[partial_table["analysis"].isin(["primary", "secondary"])]
    plot_controls(on_primary, perms, label, fig_dir / "h2_controls_primary.png")
    plot_mrm(mrm_table[mrm_table["is_primary"]], perms, label, fig_dir / "h2_mrm_primary.png")
    variants = partial_table[partial_table["analysis"].isin(["primary", "variant"]) & (partial_table["model_id"] == "primary")]
    plot_variants(variants, "partial_r", [(m, r) for m in ecfg["models"] for r in ecfg["rules"]], [m for m in MASK_COLORS if m in cfg.mask_levels],
                  primary, "Spearman partial Mantel r, text against phylogeny given ecology and length", "Primary H2 model across all 12 variants",
                  fig_dir / "h2_variants.png")

    rules = decision(head["raw_r"], head["partial_r"], head["partial_p"], alpha)
    print(f"[h2] primary: raw r = {head['raw_r']:.4f} (p = {head['raw_p']:.4f}), partial r = {head['partial_r']:.4f} (p = {head['partial_p']:.4f}), "
          f"share of raw {head['share_of_raw']:.1%}, change {head['change']:+.4f}")
    print(f"[h2] rule (a) lineage signal survives: {rules['survives']};  rule (b) signal is largely confounded: {rules['confounded']}")

    if args.skip_posterior:
        print("[h2] skipped the posterior trees")
        return
    trees, drawn = posterior_distance_matrices(cfg, order, N_TREES, seed)
    controls = [matrices[c] for c in PRIMARY_CONTROLS]
    rows = []
    for tree_matrix, tree in zip(trees, drawn.itertuples(index=False)):
        res = partial_result(text, submatrix(tree_matrix, order, species), controls, method, perms, seed)
        rows.append({"draw": tree.draw, "tree_index": tree.tree_index, "tree_file": tree.tree_file, "method": method, "controls": controls_label(PRIMARY_CONTROLS),
                     **res, **decision(res["raw_r"], res["partial_r"], res["partial_p"], alpha)})
        if tree.draw % 10 == 0:
            print(f"[h2] posterior {tree.draw}/{N_TREES} trees")
    posterior = pd.DataFrame(rows)
    write_csv(posterior, cfg.tables_dir() / "h2_posterior.csv", alternative="greater", permuted="phylogeny", n_trees=N_TREES, selection_seed=seed, **meta)
    plot_posterior(posterior, pd.Series(head), label, fig_dir / "h2_posterior.png")
    defined = posterior["share_of_raw"].notna()
    print(f"[h2] posterior: partial r median {posterior['partial_r'].median():.4f}; rule (a) holds on {int(posterior['survives'].sum())} of {len(posterior)} trees; "
          f"ratio under 50% on {int((posterior.loc[defined, 'share_of_raw'] < 0.5).sum())}, undefined on {int((~defined).sum())}")


if __name__ == "__main__":
    main()
