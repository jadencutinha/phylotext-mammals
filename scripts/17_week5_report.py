from __future__ import annotations

import os
from datetime import datetime, timezone

import pandas as pd

from lineage_charisma.config import base_parser, load_config
from lineage_charisma.h2 import CONFOUNDED_SHARE
from lineage_charisma.io_utils import atomic_write_text, read_csv, read_meta
from lineage_charisma.plotstyle import combo_label

MODEL_NAMES = {"S1": "ecology", "S2": "length", "S3": "pageviews", "primary": "ecology + length (primary)", "S4": "ecology + length + pageviews"}
PREDICTOR_NAMES = {"phylogeny": "phylogenetic distance", "ecology": "ecological distance", "length": "article-length difference", "pageviews": "pageview difference"}
TRANSFORM_NAMES = {"none": "untransformed", "log": "log(1 + distance)", "capped": "capped at 30 My"}


def md_table(df: pd.DataFrame, formats: dict[str, str] | None = None) -> str:
    formats = formats or {}
    numeric = [pd.api.types.is_numeric_dtype(df[c]) and not pd.api.types.is_bool_dtype(df[c]) for c in df.columns]
    lines = ["| " + " | ".join(str(c) for c in df.columns) + " |", "|" + "|".join("---:" if n else "---" for n in numeric) + "|"]
    for row in df.itertuples(index=False):
        cells = [format(v, formats[c]) if c in formats and pd.notna(v) else ("n/a" if pd.isna(v) else str(v)) for c, v in zip(df.columns, row)]
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def yes_no(value) -> str:
    return "not defined" if pd.isna(value) else ("yes" if bool(value) else "no")


def summary(values: pd.Series) -> dict[str, float]:
    q = values.quantile([0.025, 0.5, 0.975]).to_numpy()
    return {"median": q[1], "2.5th pct": q[0], "97.5th pct": q[2], "min": values.min(), "max": values.max()}


def main() -> None:
    args = base_parser("Week 5 report: H2 results, robustness, exploratory analyses.").parse_args()
    cfg = load_config(args.config)
    tables, fig_dir, report_path = cfg.tables_dir(), cfg.week5_figures_dir(), cfg.reports / "week5_report.md"
    scfg = cfg["stats"]
    perms, alpha, method = int(scfg["permutations"]), float(scfg["alpha"]), scfg["method"]
    primary = cfg.primary_combo
    needed = [tables / f for f in ("h2_partial_mantel.csv", "h2_mrm.csv", "exploratory_well_known.csv", "exploratory_nonlinear.csv")]
    absent = [p.name for p in needed if not p.exists()]
    if absent:
        raise SystemExit(f"missing {absent}; run make h2 and make exploratory first")

    partial, regression = read_csv(tables / "h2_partial_mantel.csv"), read_csv(tables / "h2_mrm.csv")
    for frame in (partial, regression):
        frame["model"] = frame["model"].str.split("/").str[-1]
    head = partial[partial["is_primary"]].iloc[0]
    share_fmt = {"raw r": ".4f", "raw p": ".4f", "partial r": ".4f", "p": ".4f", "share of raw r": ".0%", "raw − partial": "+.4f"}
    columns = {"raw_r": "raw r", "raw_p": "raw p", "partial_r": "partial r", "partial_p": "p", "share_of_raw": "share of raw r", "change": "raw − partial"}

    def rel(path) -> str:
        return os.path.relpath(path, report_path.parent)

    def figure(name: str, alt: str) -> list[str]:
        return [f"![{alt}]({rel(fig_dir / name)})", ""] if (fig_dir / name).exists() else []

    out = ["# Week 5 report: H2, partial Mantel tests and MRM", "",
           f"Generated {datetime.now(timezone.utc):%Y-%m-%d %H:%M} UTC by `scripts/17_week5_report.py`. Clade: {cfg.clade_rank} {cfg.clade_name}. "
           f"Primary specification: {combo_label(*primary)}, {method.capitalize()}. The H2 plan is `docs/preregistration_amendment_h2.md`, committed before any of these "
           f"numbers existed. All tests use {perms:,} permutations and seed {scfg['seed']}. Partial Mantel p-values are one-sided; MRM p-values are two-sided.", ""]

    # 1. headline
    out += ["## 1. H2 on the primary specification", "",
            f"Species: the {int(head['n_species'])} with ecology data ({int(head['n_pairs']):,} pairs). The raw correlation is recomputed on these species and is the baseline. "
            "The primary model controls for ecological distance (trait groups weighted equally) and the difference in log article length together.", "",
            md_table(pd.DataFrame([{"quantity": "raw r (no control)", "value": f"{head['raw_r']:.4f}", "one-sided p": f"{head['raw_p']:.4f}"},
                                   {"quantity": "partial r (ecology + length)", "value": f"{head['partial_r']:.4f}", "one-sided p": f"{head['partial_p']:.4f}"},
                                   {"quantity": "raw r − partial r", "value": f"{head['change']:.4f}", "one-sided p": ""},
                                   {"quantity": "partial r as a share of raw r",
                                    "value": "not defined" if pd.isna(head["share_of_raw"]) else f"{head['share_of_raw']:.1%}", "one-sided p": ""}])), "",
            "**Decision rules (fixed in the amendment, section 6):**", "",
            f"- (a) Lineage signal survives (partial r > 0 with one-sided p < {alpha}): **{yes_no(head['survives'])}**.",
            f"- (b) Signal is largely confounded (partial r under {CONFOUNDED_SHARE:.0%} of raw r): **{yes_no(head['confounded'])}**.", ""]
    out += figure("h2_controls_primary.png", "Raw and partial r for the primary specification")

    # 2. each control
    steps = partial[partial["analysis"].isin(["primary", "secondary"])].set_index("model_id").loc[["S1", "S2", "S3", "primary", "S4"]].reset_index()
    steps.insert(0, "controls for", steps["model_id"].map(MODEL_NAMES))
    out += ["## 2. How much of the raw r remains after each control", "",
            f"Raw r on these species: {head['raw_r']:.4f}. Each row removes the listed controls from both text distance and phylogenetic distance.", "",
            md_table(steps[["controls for", *columns]].rename(columns=columns).drop(columns=["raw r", "raw p"]), share_fmt), ""]
    for row in steps.itertuples(index=False):
        left = "not defined (raw r is not positive)" if pd.isna(row.share_of_raw) else f"{row.share_of_raw:.0%} of the raw r"
        out.append(f"- Controlling for {MODEL_NAMES[row.model_id].replace(' (primary)', '')}: r goes from {row.raw_r:.3f} to {row.partial_r:.3f}, which is {left}.")
    out.append("")

    # 3. MRM
    fit = regression[regression["is_primary"]]
    mrm_view = fit.assign(predictor=fit["predictor"].map(PREDICTOR_NAMES))[["predictor", "coef_std", "t", "p"]].rename(
        columns={"coef_std": "standardized coefficient", "t": "pseudo-t"})
    out += ["## 3. MRM", "",
            f"Text distance regressed on four predictors, all as ranked distances, {int(fit['n_species'].iloc[0])} species. p-values come from permuting the text matrix. "
            f"R² = {fit['r2'].iloc[0]:.4f} (permutation p = {fit['p_r2'].iloc[0]:.4f}).", "",
            md_table(mrm_view, {"standardized coefficient": ".4f", "pseudo-t": ".2f", "p": ".4f"}), "",
            "Coefficients on the unstandardized ranks are in `h2_mrm.csv` (`coef_raw`); every ranked matrix has nearly the same spread, so they match the standardized ones to four decimals.", ""]
    out += figure("h2_mrm_primary.png", "MRM standardized coefficients")

    # 4. robustness
    variants = partial[partial["analysis"].isin(["primary", "variant"]) & (partial["model_id"] == "primary")].copy()
    variants["rule (a)"], variants["rule (b)"] = variants["survives"].map(yes_no), variants["confounded"].map(yes_no)
    n_a, n_b = int(variants["survives"].sum()), int((variants["confounded"] == True).sum())  # noqa: E712
    out += ["## 4. Robustness", "", "Reported in full; none of this decides H2. No correction for multiple testing.", "",
            "### All 12 model, rule and mask variants", "",
            f"Primary H2 model (ecology + length) on each text matrix. Rule (a) holds in {n_a} of {len(variants)} and rule (b) in {n_b} of {len(variants)}.", "",
            md_table(variants[["model", "rule", "mask", *columns, "rule (a)", "rule (b)"]].rename(columns=columns), share_fmt), ""]
    out += figure("h2_variants.png", "Primary partial r across the 12 variants")

    keys = ["model", "rule", "mask"]
    fits = regression[regression["analysis"].isin(["primary", "variant"])]
    wide = variants[keys].copy()
    for name in PREDICTOR_NAMES:
        wide = wide.merge(fits.loc[fits["predictor"] == name, [*keys, "coef_std"]].rename(columns={"coef_std": f"β {name}"}), on=keys, validate="one_to_one")
    wide = wide.merge(fits.loc[fits["predictor"] == "phylogeny", [*keys, "p", "r2"]].rename(columns={"p": "p phylogeny", "r2": "R²"}), on=keys, validate="one_to_one")
    out += ["MRM standardized coefficients for the same 12 text matrices:", "",
            md_table(wide, {**{f"β {k}": ".3f" for k in PREDICTOR_NAMES}, "p phylogeny": ".4f", "R²": ".3f"}), ""]

    others = partial[partial["analysis"].isin(["primary", "unweighted_gower", "complete_traits"])].copy()
    others.insert(0, "analysis ", others["analysis"].map({"primary": "primary (for comparison)", "unweighted_gower": "unweighted Gower", "complete_traits": "complete-trait species"}))
    others["rule (a)"], others["rule (b)"] = others["survives"].map(yes_no), others["confounded"].map(yes_no)
    out += ["### Gower weighting and the complete-trait subset", "",
            md_table(others[["analysis ", "n_species", *columns, "rule (a)", "rule (b)"]].rename(columns={**columns, "n_species": "species"}), share_fmt), ""]

    posterior_path = tables / "h2_posterior.csv"
    if posterior_path.exists():
        post = read_csv(posterior_path)
        defined = post["share_of_raw"].notna()
        dist = pd.DataFrame([{"quantity": "raw r", **summary(post["raw_r"])}, {"quantity": "partial r", **summary(post["partial_r"])},
                             {"quantity": "raw r − partial r", **summary(post["change"])},
                             {"quantity": "partial r / raw r", **summary(post.loc[defined, "share_of_raw"])}])
        under = int((post.loc[defined, "share_of_raw"] < CONFOUNDED_SHARE).sum())
        out += [f"### {len(post)} posterior trees", "",
                f"The primary H2 model on the {len(post)} trees drawn in amendment A1.3, {int(post['n_species'].iloc[0])} species. MCC tree for comparison: "
                f"raw r = {head['raw_r']:.4f}, partial r = {head['partial_r']:.4f}.", "",
                md_table(dist, {c: ".4f" for c in dist.columns if c != "quantity"}), "",
                f"- Raw r is positive on {int((post['raw_r'] > 0).sum())} trees and has p < {alpha} on {int((post['raw_p'] < alpha).sum())}.",
                f"- Rule (a) (partial r > 0 with p < {alpha}) holds on **{int(post['survives'].sum())} of {len(post)}** trees.",
                f"- The partial-to-raw ratio is under {CONFOUNDED_SHARE:.0%} on **{under} of {len(post)}** trees, at or above it on {int(defined.sum()) - under}, "
                f"and not defined (raw r at or below zero) on {int((~defined).sum())}. The ratio's summary above is over the {int(defined.sum())} trees where it is defined.", ""]
        out += figure("h2_posterior.png", "Raw and partial r across posterior trees")
    else:
        out += ["### Posterior trees", "", "Not run (`make h2` was given `--skip-posterior`, or the archive is missing).", ""]

    # 5. exploratory
    well, nonlin = read_csv(tables / "exploratory_well_known.csv"), read_csv(tables / "exploratory_nonlinear.csv")
    cap_share = read_meta(tables / "exploratory_nonlinear.csv").get("share_of_pairs_at_or_below_cap")
    well_view = well[["quartile", "n_pairs", "mean_text_distance", "median_text_distance", "mean_phylo_distance", "r_text_phylo", "p"]].rename(columns={
        "n_pairs": "pairs", "mean_text_distance": "mean text distance", "median_text_distance": "median text distance",
        "mean_phylo_distance": "mean phylogenetic distance (My)", "r_text_phylo": "r text–phylogeny"})
    nonlin_view = nonlin.assign(transform=nonlin["transform"].map(TRANSFORM_NAMES))[["transform", "method", "h1_r", "h1_p", "raw_r", "partial_r", "partial_p", "share_of_raw"]].rename(
        columns={"h1_r": "H1 r (284)", "h1_p": "H1 p", "raw_r": "raw r (283)", "partial_r": "partial r", "partial_p": "partial p", "share_of_raw": "share of raw r"})
    out += ["## 5. Exploratory analyses", "",
            "**Everything in this section is exploratory.** It was planned in section 9 of the amendment but has no bearing on H1 or H2, and its p-values are descriptive.", "",
            "### Exploratory: pairs of well-known species", "",
            "All 284 species. Pairs are split into quartiles of the mean of the two species' log pageviews. The correlation is Spearman, over the pairs in the quartile, "
            "computed for the bottom and top quartiles as planned; its p-value permutes the phylogenetic matrix with the set of pairs held fixed.", "",
            md_table(well_view, {"mean text distance": ".4f", "median text distance": ".4f", "mean phylogenetic distance (My)": ".1f", "r text–phylogeny": ".4f", "p": ".4f"}).replace("n/a", ""), ""]
    out += figure("exploratory_well_known.png", "Exploratory: text distance by pageview quartile")
    out += ["### Exploratory: transformed phylogenetic distance", "",
            "H1 on 284 species and the primary H2 model on 283, with phylogenetic distance transformed. A log transform cannot change a Spearman correlation, so it is run "
            "with Pearson and compared with the untransformed Pearson row." + (f" {cap_share:.0%} of pairs are at or below the 30 My cap." if cap_share is not None else ""), "",
            md_table(nonlin_view, {"H1 r (284)": ".4f", "H1 p": ".4f", "raw r (283)": ".4f", "partial r": ".4f", "partial p": ".4f", "share of raw r": ".0%"}), ""]
    out += figure("exploratory_nonlinear.png", "Exploratory: transformed phylogenetic distance")

    bench = tables / "scaling_benchmark.csv"
    if bench.exists():
        b = read_csv(bench)
        out += ["## 6. Scaling measurements", "",
                f"Measured by `scripts/18_scaling_benchmark.py` on this machine ({read_meta(bench).get('machine', 'unknown')}). These feed the go/no-go section below.", "",
                md_table(b, {"value": ",.2f"}), ""]

    notes = cfg.reports / "week5_notes.md"
    if notes.exists():
        out += [notes.read_text().rstrip(), ""]
    atomic_write_text(report_path, "\n".join(out))
    print(f"[report] wrote {report_path.relative_to(cfg.root)}")


if __name__ == "__main__":
    main()
