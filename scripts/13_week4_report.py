from __future__ import annotations

import os
from datetime import datetime, timezone

import pandas as pd

from lineage_charisma.attention import month_range
from lineage_charisma.config import base_parser, load_config
from lineage_charisma.embed import read_species_order
from lineage_charisma.io_utils import atomic_write_text, read_csv, read_meta
from lineage_charisma.plotstyle import combo_label

MATRIX_NAMES = {
    "phylogeny": "phylogeny", "text": "text (primary)", "ecology": "ecology (groups equal)", "ecology_unweighted": "ecology (columns equal)",
    "attention_pageviews": "pageview distance", "attention_length": "length distance", "attention_combined": "combined attention",
}


def md_table(df: pd.DataFrame, formats: dict[str, str] | None = None) -> str:
    formats = formats or {}
    numeric = [pd.api.types.is_numeric_dtype(df[c]) for c in df.columns]
    lines = ["| " + " | ".join(str(c) for c in df.columns) + " |", "|" + "|".join("---:" if n else "---" for n in numeric) + "|"]
    for row in df.itertuples(index=False):
        cells = [format(v, formats[c]) if c in formats and pd.notna(v) else ("" if pd.isna(v) else str(v)) for c, v in zip(df.columns, row)]
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def p_text(p: float, perms: int) -> str:
    floor = 1 / (perms + 1)
    return f"{floor:.4f} (the smallest possible)" if p <= floor * 1.0001 else f"{p:.4f}"


def main() -> None:
    args = base_parser("Week 4 report: H1 results, robustness, trait and attention coverage, confound diagnostics.").parse_args()
    cfg = load_config(args.config)
    tables, fig_dir, report_path = cfg.tables_dir(), cfg.week4_figures_dir(), cfg.reports / "week4_report.md"
    scfg, acfg = cfg["stats"], cfg["attention"]
    perms, alpha, method = int(scfg["permutations"]), float(scfg["alpha"]), scfg["method"]
    primary = cfg.primary_combo
    needed = [tables / "h1_mantel.csv", tables / "h1_mantel_families.csv", tables / "ecology_missing.csv", tables / "confound_mantel.csv",
              cfg.processed / "attention.csv", cfg.processed / "traits.csv"]
    absent = [p.name for p in needed if not p.exists()]
    if absent:
        raise SystemExit(f"missing {absent}; run make h1, make ecology, make attention and make diagnostics first")

    order = read_species_order(cfg.species_order_path())
    h1, families = read_csv(tables / "h1_mantel.csv"), read_csv(tables / "h1_mantel_families.csv")
    h1["model"] = h1["model"].str.split("/").str[-1]
    head = h1[h1["is_primary"]].iloc[0]
    main_rows = h1[h1["method"] == method]

    def rel(path) -> str:
        return os.path.relpath(path, report_path.parent)

    out = ["# Week 4 report: H1, ecological and attention distances, confound diagnostics", "",
           f"Generated {datetime.now(timezone.utc):%Y-%m-%d %H:%M} UTC by `scripts/13_week4_report.py`. Clade: {cfg.clade_rank} {cfg.clade_name}, {len(order)} species. "
           f"Primary specification: {combo_label(*primary)}, {method.capitalize()}. The analysis plan is `docs/preregistration.md`, committed before any of these numbers existed. "
           f"All Mantel tests use {perms:,} permutations and seed {scfg['seed']}.", ""]

    out += ["## 1. H1: text distance against phylogenetic distance", "",
            f"**Primary specification: Mantel r = {head['r']:.3f}, one-sided p = {p_text(head['p'], perms)}, over {int(head['n_pairs']):,} species pairs.** "
            f"By the pre-registered rule (r > 0 and p < {alpha}), H1 is {'supported' if head['r'] > 0 and head['p'] < alpha else 'not supported'}.", "",
            f"![Text distance against phylogenetic distance, primary specification]({rel(fig_dir / 'h1_scatter_primary.png')})", "",
            f"![Permutation null distribution with the observed r]({rel(fig_dir / 'h1_null_primary.png')})", ""]

    out += ["## 2. Robustness across the 12 variants", ""]
    pearson = h1[h1["method"] != method].set_index(["model", "rule", "mask"])
    variant = main_rows[["model", "rule", "mask", "r", "p", "n_pairs", "is_primary"]].copy()
    variant[f"{pearson['method'].iloc[0]} r"] = pearson.loc[list(zip(variant["model"], variant["rule"], variant["mask"])), "r"].to_numpy()
    variant[f"{pearson['method'].iloc[0]} p"] = pearson.loc[list(zip(variant["model"], variant["rule"], variant["mask"])), "p"].to_numpy()
    variant["is_primary"] = variant["is_primary"].map({True: "yes", False: ""})
    variant = variant.rename(columns={"r": f"{method} r", "p": f"{method} p", "n_pairs": "pairs", "is_primary": "primary"})
    out += [md_table(variant, {c: ".3f" for c in variant.columns if c.endswith(" r")} | {c: ".4f" for c in variant.columns if c.endswith(" p")} | {"pairs": ","}), "",
            f"p-values are one-sided; {1 / (perms + 1):.4f} is the smallest value {perms:,} permutations can give. The full table is `{rel(tables / 'h1_mantel.csv')}`.", ""]
    n_pos, n_sig = int((main_rows["r"] > 0).sum()), int(((main_rows["r"] > 0) & (main_rows["p"] < alpha)).sum())
    n_pos_p, n_sig_p = int((pearson["r"] > 0).sum()), int(((pearson["r"] > 0) & (pearson["p"] < alpha)).sum())
    out += [f"**Direction and significance.** With {method.capitalize()}, r is positive in {n_pos} of {len(main_rows)} variants and significant at {alpha} in {n_sig}. "
            f"With {pearson['method'].iloc[0].capitalize()}, r is positive in {n_pos_p} of {len(pearson)} and significant in {n_sig_p}. No multiple-testing correction is applied.", ""]
    masks = [m for m in ("none", "taxonomy", "strict") if m in set(main_rows["mask"])]
    by_mask = main_rows.pivot(index=["model", "rule"], columns="mask", values="r")[masks].reset_index()
    by_mask.columns.name = None
    if {"none", "strict"} <= set(masks):
        by_mask["strict as share of none"] = by_mask["strict"] / by_mask["none"]
    out += [f"**{method.capitalize()} r by mask level:**", "", md_table(by_mask, {m: ".3f" for m in masks} | {"strict as share of none": ".0%"}), "",
            f"![Mantel r for all 12 variants]({rel(fig_dir / 'h1_variants.png')})", ""]

    fam = families[families["method"] == method][["family", "n_species", "n_pairs", "r", "p"]].rename(columns={"n_species": "species", "n_pairs": "pairs"})
    out += ["### Within the three largest families (descriptive)", "",
            "The primary-specification test repeated inside each family. These are descriptions, with no multiple-testing claims.", "",
            md_table(fam, {"r": ".3f", "p": ".4f", "pairs": ","}), ""]

    posterior_path = tables / "h1_posterior.csv"
    if posterior_path.exists():
        post = read_csv(posterior_path)
        q = post["r"].quantile([0.025, 0.5, 0.975]).to_numpy()
        out += ["### Tree uncertainty: posterior trees (robustness, amendment A1.3)", "",
                f"The primary-specification test repeated on {len(post)} trees drawn from the 10,000-tree posterior with the pre-registered seed, "
                "each pruned to the same 284 species. This is a robustness analysis; H1 is decided by the MCC tree above.", "",
                md_table(pd.DataFrame([{"trees": len(post), "median r": q[1], "2.5th percentile": q[0], "97.5th percentile": q[2], "minimum": post["r"].min(),
                                        "maximum": post["r"].max(), "MCC tree r": head["r"], "trees with r > 0": int((post["r"] > 0).sum()),
                                        f"trees with p < {alpha}": int((post["p"] < alpha).sum())}]),
                         {c: ".3f" for c in ("median r", "2.5th percentile", "97.5th percentile", "minimum", "maximum", "MCC tree r")}), "",
                f"{int((post['r'] < head['r']).sum())} of the {len(post)} trees give a smaller r than the MCC tree. "
                f"The posterior trees' pairwise distances have a rank correlation of {post['spearman_with_mcc'].min():.2f} to {post['spearman_with_mcc'].max():.2f} with the MCC tree's. "
                f"Per-tree results are in `{rel(posterior_path)}`.", "",
                f"![Mantel r across posterior trees]({rel(fig_dir / 'h1_posterior.png')})", ""]

    traits, rates = read_csv(cfg.processed / "traits.csv"), read_csv(tables / "ecology_missing.csv")
    name_map, residual = read_csv(cfg.interim / "traits_name_map.csv"), read_csv(cfg.interim / "traits_residual_review.csv")
    eco_species = read_species_order(cfg.matrices_dir() / "ecology_species.txt")
    out += ["## 3. Trait data and ecological distance", "",
            "Traits come from EltonTraits 1.0 (diet, foraging stratum, activity time) and COMBINE's reported values (habitat breadth, terrestrial or aquatic, body mass). "
            "Every trait, its type and its weight is in `docs/ecology_traits.md`. Ecological distance is Gower distance with the pairwise-available rule; nothing is imputed.", ""]
    match = name_map.assign(how=name_map["match_method"].fillna("no record")).groupby(["database", "how"]).size().unstack(fill_value=0).reset_index()
    match.columns.name = None
    out += ["**Name matching** (species per database, by the rule that matched them):", "", md_table(match), ""]
    if len(residual):
        out += [f"**Species without a record** ({residual['species'].nunique()}; `data/interim/traits_residual_review.csv`):", "",
                md_table(residual[["species", "common_name", "family", "database"]].rename(columns={"common_name": "common name", "database": "missing from"})), ""]
    group_rates = rates.groupby(["group", "source"], sort=False).agg(columns=("trait", "size"), missing=("n_missing", "max"), rate=("missing_rate", "max")).reset_index()
    out += ["**Missing data by trait group** (the columns in a group are missing together, except where noted in `reports/tables/ecology_missing.csv`):", "",
            md_table(group_rates.rename(columns={"missing": "species missing", "rate": "missing rate"}), {"missing rate": ".1%"}), "",
            f"{int(traits.drop(columns='species_id').notna().all(axis=1).sum())} of {len(traits)} species have every trait. "
            f"{len(eco_species)} of {len(order)} species have an ecological distance to every other species and are in `ecology.npy` and `ecology_unweighted.npy`.", ""]

    att = read_csv(cfg.processed / "attention.csv")
    months = month_range(acfg["start"], acfg["end"])
    shared = att[att.duplicated("current_title", keep=False)]
    q = att["median_monthly_pageviews"].quantile([0, 0.25, 0.5, 0.75, 1]).to_numpy()
    out += ["## 4. Attention data", "",
            f"Monthly English Wikipedia pageviews from the Wikimedia Pageviews API, {acfg['start']} to {acfg['end']} ({len(months)} months), {acfg['access']}, agent = {acfg['agent']}. "
            "Each species' count is the sum over its article and every main-namespace redirect to it, so views recorded under an earlier title are included.", "",
            f"- Coverage: {int((att['total_pageviews'] > 0).sum())} of {len(att)} species have pageviews; {int((att['months_with_views'] < len(months)).sum())} have at least one month with none.",
            f"- Median monthly pageviews: minimum {q[0]:,.0f}, lower quartile {q[1]:,.0f}, median {q[2]:,.0f}, upper quartile {q[3]:,.0f}, maximum {q[4]:,.0f}.",
            f"- Redirects: {int(att['n_redirects'].sum()):,} redirect titles in total (median {att['n_redirects'].median():.0f} per article). "
            f"They carry a median of {att['redirect_share'].median():.1%} of an article's views and at most {att['redirect_share'].max():.1%}.",
            f"- {int(att['renamed'].sum())} of the article titles recorded in Week 2 now redirect to a renamed article.",
            f"- {len(shared)} species share an article with another species" + (f" ({', '.join(shared['species_id'])})." if len(shared) else "."),
            f"- Article length: `n_tokens` from the corpus (unmasked text), minimum {att['n_tokens'].min():,}, median {att['n_tokens'].median():,.0f}, maximum {att['n_tokens'].max():,}.",
            f"- log pageviews and log length have a Spearman correlation of {att['log_pageviews'].corr(att['log_length'], method='spearman'):.2f} across species.", ""]
    top, bottom = att.nlargest(5, "median_monthly_pageviews"), att.nsmallest(5, "median_monthly_pageviews")
    ends = pd.concat([top, bottom])[["species_id", "current_title", "median_monthly_pageviews", "n_tokens"]].rename(
        columns={"species_id": "species", "current_title": "article", "median_monthly_pageviews": "median monthly pageviews", "n_tokens": "tokens"})
    out += ["Most and least viewed:", "", md_table(ends, {"median monthly pageviews": ",.0f", "tokens": ","}), "",
            "Matrices: `attention_pageviews.npy` (|difference in log(1 + median monthly pageviews)|), `attention_length.npy` (|difference in log tokens|), "
            "`attention_combined.npy` (Euclidean distance on the two standardized variables) and `attention_mean_log_pageviews.npy` "
            "(mean log pageviews of the pair, for the later \"both well-known\" check; not a distance).", ""]

    diag = read_csv(tables / "confound_mantel.csv")
    diag_meta = read_meta(tables / "confound_mantel.csv")
    shown = diag.assign(x=diag["x"].map(MATRIX_NAMES), y=diag["y"].map(MATRIX_NAMES))[["x", "y", "r", "p", "n_species"]].rename(columns={"n_species": "species"})
    out += ["## 5. Confound diagnostics", "",
            f"{method.capitalize()} Mantel correlations, {perms:,} permutations, {diag_meta.get('alternative', 'two-sided')} p-values. "
            "These are marginal correlations: nothing is partialled out, and they say nothing yet about H2.", "",
            "**Among the predictor matrices:**", "", md_table(shown[diag["role"] == "predictor overlap"], {"r": ".3f", "p": ".4f"}), "",
            f"**Text distance ({combo_label(*primary)}) against each control:**", "", md_table(shown[diag["role"] == "text against predictor"], {"r": ".3f", "p": ".4f"}), "",
            f"For comparison, text against phylogeny in the same specification is r = {head['r']:.3f}.", ""]

    notes = cfg.reports / "week4_notes.md"
    if notes.exists():
        out += [notes.read_text().rstrip(), ""]
    atomic_write_text(report_path, "\n".join(out))
    print(f"[report] wrote {report_path.relative_to(cfg.root)}")


if __name__ == "__main__":
    main()
