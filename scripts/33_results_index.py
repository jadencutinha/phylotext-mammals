from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from lineage_charisma.config import base_parser, load_config
from lineage_charisma.io_utils import atomic_write_text
from lineage_charisma.plotstyle import combo_label

# One row per pre-registered result (plus the replication checks and the disclosed post-hoc analysis), read from the
# result tables so the index cannot drift from them. Rerun after any result table changes.
ALPHA = 0.05
ROLE_ORDER = ["primary", "secondary", "robustness", "exploratory", "check", "post-hoc"]


def fmt(x, digits=3):
    return "" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{x:.{digits}f}"


def pfmt(p):
    return "" if p is None or (isinstance(p, float) and np.isnan(p)) else f"{p:.4f}"


def h1_decision(r, p):
    return "supported" if r > 0 and p < ALPHA else "not supported"


def h2_decision(row):
    a = "holds" if row["survives"] else "fails"
    conf = row["confounded"]
    b = "undefined (raw r ≤ 0)" if pd.isna(conf) else ("holds" if bool(conf) else "fails")
    return f"rule (a) {a}; rule (b) {b}"


def spec(row):
    return combo_label(row["model"], row["rule"], row["mask"])


def main() -> None:
    args = base_parser("Write reports/final_results_index.md from every result table.").parse_args()
    cfg = load_config(args.config, None)
    t = cfg.root / "reports" / "tables"
    m, s = t / "mammalia", t / "secondary"
    rel = lambda p: str(Path(p).relative_to(cfg.root))
    rows = []

    def add(group, result, role, spec_, n, stat, p, perms, decision, source):
        rows.append({"group": group, "result": result, "role": role, "spec": spec_, "n": n, "statistic": stat, "p": p, "permutations": perms,
                     "decision": decision, "source": rel(source)})

    # ---- pilot (Carnivora): docs/preregistration.md, amendment 2
    h1 = pd.read_csv(t / "h1_mantel.csv")
    for _, r in h1.iterrows():
        primary = bool(r["is_primary"])
        add("Pilot H1", "H1 Mantel r" + ("" if r["method"] == "spearman" else " (Pearson)"), "primary" if primary else "robustness", spec(r), 284,
            f"r = {fmt(r['r'])}", r["p"], 9999, h1_decision(r["r"], r["p"]) if primary else "", t / "h1_mantel.csv")
    post = pd.read_csv(t / "h1_posterior.csv")
    add("Pilot H1", f"H1 on {len(post)} posterior trees", "robustness", "primary spec", 284,
        f"median r = {fmt(post['r'].median())} [{fmt(post['r'].quantile(.025))}, {fmt(post['r'].quantile(.975))}]", None, 9999,
        f"p < 0.05 on {int((post['p'] < ALPHA).sum())}/{len(post)} trees", t / "h1_posterior.csv")
    fam = pd.read_csv(t / "h1_mantel_families.csv")
    add("Pilot H1", f"H1 within {len(fam)} families", "exploratory", "primary spec", "", f"r from {fmt(fam['r'].min())} to {fmt(fam['r'].max())}", None, 9999,
        "descriptive", t / "h1_mantel_families.csv")

    h2 = pd.read_csv(t / "h2_partial_mantel.csv")
    for _, r in h2.iterrows():
        if r["analysis"] == "primary":
            role, name = "primary", "H2 partial r (ecology + length)"
        elif r["analysis"] == "secondary":
            role, name = "secondary", f"H2 {r['model_id']} partial r ({r['controls']})"
        else:
            role, name = "robustness", f"H2 partial r, {r['analysis'].replace('_', ' ')}"
        add("Pilot H2", name, role, spec(r) + ("" if r["ecology_matrix"] == "ecology" else f", {r['ecology_matrix']}"), int(r["n_species"]),
            f"partial r = {fmt(r['partial_r'])} (raw {fmt(r['raw_r'])}, share {fmt(r['share_of_raw'], 2)})", r["partial_p"], 9999,
            h2_decision(r) if role == "primary" else "", t / "h2_partial_mantel.csv")
    mrm = pd.read_csv(t / "h2_mrm.csv")
    for _, r in mrm[mrm["is_primary"]].iterrows():
        add("Pilot H2", f"MRM coefficient: {r['predictor']}", "secondary", "primary spec, descriptive", int(r["n_species"]),
            f"β = {fmt(r['coef_std'])} (R² = {fmt(r['r2'], 4)})", r["p"], 9999, "descriptive (two-sided)", t / "h2_mrm.csv")
    hp = pd.read_csv(t / "h2_posterior.csv")
    add("Pilot H2", f"H2 on {len(hp)} posterior trees", "robustness", "primary spec", 283,
        f"median partial r = {fmt(hp['partial_r'].median())}", None, 9999,
        f"rule (a) on {int(hp['survives'].sum())}/{len(hp)}; rule (b) on {int(hp['confounded'].fillna(False).astype(bool).sum())}/{len(hp)}", t / "h2_posterior.csv")
    wk = pd.read_csv(t / "exploratory_well_known.csv").dropna(subset=["r_text_phylo"])
    for _, r in wk.iterrows():
        add("Pilot exploratory", f"text–phylogeny r, {r['quartile']} pageview quartile", "exploratory", "primary spec", int(r["n_pairs"]),
            f"r = {fmt(r['r_text_phylo'])}", r["p"], 9999, "", t / "exploratory_well_known.csv")
    for _, r in pd.read_csv(t / "exploratory_nonlinear.csv").iloc[1:].iterrows():
        add("Pilot exploratory", f"H1 with phylogeny {r['transform']} ({r['method']})", "exploratory", "primary spec", int(r["h1_n_species"]),
            f"r = {fmt(r['h1_r'])}; partial r = {fmt(r['partial_r'])}", r["h1_p"], 9999, "", t / "exploratory_nonlinear.csv")

    # ---- full scale (all mammals): docs/preregistration_amendment_fullscale.md
    for _, r in pd.read_csv(m / "replication.csv").iterrows():
        add("Full-scale checks", r["check"], "check", "primary spec, 284 pilot species", int(r["n_species"]),
            f"r = {r['r']:.6f} (pilot {r['pilot_r']:.6f})", r["p"], 9999, "passed" if r["passed"] else "FAILED", m / "replication.csv")
    fh1 = pd.read_csv(m / "h1_mantel.csv")
    for _, r in fh1.iterrows():
        primary = bool(r["is_primary"])
        add("Full-scale H1", "H1 Mantel r" + (" (Pearson)" if r["method"] == "pearson" else ""), "primary" if primary else "robustness", spec(r),
            int(r["n_species"]), f"r = {fmt(r['r'])}", r["p"], int(r["permutations"]), h1_decision(r["r"], r["p"]) if primary else "", m / "h1_mantel.csv")
    fh2 = pd.read_csv(m / "h2_partial_mantel.csv")
    for _, r in fh2.iterrows():
        a = r["analysis"]
        role = "primary" if a == "primary" else "secondary" if a.startswith("secondary") else "robustness"
        name = "H2 partial r (ecology + length)" if a == "primary" else f"H2 {a.replace('secondary ', '')} partial r ({r['controls']})" if role == "secondary" \
            else f"H2 partial r, {a.replace('_', ' ')}"
        add("Full-scale H2", name, role, spec(r) + ("" if r["ecology_matrix"] == "ecology" else f", {r['ecology_matrix']}"), int(r["n_species"]),
            f"partial r = {fmt(r['partial_r'])} (raw {fmt(r['raw_r'])})", r["partial_p"], int(r["permutations"]), h2_decision(r) if role == "primary" else "",
            m / "h2_partial_mantel.csv")
    fm = pd.read_csv(m / "h2_mrm.csv")
    for _, r in fm[fm["is_primary"]].iterrows():
        add("Full-scale H2", f"MRM coefficient: {r['predictor']}", "secondary", "primary spec, descriptive", int(r["n_species"]),
            f"β = {fmt(r['coef_std'])} (R² = {fmt(r['r2'], 4)})", r["p"], int(r["permutations"]), "descriptive (two-sided)", m / "h2_mrm.csv")
    ps = pd.read_csv(m / "posterior_summary.csv").set_index("statistic")
    meta = json.loads((m / "posterior_summary.csv.meta.json").read_text())
    add("Full-scale H1", "H1 on 20 posterior trees", "robustness", "primary spec", 4611,
        f"median r = {fmt(ps.at['h1_r', 'median'])} [{fmt(ps.at['h1_r', 'min'])}, {fmt(ps.at['h1_r', 'max'])}]", None, 999,
        f"p < 0.05 on {meta['h1_p_below_alpha']:.0%} of trees", m / "posterior_summary.csv")
    add("Full-scale H2", "H2 on 20 posterior trees", "robustness", "primary spec", 4519,
        f"median partial r = {fmt(ps.at['partial_r', 'median'])}", None, 999,
        f"rule (a) on {meta['rule_a_holds']:.0%}; rule (b) on {meta['rule_b_holds']:.0%} (raw r near 0)", m / "posterior_summary.csv")
    for _, r in pd.read_csv(m / "per_order.csv").iterrows():
        add("Full-scale per order", f"{r['order']}: H1 r / H2 partial r", "exploratory", "primary spec, descriptive", int(r["n_species"]),
            f"H1 r = {fmt(r['h1_r'])}; partial r = {fmt(r['partial_r'])}", r["h1_p"], int(r["permutations"]), "descriptive (H1 p shown)", m / "per_order.csv")
    q = pd.read_csv(m / "exploratory_pageview_quartiles.csv")
    for _, r in q[(q["group"] == "all mammals") & q["r_text_phylo"].notna()].iterrows():
        add("Full-scale exploratory", f"text–phylogeny r, {r['quartile']} pageview quartile", "exploratory", "primary spec", int(r["n_pairs"]),
            f"r = {fmt(r['r_text_phylo'])}", r["p"], 999, "", m / "exploratory_pageview_quartiles.csv")
    for _, r in pd.read_csv(m / "exploratory_nonlinear.csv").iloc[1:].iterrows():
        add("Full-scale exploratory", f"H1 with phylogeny {r['transform']} ({r['method']})", "exploratory", "primary spec", int(r["h1_n_species"]),
            f"r = {fmt(r['h1_r'])}; partial r = {fmt(r['partial_r'])}", r["h1_p"], 999, "", m / "exploratory_nonlinear.csv")
    lm = pd.read_csv(m / "posthoc_length_matched.csv").iloc[0]
    add("Full-scale post-hoc", "Length-matched within-decile r (amendment 4, D)", "post-hoc", "primary spec, 10 length deciles", int(lm["n_species"]),
        f"r_LM = {fmt(lm['r_lm'], 4)}", lm["p_two_sided"], int(lm["permutations"]), "positive, p < 0.05: consistent with length composition (post-hoc)",
        m / "posthoc_length_matched.csv")

    # ---- secondary analyses: docs/preregistration_amendment_secondary.md
    if (s / "register_comparison.csv").exists():
        rc, bs = pd.read_csv(s / "register_comparison.csv"), pd.read_csv(s / "register_bootstrap.csv").set_index("test")
        for _, r in rc.iterrows():
            add("Secondary A: register", f"{r['test']} r, {'ADW' if r['source'] == 'adw' else 'Wikipedia'}", "secondary", "bge-large / chunk / strict, pilot masks",
                int(r["n_species"]), f"r = {fmt(r['r'])}", r["p"], 9999, "", s / "register_comparison.csv")
        for test, r in bs.iterrows():
            add("Secondary A: register", f"Δ {test} (ADW − Wikipedia)", "secondary", "paired species bootstrap, 2,000 draws", int(r["n_species"]),
                f"Δ = {r['delta']:+.3f} [{r['ci_low']:+.3f}, {r['ci_high']:+.3f}]", None, "", r["reading"], s / "register_bootstrap.csv")
    if (s / "charisma_tests.csv").exists():
        for _, r in pd.read_csv(s / "charisma_tests.csv").iterrows():
            add("Secondary B: charisma", r["test"], "secondary", "claude-sonnet-5-5, retained dimensions", int(r["n_species"]), f"partial r = {fmt(r['r'])}",
                r["p"], 9999, ("explains beyond attention" if r.get("explains_beyond_attention") else "") if r["test"].startswith("B1") else "",
                s / "charisma_tests.csv")
        if (s / "charisma_mrm.csv").exists():
            for _, r in pd.read_csv(s / "charisma_mrm.csv").iterrows():
                add("Secondary B: charisma", f"B3 extended MRM coefficient: {r['predictor']}", "secondary", "claude-sonnet-5-5, descriptive",
                    int(r["n_species"]), f"β = {fmt(r['coef_std'])} (R² = {fmt(r['r2'], 4)})", r["p"], 9999, "descriptive (two-sided)", s / "charisma_mrm.csv")
    else:
        add("Secondary B: charisma", "Agreement, drop rule, B1–B3", "secondary", "claude-sonnet-5-5", "", "placeholder", None, "",
            "pending: hand labels not yet committed; no model request sent", cfg.root / "docs" / "preregistration_amendment_secondary.md")
    if (s / "zoo_summary.json").exists():
        add("Secondary C: zoo signage", "Cross-register similarity; H1 on overlap", "secondary", "descriptive", "", "see zoo_summary.json", None, "", "",
            s / "zoo_summary.json")
    else:
        # amendment 4, section C: no signage was collected, so nothing is simulated; the pipeline is tested on 3 fake rows only
        add("Secondary C: zoo signage", "Cross-register similarity; H1 on overlap", "secondary", "descriptive", "", "not run", None, "",
            "built and tested (tests/test_zoo.py, 3 fake rows), not run: no signage data collected; moved to future work",
            cfg.root / "scripts" / "28_zoo_signage.py")

    table = pd.DataFrame(rows)
    counts = table["role"].value_counts().reindex(ROLE_ORDER).fillna(0).astype(int)
    lines = ["# Final results index", "",
             "Every pre-registered result, the replication checks and the one disclosed post-hoc analysis, generated from the result tables by "
             "`scripts/33_results_index.py` (rerun it after any table changes). Roles: **primary** (the registered H1 and H2 decisions), **secondary** "
             "(registered tests that cannot change the H1/H2 conclusions), **robustness** (registered one-at-a-time variants), **exploratory** "
             "(registered, labelled exploratory or descriptive), **check** (registered replication checks), and **post-hoc** (designed after seeing "
             "the data, disclosed in amendment 4, section D). p-values are one-sided unless the decision column says otherwise; MRM p-values are two-sided.",
             "", "Counts: " + ", ".join(f"{k} {v}" for k, v in counts.items()) + ".", "",
             "**Locked 2026-11-11 (week 8).** After this date analysis code changes only for bug fixes, each reported. Every number below is "
             "read from the table in its source column; the numbers quoted in the reports, README and amendments are checked against the "
             "same tables by `scripts/34_verify_numbers.py` (`reports/number_check.md`).", ""]
    for group, sub in table.groupby("group", sort=False):
        lines += [f"## {group}", "", "| result | role | spec | n | statistic | p | perms | decision | source |", "|---|---|---|---:|---|---:|---:|---|---|"]
        for _, r in sub.iterrows():
            lines.append(f"| {r['result']} | {r['role']} | {r['spec']} | {r['n']} | {r['statistic']} | {pfmt(r['p'])} | {r['permutations']} | {r['decision']} | `{r['source']}` |")
        lines.append("")
    atomic_write_text(cfg.root / "reports" / "final_results_index.md", "\n".join(lines))
    print(f"[index] {len(table)} rows -> reports/final_results_index.md; {counts.to_dict()}")


if __name__ == "__main__":
    main()
