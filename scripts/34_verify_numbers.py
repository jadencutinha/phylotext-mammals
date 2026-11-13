from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

from lineage_charisma.config import base_parser, load_config
from lineage_charisma.distances import load_matrix, text_matrix_path
from lineage_charisma.embed import embedding_paths, read_species_order
from lineage_charisma.io_utils import atomic_write_text, read_csv, read_json

# Week 8 lock: every headline number quoted in the reports, README and amendments is re-derived from the saved result tables and
# compared at the precision the text uses. A second pass lists every number with three or more decimals in the reports that matches
# no value in any table, so numbers computed outside the tables (diagnostics) are visible. The week 6 same-order diagnostics, which
# were never saved as a table, are recomputed from the full-scale matrices. Writes reports/number_check.md and exits non-zero on any
# mismatch that is not a recorded erratum.

NUMBER = r"([−+-]?\d+\.\d+)"
SCANNED = ["reports/week4_report.md", "reports/week5_report.md", "reports/week6_report.md", "reports/week7_report.md", "README.md",
           "reports/final_facts.md", "reports/figures/final/captions.md"]


def to_float(text: str) -> float:
    return float(text.replace("−", "-").replace("+", ""))


def matches(text: str, value: float) -> bool:
    decimals = len(text.split(".")[1])
    return abs(to_float(text) - value) <= 0.5 * 10 ** -decimals + 1e-12


def claims(t: Path) -> list[tuple[str, str, str, float]]:
    """(description, file, regex with one number group, value from the tables)."""
    m, s = t / "mammalia", t / "secondary"
    h1 = pd.read_csv(t / "h1_mantel.csv").query("is_primary").iloc[0]
    h2 = pd.read_csv(t / "h2_partial_mantel.csv").set_index("analysis")
    mrm = pd.read_csv(t / "h2_mrm.csv").query("is_primary").set_index("predictor")
    fam = pd.read_csv(t / "h1_mantel_families.csv").query("method == 'spearman'").set_index("family")
    post = pd.read_csv(t / "h1_posterior.csv")
    fh1 = pd.read_csv(m / "h1_mantel.csv").set_index("analysis")
    fh1_tax = pd.read_csv(m / "h1_mantel.csv").query("mask == 'taxonomy'").iloc[0]
    fh2 = pd.read_csv(m / "h2_partial_mantel.csv").set_index("analysis")
    fmrm = pd.read_csv(m / "h2_mrm.csv").query("is_primary").set_index("predictor")
    orders = pd.read_csv(m / "per_order.csv").set_index("order")
    ps = pd.read_csv(m / "posterior_summary.csv").set_index("statistic")
    q = pd.read_csv(m / "exploratory_pageview_quartiles.csv").query("group == 'all mammals' and r_text_phylo.notna()").set_index("quartile")
    nl = pd.read_csv(m / "exploratory_nonlinear.csv").set_index(["transform", "method"])
    rep = pd.read_csv(m / "replication.csv").set_index("check")
    lm = pd.read_csv(m / "posthoc_length_matched.csv").iloc[0]
    rc = pd.read_csv(s / "register_comparison.csv").set_index(["test", "source"])
    bs = pd.read_csv(s / "register_bootstrap.csv").set_index("test")
    p = h2.loc["primary"]
    w4, w5, w6, w7, rd = "reports/week4_report.md", "reports/week5_report.md", "reports/week6_report.md", "reports/week7_report.md", "README.md"
    a3, a4 = "docs/preregistration_amendment_fullscale.md", "docs/preregistration_amendment_secondary.md"
    out = [
        ("pilot H1 r", w4, rf"Mantel r = {NUMBER}, one-sided", h1["r"]),
        ("pilot H1 p", w4, rf"Mantel r = [\d.]+, one-sided p = {NUMBER}", h1["p"]),
        ("pilot H1 r", w6, rf"\| H1, primary \| r = {NUMBER}", h1["r"]),
        ("pilot H1 r", rd, rf"\(r = {NUMBER} on 284 species", h1["r"]),
        ("pilot H1 r", a3, rf"H1 r = {NUMBER}, one-sided", h1["r"]),
        ("pilot H1 p", a3, rf"one-sided p = {NUMBER}, 284", h1["p"]),
        ("pilot H1 r", a4, rf"H1 r = {NUMBER}, p = 0", h1["r"]),
        ("pilot H1 p", a4, rf"H1 r = [\d.]+, p = {NUMBER}", h1["p"]),
        ("pilot H1 posterior median", w4, rf"\| 100 \| {NUMBER} \|", post["r"].median()),
        ("pilot H1 posterior median", a3, rf"median r = {NUMBER} over 100", post["r"].median()),
        ("pilot H1 Canidae r", w4, rf"Canidae shows a clear positive correlation \(r = {NUMBER}", fam.loc["Canidae", "r"]),
        ("pilot H1 Mustelidae r", w4, rf"show none \(r = {NUMBER}", fam.loc["Mustelidae", "r"]),
        ("pilot H2 raw r", w5, rf"raw r on 283 species is {NUMBER}", p["raw_r"]),
        ("pilot H2 partial r", w5, rf"partial r given ecology and length is {NUMBER}", p["partial_r"]),
        ("pilot H2 partial p", w5, rf"ecology and length is [\d.]+ \(p = {NUMBER}", p["partial_p"]),
        ("pilot H2 partial r", w6, rf"\| H2 partial r \(ecology \+ length\) \| {NUMBER}", p["partial_r"]),
        ("pilot H2 partial p", rd, rf"partial r = [\d.]+, p = {NUMBER}\)", p["partial_p"]),
        ("pilot H2 partial r", a4, rf"partial r \(ecology \+ length\) = {NUMBER}", p["partial_r"]),
        ("pilot H2 S4 partial r", w5, rf"lifts the partial r from [\d.]+ to {NUMBER}", h2.loc["secondary"].set_index("model_id").loc["S4", "partial_r"]),
        ("pilot MRM pageviews β", w5, rf"pageview coefficient is negative \({NUMBER}", mrm.loc["pageviews", "coef_std"]),
        ("pilot MRM phylogeny β", w6, rf"MRM phylogeny coefficient \| {NUMBER}", mrm.loc["phylogeny", "coef_std"]),
        ("replication 6a r", w7, rf"\(r = {NUMBER}, p = 0.0285\)", rep.loc["6a exact (pilot vocabulary)", "r"]),
        ("replication 6b r", w7, rf"full-vocabulary replication gave r = {NUMBER}", rep.loc["6b full-scale vocabulary", "r"]),
        ("full H1 r", w6, rf"\*\*r = {NUMBER}, p = 1", fh1.loc["primary", "r"]),
        ("full H1 p", w6, rf"\*\*r = −?[\d.]+, p = {NUMBER}\*\*", fh1.loc["primary", "p"]),
        ("full H1 r", w7, rf"\(r = {NUMBER}, p = 1.0, 4,611", fh1.loc["primary", "r"]),
        ("full H1 r", rd, rf"slightly negative \(r = {NUMBER}\)", fh1.loc["primary", "r"]),
        ("full H1 r", a4, rf"H1 r = {NUMBER}, H2 rule", fh1.loc["primary", "r"]),
        ("full H1 taxonomy-mask r", w6, rf"gives r = {NUMBER} and a partial", fh1_tax["r"]),
        ("full H1 taxonomy-mask r", w7, rf"Taxonomy mask:\*\* r = {NUMBER}", fh1_tax["r"]),
        ("full H1 taxonomy-mask r", rd, rf"unmasked it is r = {NUMBER}", fh1_tax["r"]),
        ("full H2 partial r", w6, rf"\| \*\*{NUMBER}, p = 1.000\*\* \| rule", fh2.loc["primary", "partial_r"]),
        ("full MRM phylogeny β", w6, rf"\| [\d.]+, p = 0.13 \| {NUMBER}, p", fmrm.loc["phylogeny", "coef_std"]),
        ("full MRM pageviews β", w6, rf"pageviews {NUMBER} \(p = ", fmrm.loc["pageviews", "coef_std"]),
        ("full Chiroptera H1 r", w6, rf"Chiroptera {NUMBER}", orders.loc["Chiroptera", "h1_r"]),
        ("full Lagomorpha H1 r", w6, rf"Lagomorpha {NUMBER}", orders.loc["Lagomorpha", "h1_r"]),
        ("full Primates H1 r", w6, rf"Primates {NUMBER} \(all", orders.loc["Primates", "h1_r"]),
        ("full Carnivora H1 r", w6, rf"Carnivora {NUMBER} \(p", orders.loc["Carnivora", "h1_r"]),
        ("full Chiroptera H1 r", w7, rf"bats \({NUMBER}\)", orders.loc["Chiroptera", "h1_r"]),
        ("full Artiodactyla H1 r", w7, rf"artiodactyls \({NUMBER}\)", orders.loc["Artiodactyla", "h1_r"]),
        ("full Primates H1 r", w7, rf"primates \({NUMBER}\)", orders.loc["Primates", "h1_r"]),
        ("full Lagomorpha H1 r", w7, rf"lagomorphs \({NUMBER}\)", orders.loc["Lagomorpha", "h1_r"]),
        ("full posterior median H1 r", w6, rf"H1 r median {NUMBER}, range", ps.at["h1_r", "median"]),
        ("full posterior min H1 r", w6, rf"H1 r median −?[\d.]+, range {NUMBER} to", ps.at["h1_r", "min"]),
        ("full posterior max H1 r", w6, rf"H1 r median −?[\d.]+, range −?[\d.]+ to {NUMBER}", ps.at["h1_r", "max"]),
        ("full posterior median partial r", w6, rf"partial r has median {NUMBER} \(range", ps.at["partial_r", "median"]),
        ("full posterior min partial r", w6, rf"partial r has median −?[\d.]+ \(range {NUMBER} to", ps.at["partial_r", "min"]),
        ("full posterior max partial r", w6, rf"partial r has median −?[\d.]+ \(range −?[\d.]+ to {NUMBER}\)", ps.at["partial_r", "max"]),
        ("full Q1 r", w6, rf"r is {NUMBER} \(p = [\d.]+\) in the least", q.loc["Q1 (least viewed)", "r_text_phylo"]),
        ("full Q4 r", w6, rf"and {NUMBER} \(p = [\d.]+\) in the most", q.loc["Q4 (most viewed)", "r_text_phylo"]),
        ("full 30 My cap r", w6, rf"H1 r is {NUMBER} \(Spearman", nl.loc[("capped", "spearman"), "h1_r"]),
        ("post-hoc r_LM", w7, rf"r_LM = {NUMBER}, two-sided", lm["r_lm"]),
        ("post-hoc p", w7, rf"r_LM = [+\d.]+, two-sided p = {NUMBER}", lm["p_two_sided"]),
        ("post-hoc r_LM", a4, rf"r_LM = {NUMBER}, two-sided", lm["r_lm"]),
        ("post-hoc p", a4, rf"r_LM = [\d.]+, two-sided p = {NUMBER}", lm["p_two_sided"]),
        ("post-hoc 5-bin r_LM", w7, rf"With 5 bins r_LM is {NUMBER}", lm["r_lm_5_bins"]),
        ("post-hoc 20-bin r_LM", w7, rf"with 20 bins {NUMBER}", lm["r_lm_20_bins"]),
        ("ADW H1 r, Wikipedia", w7, rf"\| H1 r \(246 species\) \| {NUMBER}", rc.loc[("H1", "wikipedia"), "r"]),
        ("ADW H1 r, ADW", w7, rf"\| H1 r \(246 species\) \| [\d.]+, p = [\d.]+ \| {NUMBER}", rc.loc[("H1", "adw"), "r"]),
        ("ADW H1 Δ", w7, rf"\*\*{NUMBER} \[\+0.019", bs.loc["H1", "delta"]),
        ("ADW H1 Δ low", w7, rf"\*\*\+0.056 \[{NUMBER},", bs.loc["H1", "ci_low"]),
        ("ADW H1 Δ high", w7, rf"\*\*\+0.056 \[\+0.019, {NUMBER}\]", bs.loc["H1", "ci_high"]),
        ("ADW H2 Δ", w7, rf"\*\*{NUMBER} \[\+0.013", bs.loc["H2 partial", "delta"]),
        ("ADW H2 partial r, ADW", w7, rf"\| primary H2 partial r \(245 species\) \| [\d.]+, p = [\d.]+ \| {NUMBER}", rc.loc[("H2 partial", "adw"), "r"]),
    ]
    return out


# quoted values found wrong by this check; the committed text is left as it was and the correction is recorded here and in the
# week 8 report
ERRATA = {("taxonomy mask, different-order mean text distance", "reports/week6_report.md"): "0.442 should read 0.441 (0.44147)",
          ("taxonomy mask, different-order mean text distance", "docs/preregistration_amendment_secondary.md"): "0.442 should read 0.441 (0.44147)"}


def diagnostic_claims(cfg) -> list[tuple[str, str, str, float]]:
    """Week 6 report, section 1 (unregistered diagnostics): mean strict- and taxonomy-mask text distance for same-order and
    different-order pairs, and the strict-mask gap among pairs of texts both below or both at/above the median token count."""
    full = load_config(None, "mammalia")
    md, op = full.matrices_dir(), full.species_order_path()
    order = read_species_order(op)
    orders = read_csv(full.matched_path()).set_index("species").loc[order, "order"].to_numpy()
    iu = np.triu_indices(len(order), 1)
    same = (orders[:, None] == orders[None, :])[iu]
    model, rule, _ = full.primary_combo
    tokens = np.array([read_json(embedding_paths(full.embeddings_dir(), model, rule, "strict").meta)["n_tokens"][s] for s in order])
    short = tokens < np.median(tokens)
    v = {}
    for mask in ("strict", "taxonomy"):
        d = load_matrix(text_matrix_path(md, model, rule, mask), op)[0][iu]
        v[mask] = (d[same].mean(), d[~same].mean())
        if mask == "strict":
            for name, group in (("short", short), ("long", ~short)):
                both = (group[:, None] & group[None, :])[iu]
                v[name] = (d[both & same].mean(), d[both & ~same].mean())
    out = []
    for file in ("reports/week6_report.md", "docs/preregistration_amendment_secondary.md"):
        out += [
            ("strict mask, same-order mean text distance", file, rf"(?:mean distance |distance ){NUMBER}\)? (?:than|and different)", v["strict"][0]),
            ("strict mask, different-order mean text distance", file, rf"(?:pairs from different orders \(|different-order pairs ){NUMBER}", v["strict"][1]),
            ("taxonomy mask, same-order mean text distance", file, rf"(?:clearly more alike \(|taxonomy masking, ){NUMBER}(?: vs| and)", v["taxonomy"][0]),
            ("taxonomy mask, different-order mean text distance", file, rf"(?:clearly more alike \(|taxonomy masking, )[\d.]+(?: vs | and ){NUMBER}", v["taxonomy"][1]),
            ("strict gap, short texts (same-order)", file, rf"{NUMBER} vs [\d.]+ (?:for|among) pairs of short", v["short"][0]),
            ("strict gap, short texts (different-order)", file, rf"[\d.]+ vs {NUMBER} (?:for|among) pairs of short", v["short"][1]),
            ("strict gap, long texts (same-order)", file, rf"{NUMBER} vs [\d.]+ (?:for|among) (?:pairs of )?long", v["long"][0]),
            ("strict gap, long texts (different-order)", file, rf"[\d.]+ vs {NUMBER} (?:for|among) (?:pairs of )?long", v["long"][1]),
        ]
    return out


def table_values(t: Path) -> np.ndarray:
    values = []
    for f in t.rglob("*.csv"):
        numeric = pd.read_csv(f).select_dtypes("number")
        values += numeric.to_numpy().ravel().tolist()
        if "posterior" in f.name:   # the reports quote percentiles of the posterior-tree tables
            values += numeric.quantile([0.025, 0.5, 0.975]).to_numpy().ravel().tolist()
    v = np.array(values, dtype=float)
    return v[np.isfinite(v)]


def main() -> None:
    args = base_parser("Re-derive every quoted headline number from the saved result tables (week 8 lock).").parse_args()
    cfg = load_config(args.config, None)
    root, t = cfg.root, cfg.root / "reports" / "tables"
    rows, bad = [], 0
    diagnostics = diagnostic_claims(cfg)
    for desc, file, pattern, value in claims(t) + diagnostics:
        text = (root / file).read_text()
        found = re.findall(pattern, text)
        status = "not found in file" if not found else ("ok" if all(matches(f, value) for f in found) else "MISMATCH")
        if status == "MISMATCH" and (desc, file) in ERRATA:
            status = f"known erratum: {ERRATA[(desc, file)]}"
        bad += status not in ("ok",) and not status.startswith("known erratum")
        rows.append(f"| {desc} | `{file}` | {', '.join(sorted(set(found))) or '—'} | {value:.6g} | {status} |")

    values = np.concatenate([table_values(t), [d[3] for d in diagnostics]])   # recomputed diagnostics count as traced
    untraced = []
    number = re.compile(r"(?<![\w.:/])[−-]?\d+\.\d{3,}(?!\d)")
    for file in SCANNED:
        if not (root / file).exists():
            continue
        for i, line in enumerate((root / file).read_text().splitlines(), 1):
            for n in number.findall(line):
                x = to_float(n)
                tol = 0.5 * 10 ** -len(n.split(".")[1]) + 1e-12
                if not ((np.abs(values - x) <= tol).any() or (np.abs(values + x) <= tol).any()):
                    untraced.append(f"| `{file}:{i}` | {n} |")

    lines = ["# Number check (week 8 lock)", "",
             "Generated by `scripts/34_verify_numbers.py`. Each headline number quoted in the reports, README and amendments is "
             "re-derived from the saved result tables in `reports/tables/` and compared at the precision the text uses. The week 6 "
             "same-order diagnostics (never saved as a table) are recomputed from the full-scale distance matrices.", "",
             f"**{sum(r.endswith('| ok |') for r in rows)} of {len(rows)} checks pass; {sum('known erratum' in r for r in rows)} known errata; "
             f"{bad} unexplained failures.**", "",
             "| quantity | file | quoted | table value | status |", "|---|---|---|---:|---|", *rows, "",
             "## Numbers with three or more decimals that match no table value", "",
             "These are not errors by themselves: they were computed outside the result tables (diagnostics, or values derived from the "
             "distance matrices). Each is listed so it can be checked by hand.", ""]
    lines += (["| location | number |", "|---|---:|", *untraced] if untraced else ["None."])
    atomic_write_text(root / "reports" / "number_check.md", "\n".join(lines) + "\n")
    print(f"[verify] {len(rows)} checks, {bad} unexplained failures; {len(untraced)} untraced numbers -> reports/number_check.md")
    if bad:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
