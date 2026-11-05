from __future__ import annotations

import platform
import time

import numpy as np
import pandas as pd

from lineage_charisma import condensed as cd
from lineage_charisma.config import base_parser, load_config
from lineage_charisma.distances import load_matrix, phylo_matrix_path, text_matrix_path
from lineage_charisma.embed import read_species_order
from lineage_charisma.h2 import CONFOUNDED_SHARE, MODELS, MRM_PREDICTORS, PRIMARY_CONTROLS, decision, share_of_raw
from lineage_charisma.io_utils import atomic_write_text, read_csv, write_csv
from lineage_charisma.plotstyle import combo_label

# docs/preregistration_amendment_fullscale.md, sections 5, 7 and 8. Tests that permute the same matrix on the same species
# with the same seed share one pass over the permutations: with a fresh generator per test their permutations are identical,
# and 999 permutations are the first 999 of 9,999 (tests/test_condensed.py), so every result equals a separate test.
ROBUST_PERMS = 999
MIN_ORDER_SIZE = 50
ECOLOGY_VARIANTS = {"ecology": "ecology.npy", "ecology_unweighted": "ecology_unweighted.npy", "ecology_no_volant": "ecology_no_volant.npy"}
ATTENTION = {"length": "attention_length.npy", "pageviews": "attention_pageviews.npy"}


def timed(label: str, fn, *args, **kw):
    started = time.time()
    out = fn(*args, **kw)
    print(f"[time] {label}: {time.time() - started:.1f}s", flush=True)
    return out


def load_all(cfg):
    out_dir, order_path = cfg.matrices_dir(), cfg.species_order_path()
    order = read_species_order(order_path)
    texts = {combo: load_matrix(text_matrix_path(out_dir, *combo), order_path)[0] for combo in cfg.combos()}
    phylo = load_matrix(phylo_matrix_path(out_dir), order_path)[0]
    species = read_species_order(out_dir / "ecology_species.txt")
    if species != [s for s in order if s in set(species)]:
        raise SystemExit("ecology_species.txt is not a subset of species_order.txt in the same order")
    controls = {name: load_matrix(out_dir / f, out_dir / "ecology_species.txt")[0] for name, f in ECOLOGY_VARIANTS.items()}
    for name, f in ATTENTION.items():
        controls[name] = cd.subset(load_matrix(out_dir / f, order_path)[0], order, species)
    return order, species, texts, phylo, controls


def h2_rows(analysis, combo, ecology, controls_used, raw, part, alpha, perms):
    return {"analysis": analysis, "model": combo[0], "rule": combo[1], "mask": combo[2], "ecology_matrix": ecology, "tree": "MCC", "method": raw.method,
            "controls": " + ".join(controls_used), "raw_r": raw.r, "raw_p": raw.p, "partial_r": part.r, "partial_p": part.p,
            "share_of_raw": share_of_raw(raw.r, part.r), "change": raw.r - part.r, "n_species": part.n, "n_pairs": part.n_pairs,
            "permutations": perms, **decision(raw.r, part.r, part.p, alpha)}


def mrm_table(fit, analysis, combo, ecology):
    return [{"analysis": analysis, "model": combo[0], "rule": combo[1], "mask": combo[2], "ecology_matrix": ecology, "predictor": name,
             "coef_std": fit.coef[j], "coef_raw": fit.coef_raw[j], "t": fit.t[j], "p": fit.p[j], "r2": fit.r2, "p_r2": fit.p_r2,
             "n_species": fit.n, "n_pairs": fit.n_pairs, "permutations": fit.permutations} for j, name in enumerate(fit.names)]


def benchmark(cfg, perms, seed, workers) -> None:
    """Section 5 / Step 3: time 100 permutations of each pass on the real matrices and project the full run."""
    order, species, texts, phylo, controls = timed("load matrices", load_all, cfg)
    primary = cfg.primary_combo
    started = time.time()
    p_phylo = cd.prepare(phylo, "spearman", "phylo")
    rank_s = time.time() - started
    p_text = cd.prepare(texts[primary], "spearman", "text")
    n_bench = 100
    rows = []
    for k in (1, 4, 9):
        started = time.time()
        cd.null_correlations(p_phylo, [p_text] * k, n_bench, seed, workers=workers)
        rows.append((k, (time.time() - started) / n_bench))
    per = {k: s for k, s in rows}
    n, m = len(order), len(species)
    scale = (m * (m - 1)) / (n * (n - 1))
    plan = [  # (pass, species set, fixed vectors, permutations)
        ("H1: phylogeny moving, 4 text matrices", n, 4, perms),
        ("H1 Pearson", n, 1, ROBUST_PERMS),
        ("H2: phylogeny moving, 4 texts + 5 controls", m, 9, perms),
        ("MRM primary: text moving, 6 predictors", m, 9, perms),
        ("MRM, 3 text variants", m, 4, 3 * ROBUST_PERMS),
        ("posterior: 20 trees x (H1 + H2)", None, None, None),
    ]
    lines = ["# Week 6 performance: full-scale permutation tests", "",
             f"Measured by `scripts/22_fullscale_tests.py --benchmark` on {platform.machine()} {platform.system()} "
             f"({platform.processor() or 'unknown CPU'}), {workers} worker threads, on the real full-scale matrices: "
             f"{n:,} species in the frozen set ({n * (n - 1) // 2:,} pairs) and {m:,} H2 species ({m * (m - 1) // 2:,} pairs).", "",
             "## Design (amendment 3, section 5)", "",
             "- Each matrix is ranked once (Spearman, average ranks) and stored as a centred unit-norm condensed vector (float64 for the observed statistic).",
             "- The permuted matrix is held as a float32 square; each permutation gathers its condensed vector through an int32 index array and multiplies it by "
             "every fixed vector at once (one float32 matrix-vector product). Permutations are drawn in chunks of 256 from one generator, identical to the "
             "unchunked draw, and run on a thread pool; each statistic is computed on its own, so results do not depend on chunk size or thread count.",
             "- Partial Mantel and MRM statistics are computed from correlations alone (the residualization is algebraic), so the controls are prepared once.",
             "- Tests that permute the same matrix on the same species share one pass; the 999-permutation robustness tests use the first 999 permutations of "
             "the 9,999-permutation pass, which is exactly the separate test with the same seed.",
             "- Patristic distances come from a vectorized version of dendropy's distance calculation (bit-identical to dendropy; 19 s for 5,583 tips instead of hours), "
             "and each posterior tree's distances are saved to disk one tree at a time.", "",
             "## Measurements", "", "| quantity | value |", "|---|---:|",
             f"| ranking one {n:,}-species matrix and standardizing it | {rank_s:.1f} s |"]
    lines += [f"| one permutation, {k} fixed vector{'s' if k > 1 else ''}, {n:,} species | {s * 1000:.1f} ms |" for k, s in rows]
    total = 0.0
    lines += ["", "## Projection for the registered run", "", "| pass | permutations | projected time |", "|---|---:|---:|"]
    for name, size, k, p in plan:
        if size is None:
            secs = 20 * ROBUST_PERMS * (per[1] + per[4] * scale)
            lines.append(f"| {name} | 20 x 2 x {ROBUST_PERMS} | {secs / 60:.1f} min |")
        else:
            secs = p * per[min(per, key=lambda c: abs(c - k))] * (scale if size == m else 1.0)
            lines.append(f"| {name} | {p:,} | {secs / 60:.1f} min |")
        total += secs
    lines += ["", f"Total projected for the pooled passes above: about {total / 60:.0f} minutes, plus per-order tests (smaller n) and the exploratory "
              f"pair-subset tests, which re-rank the subset in every permutation. No GPU is needed for the statistics. Embeddings ran on the local "
              f"Apple GPU (MPS); `scripts/slurm/embed_adroit.sh` runs the same step on Princeton Adroit if a CUDA node is preferred."]
    path = cfg.root / "reports" / "week6_performance.md"
    atomic_write_text(path, "\n".join(lines) + "\n")
    print("\n".join(lines))


def main() -> None:
    parser = base_parser("Full scale: H1, H2 (primary, secondary, MRM), robustness variants and per-order tests on the frozen species set.")
    parser.add_argument("--benchmark", action="store_true", help="time 100 permutations per pass, project the full run, write reports/week6_performance.md, and stop")
    parser.add_argument("--workers", type=int, default=None)
    args = parser.parse_args()
    cfg = load_config(args.config, args.clade or "mammalia")
    scfg = cfg["stats"]
    perms, seed, method, alpha = int(scfg["permutations"]), int(scfg["seed"]), scfg["method"], float(scfg["alpha"])
    workers = args.workers or cd.default_workers()
    if args.benchmark:
        benchmark(cfg, perms, seed, workers)
        return
    kw = dict(workers=workers)
    combos, primary = cfg.combos(), cfg.primary_combo
    variants = [c for c in combos if c != primary]
    order, species, texts, phylo, controls = timed("load matrices", load_all, cfg)
    mdd_order = read_csv(cfg.mdd_species_path()).set_index("species")["order"]
    tables = cfg.tables_dir()
    print(f"[full] {len(order)} species, {len(species)} H2 species; {perms:,} permutations (robustness {ROBUST_PERMS}), seed {seed}, {workers} threads")

    # H1 on the frozen set
    P = {c: timed(f"rank text {combo_label(*c)}", cd.prepare, texts[c], method) for c in combos}
    phylo_s = timed("rank phylogeny", cd.prepare, phylo, method)
    null = timed("H1 pass", cd.null_correlations, phylo_s, [P[c] for c in combos], perms, seed, **kw)
    h1 = []
    for j, c in enumerate(combos):
        n_used = perms if c == primary else ROBUST_PERMS
        res = cd.mantel(P[c], phylo_s, seed=seed, null=null[:n_used, j])
        h1.append({"analysis": "primary" if c == primary else "variant", "model": c[0], "rule": c[1], "mask": c[2], "method": method, "r": res.r, "p": res.p,
                   "n_species": res.n, "n_pairs": res.n_pairs, "permutations": n_used, "is_primary": c == primary})
    np.save(cfg.interim / "h1_primary_null.npy", null[:, 0])
    pearson_text, pearson_phylo = cd.prepare(texts[primary], "pearson"), cd.prepare(phylo, "pearson")
    res = timed("H1 Pearson", cd.mantel, pearson_text, pearson_phylo, ROBUST_PERMS, seed, **kw)
    h1.append({"analysis": "pearson", "model": primary[0], "rule": primary[1], "mask": primary[2], "method": "pearson", "r": res.r, "p": res.p,
               "n_species": res.n, "n_pairs": res.n_pairs, "permutations": ROBUST_PERMS, "is_primary": False})
    h1 = pd.DataFrame(h1)
    write_csv(h1, tables / "h1_mantel.csv", seed=seed, alternative="greater", permuted="phylogeny", n_species=len(order), plan="docs/preregistration_amendment_fullscale.md")
    print(h1.to_string(index=False), flush=True)
    del P, phylo_s, pearson_text, pearson_phylo

    # H2 on the ecology species: one pass with the phylogeny moving, one with the primary text moving
    text_h2 = {c: cd.subset(texts[c], order, species) for c in combos}
    Q = {c: cd.prepare(text_h2[c], method) for c in combos}
    C = {k: cd.prepare(v, method, k) for k, v in controls.items()}
    phylo_h2 = cd.prepare(cd.subset(phylo, order, species), method)
    fixed_names = [*combos, *C]
    fixed = [*(Q[c] for c in combos), *C.values()]
    null = timed("H2 pass", cd.null_correlations, phylo_h2, fixed, perms, seed, **kw)
    col = {name: j for j, name in enumerate(fixed_names)}
    partial = []

    def run(analysis, combo, ctrl_names, ecology="ecology", n_used=perms):
        names = [ecology if c == "ecology" else c for c in ctrl_names]
        raw = cd.mantel(Q[combo], phylo_h2, seed=seed, null=null[:n_used, col[combo]])
        part = cd.partial_mantel(phylo_h2, Q[combo], [C[c] for c in names], seed=seed, null_corr=null[:n_used, [col[combo], *(col[c] for c in names)]])
        row = h2_rows(analysis, combo, ecology, ctrl_names, raw, part, alpha, n_used)
        partial.append(row)
        print(f"[h2] {analysis:22s} {combo_label(*combo):42s} {ecology:18s} given {row['controls']:30s} raw {raw.r:+.4f} partial {part.r:+.4f} "
              f"p {part.p:.4f} share {row['share_of_raw']:.3f}", flush=True)
        return part

    head = run("primary", primary, MODELS["primary"])
    np.save(cfg.interim / "h2_primary_null.npy", head.null)
    for model_id in ("S1", "S2", "S3", "S4"):
        run(f"secondary {model_id}", primary, MODELS[model_id])
    for c in variants:
        run("variant", c, PRIMARY_CONTROLS, n_used=ROBUST_PERMS)
    run("unweighted_gower", primary, PRIMARY_CONTROLS, ecology="ecology_unweighted", n_used=ROBUST_PERMS)
    run("ecology_without_flight", primary, PRIMARY_CONTROLS, ecology="ecology_no_volant", n_used=ROBUST_PERMS)
    partial = pd.DataFrame(partial)
    partial["is_primary"] = partial["analysis"] == "primary"

    mrm_preds = ["phylogeny", "ecology", "length", "pageviews"]
    D_fixed = {"phylogeny": phylo_h2, **C}
    null = timed("MRM pass (primary text moving)", cd.null_correlations, Q[primary], list(D_fixed.values()), perms, seed, **kw)
    dcol = {name: j for j, name in enumerate(D_fixed)}
    regression = []
    for analysis, ecology, n_used in (("primary", "ecology", perms), ("unweighted_gower", "ecology_unweighted", ROBUST_PERMS),
                                      ("ecology_without_flight", "ecology_no_volant", ROBUST_PERMS)):
        names = ["phylogeny", ecology, "length", "pageviews"]
        fit = cd.mrm(Q[primary], {p: D_fixed[n] for p, n in zip(mrm_preds, names)}, seed=seed, null_corr=null[:n_used, [dcol[n] for n in names]])
        regression += mrm_table(fit, analysis, primary, ecology)
    for c in variants:
        fit = timed(f"MRM {combo_label(*c)}", cd.mrm, Q[c], {p: D_fixed[p] for p in mrm_preds}, ROBUST_PERMS, seed, **kw)
        regression += mrm_table(fit, "variant", c, "ecology")
    regression = pd.DataFrame(regression)
    regression["is_primary"] = regression["analysis"] == "primary"
    meta = dict(seed=seed, alpha=alpha, text_matrix=combo_label(*primary), plan="docs/preregistration_amendment_fullscale.md", confounded_share=CONFOUNDED_SHARE)
    write_csv(partial, tables / "h2_partial_mantel.csv", alternative="greater", permuted="phylogeny", **meta)
    write_csv(regression, tables / "h2_mrm.csv", alternative="two-sided on pseudo-t", permuted="text", predictors=list(MRM_PREDICTORS), **meta)
    print(regression[regression["is_primary"]].to_string(index=False), flush=True)
    del Q, C, phylo_h2, null

    # per-order breakdowns (section 7): distances from the full-scale matrices as they are
    rows = []
    in_h2 = set(species)
    counts = mdd_order.loc[order].value_counts()
    for name in counts[counts >= MIN_ORDER_SIZE].index:
        members = [s for s in order if mdd_order[s] == name]
        h2_members = [s for s in members if s in in_h2]
        t = cd.prepare(cd.subset(texts[primary], order, members), method)
        ph = cd.prepare(cd.subset(phylo, order, members), method)
        h1o = cd.mantel(t, ph, perms, seed, **kw)
        th = cd.prepare(cd.subset(texts[primary], order, h2_members), method)
        phh = cd.prepare(cd.subset(phylo, order, h2_members), method)
        ctrl = [cd.prepare(cd.subset(controls[c], species, h2_members), method) for c in PRIMARY_CONTROLS]
        nn = cd.null_correlations(phh, [th, *ctrl], perms, seed, **kw)
        raw = cd.mantel(th, phh, seed=seed, null=nn[:, 0])
        part = cd.partial_mantel(phh, th, ctrl, seed=seed, null_corr=nn)
        rows.append({"order": name, "n_species": h1o.n, "h1_r": h1o.r, "h1_p": h1o.p, "h2_n_species": part.n, "raw_r": raw.r, "raw_p": raw.p,
                     "partial_r": part.r, "partial_p": part.p, "share_of_raw": share_of_raw(raw.r, part.r), "permutations": perms})
        print(f"[order] {name:16s} n={h1o.n:5d} H1 r {h1o.r:+.4f} (p {h1o.p:.4f})  H2 n={part.n:5d} raw {raw.r:+.4f} partial {part.r:+.4f} (p {part.p:.4f})", flush=True)
    write_csv(pd.DataFrame(rows), tables / "per_order.csv", status="descriptive", seed=seed, alternative="greater", min_order_size=MIN_ORDER_SIZE,
              controls=list(PRIMARY_CONTROLS), text_matrix=combo_label(*primary))

    h = partial.set_index("analysis").loc["primary"]
    rules = decision(h["raw_r"], h["partial_r"], h["partial_p"], alpha)
    top = h1[h1["is_primary"]].iloc[0]
    print(f"[full] H1: r = {top['r']:.4f}, p = {top['p']:.4f} ({'supported' if top['r'] > 0 and top['p'] < alpha else 'not supported'})")
    print(f"[full] H2: raw r = {h['raw_r']:.4f}, partial r = {h['partial_r']:.4f} (p = {h['partial_p']:.4f}), share {h['share_of_raw']:.1%}; "
          f"rule (a) {rules['survives']}, rule (b) {rules['confounded']}")


if __name__ == "__main__":
    main()
