from __future__ import annotations

import numpy as np
import pandas as pd

from lineage_charisma import condensed as cd
from lineage_charisma.config import base_parser, load_config
from lineage_charisma.corpus import read_corpus
from lineage_charisma.distances import load_matrix, phylo_matrix_path, text_matrix_path
from lineage_charisma.embed import read_species_order
from lineage_charisma.io_utils import read_csv, write_csv
from lineage_charisma.plotstyle import combo_label

# docs/preregistration_amendment_secondary.md, section D: a disclosed POST-HOC analysis, not a registered H1/H2 outcome
N_BINS = 10
CHECK_BINS = (5, 20)


def length_bins(n_tokens: pd.Series, n_bins: int) -> np.ndarray:
    """Equal-size bins by token count; ties broken by position in species_order.txt (n_tokens is indexed in that order)."""
    return pd.qcut(n_tokens.rank(method="first"), n_bins, labels=False).to_numpy()


def bin_members(bins: np.ndarray, n_bins: int) -> list[np.ndarray]:
    return [np.flatnonzero(bins == b) for b in range(n_bins)]


def observed(text, phylo, members):
    rs, weights = [], []
    for idx in members:
        sub = np.ix_(idx, idx)
        rs.append(cd.correlation(cd.prepare(text[sub], "spearman", check=False), cd.prepare(phylo[sub], "spearman", check=False)))
        weights.append(len(idx) * (len(idx) - 1) // 2)
    w = np.array(weights) / sum(weights)
    return float(np.dot(w, rs)), np.array(rs), np.array(weights)


def main() -> None:
    parser = base_parser("Post-hoc (amendment 4, section D): length-matched text-phylogeny correlation on the full-scale frozen set.")
    args = parser.parse_args()
    cfg = load_config(args.config, args.clade or "mammalia")
    perms, seed = int(cfg["stats"]["permutations"]), int(cfg["stats"]["seed"])
    out_dir, order_path = cfg.matrices_dir(), cfg.species_order_path()
    primary = cfg.primary_combo
    order = read_species_order(order_path)
    text = load_matrix(text_matrix_path(out_dir, *primary), order_path)[0]
    phylo = load_matrix(phylo_matrix_path(out_dir), order_path)[0]
    n_tokens = read_corpus(cfg.corpus_path()).set_index("species_id").loc[order, "n_tokens"].astype(int)
    mdd_order = read_csv(cfg.mdd_species_path()).set_index("species").loc[order, "order"].to_numpy()

    bins = length_bins(n_tokens, N_BINS)
    members = bin_members(bins, N_BINS)
    r_lm, r_b, pairs = observed(text, phylo, members)
    w = pairs / pairs.sum()

    # null: species permuted within bins; one generator, bins permuted in order 0..9 within each permutation
    rng = np.random.default_rng(seed)
    within = [np.empty((perms, len(idx)), dtype=np.int32) for idx in members]
    for k in range(perms):
        for b, idx in enumerate(members):
            within[b][k] = rng.permutation(len(idx))
    null = np.zeros(perms)
    for b, idx in enumerate(members):
        sub = np.ix_(idx, idx)
        t, ph = cd.prepare(text[sub], "spearman", check=False), cd.prepare(phylo[sub], "spearman", check=False)
        chunks = (within[b][s:s + 512] for s in range(0, perms, 512))
        null += w[b] * cd.correlations_for_permutations(ph, [t], chunks, perms)[:, 0]
    p = (1 + int((np.abs(null) >= abs(r_lm) - 1e-12).sum())) / (1 + perms)

    rows = []
    for b, idx in enumerate(members):
        sub = np.ix_(idx, idx)
        tv = cd.condensed(text[sub])
        o = mdd_order[idx]
        rr, cc = np.triu_indices(len(idx), k=1)
        same = o[rr] == o[cc]
        rows.append({"status": "post-hoc", "bin": b, "n_species": len(idx), "n_pairs": int(pairs[b]), "weight": w[b],
                     "min_tokens": int(n_tokens.iloc[idx].min()), "max_tokens": int(n_tokens.iloc[idx].max()), "r_spearman": r_b[b],
                     "n_same_order_pairs": int(same.sum()), "mean_text_same_order": tv[same].mean() if same.any() else np.nan,
                     "mean_text_different_order": tv[~same].mean(), "different_minus_same": tv[~same].mean() - tv[same].mean() if same.any() else np.nan})
    per_bin = pd.DataFrame(rows)
    gap = per_bin.dropna(subset=["different_minus_same"])
    weighted_gap = float(np.average(gap["different_minus_same"], weights=gap["n_pairs"]))
    checks = {f"r_lm_{k}_bins": observed(text, phylo, bin_members(length_bins(n_tokens, k), k))[0] for k in CHECK_BINS}

    summary = pd.DataFrame([{"status": "post-hoc", "statistic": "pair-weighted mean within-decile Spearman r", "n_bins": N_BINS, "r_lm": r_lm,
                             "p_two_sided": p, "permutations": perms, "seed": seed, "n_species": len(order), "n_within_bin_pairs": int(pairs.sum()),
                             "null_mean": null.mean(), "null_sd": null.std(), "unmatched_h1_r": float(read_csv(cfg.tables_dir() / "h1_mantel.csv").query("is_primary")["r"].iloc[0]),
                             "weighted_different_minus_same": weighted_gap, **checks}])
    meta = dict(status="post-hoc, designed after the Week 6 results (docs/preregistration_amendment_secondary.md, section D)", text_matrix=combo_label(*primary),
                permuted="phylogeny, species permuted within length deciles", alternative="two-sided", seed=seed)
    write_csv(summary, cfg.tables_dir() / "posthoc_length_matched.csv", **meta)
    write_csv(per_bin, cfg.tables_dir() / "posthoc_length_matched_bins.csv", **meta)
    print(per_bin[["bin", "n_species", "min_tokens", "max_tokens", "r_spearman", "mean_text_same_order", "mean_text_different_order", "different_minus_same"]].to_string(index=False))
    print(f"[post-hoc D] r_LM = {r_lm:+.4f}, two-sided p = {p:.4f} ({perms:,} within-bin permutations; null mean {null.mean():+.4f}, sd {null.std():.4f}); "
          f"5 bins {checks['r_lm_5_bins']:+.4f}, 20 bins {checks['r_lm_20_bins']:+.4f}; weighted different-minus-same text distance {weighted_gap:+.4f}")


if __name__ == "__main__":
    main()
