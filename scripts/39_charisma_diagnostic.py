from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from lineage_charisma.charisma import charisma_distance
from lineage_charisma.config import base_parser, load_config
from lineage_charisma.corpus import read_corpus
from lineage_charisma.distances import load_matrix, text_matrix_path, validate_distance_matrix
from lineage_charisma.embed import read_species_order
from lineage_charisma.io_utils import read_csv, read_meta, write_csv
from lineage_charisma.stats import submatrix

# Diagnostic of the B1 sign (amendment 4, section B), requested after analysis B ran. Not pre-registered: correlations and
# means only, no p-values, nothing here changes B1-B3. Pairs are the 283 H2 species used by B1.
STATUS = "diagnostic, not pre-registered (requested after analysis B ran)"
STUB_LIKE = 100       # config.yaml matrices.sensitivity_min_tokens: the registered "short text" sensitivity threshold
N_BINS = 10


def partial_spearman(x: np.ndarray, y: np.ndarray, controls: list[np.ndarray]) -> float:
    """Spearman partial correlation of x and y given the controls (ranks, then linear residuals)."""
    def ranks(v):
        r = pd.Series(v).rank().to_numpy()
        return r - r.mean()
    Z = np.column_stack([np.ones(len(x)), *[ranks(c) for c in controls]])
    res = [ranks(v) - Z @ np.linalg.lstsq(Z, ranks(v), rcond=None)[0] for v in (x, y)]
    return float(res[0] @ res[1] / np.sqrt((res[0] @ res[0]) * (res[1] @ res[1])))


def main() -> None:
    args = base_parser("Diagnostic of the sign of B1 (charisma content vs text distance); not pre-registered.").parse_args()
    cfg = load_config(args.config, None)
    md, op = cfg.matrices_dir(), cfg.species_order_path()
    order = read_species_order(op)
    species = read_species_order(md / "ecology_species.txt")
    path = md / "charisma_content.npy"
    retained = read_meta(path)["dimensions"]
    ann = cfg.interim / "annotation"
    r1, r2 = (read_csv(ann / f"scores_run{k}.csv").set_index("species_id") for k in (1, 2))
    mean = ((r1 + r2) / 2).loc[order, retained]

    # 1. the saved matrix is a distance, and it is the registered one
    saved = load_matrix(path, op)[0]
    validate_distance_matrix(saved, path.name)
    rebuilt = charisma_distance(mean.to_numpy())
    raw_scores = mean.to_numpy()
    same = (raw_scores[:, None, :] == raw_scores[None, :, :]).all(-1)
    iu = np.triu_indices(len(order), 1)
    checks = {"max |saved - rebuilt|": float(np.abs(saved - rebuilt).max()),
              "max distance between species with identical scores": float(saved[same].max()),
              "Spearman r with unstandardized score difference": float(spearmanr(saved[iu], np.sqrt(((raw_scores[:, None] - raw_scores[None]) ** 2).sum(-1))[iu])[0])}

    corpus = read_corpus(cfg.corpus_path()).set_index("species_id")
    tokens = corpus.loc[species, "n_tokens"].to_numpy(dtype=float)
    text = submatrix(load_matrix(text_matrix_path(md, *cfg.primary_combo), op)[0], order, species)
    charisma = submatrix(saved, order, species)
    length = submatrix(load_matrix(md / "attention_length.npy", op)[0], order, species)
    pageviews = submatrix(load_matrix(md / "attention_pageviews.npy", op)[0], order, species)
    zero = (mean.loc[species].to_numpy() == 0).all(axis=1)
    i, j = np.triu_indices(len(species), 1)
    t, c, ln, pv = text[i, j], charisma[i, j], length[i, j], pageviews[i, j]
    level = (np.log(tokens[i]) + np.log(tokens[j])) / 2          # pair length level: mean log token count (not a registered control)

    # 2. pair classes
    classes = {"both zero on every retained dimension": zero[i] & zero[j], "one zero, one not": zero[i] != zero[j],
               "both non-zero": ~zero[i] & ~zero[j]}
    rows = [{"status": STATUS, "section": "pair class", "item": name, "n_pairs": int(mask.sum()), "mean_text_distance": t[mask].mean(),
             "mean_charisma_distance": c[mask].mean(), "median_tokens_of_pair_members": float(np.median(np.r_[tokens[i[mask]], tokens[j[mask]]]))}
            for name, mask in classes.items()]
    rows.append({"status": STATUS, "section": "species", "item": "all-zero species", "n_species": int(zero.sum()),
                 "median_tokens": float(np.median(tokens[zero]))})
    rows.append({"status": STATUS, "section": "species", "item": "other species", "n_species": int((~zero).sum()),
                 "median_tokens": float(np.median(tokens[~zero]))})

    # 3. what the correlation looks like once length level is held fixed or stub-like texts are removed
    def corr(name, mask=None, controls=None, **extra):
        m = np.ones(len(t), bool) if mask is None else mask
        r = spearmanr(t[m], c[m])[0] if not controls else partial_spearman(t[m], c[m], [x[m] for x in controls])
        rows.append({"status": STATUS, "section": "correlation", "item": name, "n_pairs": int(m.sum()), "spearman_r": float(r), **extra})

    corr("raw r(text, charisma), as registered (all 283 H2 species)")
    corr("partial | length difference + pageviews difference (B1 as registered)", controls=[ln, pv])
    for name, v in (("text", t), ("charisma", c)):
        rows.append({"status": STATUS, "section": "correlation", "item": f"r({name} distance, pair length level)", "n_pairs": len(t),
                     "spearman_r": float(spearmanr(v, level)[0])})
    corr("partial | length difference + pageviews difference + pair length level", controls=[ln, pv, level])
    long_enough = tokens >= STUB_LIKE
    keep = long_enough[i] & long_enough[j]
    corr(f"raw r, species with >= {STUB_LIKE} tokens only ({int(long_enough.sum())} species)", mask=keep)
    corr(f"partial | length + pageviews difference, species with >= {STUB_LIKE} tokens only", mask=keep, controls=[ln, pv])
    corr("raw r, pairs where neither species is all-zero", mask=classes["both non-zero"])

    bins = pd.qcut(pd.Series(tokens).rank(method="first"), N_BINS, labels=False).to_numpy()
    within = bins[i] == bins[j]
    rb, wb = [], []
    for b in range(N_BINS):
        m = within & (bins[i] == b)
        rb.append(spearmanr(t[m], c[m])[0])
        wb.append(m.sum())
    rows.append({"status": STATUS, "section": "correlation", "item": f"pair-weighted within-length-decile r ({N_BINS} equal bins, as post-hoc D)",
                 "n_pairs": int(within.sum()), "spearman_r": float(np.average(rb, weights=wb))})

    table = pd.DataFrame(rows)
    write_csv(table, cfg.root / "reports" / "tables" / "secondary" / "charisma_diagnostic.csv", status=STATUS, retained=retained,
              matrix_checks=checks, note="correlations and means only; no p-values; does not change B1-B3")
    print("matrix checks:", checks)
    print(table.drop(columns="status").to_string(index=False, float_format=lambda x: f"{x:.4f}"))


if __name__ == "__main__":
    main()
