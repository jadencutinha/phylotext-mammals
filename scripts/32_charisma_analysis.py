from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from lineage_charisma.charisma import AGREEMENT_FLOOR, DIMENSIONS, charisma_distance
from lineage_charisma.config import base_parser, load_config
from lineage_charisma.distances import load_matrix, phylo_matrix_path, save_matrix, text_matrix_path
from lineage_charisma.embed import read_species_order
from lineage_charisma.io_utils import read_csv, write_csv, write_json
from lineage_charisma.stats import mantel, mrm, partial_mantel, submatrix

# docs/preregistration_amendment_secondary.md, section B (secondary)


def agreement(a: pd.Series, b: pd.Series) -> float:
    if a.nunique() < 2 or b.nunique() < 2:
        return float("nan")      # undefined when a rater gives every species the same score: the dimension is dropped
    return float(spearmanr(a, b).statistic)


def main() -> None:
    args = base_parser("Charisma content (amendment 4, section B): agreement, drop rule, distance matrix and tests B1-B3.").parse_args()
    cfg = load_config(args.config, None)
    scfg = cfg["stats"]
    perms, seed, method = int(scfg["permutations"]), int(scfg["seed"]), scfg["method"]
    order = read_species_order(cfg.species_order_path())
    ann = cfg.interim / "annotation"
    r1, r2 = (read_csv(ann / f"scores_run{k}.csv").set_index("species_id") for k in (1, 2))
    valid = [s for s in order if r1.loc[s].notna().all() and r2.loc[s].notna().all()]
    mean = (r1.loc[valid, list(DIMENSIONS)] + r2.loc[valid, list(DIMENSIONS)]) / 2
    human = read_csv(cfg.interim / "annotation_validation.csv").merge(read_csv(cfg.interim / "annotation_validation_key.csv"), on="item").set_index("species_id")
    both = [s for s in human.index if s in mean.index]

    rows = []
    for d in DIMENSIONS:
        hm = agreement(human.loc[both, d], mean.loc[both, d])
        mm = agreement(r1.loc[valid, d], r2.loc[valid, d])
        keep = bool(np.isfinite(hm) and np.isfinite(mm) and hm >= AGREEMENT_FLOOR and mm >= AGREEMENT_FLOOR)
        dist = mean[d].round().value_counts().sort_index().to_dict()
        rows.append({"status": "secondary", "dimension": d, "author_vs_model": hm, "n_validation": len(both), "run1_vs_run2": mm, "n_species": len(valid),
                     "retained": keep, "mean_score": mean[d].mean(), "score_distribution_rounded_mean": str(dist)})
    agree = pd.DataFrame(rows)
    out = cfg.root / "reports" / "tables" / "secondary"
    write_csv(agree, out / "charisma_agreement.csv", floor=AGREEMENT_FLOOR, dropped_species=[s for s in order if s not in valid])
    print(agree.to_string(index=False))
    retained = agree.loc[agree["retained"], "dimension"].tolist()
    if not retained:
        write_json(out / "charisma_stopped.json", {"reason": "every dimension fell below the agreement floor", "floor": AGREEMENT_FLOOR})
        raise SystemExit("[charisma] all four dimensions dropped: analysis B stops with the agreement table (amendment 4, section B)")
    if len(valid) != len(order):
        raise SystemExit(f"[charisma] {len(order) - len(valid)} species have no valid scores; the registered matrix needs every pilot species or a written exclusion")

    dist = charisma_distance(mean.loc[order, retained].to_numpy())
    path = cfg.matrices_dir() / "charisma_content.npy"
    save_matrix(path, dist, order, dimensions=retained, metric="Euclidean on standardized mean scores", model="claude-sonnet-5-5")

    md, op = cfg.matrices_dir(), cfg.species_order_path()
    species = read_species_order(md / "ecology_species.txt")
    m = {"text": submatrix(load_matrix(text_matrix_path(md, *cfg.primary_combo), op)[0], order, species),
         "phylogeny": submatrix(load_matrix(phylo_matrix_path(md), op)[0], order, species),
         "ecology": load_matrix(md / "ecology.npy", md / "ecology_species.txt")[0],
         "length": submatrix(load_matrix(md / "attention_length.npy", op)[0], order, species),
         "pageviews": submatrix(load_matrix(md / "attention_pageviews.npy", op)[0], order, species),
         "charisma": submatrix(load_matrix(path, op)[0], order, species)}

    marg = [{"with": k, **{f: getattr(res, f) for f in ("r", "p")}} for k in ("text", "phylogeny", "length", "pageviews")
            for res in [mantel(m["charisma"], m[k], method=method, permutations=perms, seed=seed, alternative="two-sided")]]
    b1 = partial_mantel(m["charisma"], m["text"], [m["length"], m["pageviews"]], method=method, permutations=perms, seed=seed)
    b2 = partial_mantel(m["phylogeny"], m["text"], [m["ecology"], m["length"], m["charisma"]], method=method, permutations=perms, seed=seed)
    primary = partial_mantel(m["phylogeny"], m["text"], [m["ecology"], m["length"]], method=method, permutations=perms, seed=seed)
    b3 = mrm(m["text"], {k: m[k] for k in ("phylogeny", "ecology", "length", "pageviews", "charisma")}, method=method, permutations=perms, seed=seed)

    tests = pd.DataFrame([
        {"status": "secondary", "test": "B1 text ~ charisma | length + pageviews", "r": b1.r, "p": b1.p, "raw_r": b1.r_marginal, "n_species": b1.n,
         "explains_beyond_attention": bool(b1.r > 0 and b1.p < 0.05)},
        {"status": "secondary", "test": "B2 text ~ phylogeny | ecology + length + charisma", "r": b2.r, "p": b2.p, "raw_r": b2.r_marginal, "n_species": b2.n,
         "primary_h2_partial_r_same_run": primary.r, "change_from_primary": b2.r - primary.r},
    ])
    mrm_rows = pd.DataFrame([{"status": "secondary (descriptive)", "predictor": n, "coef_std": b3.coef[j], "t": b3.t[j], "p": b3.p[j], "r2": b3.r2,
                              "p_r2": b3.p_r2, "n_species": b3.n} for j, n in enumerate(b3.names)])
    meta = dict(status="secondary (amendment 4, section B)", retained=retained, permutations=perms, seed=seed, model="claude-sonnet-5-5")
    write_csv(tests, out / "charisma_tests.csv", **meta)
    write_csv(mrm_rows, out / "charisma_mrm.csv", **meta)
    write_csv(pd.DataFrame(marg), out / "charisma_marginal.csv", **meta, alternative="two-sided (descriptive)")
    print(tests.to_string(index=False))
    print(mrm_rows.to_string(index=False))


if __name__ == "__main__":
    main()
