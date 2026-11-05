from __future__ import annotations

import zipfile

import dendropy
import numpy as np
import pandas as pd

from lineage_charisma import condensed as cd
from lineage_charisma.config import base_parser, load_config
from lineage_charisma.distances import load_matrix, order_hash, phylo_distance_matrix, phylo_matrix_path, text_matrix_path
from lineage_charisma.embed import read_species_order
from lineage_charisma.h2 import PRIMARY_CONTROLS, decision, share_of_raw
from lineage_charisma.io_utils import read_csv, read_json, slugify_species, write_csv, write_json
from lineage_charisma.phylo import prune_tree, relabel_tips
from lineage_charisma.posterior import N_TREES, posterior_label, select_trees

# docs/preregistration_amendment_fullscale.md, section 8: the first 20 of the 100 trees drawn in A1.3, 999 permutations
N_FULLSCALE_TREES = 20
PERMS = 999
CONTROL_FILES = {"ecology": "ecology.npy", "length": "attention_length.npy"}


def tree_distances(cfg, zf, name, order, tips, relabel, cache_dir):
    """Condensed patristic distances of one posterior tree over `order`, computed once and saved to disk."""
    path = cache_dir / f"{name.rsplit('/', 1)[-1].removesuffix('.tre')}.npy"
    key_path = path.with_suffix(".json")
    if path.exists() and key_path.exists() and read_json(key_path).get("species_order_sha256") == order_hash(order):
        return np.load(path)
    tree = dendropy.Tree.get(data=zf.read(name).decode(), schema="newick", preserve_underscores=True, rooting="default-rooted")
    matrix = phylo_distance_matrix(relabel_tips(prune_tree(tree, tips), relabel), order)
    values = cd.condensed(matrix)
    np.save(path, values)
    write_json(key_path, {"species_order_sha256": order_hash(order), "tree_file": name, "n_species": len(order), "max_distance_my": float(values.max())})
    return values


def main() -> None:
    parser = base_parser("Full scale: H1 and the primary H2 model on 20 posterior trees (999 permutations each).")
    parser.add_argument("--trees-only", action="store_true", help="only prune the trees and save their distances to disk")
    args = parser.parse_args()
    cfg = load_config(args.config, args.clade or "mammalia")
    scfg = cfg["stats"]
    seed, method, alpha = int(scfg["seed"]), scfg["method"], float(scfg["alpha"])
    out_dir, order_path = cfg.matrices_dir(), cfg.species_order_path()
    order = read_species_order(order_path)
    matched = read_csv(cfg.matched_path()).set_index("species").loc[order]
    tips = matched["tree_tip"].map(posterior_label)
    relabel = dict(zip(tips, map(slugify_species, order)))
    archive = cfg.tree_path().parent / cfg["sources"]["tree"]["posterior_filename"]
    cache_dir = cfg.interim / "posterior"
    cache_dir.mkdir(parents=True, exist_ok=True)
    if args.trees_only:
        with zipfile.ZipFile(archive) as zf:
            names = sorted(x for x in zf.namelist() if x.endswith(".tre"))
            for draw, index in enumerate(select_trees(names, N_TREES, seed)[:N_FULLSCALE_TREES], 1):
                tree_distances(cfg, zf, names[index], order, tips, relabel, cache_dir)
                print(f"[posterior] saved tree {draw}/{N_FULLSCALE_TREES} ({index})", flush=True)
        return
    species = read_species_order(out_dir / "ecology_species.txt")
    primary = cfg.primary_combo
    text = load_matrix(text_matrix_path(out_dir, *primary), order_path)[0]
    mcc = load_matrix(phylo_matrix_path(out_dir), order_path)[0]
    p_text = cd.prepare(text, method)
    p_text_h2 = cd.prepare(cd.subset(text, order, species), method)
    controls = [cd.prepare(load_matrix(out_dir / CONTROL_FILES["ecology"], out_dir / "ecology_species.txt")[0], method),
                cd.prepare(cd.subset(load_matrix(out_dir / CONTROL_FILES["length"], order_path)[0], order, species), method)]
    mcc_vector = cd.condensed(mcc)
    h2_index = np.array([i for i, s in enumerate(order) if s in set(species)])
    h2_rows, h2_cols = np.triu_indices(len(species), k=1)
    del text, mcc

    n = len(order)
    rows = []
    with zipfile.ZipFile(archive) as zf:
        names = sorted(x for x in zf.namelist() if x.endswith(".tre"))
        chosen = select_trees(names, N_TREES, seed)[:N_FULLSCALE_TREES]
        for draw, index in enumerate(chosen, 1):
            values = tree_distances(cfg, zf, names[index], order, tips, relabel, cache_dir)
            p_phylo = cd.prepare_vector(values, n, method)
            h1 = cd.mantel(p_text, p_phylo, PERMS, seed)
            square = np.zeros((n, n))
            square[np.triu_indices(n, k=1)] = values
            sub = square[np.ix_(h2_index, h2_index)][h2_rows, h2_cols]
            del square
            p_phylo_h2 = cd.prepare_vector(sub, len(species), method)
            null = cd.null_correlations(p_phylo_h2, [p_text_h2, *controls], PERMS, seed)
            raw = cd.mantel(p_text_h2, p_phylo_h2, seed=seed, null=null[:, 0])
            part = cd.partial_mantel(p_phylo_h2, p_text_h2, controls, seed=seed, null_corr=null)
            rows.append({"draw": draw, "tree_index": index, "tree_file": names[index].rsplit("/", 1)[-1], "method": method, "h1_r": h1.r, "h1_p": h1.p,
                         "raw_r": raw.r, "raw_p": raw.p, "partial_r": part.r, "partial_p": part.p, "share_of_raw": share_of_raw(raw.r, part.r),
                         "change": raw.r - part.r, "n_species": n, "h2_n_species": len(species), "max_distance_my": float(values.max()),
                         "spearman_with_mcc": cd.correlation(cd.prepare_vector(mcc_vector, n), cd.prepare_vector(values, n)),
                         **decision(raw.r, part.r, part.p, alpha)})
            print(f"[posterior] tree {draw}/{N_FULLSCALE_TREES} ({index}): H1 r {h1.r:+.4f} (p {h1.p:.4f})  raw {raw.r:+.4f}  partial {part.r:+.4f} "
                  f"(p {part.p:.4f})", flush=True)
    table = pd.DataFrame(rows)
    write_csv(table, cfg.tables_dir() / "posterior.csv", permutations=PERMS, seed=seed, alternative="greater", n_trees=N_FULLSCALE_TREES,
              selection="first 20 of numpy.random.default_rng(20261004).choice(N, size=100, replace=False)", controls=list(PRIMARY_CONTROLS))

    summary = []
    for column in ("h1_r", "raw_r", "partial_r", "share_of_raw", "change"):
        v = table[column].dropna()
        summary.append({"statistic": column, "median": v.median(), "p2.5": v.quantile(0.025), "p97.5": v.quantile(0.975), "min": v.min(), "max": v.max(), "n_trees": len(v)})
    summary = pd.DataFrame(summary)
    shares = {"h1_p_below_alpha": float((table["h1_p"] < alpha).mean()), "rule_a_holds": float(table["survives"].mean()),
              "rule_b_holds": float(table["confounded"].fillna(False).astype(bool).mean())}
    write_csv(summary, cfg.tables_dir() / "posterior_summary.csv", n_trees=N_FULLSCALE_TREES, permutations=PERMS, **shares,
              note="with 20 trees the 2.5th and 97.5th percentiles are interpolated between the two most extreme trees")
    print(summary.to_string(index=False))
    print(f"[posterior] share of trees: H1 p < {alpha} {shares['h1_p_below_alpha']:.0%}, rule (a) {shares['rule_a_holds']:.0%}, rule (b) {shares['rule_b_holds']:.0%}")


if __name__ == "__main__":
    main()
