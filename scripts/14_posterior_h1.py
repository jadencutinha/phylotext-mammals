from __future__ import annotations

import zipfile

import dendropy
import numpy as np
import pandas as pd

from lineage_charisma.config import base_parser, load_config
from lineage_charisma.distances import load_matrix, phylo_distance_matrix, phylo_matrix_path, text_matrix_path, validate_distance_matrix
from lineage_charisma.io_utils import is_cached, read_csv, slugify_species, write_csv
from lineage_charisma.phylo import prune_tree, relabel_tips
from lineage_charisma.plotstyle import BLUE, HALO, INK, SURFACE, apply_style, combo_label, plt, save
from lineage_charisma.stats import mantel, upper_triangle

N_TREES = 100   # fixed in docs/preregistration.md, amendment A1.3


def posterior_label(mcc_tip: str) -> str:
    """Posterior trees label tips Genus_species; the MCC tree appends _FAMILY_ORDER to the same names."""
    return "_".join(mcc_tip.split("_")[:2])


def select_trees(names: list[str], n_trees: int, seed: int) -> list[int]:
    """The pre-registered draw: indices into the name-sorted tree files, without replacement."""
    return [int(i) for i in np.random.default_rng(seed).choice(len(names), size=n_trees, replace=False)]


def plot_posterior(table: pd.DataFrame, mcc_r: float, label: str, path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 3.8), constrained_layout=True)
    ax.hist(table["r"], bins=20, color=BLUE, edgecolor=SURFACE, linewidth=0.4)
    ax.axvline(mcc_r, color=INK, linewidth=1.6)
    right = mcc_r < table["r"].median()
    ax.annotate(f"MCC tree, r = {mcc_r:.3f}", xy=(mcc_r, 0.92), xycoords=("data", "axes fraction"), xytext=(6 if right else -6, 0),
                textcoords="offset points", ha="left" if right else "right", va="center", color=INK, path_effects=HALO)
    ax.axvline(0, color=INK, linewidth=0.6, linestyle=":")
    ax.set_xlabel("Spearman Mantel r, text distance against phylogenetic distance")
    ax.set_ylabel("posterior trees")
    ax.set_title(f"H1 on {len(table)} posterior trees ({label})")
    ax.grid(True, axis="y")
    ax.set_axisbelow(True)
    save(fig, path)


def main() -> None:
    args = base_parser("Robustness: repeat the primary H1 test on a pre-registered sample of posterior trees.").parse_args()
    cfg = load_config(args.config)
    apply_style()
    scfg, src = cfg["stats"], cfg["sources"]["tree"]
    perms, seed, method = int(scfg["permutations"]), int(scfg["seed"]), scfg["method"]
    archive = cfg.tree_path().parent / src["posterior_filename"]
    if not is_cached(archive):
        raise SystemExit(f"{archive.name} is not downloaded; run make posterior (about 1.2 GB) first")

    primary, order_path = cfg.primary_combo, cfg.species_order_path()
    text, order = load_matrix(text_matrix_path(cfg.matrices_dir(), *primary), order_path)
    mcc, _ = load_matrix(phylo_matrix_path(cfg.matrices_dir()), order_path)
    matched = read_csv(cfg.matched_path()).set_index("species").loc[order]
    tips = matched["tree_tip"].map(posterior_label)
    relabel = dict(zip(tips, map(slugify_species, order)))
    mcc_result = mantel(text, mcc, method=method, permutations=perms, seed=seed)

    rows = []
    with zipfile.ZipFile(archive) as zf:
        names = sorted(n for n in zf.namelist() if n.endswith(".tre"))
        chosen = select_trees(names, N_TREES, seed)
        print(f"[posterior] {len(names)} trees in {archive.name}; {N_TREES} drawn with seed {seed}")
        for k, index in enumerate(chosen, 1):
            tree = dendropy.Tree.get(data=zf.read(names[index]).decode(), schema="newick", preserve_underscores=True, rooting="default-rooted")
            phylo = phylo_distance_matrix(relabel_tips(prune_tree(tree, tips), relabel), order)
            validate_distance_matrix(phylo, names[index])
            res = mantel(text, phylo, method=method, permutations=perms, seed=seed)
            agreement = mantel(phylo, mcc, method="spearman", permutations=0).r
            rows.append({"draw": k, "tree_index": index, "tree_file": names[index].rsplit("/", 1)[-1], "method": method, "r": res.r, "p": res.p,
                         "n_pairs": res.n_pairs, "max_distance_my": float(upper_triangle(phylo).max()), "spearman_with_mcc": agreement})
            if k % 10 == 0:
                print(f"[posterior] {k}/{N_TREES} trees, r so far {min(r['r'] for r in rows):+.4f} to {max(r['r'] for r in rows):+.4f}")

    table = pd.DataFrame(rows)
    write_csv(table, cfg.tables_dir() / "h1_posterior.csv", permutations=perms, seed=seed, selection_seed=seed, n_trees=N_TREES, alternative="greater",
              archive=archive.name, text_matrix=combo_label(*primary), mcc_r=mcc_result.r, mcc_p=mcc_result.p)
    plot_posterior(table, mcc_result.r, combo_label(*primary), cfg.week4_figures_dir() / "h1_posterior.png")
    q = table["r"].quantile([0.025, 0.5, 0.975]).to_numpy()
    print(f"[posterior] r: median {q[1]:.4f}, 2.5th to 97.5th percentile {q[0]:.4f} to {q[2]:.4f}, min {table['r'].min():.4f}, max {table['r'].max():.4f}; "
          f"MCC tree r = {mcc_result.r:.4f}; {int((table['p'] < float(scfg['alpha'])).sum())} of {len(table)} trees have p < {scfg['alpha']}")


if __name__ == "__main__":
    main()
