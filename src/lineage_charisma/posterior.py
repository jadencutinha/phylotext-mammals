from __future__ import annotations

import zipfile
from typing import Sequence

import dendropy
import numpy as np
import pandas as pd

from lineage_charisma.config import Config
from lineage_charisma.distances import order_hash, phylo_distance_matrix, validate_distance_matrix
from lineage_charisma.io_utils import is_cached, read_csv, slugify_species
from lineage_charisma.phylo import prune_tree, relabel_tips

N_TREES = 100   # fixed in docs/preregistration.md, amendment A1.3


def posterior_label(mcc_tip: str) -> str:
    """Posterior trees label tips Genus_species; the MCC tree appends _FAMILY_ORDER to the same names."""
    return "_".join(mcc_tip.split("_")[:2])


def select_trees(names: list[str], n_trees: int, seed: int) -> list[int]:
    """The pre-registered draw: indices into the name-sorted tree files, without replacement."""
    return [int(i) for i in np.random.default_rng(seed).choice(len(names), size=n_trees, replace=False)]


def posterior_distance_matrices(cfg: Config, order: Sequence[str], n_trees: int = N_TREES, seed: int | None = None,
                                force: bool = False) -> tuple[np.ndarray, pd.DataFrame]:
    """Patristic distance matrices of the pre-registered posterior trees, pruned to `order`: (trees x n x n, one row per tree).

    Reading and pruning 100 trees of 5,911 tips takes several minutes, so the result is cached in
    data/interim and reused while the species order, tree count and seed are unchanged.
    """
    seed = int(cfg["stats"]["seed"]) if seed is None else seed
    cache = cfg.interim / "posterior_phylo.npz"
    key = f"{order_hash(order)}:{n_trees}:{seed}"
    if cache.exists() and not force:
        with np.load(cache, allow_pickle=False) as stored:
            if str(stored["key"]) == key:
                table = pd.DataFrame({"draw": np.arange(1, n_trees + 1), "tree_index": stored["tree_index"], "tree_file": stored["tree_file"]})
                return stored["matrices"], table

    archive = cfg.tree_path().parent / cfg["sources"]["tree"]["posterior_filename"]
    if not is_cached(archive):
        raise SystemExit(f"{archive.name} is not downloaded; run make posterior (about 1.2 GB) first")
    matched = read_csv(cfg.matched_path()).set_index("species").loc[list(order)]
    tips = matched["tree_tip"].map(posterior_label)
    relabel = dict(zip(tips, map(slugify_species, order)))
    matrices, rows = [], []
    with zipfile.ZipFile(archive) as zf:
        names = sorted(n for n in zf.namelist() if n.endswith(".tre"))
        chosen = select_trees(names, n_trees, seed)
        print(f"[posterior] {len(names)} trees in {archive.name}; {n_trees} drawn with seed {seed}")
        for k, index in enumerate(chosen, 1):
            tree = dendropy.Tree.get(data=zf.read(names[index]).decode(), schema="newick", preserve_underscores=True, rooting="default-rooted")
            phylo = phylo_distance_matrix(relabel_tips(prune_tree(tree, tips), relabel), order)
            validate_distance_matrix(phylo, names[index])
            matrices.append(phylo)
            rows.append({"draw": k, "tree_index": index, "tree_file": names[index].rsplit("/", 1)[-1]})
            if k % 10 == 0:
                print(f"[posterior] pruned {k}/{n_trees} trees")
    stack, table = np.stack(matrices), pd.DataFrame(rows)
    np.savez_compressed(cache, matrices=stack, key=np.array(key), tree_index=table["tree_index"].to_numpy(), tree_file=table["tree_file"].to_numpy(dtype=str))
    return stack, table
