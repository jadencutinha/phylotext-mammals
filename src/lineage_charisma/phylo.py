from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable, Mapping

import dendropy
import numpy as np
import pandas as pd

TIP_RE = re.compile(r"^(?P<genus>[A-Z][A-Za-z-]+)_(?P<epithet>[a-z][a-z-]*)(?:_(?P<family>[A-Za-z]+))?(?:_(?P<order>[A-Za-z]+))?$")


def parse_tip_label(label: str) -> dict[str, str | None] | None:
    text = label.strip().strip("'\"").replace(" ", "_")
    m = TIP_RE.match(text)
    if not m:
        return None
    d = m.groupdict()
    return {
        "tip_label": text,
        "binomial": f"{d['genus']} {d['epithet']}",
        "genus": d["genus"],
        "epithet": d["epithet"],
        "family": d["family"].upper() if d["family"] else None,
        "order": d["order"].upper() if d["order"] else None,
    }


def load_tree(path: Path, schema: str | None = None) -> dendropy.Tree:
    path = Path(path)
    if schema is None:
        with open(path) as fh:
            schema = "nexus" if fh.read(64).lstrip().upper().startswith("#NEXUS") else "newick"
    return dendropy.Tree.get(path=str(path), schema=schema, preserve_underscores=True, rooting="default-rooted")


def tip_table(tree: dendropy.Tree) -> pd.DataFrame:
    rows = []
    for leaf in tree.leaf_node_iter():
        parsed = parse_tip_label(leaf.taxon.label)
        if parsed is not None:
            rows.append(parsed)
    return pd.DataFrame(rows, columns=["tip_label", "binomial", "genus", "epithet", "family", "order"])


def filter_tips(tips: pd.DataFrame, orders: Iterable[str] | None = None, families: Iterable[str] | None = None) -> pd.DataFrame:
    out = tips
    if orders:
        out = out[out["order"].isin({o.upper() for o in orders})]
    if families:
        out = out[out["family"].isin({f.upper() for f in families})]
    return out.reset_index(drop=True)


def prune_tree(tree: dendropy.Tree, keep_labels: Iterable[str]) -> dendropy.Tree:
    keep = set(keep_labels)
    present = {t.label for t in tree.taxon_namespace}
    missing = keep - present
    if missing:
        raise KeyError(f"{len(missing)} labels not in tree, e.g. {sorted(missing)[:5]}")
    pruned = tree.extract_tree_with_taxa_labels(labels=keep)
    newick = pruned.as_string(schema="newick", suppress_rooting=False, suppress_annotations=True, unquoted_underscores=True)
    return dendropy.Tree.get(data=newick, schema="newick", preserve_underscores=True, rooting="default-rooted")


def relabel_tips(tree: dendropy.Tree, mapping: Mapping[str, str]) -> dendropy.Tree:
    out = tree.clone(depth=2)
    for leaf in out.leaf_node_iter():
        if leaf.taxon.label in mapping:
            leaf.taxon.label = mapping[leaf.taxon.label]
    return out


def cophenetic_matrix(tree: dendropy.Tree) -> pd.DataFrame:
    pdm = tree.phylogenetic_distance_matrix()
    taxa = sorted((leaf.taxon for leaf in tree.leaf_node_iter()), key=lambda t: t.label)
    labels = [t.label for t in taxa]
    n = len(taxa)
    mat = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = pdm.patristic_distance(taxa[i], taxa[j])
            mat[i, j] = mat[j, i] = d
    return pd.DataFrame(mat, index=labels, columns=labels)


def group_separation(dist: pd.DataFrame, members: Iterable[str]) -> dict[str, float | int | bool]:
    inside = [m for m in members if m in dist.index]
    outside = [x for x in dist.index if x not in set(inside)]
    if len(inside) < 2 or not outside:
        return {"n": len(inside), "max_within": float("nan"), "min_to_outside": float("nan"), "mean_within": float("nan"), "mean_to_outside": float("nan"), "clusters": False}
    within = dist.loc[inside, inside].to_numpy()
    iu = np.triu_indices(len(inside), k=1)
    between = dist.loc[inside, outside].to_numpy()
    max_within = float(within[iu].max())
    min_out = float(between.min())
    return {
        "n": len(inside),
        "max_within": round(max_within, 3),
        "min_to_outside": round(min_out, 3),
        "mean_within": round(float(within[iu].mean()), 3),
        "mean_to_outside": round(float(between.mean()), 3),
        "clusters": max_within < min_out,
    }


def nearest_neighbors(dist: pd.DataFrame, label: str, k: int = 5) -> list[tuple[str, float]]:
    row = dist.loc[label].drop(label).sort_values()
    return [(idx, round(float(v), 3)) for idx, v in row.head(k).items()]


def is_ultrametric(tree: dendropy.Tree, rel_tol: float = 1e-3) -> bool:
    depths = [leaf.distance_from_root() for leaf in tree.leaf_node_iter()]
    return (max(depths) - min(depths)) <= rel_tol * max(depths)
