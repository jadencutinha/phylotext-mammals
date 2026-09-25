from __future__ import annotations

from lineage_charisma.config import base_parser, load_config
from lineage_charisma.io_utils import atomic_write_text, read_csv, skip_if_cached, slugify_species, write_csv, write_meta
from lineage_charisma.phylo import cophenetic_matrix, group_separation, is_ultrametric, load_tree, nearest_neighbors, prune_tree, relabel_tips


def sanity_report(dist, matched, groups: dict, examples: list[str]) -> tuple[list[str], bool]:
    lines, ok = [], True
    for name, spec in groups.items():
        col = spec["rank"]
        members = matched.loc[matched[col].str.lower() == str(spec["name"]).lower(), "slug"].tolist()
        stats = group_separation(dist, members)
        passed = bool(stats["clusters"]) or stats["n"] < 2
        ok &= passed
        lines.append(
            f"{'PASS' if passed else 'FAIL'} {name} ({col}={spec['name']}, n={stats['n']}): "
            f"max within={stats['max_within']} < min to outside={stats['min_to_outside']} "
            f"(mean within={stats['mean_within']}, mean to outside={stats['mean_to_outside']})"
        )
    fam = matched.set_index("slug")["family"]
    same = sum(fam[nearest_neighbors(dist, s, k=1)[0][0]] == fam[s] for s in dist.index)
    rate = same / len(dist)
    ok &= rate >= 0.95
    lines.append(f"{'PASS' if rate >= 0.95 else 'FAIL'} nearest phylogenetic neighbour shares family for {same}/{len(dist)} species ({rate:.1%})")
    for example in [slugify_species(s) for s in examples if slugify_species(s) in dist.index]:
        nn = ", ".join(f"{n} ({d:.1f})" for n, d in nearest_neighbors(dist, example, k=4))
        lines.append(f"     nearest to {example}: {nn}")
    return lines, ok


def main() -> None:
    args = base_parser("Prune the tree to matched species and compute cophenetic distances.").parse_args()
    cfg = load_config(args.config)
    out_tree = cfg.interim / "tree_pruned.nwk"
    out_tree_mdd = cfg.interim / "tree_pruned_mdd_names.nwk"
    out_dist = cfg.interim / "cophenetic_distance.csv"
    out_report = cfg.interim / "sanity_check.txt"
    if skip_if_cached([out_tree, out_tree_mdd, out_dist, out_report], args.force):
        return

    matched = read_csv(cfg.matched_path())
    matched["slug"] = matched["species"].map(slugify_species)
    print(f"[tree] loading {cfg.tree_path().name}")
    tree = load_tree(cfg.tree_path(), schema=cfg["sources"]["tree"].get("schema"))
    pruned = prune_tree(tree, matched["tree_tip"])
    renamed = relabel_tips(pruned, dict(zip(matched["tree_tip"], matched["slug"])))
    n = len(list(pruned.leaf_node_iter()))
    print(f"[tree] pruned to {n} tips (ultrametric={is_ultrametric(pruned)})")

    atomic_write_text(out_tree, pruned.as_string(schema="newick", unquoted_underscores=True, suppress_rooting=False))
    atomic_write_text(out_tree_mdd, renamed.as_string(schema="newick", unquoted_underscores=True, suppress_rooting=False))
    write_meta(out_tree, tips=n, source_tree=cfg.tree_path().name)
    write_meta(out_tree_mdd, tips=n, source_tree=cfg.tree_path().name, labels="MDD species (Genus_species)")

    dist = cophenetic_matrix(renamed)
    atomic_write_text(out_dist, dist.round(6).to_csv(index_label="species"))
    write_meta(out_dist, shape=list(dist.shape), units="million years (sum of branch lengths)")
    print(f"[dist] {dist.shape[0]}x{dist.shape[1]} cophenetic matrix, max={dist.values.max():.1f}")

    sc = cfg.get("sanity_check", {})
    lines, ok = sanity_report(dist, matched, sc.get("groups", {}), sc.get("example_species", []))
    lines.insert(0, f"tips={n} ultrametric={is_ultrametric(pruned)} all_matched_present={set(dist.index) == set(matched['slug'])}")
    atomic_write_text(out_report, "\n".join(lines) + "\n")
    print("\n".join(lines))
    if not ok:
        raise SystemExit("[sanity] FAILED - see sanity_check.txt")


if __name__ == "__main__":
    main()
