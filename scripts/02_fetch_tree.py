from __future__ import annotations

from lineage_charisma.config import base_parser, load_config
from lineage_charisma.io_utils import download, is_cached, read_csv, skip_if_cached, write_csv
from lineage_charisma.phylo import filter_tips, load_tree, tip_table


def main() -> None:
    parser = base_parser("Download the Upham et al. 2019 MCC tree and extract target-clade tip labels.")
    parser.add_argument("--posterior", action="store_true", help="also download the full 10k-tree posterior zip (~1.2 GB)")
    args = parser.parse_args()
    cfg = load_config(args.config)
    src = cfg["sources"]["tree"]
    tree_path = cfg.tree_path()
    tree_path.parent.mkdir(parents=True, exist_ok=True)

    if args.force or not is_cached(tree_path):
        print(f"[download] {src['mcc_url']}")
        download(src["mcc_url"], tree_path, user_agent=cfg["wikipedia"]["user_agent"])
    else:
        print(f"[cache] {tree_path} exists")

    if args.posterior:
        post = tree_path.parent / src["posterior_filename"]
        if args.force or not is_cached(post):
            print(f"[download] {src['posterior_url']} (large)")
            download(src["posterior_url"], post, user_agent=cfg["wikipedia"]["user_agent"], timeout=600)

    all_tips_path = cfg.interim / "tree_tips_all.csv"
    clade_tips_path = cfg.tips_path()
    if skip_if_cached([all_tips_path, clade_tips_path], args.force):
        return

    print(f"[tree] loading {tree_path.name}")
    tree = load_tree(tree_path, schema=src.get("schema"))
    tips = tip_table(tree)
    n_leaves = len(list(tree.leaf_node_iter()))
    print(f"[tree] {n_leaves} leaves, {len(tips)} parsed as Genus_species_FAMILY_ORDER")

    species = read_csv(cfg.mdd_species_path())
    orders = sorted(species["order"].dropna().str.upper().unique())
    clade_tips = filter_tips(tips, orders=orders)
    write_csv(tips, all_tips_path, tree_file=tree_path.name)
    write_csv(clade_tips, clade_tips_path, tree_file=tree_path.name, orders=orders)
    print(f"[tree] {len(clade_tips)} tips in {orders} -> {clade_tips_path.name}")


if __name__ == "__main__":
    main()
