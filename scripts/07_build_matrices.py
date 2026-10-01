from __future__ import annotations

from lineage_charisma.config import base_parser, load_config
from lineage_charisma.corpus import read_corpus
from lineage_charisma.distances import (
    assert_same_order,
    canonical_order,
    cosine_distance_matrix,
    load_matrix,
    phylo_distance_matrix,
    phylo_matrix_path,
    reindex_rows,
    save_matrix,
    text_matrix_path,
)
from lineage_charisma.embed import embedding_paths, load_embeddings, read_species_order, write_species_order
from lineage_charisma.io_utils import is_cached, read_meta, sha256, slugify_species
from lineage_charisma.phylo import load_tree


def main() -> None:
    parser = base_parser("Build aligned text-distance and phylogenetic-distance matrices.")
    parser.add_argument("--sensitivity", action="store_true",
                        help="drop species under matrices.sensitivity_min_tokens instead of the stubs, and write to a separate matrices_min<N>/ directory")
    args = parser.parse_args()
    cfg = load_config(args.config)
    mcfg = cfg.get("matrices", {})
    threshold = int(mcfg["sensitivity_min_tokens"]) if args.sensitivity else None
    out_dir, order_path = cfg.matrices_dir(threshold), cfg.species_order_path(threshold)
    combos = cfg.combos()

    corpus = read_corpus(cfg.corpus_path())
    excluding = args.sensitivity or mcfg.get("exclude_stubs", False)
    if excluding:
        # both rules use n_tokens of the unmasked text, so the species set is the same at every mask level
        drop = corpus["n_tokens"] < threshold if args.sensitivity else corpus["is_stub"]
        dropped = corpus.loc[drop, "species_id"].tolist()
        corpus = corpus[~drop]
        rule = f"n_tokens < {threshold}" if args.sensitivity else f"is_stub (n_tokens < {cfg['corpus']['min_tokens']})"
        print(f"[matrices] excluding {len(dropped)} species with {rule}: {dropped}")
    order = canonical_order(corpus["species_id"])

    tree_path = cfg.pruned_tree_path()
    emb_paths = {combo: embedding_paths(cfg.embeddings_dir(), *combo) for combo in combos}
    absent = [p.array.name for p in emb_paths.values() if not is_cached(p.all())]
    if absent:
        raise SystemExit(f"missing embeddings {absent}; run scripts/06_embed.py (make embed) first")
    sources = {"phylo": sha256(tree_path), **{"__".join(combo): sha256(p.array) for combo, p in emb_paths.items()}}

    outputs = {"phylo": phylo_matrix_path(out_dir), **{"__".join(combo): text_matrix_path(out_dir, *combo) for combo in combos}}
    fresh = (
        not args.force
        and is_cached([order_path, *outputs.values()])
        and read_species_order(order_path) == order
        and all(read_meta(path).get("source_sha256") == sources[key] for key, path in outputs.items())
    )
    if fresh:
        print(f"[cache] {len(outputs)} matrices in {out_dir.relative_to(cfg.root)} are up to date, skipping (use --force to rebuild)")
        return

    write_species_order(order_path, order)
    print(f"[matrices] canonical order: {len(order)} species -> {order_path.relative_to(cfg.root)}")

    tree = load_tree(tree_path)
    tips = {leaf.taxon.label for leaf in tree.leaf_node_iter()}
    if not excluding and tips != {slugify_species(s) for s in order}:
        raise SystemExit("[matrices] the pruned tree's tips are not exactly the corpus species; rerun make prune and make corpus")
    phylo = phylo_distance_matrix(tree, order)
    save_matrix(outputs["phylo"], phylo, order, source=tree_path.name, source_sha256=sources["phylo"],
                units="million years of branch length between tips (patristic; twice the divergence time on this ultrametric tree)")
    print(f"[matrices] {outputs['phylo'].name}: {phylo.shape}, max={phylo.max():.2f} My")

    for combo in combos:
        model_name, rule, mask = combo
        key = "__".join(combo)
        embeddings, emb_species = load_embeddings(emb_paths[combo])
        embeddings = reindex_rows(embeddings, emb_species, order, emb_paths[combo].array.name)
        text = cosine_distance_matrix(embeddings)
        save_matrix(outputs[key], text, order, model=model_name, rule=rule, mask_level=mask, metric="1 - cosine similarity",
                    source=emb_paths[combo].array.name, source_sha256=sources[key])
        print(f"[matrices] {outputs[key].name}: {text.shape}, mean={text.mean():.3f}, max={text.max():.3f}")

    for key, path in outputs.items():
        matrix, loaded_order = load_matrix(path, order_path)
        assert_same_order(loaded_order, order, path.name)
        assert matrix.shape == phylo.shape, f"{path.name}: shape {matrix.shape} != phylo {phylo.shape}"
    print(f"[matrices] verified: {len(outputs)} matrices share one species order ({len(order)} species)")


if __name__ == "__main__":
    main()
