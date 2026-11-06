from __future__ import annotations

from collections import Counter

import numpy as np
import pandas as pd

from lineage_charisma import condensed as cd
from lineage_charisma.config import base_parser, load_config
from lineage_charisma.corpus import load_texts, read_corpus
from lineage_charisma.distances import cosine_distance_matrix, load_matrix, phylo_distance_matrix, phylo_matrix_path, text_matrix_path
from lineage_charisma.embed import SentenceTransformerEncoder, embed_texts, embedding_paths, load_embeddings, read_species_order
from lineage_charisma.io_utils import read_csv, write_csv, write_json
from lineage_charisma.masking import PILOT_NAME_RANKS, Masker, build_vocabulary, text_column
from lineage_charisma.phylo import load_tree

THRESHOLD = 0.01   # docs/preregistration_amendment_fullscale.md, section 6b


def vocabulary(cfg, species_path, synonyms_path):
    m = cfg.get("masking") or {}
    return build_vocabulary(read_csv(species_path), read_csv(synonyms_path), extra_terms=m.get("extra_taxonomy_terms", []),
                            keep_terms=m.get("keep_terms", []), keep_when_followed_by=m.get("keep_when_followed_by", {}),
                            head_noun_rule=m.get("head_noun_rule", "pilot"), epithets_anywhere=m.get("epithets_anywhere", True),
                            name_ranks=m.get("name_ranks", PILOT_NAME_RANKS))


def term_counts(masker: Masker, text: str, level: str) -> Counter:
    counts: Counter = Counter()
    masker.mask(text, level, counts)
    return Counter({surface.lower(): n for (_, surface), n in counts.items()})


def main() -> None:
    parser = base_parser("Full scale, section 6: exact replication of the pilot H1 (6a) and replication with the full-scale vocabulary (6b).")
    args = parser.parse_args()
    full = load_config(args.config, args.clade or "mammalia")
    pilot = load_config(args.config, None)
    scfg = full["stats"]
    perms, seed, method = int(scfg["permutations"]), int(scfg["seed"]), scfg["method"]
    combo = full.primary_combo
    model, rule, level = combo
    ecfg = full["embedding"]

    pilot_order = read_species_order(pilot.species_order_path())
    frozen = read_species_order(full.species_order_path())
    reference = read_csv(pilot.tables_dir() / "h1_mantel.csv")
    reference = reference[reference["is_primary"]].iloc[0]
    rows, notes = [], {}

    # 6a: the full-scale data path on the pilot species, masked with the pilot's own vocabulary
    matched = read_csv(full.matched_path()).set_index("species")
    pilot_matched = read_csv(pilot.matched_path()).set_index("species")
    missing = [s for s in pilot_order if s not in matched.index]
    changed = [s for s in pilot_order if s in matched.index and (matched.at[s, "tree_tip"], matched.at[s, "wiki_title"]) != (pilot_matched.at[s, "tree_tip"], pilot_matched.at[s, "wiki_title"])]
    texts = load_texts(pilot_order, full.wiki_dir() / "clean")
    pilot_masker = Masker(vocabulary(pilot, pilot.mdd_species_path(), pilot.mdd_synonyms_path()), (pilot.get("masking") or {}).get("token", "[TAXON]"))
    masked = [pilot_masker.mask(texts[s], level) for s in pilot_order]
    pilot_corpus = read_corpus(pilot.corpus_path()).set_index("species_id")
    same_text = sum(t == pilot_corpus.at[s, text_column(level)] for s, t in zip(pilot_order, masked))

    cache = full.interim / "replication" / f"exact_{'__'.join(combo).replace('/', '_')}.npy"
    if cache.exists() and not args.force:
        emb = np.load(cache)
    else:
        encoder = SentenceTransformerEncoder(model, device=ecfg.get("device", "auto"), max_tokens=ecfg.get("max_tokens"))
        emb = embed_texts(masked, encoder, rule, overlap=ecfg["chunk_overlap"], batch_size=ecfg["batch_size"]).embeddings
        cache.parent.mkdir(parents=True, exist_ok=True)
        np.save(cache, emb)
    pilot_emb, pilot_emb_species = load_embeddings(embedding_paths(pilot.embeddings_dir(), *combo))
    pilot_emb = pilot_emb[[pilot_emb_species.index(s) for s in pilot_order]]
    emb_diff = float(np.abs(emb - pilot_emb).max())

    tree = load_tree(full.pruned_tree_path())
    phylo_6a = phylo_distance_matrix(tree, pilot_order)
    pilot_phylo, _ = load_matrix(phylo_matrix_path(pilot.matrices_dir()), pilot.species_order_path())
    text_6a = cosine_distance_matrix(emb)
    res = cd.mantel(cd.prepare(text_6a, method), cd.prepare(phylo_6a, method), permutations=perms, seed=seed)
    exact = round(res.r, 3) == round(float(reference["r"]), 3) and res.p == float(reference["p"])
    rows.append({"check": "6a exact (pilot vocabulary)", "n_species": res.n, "n_pairs": res.n_pairs, "pilot_r": reference["r"], "pilot_p": reference["p"],
                 "r": res.r, "p": res.p, "difference": res.r - reference["r"], "passed": bool(exact)})
    notes["6a"] = {"pilot_species_missing_from_full_scale_matches": missing, "pilot_species_with_a_different_tip_or_article": changed,
                   "masked_texts_identical_to_pilot_corpus": f"{same_text}/{len(pilot_order)}", "max_abs_embedding_difference": emb_diff,
                   "phylo_identical_to_pilot": bool(np.array_equal(phylo_6a, pilot_phylo)),
                   "r_exact_difference": res.r - float(reference["r"])}
    print(f"[6a] r = {res.r:.6f} (pilot {reference['r']:.6f}), p = {res.p:.4f} (pilot {reference['p']:.4f}); "
          f"texts identical {same_text}/{len(pilot_order)}; max embedding difference {emb_diff:.2e}; phylo identical {notes['6a']['phylo_identical_to_pilot']}")
    if not exact:
        write_csv(pd.DataFrame(rows), full.tables_dir() / "replication.csv")
        raise SystemExit("[6a] FAILED: the full-scale path does not reproduce the pilot H1; this is a bug to fix before going on")

    # 6b: the full-scale primary matrices cut to the pilot species
    shared = [s for s in pilot_order if s in set(frozen)]
    text_full, _ = load_matrix(text_matrix_path(full.matrices_dir(), *combo), full.species_order_path())
    phylo_full, _ = load_matrix(phylo_matrix_path(full.matrices_dir()), full.species_order_path())
    text_b, phylo_b = cd.subset(text_full, frozen, shared), cd.subset(phylo_full, frozen, shared)
    phylo_ref = cd.subset(pilot_phylo, pilot_order, shared)
    phylo_diff = float(np.abs(phylo_b - phylo_ref).max())
    if phylo_diff != 0:
        raise SystemExit(f"[6b] phylogenetic distances differ from the pilot by up to {phylo_diff}; this is a bug")
    res_b = cd.mantel(cd.prepare(text_b, method), cd.prepare(phylo_b, method), permutations=perms, seed=seed)
    pilot_r_shared = reference["r"]
    if len(shared) != len(pilot_order):
        pilot_text, _ = load_matrix(text_matrix_path(pilot.matrices_dir(), *combo), pilot.species_order_path())
        pilot_r_shared = cd.correlation(cd.prepare(cd.subset(pilot_text, pilot_order, shared), method), cd.prepare(phylo_ref, method))
    gap = res_b.r - pilot_r_shared
    rows.append({"check": "6b full-scale vocabulary", "n_species": res_b.n, "n_pairs": res_b.n_pairs, "pilot_r": pilot_r_shared, "pilot_p": reference["p"],
                 "r": res_b.r, "p": res_b.p, "difference": gap, "passed": bool(abs(gap) <= THRESHOLD)})

    full_masker = Masker(vocabulary(full, full.mdd_species_path(), full.mdd_synonyms_path()), (full.get("masking") or {}).get("token", "[TAXON]"))
    full_corpus = read_corpus(full.corpus_path()).set_index("species_id")
    term_rows, n_changed = [], 0
    for s in shared:
        if full_corpus.at[s, text_column(level)] == pilot_corpus.at[s, text_column(level)]:
            continue
        n_changed += 1
        before, after = term_counts(pilot_masker, texts[s], level), term_counts(full_masker, texts[s], level)
        for term in set(before) | set(after):
            if before[term] != after[term]:
                term_rows.append({"species": s, "term": term, "pilot_masks": before[term], "full_scale_masks": after[term]})
    terms = pd.DataFrame(term_rows)
    summary = (terms.assign(direction=np.where(terms["full_scale_masks"] > terms["pilot_masks"], "newly masked", "unmasked"))
               .groupby(["term", "direction"]).agg(species=("species", "nunique"), occurrences_pilot=("pilot_masks", "sum"),
                                                   occurrences_full=("full_scale_masks", "sum")).reset_index()
               .sort_values("species", ascending=False))
    notes["6b"] = {"pilot_species_missing_from_frozen_set": [s for s in pilot_order if s not in set(frozen)], "phylo_max_abs_difference": phylo_diff,
                   "pilot_texts_changed_by_full_scale_masks": n_changed, "threshold": THRESHOLD}
    print(f"[6b] r = {res_b.r:.6f} on {res_b.n} species (pilot {pilot_r_shared:.6f}), difference {gap:+.4f}, p = {res_b.p:.4f}; "
          f"{n_changed} pilot texts masked differently; {'within' if abs(gap) <= THRESHOLD else 'OUTSIDE'} the {THRESHOLD} threshold")

    tables = full.tables_dir()
    write_csv(pd.DataFrame(rows), tables / "replication.csv", permutations=perms, seed=seed, method=method, spec=list(combo))
    write_csv(summary, tables / "replication_mask_terms.csv", level=level, n_texts_changed=n_changed)
    write_json(tables / "replication_notes.json", notes)
    print(f"[replication] wrote {tables.relative_to(full.root)}/replication.csv, replication_mask_terms.csv, replication_notes.json")


if __name__ == "__main__":
    main()
