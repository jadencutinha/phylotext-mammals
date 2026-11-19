from __future__ import annotations

import hashlib
import json
from collections import Counter

import numpy as np
import pandas as pd

from lineage_charisma import condensed as cd
from lineage_charisma.config import base_parser, load_config
from lineage_charisma.corpus import load_token_counter
from lineage_charisma.distances import load_matrix, text_matrix_path
from lineage_charisma.embed import (
    SentenceTransformerEncoder,
    cache_is_current,
    embed_texts,
    embedding_paths,
    load_embeddings,
    read_species_order,
    save_embeddings,
)
from lineage_charisma.extensions import (
    LEAD,
    MIN_GROUP_SPECIES,
    MIN_TOKENS,
    REPORTED,
    STATUS,
    amendment_is_committed,
    assign_groups,
    group_texts,
    h1_test,
    scale_masker,
    tables_dir,
    text_distance,
    work_dir,
)
from lineage_charisma.io_utils import slugify_species, write_csv
from lineage_charisma.masking import term_leaks
from lineage_charisma.wiki import page_from_response, split_sections

# docs/preregistration_amendment_extensions.md, section C (secondary, descriptive, designed after the 2026-11-11 lock)
SCALES = {"pilot": None, "full": "mammalia"}
NOTES = {"other": "mixes taxonomy, conservation and human relations",
         LEAD: "reference row: the lead alone, which is part of the corpus text"}


def extracts(cfg, species) -> dict[str, str]:
    raw = cfg.wiki_dir() / "raw"
    return {s: page_from_response(json.loads((raw / f"{slugify_species(s)}.json").read_text()))["extract"] for s in species}


def heading_map(cfg) -> pd.DataFrame:
    """Every heading over the frozen full-scale species, with its count and the group and rule the amendment assigns it.
    Inherited assignments depend on the enclosing section, so a heading can appear under more than one group."""
    counts: Counter = Counter()
    for text in extracts(cfg, read_species_order(cfg.species_order_path())).values():
        sections = split_sections(text)
        for sec, (group, rule) in zip(sections[1:], assign_groups(sections)):
            counts[(sec.heading.strip().lower(), sec.level, group, rule)] += 1
    rows = [{"heading": h, "level": lvl, "group": g, "rule": r, "count": n} for (h, lvl, g, r), n in counts.items()]
    return pd.DataFrame(rows).sort_values(["count", "heading"], ascending=[False, True]).reset_index(drop=True)


def main() -> None:
    parser = base_parser("Amendment 5, section C: H1 per Wikipedia section group (descriptive).")
    parser.add_argument("--scale", action="append", choices=list(SCALES), help="default: pilot then full")
    parser.add_argument("--map-only", action="store_true", help="write the heading-to-group table only (no text distance or statistic)")
    args = parser.parse_args()
    full = load_config(args.config, "mammalia")
    table = heading_map(full)
    write_csv(table, tables_dir(load_config(args.config, None)) / "section_heading_map.csv", status=STATUS,
              plan="docs/preregistration_amendment_extensions.md, section C", species="the 4,611 frozen full-scale species")
    print(f"[sections] heading map: {len(table):,} heading/level/group rows; sections by group: "
          f"{table.groupby('group')['count'].sum().to_dict()}")
    if args.map_only:
        return
    ok, why = amendment_is_committed(full.root)
    if not ok:
        raise SystemExit(f"[sections] {why}")

    for scale in args.scale or list(SCALES):
        cfg = load_config(args.config, SCALES[scale])
        s, ecfg = cfg["stats"], cfg["embedding"]
        perms, seed, method = int(s["permutations"]), int(s["seed"]), s["method"]
        model, rule, level = cfg.primary_combo
        order_path = cfg.species_order_path()
        order = read_species_order(order_path)
        phylo = load_matrix(cfg.matrices_dir() / "phylo.npy", order_path)[0]
        primary = load_matrix(text_matrix_path(cfg.matrices_dir(), *cfg.primary_combo), order_path)[0]
        groups = {sp: group_texts(text) for sp, text in extracts(cfg, order).items()}
        count = load_token_counter(cfg["corpus"]["tokenizer"])
        masker = scale_masker(cfg)
        encoder = None
        rows = []
        for group in REPORTED:
            have = [sp for sp in order if group in groups[sp]]
            tokens = dict(zip(have, count([groups[sp][group] for sp in have])))
            species = [sp for sp in have if tokens[sp] >= MIN_TOKENS]
            row = {"status": STATUS, "scale": scale, "group": group, "n_with_section": len(have), "n_species": len(species),
                   "median_tokens": float(np.median([tokens[sp] for sp in species])) if species else float("nan")}
            if len(species) < MIN_GROUP_SPECIES:
                rows.append({**row, "note": f"fewer than {MIN_GROUP_SPECIES} species; not tested"})
                print(f"[sections] {scale} {group}: {len(species)} species, not tested")
                continue
            masked = [masker.mask(groups[sp][group], level) for sp in species]
            leaks = [sp for sp, t in zip(species, masked) if term_leaks(t, sp.split()[0])]
            if leaks:
                raise SystemExit(f"[sections] {scale} {group}: masked text still names its own genus for {leaks[:5]}")
            paths = embedding_paths(work_dir(cfg), model, rule, f"{level}_section_{group}")
            params = {"model": model, "rule": rule, "mask_level": level, "chunk_overlap": ecfg["chunk_overlap"],
                      "texts_sha256": hashlib.sha256("\x00".join([*species, *masked]).encode()).hexdigest()}
            if not cache_is_current(paths, params)[0] or args.force:
                encoder = encoder or SentenceTransformerEncoder(model, device=ecfg.get("device", "auto"), max_tokens=ecfg.get("max_tokens"))
                result = embed_texts(masked, encoder, rule, overlap=ecfg["chunk_overlap"], batch_size=ecfg["batch_size"])
                save_embeddings(paths, result.embeddings, species, **params, status=STATUS)
            emb, emb_species = load_embeddings(paths)
            assert emb_species == species
            res = h1_test(text_distance(emb), cd.subset(phylo, order, species), perms, seed, method)
            ref = h1_test(cd.subset(primary, order, species), cd.subset(phylo, order, species), perms, seed, method)
            rows.append({**row, "r": res["r"], "p": res["p"], "n_pairs": res["n_pairs"], "reference_r": ref["r"], "reference_p": ref["p"],
                         "permutations": perms, "note": NOTES.get(group, "")})
            print(f"[sections] {scale} {group:22s} {len(species):5d} species  r = {res['r']:+.4f} (p = {res['p']:.4f})  "
                  f"corpus text on the same species r = {ref['r']:+.4f}")
        write_csv(pd.DataFrame(rows), tables_dir(cfg) / "section_groups_h1.csv", status=STATUS,
                  plan="docs/preregistration_amendment_extensions.md, section C", reading="descriptive only; p-values for completeness, no multiple-testing claim",
                  permutations=perms, seed=seed, alternative="greater", permuted="phylogeny", min_tokens=MIN_TOKENS, min_species=MIN_GROUP_SPECIES)


if __name__ == "__main__":
    main()
