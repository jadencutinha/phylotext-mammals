from __future__ import annotations

from collections import Counter

import pandas as pd

from lineage_charisma.config import base_parser, load_config
from lineage_charisma.corpus import build_corpus, length_summary, load_texts, load_token_counter, stub_table, write_corpus
from lineage_charisma.io_utils import read_csv, skip_if_cached, write_csv, write_meta
from lineage_charisma.masking import MASK_LEVELS, Masker, build_vocabulary, leak_counts, term_leaks, text_column


def main() -> None:
    args = base_parser("Assemble the cleaned descriptions into one corpus table, flag stubs, and build the masked variants.").parse_args()
    cfg = load_config(args.config)
    out_corpus, out_stubs, out_terms = cfg.corpus_path(), cfg.stub_species_path(), cfg.mask_terms_path()
    if skip_if_cached([out_corpus, out_stubs, out_terms], args.force):
        return

    ccfg, mcfg = cfg["corpus"], cfg.get("masking") or {}
    levels = cfg.mask_levels
    unknown = [k for k in levels if k not in MASK_LEVELS]
    if unknown:
        raise SystemExit(f"unknown mask levels in config: {unknown}; expected a subset of {list(MASK_LEVELS)}")
    matched = read_csv(cfg.matched_path())
    texts = load_texts(matched["species"], cfg.wiki_dir() / "clean")
    empty = [s for s, t in texts.items() if not t]
    if empty:
        raise SystemExit(f"{len(empty)} species have empty cleaned text: {empty[:10]}")

    vocabulary = build_vocabulary(
        read_csv(cfg.mdd_species_path()), read_csv(cfg.mdd_synonyms_path()),
        extra_terms=mcfg.get("extra_taxonomy_terms", []), keep_terms=mcfg.get("keep_terms", []),
        keep_when_followed_by=mcfg.get("keep_when_followed_by", {}),
    )
    masker = Masker(vocabulary, mcfg.get("token", "[TAXON]"))
    corpus = build_corpus(matched, texts, load_token_counter(ccfg["tokenizer"]), ccfg["min_tokens"], masker, levels)

    genus = [s.split()[0] for s in corpus["species_id"]]
    epithet = [s.split()[1] for s in corpus["species_id"]]
    term_rows, leaks = [], {}
    for level in levels:
        col = text_column(level)
        leaks[level] = {
            "own_genus": leak_counts(corpus[col], genus),
            "own_epithet": leak_counts(corpus[col], epithet),
            "own_family": leak_counts(corpus[col], corpus["family"]),
        }
        if level == "none":
            continue
        counts: Counter = Counter()
        for text in corpus["text"]:
            masker.mask(text, level, counts)
        term_rows += [{"mask_level": level, "kind": kind, "surface": surface, "count": n} for (kind, surface), n in counts.most_common()]
        if leaks[level]["own_genus"]:
            bad = [s for s, t, g in zip(corpus["species_id"], corpus[col], genus) if term_leaks(t, g)]
            raise SystemExit(f"[corpus] mask level {level}: {len(bad)} texts still contain their own genus name: {bad[:10]}")

    write_corpus(corpus, out_corpus)
    summary = length_summary(corpus)
    write_meta(out_corpus, rows=len(corpus), columns=list(corpus.columns), tokenizer=ccfg["tokenizer"], min_tokens=ccfg["min_tokens"],
               mask_levels=levels, mask_token=masker.token, own_name_leaks=leaks, **summary)
    wiki = read_csv(cfg.interim / "wiki_texts.csv").rename(columns={"species": "species_id"})
    stubs = stub_table(corpus, wiki[["species_id", "wiki_title", "has_description"]])
    write_csv(stubs, out_stubs, min_tokens=ccfg["min_tokens"], tokenizer=ccfg["tokenizer"])
    write_csv(pd.DataFrame(term_rows, columns=["mask_level", "kind", "surface", "count"]), out_terms, mask_token=masker.token)

    print(f"[corpus] {len(corpus)} species -> {out_corpus.relative_to(cfg.root)}")
    print(f"[corpus] tokens ({ccfg['tokenizer']}): min={summary['min']} median={summary['median']} max={summary['max']}; "
          f"over 512: {summary['over_512']}; over 384: {summary['over_384']}")
    print(f"[corpus] {len(stubs)} species below min_tokens={ccfg['min_tokens']} (flagged, not dropped) -> {out_stubs.relative_to(cfg.root)}")
    for rec in stubs.itertuples(index=False):
        print(f"         {rec.species_id} ({rec.wiki_title}): {rec.n_tokens} tokens")
    for level in levels:
        note = "" if level == "none" else f"; median {int(corpus[f'n_masked_{level}'].median())} masks per text"
        print(f"[corpus] mask level {level}: texts still naming their own genus {leaks[level]['own_genus']}, "
              f"epithet {leaks[level]['own_epithet']}, family {leaks[level]['own_family']}{note}")
    print(f"[corpus] every replaced term -> {out_terms.relative_to(cfg.root)}")


if __name__ == "__main__":
    main()
