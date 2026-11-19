from __future__ import annotations

import pandas as pd

from lineage_charisma.condensed import subset
from lineage_charisma.config import base_parser, load_config
from lineage_charisma.corpus import read_corpus
from lineage_charisma.distances import load_matrix, save_matrix, text_matrix_path
from lineage_charisma.embed import read_species_order
from lineage_charisma.extensions import (
    STATUS,
    TFIDF_SETTINGS,
    amendment_is_committed,
    load_registered,
    matrix_agreement,
    registered_tests,
    tables_dir,
    tfidf_distance,
    work_dir,
)
from lineage_charisma.io_utils import write_csv

# docs/preregistration_amendment_extensions.md, section B (secondary, designed after the 2026-11-11 lock)
SCALES = {"pilot": None, "full": "mammalia"}


def main() -> None:
    parser = base_parser("Amendment 5, section B: TF-IDF baseline on the strict-masked texts, registered H1 and H2 tests.")
    parser.add_argument("--scale", action="append", choices=list(SCALES), help="default: pilot then full")
    args = parser.parse_args()
    ok, why = amendment_is_committed(load_config(args.config, None).root)
    if not ok:
        raise SystemExit(f"[tfidf] {why}")
    for scale in args.scale or list(SCALES):
        cfg = load_config(args.config, SCALES[scale])
        s = cfg["stats"]
        perms, seed, alpha, method = int(s["permutations"]), int(s["seed"]), float(s["alpha"]), s["method"]
        order = read_species_order(cfg.species_order_path())
        corpus = read_corpus(cfg.corpus_path()).set_index("species_id")
        result = tfidf_distance(corpus.loc[order, "text_strict"].tolist())
        dropped = [order[i] for i in result.dropped]
        kept = [order[i] for i in result.kept]
        reg = load_registered(cfg)
        if dropped:
            print(f"[tfidf] {scale}: {len(dropped)} species have an all-zero vector and are dropped: {dropped}")
            reg = reg.restrict(kept)
        save_matrix(work_dir(cfg) / "text__tfidf__strict.npy", result.distance, kept, metric="1 - cosine similarity of TF-IDF vectors",
                    status=STATUS, **{k: str(v) for k, v in TFIDF_SETTINGS.items()})
        print(f"[tfidf] {scale}: {len(kept)} species, vocabulary {result.vocabulary_size:,} terms, median {result.median_terms:.0f} terms per text")

        rows = [{"status": STATUS, "scale": scale, "representation": "tfidf", **row} for row in registered_tests(result.distance, reg, perms, seed, alpha, method)]
        for row in rows:
            print(f"[tfidf] {scale} {row['test']:30s} r = {row['r']:+.4f} p = {row['p']:.4f}")
        primary = load_matrix(text_matrix_path(cfg.matrices_dir(), *cfg.primary_combo), cfg.species_order_path())[0]
        agreement = matrix_agreement({"bge-large (primary)": subset(primary, order, kept), "tfidf": result.distance}, method)
        meta = dict(status=STATUS, plan="docs/preregistration_amendment_extensions.md, section B", permutations=perms, seed=seed, alternative="greater",
                    permuted="phylogeny", vocabulary_size=result.vocabulary_size, median_terms=result.median_terms, dropped_species=dropped,
                    tfidf={k: str(v) for k, v in TFIDF_SETTINGS.items()})
        write_csv(pd.DataFrame(rows), tables_dir(cfg) / "tfidf_tests.csv", **meta)
        write_csv(pd.DataFrame(agreement), tables_dir(cfg) / "tfidf_matrix_agreement.csv", status=STATUS,
                  note="descriptive Spearman correlation of upper triangles, no p-values")


if __name__ == "__main__":
    main()
