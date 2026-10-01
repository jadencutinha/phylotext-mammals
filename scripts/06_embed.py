from __future__ import annotations

import time

from lineage_charisma.config import base_parser, load_config
from lineage_charisma.corpus import read_corpus
from lineage_charisma.embed import RULES, SentenceTransformerEncoder, cache_is_current, embed_texts, embedding_paths, save_embeddings
from lineage_charisma.io_utils import sha256
from lineage_charisma.masking import text_column


def main() -> None:
    parser = base_parser("Embed every species description with each configured model, long-text rule, and mask level.")
    parser.add_argument("--model", action="append", help="only this model (repeatable); default: all in config")
    parser.add_argument("--rule", action="append", choices=RULES, help="only this rule (repeatable); default: all in config")
    parser.add_argument("--mask", action="append", help="only this mask level (repeatable); default: all in config")
    args = parser.parse_args()
    cfg = load_config(args.config)
    ecfg = cfg["embedding"]
    models = args.model or ecfg["models"]
    rules = args.rule or ecfg["rules"]
    masks = args.mask or cfg.mask_levels
    unknown = [r for r in rules if r not in RULES]
    if unknown:
        raise SystemExit(f"unknown rules in config: {unknown}; expected a subset of {list(RULES)}")

    corpus_path = cfg.corpus_path()
    if not corpus_path.exists():
        raise SystemExit(f"{corpus_path} not found; run scripts/05b_build_corpus.py (make corpus) first")
    corpus = read_corpus(corpus_path)
    absent = [k for k in masks if text_column(k) not in corpus.columns]
    if absent:
        raise SystemExit(f"{corpus_path.name} has no text for mask levels {absent}; rerun scripts/05b_build_corpus.py --force")
    species = corpus["species_id"].tolist()
    corpus_sha = sha256(corpus_path)
    print(f"[embed] corpus: {len(species)} species ({int(corpus['is_stub'].sum())} flagged as stubs, still embedded)")

    for model_name in models:
        todo = []
        for mask, rule in [(k, r) for k in masks for r in rules]:
            paths = embedding_paths(cfg.embeddings_dir(), model_name, rule, mask)
            params = {"model": model_name, "rule": rule, "mask_level": mask, "max_tokens_config": ecfg.get("max_tokens"), "corpus_sha256": corpus_sha}
            if rule == "chunk":
                params["chunk_overlap"] = ecfg["chunk_overlap"]
            current, why = cache_is_current(paths, params)
            if current and not args.force:
                print(f"[cache] {paths.array.name} exists, skipping (use --force to rebuild)")
            else:
                if not args.force and why != "not cached":
                    print(f"[embed] {paths.array.name} is stale: {why}; rebuilding")
                todo.append((rule, mask, paths, params))
        if not todo:
            continue

        encoder = SentenceTransformerEncoder(model_name, device=ecfg.get("device", "auto"), max_tokens=ecfg.get("max_tokens"))
        print(f"[embed] {model_name} on {encoder.device}, max_seq_length={encoder.max_seq_length}")
        for rule, mask, paths, params in todo:
            started = time.time()
            texts = corpus[text_column(mask)].tolist()
            result = embed_texts(texts, encoder, rule, overlap=ecfg["chunk_overlap"], batch_size=ecfg["batch_size"])
            save_embeddings(
                paths, result.embeddings, species, **params,
                max_seq_length=encoder.max_seq_length,
                device=encoder.device,
                n_passages_total=sum(result.n_passages),
                n_species_multi_passage=sum(n > 1 for n in result.n_passages),
                n_species_truncated=sum(used < total for used, total in zip(result.n_tokens_used, result.n_tokens)),
                n_passages=dict(zip(species, result.n_passages)),
                n_tokens=dict(zip(species, result.n_tokens)),
            )
            print(f"[embed] {paths.array.name}: {result.embeddings.shape}, {sum(result.n_passages)} passages, "
                  f"{sum(n > 1 for n in result.n_passages)} species with >1 passage, {time.time() - started:.0f}s")


if __name__ == "__main__":
    main()
