from __future__ import annotations

import os

import pandas as pd

from lineage_charisma.config import load_config
from lineage_charisma.corpus import read_corpus
from lineage_charisma.distances import load_matrix, text_matrix_path
from lineage_charisma.embed import cache_is_current, embedding_paths, load_embeddings, read_species_order, save_embeddings
from lineage_charisma.extensions import (
    COST_STOP_USD,
    PROVIDERS,
    STATUS,
    ApiEmbedder,
    amendment_is_committed,
    cost_estimate,
    embed_with_provider,
    load_registered,
    matrix_agreement,
    provider_tokenizer,
    registered_tests,
    tables_dir,
    text_distance,
    work_dir,
)
from lineage_charisma.io_utils import sha256, write_csv

# docs/preregistration_amendment_extensions.md, section A (secondary, designed after the 2026-11-11 lock).
# Default action: the cost estimate, which reads local files only. --embed sends requests; --test computes statistics from the cache.
SCALES = {"pilot": None, "full": "mammalia"}
OVERLAP = 32
MPNET = ("sentence-transformers/all-mpnet-base-v2", "chunk", "strict")


def api_paths(cfg, provider):
    return embedding_paths(work_dir(cfg), provider.model, "api_chunk", "strict")     # the passage limit is in the .meta.json


def strict_texts(cfg) -> tuple[list[str], list[str]]:
    order = read_species_order(cfg.species_order_path())
    corpus = read_corpus(cfg.corpus_path()).set_index("species_id")
    return order, corpus.loc[order, "text_strict"].tolist()


def estimate(cfgs) -> None:
    import tiktoken

    enc = tiktoken.get_encoding("cl100k_base")
    tokens = {scale: sum(len(enc.encode(t)) for t in strict_texts(cfg)[1]) for scale, cfg in cfgs.items()}
    rows = pd.DataFrame([r for p in PROVIDERS.values() for r in cost_estimate(tokens, p)])
    print(rows.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    print(f"total: US${rows['usd'].sum():.2f} at list prices, US${rows['usd_pessimistic'].sum():.2f} pessimistic "
          f"(exact cl100k_base token counts)")
    for p in PROVIDERS.values():
        print(f"{p.key_env}: {'set' if os.environ.get(p.key_env) else 'not set'}")


def embed(cfgs, providers, force: bool) -> None:
    for name in providers:
        provider = PROVIDERS[name]
        if not os.environ.get(provider.key_env):
            raise SystemExit(f"[api] {provider.key_env} is not set; nothing sent")
        tokenizer, limit = provider_tokenizer(provider)
        print(f"[api] {provider.name} {provider.model}: tokenizer {tokenizer.name}, passage limit {limit}")
        spent = 0.0
        for scale, cfg in cfgs.items():
            order, texts = strict_texts(cfg)
            paths = api_paths(cfg, provider)
            params = {"model": provider.model, "provider": provider.name, "corpus_sha256": sha256(cfg.corpus_path()), "limit": limit, "overlap": OVERLAP}
            if cache_is_current(paths, params)[0] and not force:
                print(f"[cache] {paths.array.name} exists, skipping")
                continue
            embedder = ApiEmbedder(provider, cfg.path("interim", shared=True) / "api_embeddings")
            result = embed_with_provider(texts, embedder, tokenizer, limit, OVERLAP)
            cost = embedder.tokens_billed / 1e6 * provider.usd_per_mtok
            spent += cost
            save_embeddings(paths, result.embeddings, order, **params, tokenizer=tokenizer.name, status=STATUS,
                            n_species_multi_passage=sum(n > 1 for n in result.n_passages), requests_sent=embedder.requests_sent,
                            tokens_billed=embedder.tokens_billed, usd=cost)
            print(f"[api] {scale}: {result.embeddings.shape}, {embedder.requests_sent} requests, {embedder.tokens_billed:,} tokens billed, "
                  f"US${cost:.3f}; {sum(n > 1 for n in result.n_passages)} texts needed more than one passage")
            if scale == "pilot" and "full" in cfgs and cost > 0:
                full_tokens = sum(len(tokenizer.encode(t)) for t in strict_texts(cfgs["full"])[1])
                pilot_tokens = sum(len(tokenizer.encode(t)) for t in texts)
                projected = spent + cost * full_tokens / pilot_tokens
                if projected > COST_STOP_USD:
                    raise SystemExit(f"[api] projected cost US${projected:.2f} is above the US${COST_STOP_USD:.0f} stop; full scale not sent")


def test(cfgs) -> None:
    for scale, cfg in cfgs.items():
        s = cfg["stats"]
        perms, seed, alpha, method = int(s["permutations"]), int(s["seed"]), float(s["alpha"]), s["method"]
        reg = load_registered(cfg)
        rows, matrices = [], {}
        order_path = cfg.species_order_path()
        matrices["bge-large (primary)"] = load_matrix(text_matrix_path(cfg.matrices_dir(), *cfg.primary_combo), order_path)[0]
        matrices["mpnet"] = load_matrix(text_matrix_path(cfg.matrices_dir(), *MPNET), order_path)[0]
        for provider in PROVIDERS.values():
            paths = api_paths(cfg, provider)
            if not paths.array.exists():
                print(f"[api] {scale}: no {provider.name} embeddings yet; run --embed first")
                continue
            emb, species = load_embeddings(paths)
            if species != reg.order:
                raise SystemExit(f"[api] {scale}: {provider.name} embeddings are not in the frozen species order")
            text = text_distance(emb)
            matrices[provider.model] = text
            for row in registered_tests(text, reg, perms, seed, alpha, method):
                rows.append({"status": STATUS, "scale": scale, "provider": provider.name, "model": provider.model, **row})
                print(f"[api] {scale} {provider.model:24s} {row['test']:30s} r = {row['r']:+.4f} p = {row['p']:.4f}")
        if not rows:
            continue
        meta = dict(status=STATUS, plan="docs/preregistration_amendment_extensions.md, section A", permutations=perms, seed=seed,
                    alternative="greater", permuted="phylogeny", texts="text_strict of the scale's corpus")
        write_csv(pd.DataFrame(rows), tables_dir(cfg) / "api_embeddings_tests.csv", **meta)
        write_csv(pd.DataFrame(matrix_agreement(matrices, method)), tables_dir(cfg) / "api_embeddings_matrix_agreement.csv",
                  status=STATUS, note="descriptive Spearman correlation of upper triangles, no p-values")


def main() -> None:
    from lineage_charisma.config import base_parser

    parser = base_parser("Amendment 5, section A (and its Change 1): commercial embeddings with OpenAI text-embedding-3-large. Default: print the cost estimate only.")
    parser.add_argument("--provider", action="append", choices=list(PROVIDERS), help="default: openai, the only provider (amendment 5, Change 1)")
    parser.add_argument("--scale", action="append", choices=list(SCALES), help="default: pilot then full")
    parser.add_argument("--embed", action="store_true", help="send the embedding requests (needs the API keys and the committed amendment)")
    parser.add_argument("--test", action="store_true", help="run the registered tests on the cached embeddings")
    args = parser.parse_args()
    cfgs = {s: load_config(args.config, SCALES[s]) for s in (args.scale or list(SCALES))}
    root = next(iter(cfgs.values())).root
    if not (args.embed or args.test):
        estimate(cfgs)
        return
    ok, why = amendment_is_committed(root)
    if not ok:
        raise SystemExit(f"[api] {why}")
    if args.embed:
        embed(cfgs, args.provider or list(PROVIDERS), args.force)
    if args.test:
        test(cfgs)


if __name__ == "__main__":
    main()
