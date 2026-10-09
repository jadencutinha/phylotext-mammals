from __future__ import annotations

import platform
import time

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from lineage_charisma.config import base_parser, load_config
from lineage_charisma.corpus import read_corpus
from lineage_charisma.io_utils import write_csv
from lineage_charisma.stats import _standardized_square, draw_permutations, permuted_correlations, permuted_statistics, upper_triangle

# sizing for Week 6: nothing here touches the project's results, it only times the code on matrices of the full-mammal size
N_SCALE = 5911          # tips in the Upham et al. tree
EMBED_PASSAGES = 64


def random_square(n: int, rng: np.random.Generator) -> np.ndarray:
    values = rng.random(n * (n - 1) // 2)
    return _standardized_square(values, n, "benchmark")


def seconds_per_permutation(squares: list[np.ndarray], n: int, repeats: int, batch_size: int) -> float:
    perms = draw_permutations(n, repeats, 0)
    start = time.perf_counter()
    if len(squares) == 2:
        permuted_statistics(squares[0], squares[1], perms, batch_size)
    else:
        permuted_correlations(squares[:-1], squares[-1], perms, batch_size)
    return (time.perf_counter() - start) / repeats


def main() -> None:
    parser = base_parser("Time the Mantel machinery and the embedding model at full-mammal size, for the Week 6 scaling decision.")
    parser.add_argument("--skip-embedding", action="store_true", help="do not load the embedding model")
    args = parser.parse_args()
    cfg = load_config(args.config)
    perms, ecfg = int(cfg["stats"]["permutations"]), cfg["embedding"]
    rng = np.random.default_rng(0)
    rows = []

    def add(quantity: str, value: float, unit: str) -> None:
        rows.append({"quantity": quantity, "value": value, "unit": unit})
        print(f"[benchmark] {quantity}: {value:,.2f} {unit}")

    pairs = N_SCALE * (N_SCALE - 1) // 2
    add(f"species pairs at n = {N_SCALE}", pairs / 1e6, "million")
    add("one distance matrix, float64", N_SCALE**2 * 8 / 1e9, "GB")
    add("one distance matrix, float32", N_SCALE**2 * 4 / 1e9, "GB")
    add("permutation batch of 128 at float64 (current default)", 128 * N_SCALE**2 * 8 / 1e9, "GB")

    small = [random_square(284, rng) for _ in range(4)]
    add("Mantel test at n = 284, 9,999 permutations", seconds_per_permutation(small[:2], 284, 2000, 128) * perms, "seconds")
    add("partial Mantel (2 controls) at n = 284, 9,999 permutations", seconds_per_permutation(small, 284, 2000, 128) * perms, "seconds")

    big = [random_square(N_SCALE, rng) for _ in range(4)]
    start = time.perf_counter()
    rankdata(upper_triangle(big[0]))
    add(f"ranking one matrix at n = {N_SCALE}", time.perf_counter() - start, "seconds")
    one = seconds_per_permutation(big[:2], N_SCALE, 6, 2)
    add(f"Mantel test at n = {N_SCALE}, per permutation", one, "seconds")
    add(f"Mantel test at n = {N_SCALE}, 9,999 permutations", one * perms / 3600, "hours")
    add(f"Mantel test at n = {N_SCALE}, 999 permutations", one * 999 / 60, "minutes")
    three = seconds_per_permutation(big, N_SCALE, 6, 2)
    add(f"partial Mantel (2 controls) at n = {N_SCALE}, 9,999 permutations", three * perms / 3600, "hours")
    del big

    if not args.skip_embedding:
        from lineage_charisma.embed import SentenceTransformerEncoder

        corpus = read_corpus(cfg.corpus_path())
        add("passages per species in the Carnivora corpus (chunk rule, 512 tokens)", float(np.ceil(corpus["n_tokens"] / 478).clip(lower=1).mean()), "passages")
        for model in ecfg["models"]:
            encoder = SentenceTransformerEncoder(model, ecfg.get("device", "auto"), ecfg.get("max_tokens"))
            ids = encoder.tokenize(" ".join(corpus["text"].head(40)))
            width = encoder.max_content_tokens
            passages = [ids[i * width:(i + 1) * width] for i in range(EMBED_PASSAGES)]
            encoder.encode_ids(passages[:8], int(ecfg["batch_size"]))      # warm-up
            start = time.perf_counter()
            encoder.encode_ids(passages, int(ecfg["batch_size"]))
            add(f"{model.split('/')[-1]}: full-length passages per second on {encoder.device}", EMBED_PASSAGES / (time.perf_counter() - start), "passages/s")

    write_csv(pd.DataFrame(rows), cfg.tables_dir() / "scaling_benchmark.csv", machine=f"{platform.machine()} {platform.system()}, {platform.processor() or 'unknown CPU'}",
              n_scale=N_SCALE, note="timings on random matrices of the stated size; nothing here uses project results")


if __name__ == "__main__":
    main()
