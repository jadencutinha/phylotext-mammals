# Week 6 performance: full-scale permutation tests

Measured by `scripts/22_fullscale_tests.py --benchmark` on arm64 Darwin (arm), 8 worker threads, on the real full-scale matrices: 4,611 species in the frozen set (10,628,355 pairs) and 4,519 H2 species (10,208,421 pairs).

## Design (amendment 3, section 5)

- Each matrix is ranked once (Spearman, average ranks) and stored as a centred unit-norm condensed vector (float64 for the observed statistic).
- The permuted matrix is held as a float32 square; each permutation gathers its condensed vector through an int32 index array and multiplies it by every fixed vector at once (one float32 matrix-vector product). Permutations are drawn in chunks of 256 from one generator, identical to the unchunked draw, and run on a thread pool; each statistic is computed on its own, so results do not depend on chunk size or thread count.
- Partial Mantel and MRM statistics are computed from correlations alone (the residualization is algebraic), so the controls are prepared once.
- Tests that permute the same matrix on the same species share one pass; the 999-permutation robustness tests use the first 999 permutations of the 9,999-permutation pass, which is exactly the separate test with the same seed.
- Patristic distances come from a vectorized version of dendropy's distance calculation (bit-identical to dendropy; 19 s for 5,583 tips instead of hours), and each posterior tree's distances are saved to disk one tree at a time.

## Measurements

| quantity | value |
|---|---:|
| ranking one 4,611-species matrix and standardizing it | 1.8 s |
| one permutation, 1 fixed vector, 4,611 species | 25.1 ms |
| one permutation, 4 fixed vectors, 4,611 species | 32.1 ms |
| one permutation, 9 fixed vectors, 4,611 species | 38.9 ms |

## Projection for the registered run

| pass | permutations | projected time |
|---|---:|---:|
| H1: phylogeny moving, 4 text matrices | 9,999 | 5.4 min |
| H1 Pearson | 999 | 0.4 min |
| H2: phylogeny moving, 4 texts + 5 controls | 9,999 | 6.2 min |
| MRM primary: text moving, 6 predictors | 9,999 | 6.2 min |
| MRM, 3 text variants | 2,997 | 1.5 min |
| posterior: 20 trees x (H1 + H2) | 20 x 2 x 999 | 18.6 min |

Total projected for the pooled passes above: about 38 minutes, plus per-order tests (smaller n) and the exploratory pair-subset tests, which re-rank the subset in every permutation. No GPU is needed for the statistics. Embeddings ran on the local Apple GPU (MPS); `scripts/slurm/embed_adroit.sh` runs the same step on Princeton Adroit if a CUDA node is preferred.

## Measured on the registered run (2026-11-05)

| step | time |
|---|---:|
| embeddings, bge-large / chunk / strict (5,583 texts, 6,397 passages, MPS) | 206 s |
| embeddings, bge-large / truncate / strict (5,583 passages) | 166 s |
| embeddings, bge-large / chunk / taxonomy (6,350 passages) | 202 s |
| embeddings, all-mpnet-base-v2 / chunk / strict (6,980 passages) | 84 s |
| pruning the MCC tree to 5,583 tips and all patristic distances | 19 s |
| pruning 20 posterior trees and saving their distances | under 60 s |
| building the five frozen-set matrices | 5 s |
| exact and full-vocabulary replication (including re-embedding 284 pilot texts) | 42 s |
| H1 pass, 9,999 permutations, 4 text matrices | 224 s |
| H2 pass, 9,999 permutations, 9 fixed vectors | 293 s |
| MRM pass, 9,999 permutations | 247 s |

The measured passes took 20 to 35% less time than projected from the 100-permutation benchmark.
