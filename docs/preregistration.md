# Pre-registration: lineage or charisma, Carnivora pilot

Written 2026-10-04 21:09 UTC, before any Mantel test or other correlation between text distance and phylogenetic distance was computed. The commit that adds this file is the timestamp of record. Changes after that commit go in the "Amendments" section at the end, with a date and a reason, and never by editing the sections above it.

## 1. What had been seen before this was written

The analysis plan is fixed before any test, but not before any look at the data. The Week 3 sanity report (`reports/week3_sanity.md`, `reports/week3_notes.md`) is descriptive and was read before this plan was written. It showed:

- scatter plots of text distance against phylogenetic distance, and median text distance by divergence-time bin, at the taxonomy mask level;
- within-family against between-family text distance at all three mask levels, for all four model and rule combinations;
- nearest-neighbour lists for ten showcase species.

No correlation, Mantel statistic or p-value between text distance and phylogenetic distance has been computed. The Week 3 report computed rank correlations only between pairs of text matrices.

**The primary mask level was changed from taxonomy to strict after reading that report.** The reason is that taxonomy masking behaved like no masking (family structure unchanged), while strict masking removed most of the family structure that shared common names had carried. Strict is therefore the harder test of H1, and the choice was made without knowing the H1 statistic at any mask level.

## 2. Data

- **Species**: the 284 Carnivora species in `data/processed/matrices/species_order.txt` (285 matched species less one stub under 50 tokens). Every matrix uses this order.
- **Text distance**: 1 − cosine similarity between embeddings of each species' English Wikipedia description text.
- **Phylogenetic distance**: patristic distance in million years on the pruned MCC tree of Upham et al. (2019), `data/processed/matrices/phylo.npy`.

## 3. Primary specification

| setting | value |
|---|---|
| embedding model | `BAAI/bge-large-en-v1.5` |
| long-text rule | chunk |
| mask level | strict (scientific names and common-name head nouns replaced by `[TAXON]`) |
| correlation | Spearman |
| species | all 284 |

## 4. H1

**H1: text distance is positively associated with phylogenetic distance.** Mantel r > 0 between the text-distance and phylogenetic-distance matrices.

- **Statistic**: Spearman correlation between the upper triangles of the two matrices (284 × 283 / 2 = 40,186 pairs), with average ranks for ties.
- **Null distribution**: 9,999 permutations. Each permutation reorders the rows and columns of the phylogenetic matrix together with one random permutation of the species, and the statistic is recomputed.
- **p-value**: one-sided, p = (1 + number of permuted r ≥ observed r) / (1 + 9,999). The smallest possible p is 0.0001.
- **Alpha**: 0.05.
- **Seed**: 20261004 (`stats.seed` in `config.yaml`). Every test starts a fresh generator from this seed, so each result can be reproduced on its own and does not depend on the order tests are run in.
- **Decision rule**: H1 is supported if the primary specification gives r > 0 with p < 0.05. Only the primary specification decides H1. The size of r is reported and interpreted alongside the p-value, since with 40,186 pairs a small r can be significant.

## 5. Robustness variants

The same test is run on all 12 combinations of model, rule and mask level. The 11 that are not the primary specification are robustness variants:

| model | rule | mask levels |
|---|---|---|
| bge-large-en-v1.5 | chunk | none, taxonomy (strict is primary) |
| bge-large-en-v1.5 | truncate | none, taxonomy, strict |
| all-mpnet-base-v2 | chunk | none, taxonomy, strict |
| all-mpnet-base-v2 | truncate | none, taxonomy, strict |

Pearson correlation is reported next to Spearman for all 12 as a further variant, with the same permutations and seed.

The variants are reported in full, whatever they show. No multiple-testing correction is applied and no variant is used to claim H1; the report states how many agree with the primary result in sign and in significance at 0.05, and how r changes from none to taxonomy to strict masking.

## 6. Descriptive checks (no confirmatory claims)

- **Within families**: the primary-specification Mantel test repeated inside each of the three largest families in the 284-species set, fixed here by species count: Mustelidae (58), Felidae (38), Canidae (37). These are reported as description only, with no multiple-testing claims.
- **Confound diagnostics**: Spearman Mantel correlations (9,999 permutations, same seed) between phylogeny and ecology, phylogeny and each attention matrix, and ecology and each attention matrix; and between primary-specification text distance and ecology and each attention matrix. These are marginal correlations. Nothing is partialled in Week 4.

## 7. Planned H2 controls

H2 asks whether the text–phylogeny association survives controlling for ecology and for human attention. The control matrices are specified here and built in Week 4. The H2 tests themselves (partial Mantel and MRM) are Week 5; their exact model will be added as a dated amendment before any of them is run.

### Ecological distance

Gower distance over these traits and no others:

| group | traits | source | type |
|---|---|---|---|
| diet | the ten `Diet-*` percentage columns | EltonTraits 1.0 (Wilman et al. 2014) | continuous, 0–100 |
| foraging stratum | `ForStrat-Value` | EltonTraits 1.0 | categorical |
| activity time | nocturnal, crepuscular, diurnal flags | EltonTraits 1.0 | binary, not exclusive |
| habitat breadth | habitat breadth | COMBINE (Soria et al. 2021) | continuous |
| terrestrial or aquatic | terrestrial, marine and freshwater flags | COMBINE | binary |
| body mass | log10 adult body mass | COMBINE | continuous |

- Continuous traits contribute |difference| / range over the 284 species; categorical and binary traits contribute 0 for a match and 1 for a mismatch.
- **The six groups are weighted equally.** Columns inside a group share that group's weight, so the ten diet columns together count as much as body mass alone. Without this, diet would carry ten of roughly eighteen columns.
- **Missing values** are handled by Gower's pairwise-available rule: a pair's distance is the weighted mean over the traits both species have. Nothing is imputed. Missing-data rates are reported per trait.
- A species with no trait data at all, or a pair with no trait in common, cannot be given an ecological distance. Such species are listed, and analyses that use ecology are run on the species that remain. H1 always uses all 284.
- Exact column names are confirmed against the downloaded files and recorded in `docs/ecology_traits.md`. If a listed trait does not exist in the source as described, that is recorded as an amendment rather than replaced silently.

### Attention distance

Two variables per species:

- **log(1 + median monthly pageviews)** of the species' English Wikipedia article, from the Wikimedia Pageviews API: 2023-01 through 2025-12, all access methods, agent = user, with pageviews of redirects to the article included.
- **log(article length in tokens)**: `n_tokens` in the corpus (unmasked text, bge-large tokenizer).

Three matrices, all kept:

- pageview distance: |difference in log pageviews|;
- length distance: |difference in log length|;
- combined attention distance: Euclidean distance on the two variables after each is standardised to mean 0 and standard deviation 1.

A fourth pairwise quantity, the mean log pageviews of the two species, is built to examine later whether pairs of well-known species are described alike. It is a similarity-type covariate, not a distance, and has no confirmatory role.

## 8. Not covered by this plan

- The 277-species matrices that drop texts under 100 tokens (`matrices_min100/`) are not part of H1 or its robustness variants. Any analysis of them is exploratory and will be labelled so.
- The MCC tree is used alone; uncertainty across the posterior sample of trees is not analysed here.
- Anything else run that is not listed above is exploratory and will be labelled so in the report.

## 9. Amendments

None.
