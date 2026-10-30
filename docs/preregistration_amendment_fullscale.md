# Amendment 3 to the pre-registration: full-scale analysis, all mammals (Week 6)

Written 2026-10-29, before any full-scale species set, text, embedding, distance matrix or test existed. The commit that adds this file is the timestamp of record. It amends `docs/preregistration.md` (sections 1 to 8, Amendment 1) and `docs/preregistration_amendment_h2.md` (Amendment 2). Both stand for the Carnivora pilot, and they apply to the full-scale analysis unless an item below says otherwise. Later changes go in the "Changes" section at the end of this file, with a date and a reason.

## 1. What had been seen before this was written

- Every pilot result in `reports/week4_report.md` and `reports/week5_report.md`. On the primary specification these were: H1 r = 0.033, one-sided p = 0.0285, 284 species, with median r = 0.080 over 100 posterior trees. For H2, raw r = 0.035 and partial r (ecology + length) = 0.020, p = 0.127, which is 57.5% of the raw r, so neither rule (a) nor rule (b) holds. In the MRM, ecology and length were significant and phylogeny was not. Exploratory results: pairs in the top pageview quartile are described more alike, and the text–phylogeny correlation is 0.133 in the bottom quartile and 0.035 in the top.
- The Week 5 scaling benchmark (`scripts/18_scaling_benchmark.py`, section 8 of the Week 5 report), which timed the permutation code on synthetic matrices of full size.
- The size of the source tables: MDD v2.5 lists 6,904 species, the Upham tree has 5,911 tips, and 320 MDD species are in Carnivora.

No full-scale name reconciliation, Wikipedia text, stub count, trait coverage, pageview count, embedding or distance matrix exists. No correlation of any kind has been computed on species outside the pilot set.

## 2. Role of the full-scale analysis

The full-scale analysis is a separate pre-registered test of H1 and H2, run on a larger and different species set. The pilot results stand as recorded and are not re-decided. The report gives the two side by side. Neither result overrides the other, and the report says so where they disagree.

## 3. Species

The full-scale species set is every species that passes all of these steps:

1. **MDD.** Listed in MDD v2.5 (the file the pilot used).
2. **Tree.** Matched to an Upham MCC tip with the existing reconciliation engine (`scripts/03_reconcile_names.py`, `matching` settings in `config.yaml` unchanged). Exact and synonym matches are accepted automatically. Fuzzy matches at or above `fuzzy_accept_threshold` (92) with the required margin are accepted as in the pilot. Matches between the review and accept thresholds, and species with no match, go to `data/interim/mammalia/residual_review.csv`. A species there enters only if it is resolved by a rule written in the Changes section before the set is frozen. Otherwise it is excluded, and the exclusion is listed by order. Extinct and domestic species without their own tip are excluded, as in the pilot.
3. **Wikipedia.** Has an English Wikipedia article reached from the matched name. The text is built by the same rule as the pilot: the cleaned lead plus the first section whose heading is in `wikipedia.description_headings`, joined.
4. **Stub rule.** The unmasked text has at least 50 tokens (bge-large tokenizer, special tokens not counted). Species under 50 are stubs and are dropped from every matrix. The stub count is reported by order.

**Freezing.** The final list is written to `data/processed/mammalia/matrices/species_order.txt` and committed after these steps and before any test is run. Every full-scale matrix uses this order. The species count and the list of orders with at least 50 species (section 7) are then added to the Changes section of this file in the same commit, without changes to anything else.

**H2 species.** These are the species from the frozen list that have an ecological distance to every other species. They are found by the existing rule (`traits.complete_species`: drop the species with the most undefined pairs, repeat) and written to `ecology_species.txt`. As in Amendment 2, the raw r is recomputed on this set and is the H2 baseline.

## 4. Matrices

These are the same definitions as the pilot (Amendment 2, section 4), built on the full species set:

- **Text distance**: 1 − cosine similarity of the primary embeddings.
- **Phylogenetic distance**: patristic distance in My on the MCC tree pruned to the frozen set.
- **Ecological distance**: Gower over the same traits and groups (`docs/ecology_traits.md`), groups weighted equally, pairwise-available missing-data rule, no imputation. Continuous traits are scaled by their range over the full-scale species. A trait column that is constant within Carnivora but varies across mammals is used. In particular, `terrestrial_volant` (bats) joins the terrestrial-or-aquatic group if COMBINE has values for it. That column is listed in `docs/ecology_traits.md` when the matrix is built. Because flight status nearly coincides with Chiroptera, this column carries lineage information as well as ecology (section 10). The primary H2 partial r is therefore also reported with an ecology matrix built without `terrestrial_volant`, as a declared robustness variant (section 8). If COMBINE has no values for the column, neither matrix uses it and the variant is not run.
- **Length distance** and **pageview distance**: as before. Pageviews cover the same window (2023-01 to 2025-12, all access, agent = user, redirects included).

**Masking.** The taxonomy and strict mask vocabularies are rebuilt from MDD for all mammals by the same procedure. The strict level adds the common-name head nouns that the larger set needs (for example mouse, rat, bat, shrew, monkey, deer, whale), drawn from MDD common names by the existing rule. The check that fails if any masked text still contains its own genus name stays. `extra_taxonomy_terms`, `keep_terms` and `keep_when_followed_by` are extended for the new orders before any embedding is computed, and the final lists are recorded in `config.yaml` and `data/interim/mammalia/mask_terms.csv`.

## 5. H1 and H2: unchanged settings

| setting | value |
|---|---|
| primary specification | `BAAI/bge-large-en-v1.5` / chunk / strict mask |
| correlation | Spearman on the upper triangles, average ranks for ties, untransformed distances |
| permutations | 9,999, rows and columns of the phylogenetic matrix permuted together (text matrix for MRM) |
| seed | 20261004, a fresh generator per test |
| p-value | one-sided (1 + count ≥ observed) / 10,000; MRM two-sided on pseudo-t |
| alpha | 0.05 |
| H1 decision | supported if r > 0 with p < 0.05 on the primary specification and the MCC tree |
| primary H2 model | partial Mantel, text on phylogeny, controlling for ecology + length jointly |
| H2 rule (a) | partial r > 0 with one-sided p < 0.05 |
| H2 rule (b) | partial r < 50% of the raw r on the H2 species, MCC tree; not defined if raw r ≤ 0 |
| also reported for H2 | raw − partial r, and the partial-to-raw ratio across posterior trees (Amendment 2, Change 1) |
| secondary H2 models | S1 ecology, S2 length, S3 pageviews, S4 ecology + length + pageviews (descriptive) |
| MRM | text ~ phylogeny + ecology + length + pageviews, ranked and standardized (descriptive) |

Spearman on untransformed distance stays primary. The Week 5 nonlinear result is exploratory and does not change the primary correlation.

**Implementation is not a change of method.** At this scale the code works on condensed upper-triangle vectors in float32, ranks each matrix once, and permutes through index arrays in chunks. Before any full-scale test is run, tests check that this code gives the same statistic as the square-matrix pilot code on small inputs, and the same permutation distribution with the same seed for any chunk size. If the numbers differ, the pilot code is the reference.

## 6. Replication check (before any full-scale result is read)

There are two checks, run in this order. Both use the pilot settings: Spearman, 9,999 permutations, seed 20261004.

**6a. Exact replication with the pilot vocabulary.** The full-scale pipeline (reconciliation, text building, masking, embedding, distance matrices and the condensed-vector Mantel code) is run on the 284 pilot species (`data/processed/matrices/species_order.txt`), using the pilot's own mask vocabulary (`data/interim/mask_terms.csv`) in place of the full-scale one. This must reproduce the pilot H1 r = 0.033 and its p-value exactly. Any difference is a bug in the full-scale code or data path and is fixed before going on. Embeddings may differ in the last bits of floating-point precision if they are computed on different hardware. If so, the check requires the Spearman r to agree to the reported precision and the p-value to be identical, and any remaining floating-point difference is reported.

**6b. Replication with the full-scale vocabulary.** The full-scale primary text and phylogenetic matrices are cut to the 284 pilot species, H1 is recomputed, and the result is compared with the pilot r = 0.033. Because 6a has shown the pipeline is the same, any gap here is attributed to the change in mask vocabulary.

- The phylogenetic distances should be identical: patristic distance between two tips does not depend on which other tips are kept. Any non-zero difference is a bug and is fixed before going on.
- Text distances can differ because the strict mask vocabulary is larger. A term added for another order can now mask a word in a carnivore's text. The report gives the pilot r, the replication r, the difference, the number of pilot species whose masked text changed, and the terms responsible.
- A pilot species that is missing from the full-scale set (for example, a changed match) is listed. The check is then run on the species the two sets share, and the pilot r is recomputed on that same subset for comparison.
- If |replication r − pilot r| > 0.01, the cause is identified and written up before the full-scale H1 is read. The full-scale analysis then goes ahead on the full-scale masks as specified here, and the pilot masks are not used.

## 7. Per-order breakdowns (descriptive)

H1 and the primary H2 model (raw and partial r, both with the settings in section 5) are run within each mammalian order that has at least 50 species in the frozen set. The list of orders is fixed in the Changes section when the set is frozen, using MDD's `order` column. Each order uses its own species from the frozen list (and from the H2 species list for H2), with distances taken from the full-scale matrices as they are. Nothing is re-masked or re-scaled within an order. Reported per order: n, raw r and one-sided p, partial r and one-sided p, and the partial-to-raw ratio. These are descriptions. There is no multiple-testing correction and no claim that any order supports or contradicts H1 or H2. Carnivora appears in this table. Its value differs from the pilot because the masks and Gower ranges are full-scale, and the replication check in section 6 accounts for the first difference.

## 8. Robustness (reported in full, not decisive)

These are one-at-a-time changes from the primary specification. For each one: H1 r and p on the frozen set, and the raw r, the primary partial r and the MRM on the H2 species. There is no correction for multiple testing, and none of these changes the outcome of H1 or of rules (a) and (b).

| variant | change |
|---|---|
| second model | `sentence-transformers/all-mpnet-base-v2` in place of bge-large (chunk, strict) |
| truncate rule | truncate in place of chunk (bge-large, strict) |
| taxonomy mask | taxonomy in place of strict (bge-large, chunk) |
| unweighted Gower | `ecology_unweighted.npy` in place of group-weighted (H2 only), per A1.1 |
| ecology without flight | group-weighted Gower built without `terrestrial_volant` (H2 only): raw r, primary partial r and its share of raw r, and MRM |
| Pearson | Pearson next to Spearman for the primary H1 test |

The full 12-variant grid is not run at full scale.

**Posterior trees.** 20 trees, the first 20 entries of the 100 indices drawn in A1.3: `numpy.random.default_rng(20261004).choice(N, size=100, replace=False)[:20]`. These are a subset of the pilot's trees. For each tree, the tree is pruned to the frozen set and patristic distances are computed and saved to disk one tree at a time. Then: H1 raw r on the frozen set, and raw r and primary partial r on the H2 species. Reported as distributions (median, 2.5th and 97.5th percentiles, minimum, maximum) next to the MCC values, with the share of trees meeting p < 0.05 for H1, rule (a) and rule (b), and the partial-to-raw ratio summarised as in Amendment 2, Change 1. With 20 trees the percentiles are coarse, and the report says so.

**Permutations.** The primary H1 test, the primary H2 model, the secondary H2 models, the MRM, the per-order tests and the replication check use 9,999 permutations. The robustness variants and the posterior trees use 999, with the same seed. With 999 the smallest possible p is 0.001.

## 9. Exploratory analyses

These are labelled exploratory in tables, figures and the report. None of them bears on H1 or H2.

- **Pageview quartiles.** The same analysis as the pilot (Amendment 2, section 9, "Both well-known") on the frozen set: pairs split into quartiles of mean log pageviews, with mean and median text distance per quartile and the text–phylogeny Spearman correlation within the top and bottom quartiles. The permutation p reorders the phylogenetic matrix with the set of pairs held fixed, using 999 permutations. The same analysis is also reported within each order from section 7.
- **30 My cap.** The phylogenetic matrix capped at 30 My, as in the pilot, with H1 and the primary H2 model rerun with Spearman and Pearson, 999 permutations. The share of pairs at or below 30 My is reported. Across all mammals it will be much smaller than in Carnivora (about 13%), so many more pairs tie at the cap. The log-transform check is run with Pearson only, as in the pilot.

## 10. Known limitations, stated in advance

- **Everything in Amendment 2, section 10, still applies**, including inflated partial Mantel type I error under autocorrelation. With about 17 million pairs, very small correlations will be significant. Effect sizes are reported next to every p-value and are interpreted on their own terms.
- **Pairs are not independent.** The Mantel permutation accounts for this. The number of pairs still overstates how much information there is: the information is closer to the number of species than to the number of pairs.
- **Stubs are not random.** Wikipedia coverage is expected to be thinner for rodents, bats and shrews. The stub rule therefore removes species non-randomly with respect to both phylogeny and attention. The frozen set over-represents well-documented species, and the stub table by order is reported as a result in its own right.
- **Name matching at scale.** Residual cases are resolved by written rules, not case-by-case review, and some species will be excluded who could have been matched by hand. Exclusions are counted by order.
- **The ecology control absorbs some lineage signal for bats.** Flight status (`terrestrial_volant`) nearly coincides with Chiroptera, so in the primary H2 model, controlling for ecology partly controls for being a bat. This lowers the partial r for reasons that are phylogenetic, not ecological, and makes rule (b) easier to meet and rule (a) harder. The primary H2 decision uses the matrix with `terrestrial_volant` as specified. The partial r without it is reported next to it, and the report reads the difference between the two as a measure of this overlap, not as evidence for or against H2.
- **Order-level structure dominates.** Most of the variance in patristic distance across mammals lies between orders, so the full-scale r may mostly reflect whether texts separate orders. The per-order results show what happens within orders.

## 11. Changes

### Change 1 (2026-10-30, before the species set was frozen and before any embedding or test): residual rule for tree matches

**What had been seen:** the reconciliation output (`matched.csv`, `residual_review.csv`, `excluded_species.csv`) and the Wikipedia texts fetched so far. No embedding, distance matrix or test existed.

**Why:** at full scale the engine's epithet-within-family rule makes wrong matches across genera within large families (for example *Calomys miurus* to the tip *Microtus miurus*, *Acomys selousi* to *Grammomys selousi*, *Sorex obscurus* to *Cryptotis obscura*). This replaces the fuzzy-acceptance sentence in section 3, step 2.

**Rule** (`scripts/03b_residual_rules.py`): a tree match is accepted without review only if the engine made it by exact name, MSW3/CMW legacy name, or MDD synonym. Every other automatic match (fuzzy, epithet within family, Wikipedia bridge) goes to `residual_review.csv` and is excluded under section 3, unless the pilot hand-reviewed the same match to the same tip. This exception applies to two species, *Canis lupaster* and *Otaria flavescens*, and keeps all 284 pilot species for the replication check. This moved 16 species: 12 epithet within family, 3 fuzzy and 1 Wikipedia bridge. A random sample of 50 accepted matches, stratified by method (seed 20261004), is in `data/interim/mammalia/match_spotcheck_sample.csv` for manual checking.

### Change 2 (2026-10-30, before the species set was frozen and before any embedding or test): masking rules at full scale

**What had been seen:** the terms that the pilot rules mask in the Wikipedia texts fetched so far (about 3,400 species), counted by term. No embedding, distance matrix or test existed.

**Why:** with the all-mammal vocabulary, the pilot rules mask common English words. Examples: "tail" (about 3,000 times), "head", "two", "top" and "grass", which come from rare alternative common names such as "Cameroon Scaly-tail", "Hard-head", "One-Two", "Steeple-top" and "Cutting-grass". They also mask epithets used as ordinary words ("major", "minor", "comes", "niger", "papua"). This replaces "by the same procedure" in section 4. The pilot's own masking is unchanged (`config.yaml` keeps the pilot rules as defaults), so the exact replication in section 6a is unaffected. The difference the full-scale rules make on the pilot texts is part of the vocabulary effect measured in section 6b.

**Rules** (`clades.mammalia.masking` in `config.yaml`, `src/lineage_charisma/masking.py`):

1. **Head nouns** (strict level): the whole final word of each MDD common name, in any script. A hyphenated final word is kept whole ("scaly-tail", not "tail"). Names containing a preposition (de, del, do, da, dos, das, of, des, du, la, le) give no head noun, because their final word is not the head noun ("Tuco-tuco of the Dunes", "Hocicudo de Azara"). Modifiers ("long-tailed", "greater") are never masked as such.
2. **Epithets**: masked only inside a binomial or an abbreviated binomial ("Panthera leo", "P. leo"), not wherever they occur. The pilot masked current epithets anywhere.
3. **Higher-rank names**: names from every MDD rank column (subclass, infraclass, magnorder, superorder, order, suborder, infraorder, parvorder, superfamily, family, subfamily, tribe), with derived forms (-idae to -id and -ine, -inae to -ine, -oidea to -oid). In the pilot, these names were hand-listed for Carnivora only.
4. **Vernacular names of orders and higher clades**: an explicit list, the full-scale counterpart of the pilot's "carnivoran" (marsupial, primate, rodent, cetacean, lagomorph and others; listed in `config.yaml`).
5. **`keep_terms`**: an explicit list for words that rules 1 to 4 still mask but are not taxon names in these texts. These are body-part and generic words from one-off alternative names (tail, head, sides, human); sex terms used for many species (bull, cow, buck, stag); a few others (ounce, jack, digger, runner, canary, pinyin); the historical genus *Fiber*; "boss" (Bos plus a plural s); and "fossorial", which contains a genus name. Plural forms are listed with each.

**Effect on the texts fetched so far** (3,410 species, strict level): rules 1 and 2 unmask 284 terms (2,899 occurrences); rules 3 and 4 newly mask 39 terms (485 occurrences); `keep_terms` unmasks 27 terms (5,782 occurrences). The guard that fails if any masked text still contains its own genus name is unchanged. These counts are updated on the full corpus when the species set is frozen.

### Frozen species set

To be added when `data/processed/mammalia/matrices/species_order.txt` is written: the species count and the orders with at least 50 species.
