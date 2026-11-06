# Week 6 report: the full-scale analysis, all mammals

Run 2026-11-05 under pre-registration Amendment 3 (`docs/preregistration_amendment_fullscale.md`, committed before any full-scale species set existed, with Changes 1 to 3 and the frozen species set committed before any embedding or test). Every number below comes from a file in `reports/tables/mammalia/` or `data/processed/mammalia/`. Settings are the pilot's: bge-large-en-v1.5 / chunk / strict as the primary specification, Spearman on untransformed distances, 9,999 permutations (999 for robustness variants, posterior trees and exploratory tests), seed 20261004, one-sided p.

## 1. Headline

| test | pilot (Carnivora) | full scale (all mammals) | registered outcome at full scale |
|---|---|---|---|
| H1, primary | r = 0.033, p = 0.0285, 284 species | **r = −0.028, p = 1.000**, 4,611 species | **not supported** |
| H2 raw r (H2 species) | 0.035, 283 species | −0.029, 4,519 species | |
| H2 partial r (ecology + length) | 0.020, p = 0.127 | **−0.034, p = 1.000** | rule (a) **fails**; rule (b) **undefined** (raw r ≤ 0) |
| MRM phylogeny coefficient | 0.027, p = 0.13 | −0.023, p = 0.006 (two-sided) | descriptive |

Section 2 of the amendment applies: the pilot results stand as recorded and are not re-decided. The two analyses disagree in sign, so neither supports H1 as a general claim about mammal descriptions.

**The full-scale primary r is negative because strict masking removes almost everything that separates orders.** Under strict masking, pairs from the same order are on average slightly *less* alike in text (mean distance 0.412) than pairs from different orders (0.403). Under taxonomy masking, which leaves common-name head nouns such as "bat", "mouse" and "whale" in place, same-order pairs are clearly more alike (0.401 vs 0.442) and H1 r is +0.186. Splitting by text length mostly removes the strict-mask gap: 0.441 vs 0.444 for pairs of short texts (below the median token count) and 0.348 vs 0.350 for long ones. Same-order pairs are dominated by rodent–rodent and bat–bat pairs, and these orders have short, sparse texts that sit far apart; this composition, not a text–phylogeny relation, drives the sign. These checks are diagnostic, not pre-registered.

**Within orders the signal is positive where the texts are informative.** Descriptive per-order H1 r: Artiodactyla 0.283, Chiroptera 0.214, Lagomorpha 0.143, Primates 0.139 (all p ≤ 0.0011); Carnivora 0.029 (p = 0.048); Rodentia −0.006 and Eulipotyphla −0.008 (n.s.).

**The attention pattern from the pilot replicates across mammals and within most orders (exploratory).** Pairs of popular species are described more alike, and the text–phylogeny correlation is stronger among obscure species: across all mammals, r = 0.043 in the least-viewed quartile of pairs and 0.024 in the most-viewed; within Rodentia, 0.151 vs 0.010; within Primates, 0.251 vs 0.054.

## 2. Species set (Step 2 and the freeze)

| step | species |
|---|---:|
| MDD v2.5 | 6,904 |
| excluded: no Upham tip | 1,190 |
| excluded: extinct or domestic without a tip | 43 |
| residual (not resolved by a written rule) | 88 |
| matched to a tip and a Wikipedia article | 5,583 |
| stubs (< 50 tokens) | 972 |
| **frozen set** (`species_order.txt`, SHA-256 `4730f407…763f0e`) | **4,611** |
| H2 species (ecological distance to every other species) | 4,519 |

Residual pile (88): no validated Wikipedia article 60, ambiguous tree match 11, and 17 matches the engine had accepted that Change 1 sends to review (12 epithet within family, 3 fuzzy, 1 Wikipedia bridge, and *Homo sapiens*, whose article had been found through the synonym *Homo daliensis*, "Dali Man"). All 284 pilot species are in the frozen set with the same tip and article.

Orders with at least 50 species: Rodentia 1,722; Chiroptera 989; Primates 406; Artiodactyla 329; Eulipotyphla 325; Carnivora 284; Diprotodontia 139; Lagomorpha 85; Didelphimorphia 79; Dasyuromorphia 69. Stub rates are highest in Eulipotyphla (31%), Rodentia (23%) and Chiroptera (19%), against 0.4% in Carnivora (`reports/tables/mammalia/stubs_by_order.csv`, `reports/figures/week6/stubs_by_order.png`). Microbiotheria loses its only matched species to the stub rule.

Coverage of the frozen set (`coverage_by_order.csv`): 4,191 species have an EltonTraits record and 4,529 a COMBINE record; 4,025 have every trait and 64 none. Every species has pageviews in the 2023–2025 window, and 4,585 have views in all 36 months. The median of median monthly pageviews is 141 overall, 59 to 80 in shrews, bats and rodents, and 3,590 in Carnivora.

**Known edge cases in the matching (accepted as MDD synonymy).** From the 50-species stratified spot-check (`data/interim/mammalia/match_spotcheck_sample.csv`): *Cervus canadensis* is matched to the tip `Cervus_wallichi` (article "Elk"), because the Upham tree has no *canadensis* tip and MDD lists *wallichi* as a synonym; *Pseudomys mimulus* is matched to `Pseudomys_pilligaensis` (article "Pilliga mouse"), a name MDD treats as a synonym; and a few accepted titles are genus articles or redirects rather than species articles (for example "Setirostris" for the monotypic *Setirostris eleryi*, and "Crocidura kegoensis"). No wrong cross-genus matches like the earlier *Calomys miurus* appeared in the sample.

## 3. Performance (Step 3)

Details in `reports/week6_performance.md`. The statistics ran on the CPU of an Apple M3 Pro (8 threads), with no GPU needed: about 25 ms per permutation at 10.6 million pairs, 224 s for the 9,999-permutation H1 pass, 293 s for the H2 pass and 247 s for the MRM pass. Embeddings ran on the local Apple GPU (84 to 206 s per specification). Two implementation changes made this possible, both checked against the pilot code: condensed-vector permutation tests (`src/lineage_charisma/condensed.py`, `tests/test_condensed.py`), and a vectorized patristic distance that is bit-identical to dendropy's (19 s for 5,583 tips instead of hours; `tests/test_phylo.py`).

## 4. Replication (Step 4, before any full-scale result was read)

| check | species | pilot r | replication r | difference | p | outcome |
|---|---:|---:|---:|---:|---:|---|
| 6a exact, pilot vocabulary | 284 | 0.032706 | 0.032706 | 7 × 10⁻¹⁷ | 0.0285 (pilot 0.0285) | **passed** |
| 6b full-scale vocabulary | 284 | 0.032706 | 0.028857 | −0.0038 | 0.0480 | **passed** (within 0.01) |

In 6a the full-scale path re-masked all 284 pilot texts identically to the pilot corpus, re-embedded them to bit-identical vectors, and reproduced the pilot's phylogenetic distances exactly. In 6b the phylogenetic distances were again identical (maximum difference 0), and 94 of the 284 pilot texts were masked differently. Most changes newly mask other mammals that carnivore texts mention as prey or competitors: "rodents" (32 texts), "rabbits" (9), "primates" (7), "hog" (7), "agouti" (5), "deer" (4) and others, 77 terms in all; 19 terms were unmasked, mostly epithets that now appear outside a binomial (`replication_mask_terms.csv`). The gap is attributed to this vocabulary change, as the amendment specifies.

## 5. H1 and robustness

| specification | r | p | permutations |
|---|---:|---:|---:|
| **bge-large / chunk / strict (primary)** | **−0.0282** | **1.000** | 9,999 |
| all-mpnet-base-v2 / chunk / strict | 0.0006 | 0.444 | 999 |
| bge-large / truncate / strict | −0.0251 | 0.999 | 999 |
| bge-large / chunk / taxonomy | 0.1861 | 0.001 | 999 |
| Pearson, primary | −0.0221 | 0.984 | 999 |

10,628,355 pairs. Every strict-mask specification is at or below zero; only the taxonomy mask, which leaves the common-name head nouns in place, gives a positive correlation.

**Posterior trees** (20 trees, 999 permutations each; `posterior.csv`, `posterior_summary.csv`, `reports/figures/week6/posterior.png`): H1 r median −0.015, range −0.060 to 0.011; no tree gives p < 0.05. With 20 trees the percentiles are coarse. The MCC value (−0.028) lies inside the posterior range, so the negative result is not specific to the summary tree.

## 6. H2

| model | controls | partial r | p | share of raw |
|---|---|---:|---:|---:|
| raw (H2 species) | none | −0.0289 | 0.9999 | |
| **primary** | **ecology + length** | **−0.0337** | **1.000** | undefined |
| S1 | ecology | −0.0330 | 1.000 | undefined |
| S2 | length | −0.0299 | 0.9999 | undefined |
| S3 | pageviews | −0.0154 | 0.972 | undefined |
| S4 | ecology + length + pageviews | −0.0233 | 0.998 | undefined |
| robustness: unweighted Gower | ecology (columns equal) + length | −0.0377 | 1.000 | undefined |
| **robustness: ecology without flight** | ecology without `terrestrial_volant` + length | **−0.0311** | 1.000 | undefined |
| robustness: all-mpnet-base-v2 | ecology + length | −0.0120 (raw 0.0001) | 0.943 | not meaningful (raw r ≈ 0) |
| robustness: truncate | ecology + length | −0.0355 (raw −0.0258) | 1.000 | undefined |
| robustness: taxonomy mask | ecology + length | 0.1625 (raw 0.1865) | 0.001 | 0.87 |

**Rule (a) fails** (partial r ≤ 0). **Rule (b) is undefined**, because the raw r on the H2 species is not positive. The amendment states this in advance (section 5) rather than treating it as "confounded".

Flight status (Change 3: `terrestrial_volant` coincides exactly with Chiroptera) makes no material difference here: the partial r is −0.034 with it and −0.031 without it. With the raw r already negative, there is no lineage signal for the column to absorb.

Across the 20 posterior trees the partial r has median −0.022 (range −0.067 to 0.004). Rule (a) holds on no tree. On the 4 trees where the raw r is positive, the partial r is below half of it (rule (b) "holds" on 20% of trees). Every one of those raw r values is under 0.011, so the ratio describes noise around zero.

**MRM** (text ~ phylogeny + ecology + length + pageviews, 4,519 species, 9,999 permutations; `h2_mrm.csv`, `reports/figures/week6/h2_mrm.png`): pageviews −0.190 (p = 0.0001), length 0.085 (0.0001), ecology 0.058 (0.0001), phylogeny −0.023 (0.006); R² = 0.034. The large negative pageview coefficient alongside a positive length coefficient repeats the pilot's suppressor pattern: the two attention measures overlap, and the pageview term should not be read alone as "popular pairs differ less". Under the taxonomy mask the phylogeny coefficient is 0.165, about as large as ecology (0.183).

## 7. Per-order breakdowns (descriptive)

| order | n | H1 r | p | H2 n | raw r | partial r | p | share |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Rodentia | 1,722 | −0.006 | 0.689 | 1,685 | −0.003 | −0.020 | 0.950 | |
| Chiroptera | 989 | 0.214 | 0.0001 | 968 | 0.215 | 0.204 | 0.0001 | 0.95 |
| Primates | 406 | 0.139 | 0.0001 | 398 | 0.134 | 0.159 | 0.0001 | 1.18 |
| Artiodactyla | 329 | 0.283 | 0.0001 | 316 | 0.316 | 0.175 | 0.0001 | 0.55 |
| Eulipotyphla | 325 | −0.008 | 0.615 | 322 | −0.018 | −0.027 | 0.816 | |
| Carnivora | 284 | 0.029 | 0.048 | 283 | 0.031 | 0.017 | 0.175 | 0.53 |
| Diprotodontia | 139 | 0.063 | 0.102 | 139 | 0.063 | 0.027 | 0.289 | 0.42 |
| Lagomorpha | 85 | 0.143 | 0.0011 | 85 | 0.143 | 0.063 | 0.082 | 0.44 |
| Didelphimorphia | 79 | 0.021 | 0.351 | 77 | 0.020 | 0.023 | 0.338 | 1.15 |
| Dasyuromorphia | 69 | −0.033 | 0.735 | 65 | −0.028 | −0.034 | 0.737 | |

9,999 permutations, no correction for multiple tests, no claim per order (`per_order.csv`, `reports/figures/week6/per_order.png`). Carnivora here is the 6b value (full-scale masks; Gower ranges from all mammals). The two orders with the most stubs and the least-viewed articles, rodents and shrews, show no signal; bats, primates and artiodactyls show a clear one, and in bats and primates the ecology and length controls barely change it.

## 8. Exploratory analyses

**Pageview quartiles** (`exploratory_pageview_quartiles.csv`, `reports/figures/week6/pageview_quartiles.png`; 999 permutations, pairs held fixed). Across all mammals, mean text distance falls steadily from the least- to the most-viewed quartile of pairs (0.439, 0.422, 0.401, 0.357). The text–phylogeny r is 0.043 (p = 0.001) in the least-viewed quartile and 0.024 (p = 0.006) in the most-viewed. Within orders, the least-viewed quartile has the larger r in 9 of 10 orders (Lagomorpha is the exception, 0.075 vs 0.098). The gap is large in Rodentia (0.151 vs 0.010), Primates (0.251 vs 0.054), Eulipotyphla (0.173 vs −0.018) and Carnivora (0.123 vs 0.031). Both pilot exploratory findings, that popular pairs are described more alike and that the tree signal is stronger for obscure species, hold across mammals.

**Transformed phylogenetic distance** (`exploratory_nonlinear.csv`). Only 0.83% of pairs are within 30 My of each other (13% in the Carnivora pilot). With distance capped at 30 My, H1 r is 0.074 (Spearman, p = 0.001) and 0.086 (Pearson), and the primary partial r is 0.071 and 0.083. The log transform gives Pearson r = 0.018 (p = 0.019). The text tracks recent divergences, among close relatives, and loses or reverses the relation over deep divergences between orders; that is where the strict mask has removed the cues. With so many pairs tied at the cap, this is a description, not a test.

## 9. Interpretation and limits

- **The registered answer at full scale is no.** On the primary specification, text distance between mammal descriptions does not track phylogenetic distance, before or after the controls. This holds on every posterior tree.
- **What the primary specification measures.** Strict masking was chosen to stop the embedding from matching common-name words. At full scale, those words (bat, mouse, shrew, whale) are most of what separates orders in these short texts. The taxonomy-mask variant, which keeps them, gives r = 0.19 and a partial r of 0.16. The difference between the two specifications is a measure of how much of the cross-order signal is carried by those words. It is not evidence that ecology or attention explains the signal.
- **Attention organizes the texts.** Popular species pairs are described more alike in every order, and the tree signal is mostly found among obscure species. Pageviews have the largest MRM coefficient. Article length and stub status are strongly patterned by order (section 2), which is part of why the cross-order relation is negative.
- **Limitations stated in advance still apply** (amendment section 10): pairs are not independent, partial Mantel tests have inflated type I error under autocorrelation, the stub rule removes species non-randomly, and with 10 million pairs very small correlations reach significance. Effect sizes here are small throughout, except within bats, primates and artiodactyls, and under the taxonomy mask.

## 10. Files

- Tables: `reports/tables/mammalia/` — `h1_mantel.csv`, `h2_partial_mantel.csv`, `h2_mrm.csv`, `per_order.csv`, `posterior.csv`, `posterior_summary.csv`, `replication.csv`, `replication_mask_terms.csv`, `replication_notes.json`, `exploratory_pageview_quartiles.csv`, `exploratory_nonlinear.csv`, `species_by_order.csv`, `stubs_by_order.csv`, `coverage_by_order.csv`, `trait_missing_rate_by_order.csv`, `ecology_missing.csv`.
- Figures: `reports/figures/week6/` — `h2_controls.png`, `h2_mrm.png`, `posterior.png`, `per_order.png`, `pageview_quartiles.png`, `stubs_by_order.png`.
- Scripts: `19_freeze_species.py`, `20_coverage_by_order.py`, `21_replication.py`, `22_fullscale_tests.py`, `23_posterior_fullscale.py`, `24_exploratory_fullscale.py`, `25_week6_figures.py`; optional `scripts/slurm/embed_adroit.sh`.
- Performance: `reports/week6_performance.md`.
