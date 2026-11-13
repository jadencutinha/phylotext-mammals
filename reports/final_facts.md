# Final facts sheet (supporting material for the report)

Verified facts, numbers and decisions, grouped by paper section, each with its source. There is no prose here: the report text is the author's.

- Every number in the results section is read from a table in `reports/tables/`, and the headline numbers are re-derived by `scripts/34_verify_numbers.py` (`reports/number_check.md`).
- "pilot" means 284 Carnivora and "full scale" means 4,611 mammals.
- p-values are one-sided permutation p unless marked; MRM p-values are two-sided.

## Pre-registration history

| document | committed | commit | written before | what it fixed |
|---|---|---|---|---|
| Pre-registration (`docs/preregistration.md`, sections 1–8) | 2026-10-04 21:09 UTC | `3f66a97` | any text–phylogeny correlation | Primary spec (bge-large-en-v1.5 / chunk / strict). H1 Mantel test (Spearman, 9,999 permutations, seed 20261004, one-sided, α = 0.05) and decision rule. 11 robustness variants plus Pearson. Within-family checks. Ecology traits and Gower weighting; attention matrices. Primary mask changed from taxonomy to strict after the descriptive Week 3 report and before any H1 statistic (section 1). |
| Amendment 1 (`docs/preregistration.md` §9) | 2026-10-04 | `b64da42` | any results | A1.1 unweighted Gower variant. A1.2 raw r recomputed on the H2 species. A1.3 100 posterior trees as robustness, with the selection rule fixed. |
| Amendment 2 (`docs/preregistration_amendment_h2.md`) | 2026-10-09 18:30 UTC | `cb0211f` | any partial Mantel or MRM | The H2 statement. Primary partial Mantel (ecology + length) and decision rules (a) and (b). Models S1–S4 and the MRM. Robustness: 12 variants, unweighted Gower, 272 complete-trait species, 100 trees. Exploratory analyses: well-known pairs; log and 30 My cap. Limitations. |
| Amendment 2, Change 1 | 2026-10-09 | `329d1b7` | any H2 results | Absolute change (raw − partial r) and the partial/raw ratio over posterior trees are reported next to rule (b). Neither is decisive. |
| Amendment 3 (`docs/preregistration_amendment_fullscale.md`) | 2026-10-29 | `92484f9` | any full-scale species set, matrix or test | Full-scale species pipeline and freezing. Same H1/H2 settings. Replication checks 6a (exact) and 6b (threshold 0.01). Per-order breakdowns (≥ 50 species). Robustness (3 variants, unweighted Gower, ecology without flight, Pearson, 20 trees at 999 permutations). Exploratory pageview quartiles and 30 My cap. Limitations. |
| Amendment 3, Changes 1–2 | 2026-10-30 | `0384027` | the freeze and any embedding or test | Change 1: residual rule (only exact, legacy-name and synonym matches accepted automatically); synonym-only Wikipedia titles go to review; trait matches follow the same rule. Change 2: full-scale masking rules. |
| Amendment 3, Change 3 and frozen set | 2026-11-05 | `3c171ca` | any embedding or test | `terrestrial_volant` coincides exactly with Chiroptera. Frozen set of 4,611 species (SHA-256 `4730f407…763f0e`) and the 10 orders with ≥ 50 species. |
| Amendment 4 (`docs/preregistration_amendment_secondary.md`) | 2026-11-05 | `936ef8b` | any ADW page, annotation, hand label or signage | Secondary A (ADW register comparison, paired bootstrap), B (LLM charisma: rubric, exact prompt, validation, 0.5 drop rule, B1–B3) and C (zoo signage, descriptive). Disclosed post-hoc D (length-matched, two-sided, within-decile permutations). |
| Amendment 4, Change 1 | 2026-11-05 | `0b886b1` | any ADW request, LLM request or hand label | ADW: 8 s crawl delay, contact User-Agent, caching, CC BY-NC-SA 3.0 handling. LLM: claude-sonnet-5-5, default sampling (temperature 0 is rejected), thinking off, JSON-schema output, max_tokens 256, no fallback model, US$10 stop. Zoo: implementation details. |

## Introduction

- **Question:** does text-embedding similarity between Wikipedia descriptions of mammal species recover phylogeny (lineage), or does it mostly reflect human attention (charisma)? Source: `README.md`, `docs/preregistration.md` §4, `docs/preregistration_amendment_h2.md` §2.
- **H1:** text distance is positively associated with phylogenetic distance (Mantel r > 0). Source: `docs/preregistration.md` §4.
- **H2:** text distance remains positively associated with phylogenetic distance after ecological distance and article-length distance are controlled for. Source: `docs/preregistration_amendment_h2.md` §2.
- **Design:** a pilot on Carnivora (284 species), then a separately pre-registered full-scale test on all mammals (4,611 species). The pilot is not re-decided by the full-scale result. Source: `docs/preregistration_amendment_fullscale.md` §2.

## Related work (references the repository uses)

The repository holds only the references below; a literature review is not part of it.

- **Taxonomy:** Mammal Diversity Database (MDD) v2.5. Source: `config.yaml` `sources.mdd`, `docs/preregistration_amendment_fullscale.md` §1.
- **Phylogeny:** Upham, Esselstyn and Jetz (2019) mammal tree. MCC tree `MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre`; posterior `Completed_5911sp_topoCons_NDexp.zip`, 10,000 trees. Source: `config.yaml` `sources.tree`.
- **Traits:** EltonTraits 1.0, Wilman et al. 2014, *Ecology* 95:2027, doi:10.6084/m9.figshare.3559887; COMBINE, Soria et al. 2021, *Ecology* 102:e03344, doi:10.6084/m9.figshare.13028255 (reported, not imputed, values). Source: `docs/ecology_traits.md`, `config.yaml`.
- **Partial Mantel type I inflation:** Guillot and Rousset 2013, *Methods in Ecology and Evolution* 4:336–344. Source: `docs/preregistration_amendment_h2.md` §10.
- **MRM:** Lichstein 2007. Source: `docs/preregistration_amendment_h2.md` §7.
- **Animal Diversity Web:** text under CC BY-NC-SA 3.0. Source: `docs/preregistration_amendment_secondary.md` Change 1, `README.md`.
- **Embedding models:** `BAAI/bge-large-en-v1.5` and `sentence-transformers/all-mpnet-base-v2`. Source: `config.yaml` `embedding.models`.
- **Colour-vision simulation for the figures:** Machado et al. 2009. Source: `reports/week8_figure_plan.md`.

## Data

**Pilot species** (Carnivora):
- 320 MDD species; 35 excluded (32 with no tree tip, 2 domestic and 1 extinct without a tip).
- 285 matched to a tip and a Wikipedia article. Match methods: 248 exact, 27 legacy name, 8 synonym, 1 Wikipedia bridge, 1 fuzzy.
- 1 stub removed, leaving 284.
- 283 H2 species: *Felis catus* shares no trait with 8 species.
- 272 have complete traits.

Source: `data/interim/excluded_species.csv`, `data/interim/matched.csv`, `docs/preregistration_amendment_h2.md` §3 and §8.

**Full-scale species:**
- 6,904 MDD v2.5 species.
- Excluded: 1,190 with no Upham tip, 43 extinct or domestic without a tip, and 88 residual cases not resolved by a written rule.
- 5,583 matched to a tip and an article; 972 stubs removed.
- **4,611 in the frozen set**, of which **4,519 are H2 species**.

Source: `reports/week6_report.md` §2, `reports/tables/mammalia/species_by_order.csv`.

**Stubs and coverage:**
- The stub rule is fewer than 50 tokens (bge-large tokenizer, unmasked text, special tokens not counted).
- Stub rate by order: Eulipotyphla 31%, Rodentia 23%, Chiroptera 19%, Carnivora 0.4%.
- Microbiotheria loses its only matched species.

Source: `reports/tables/mammalia/stubs_by_order.csv`.

**Text:**
- Each species' text is the cleaned article lead plus the first section whose heading is a description heading (description, physical description, characteristics, appearance, anatomy, morphology, and similar).
- In the pilot, the median is 13 masked terms per text under strict masking.

Source: `config.yaml` `wikipedia.description_headings`, `docs/preregistration_amendment_fullscale.md` §3, `reports/week3_notes.md`.

**Pageviews:**
- Wikimedia Pageviews API, 2023-01 to 2025-12, all access methods, agent = user, redirects included. Every full-scale species has pageviews.
- Median of the median monthly views: 141 overall, 59–80 in shrews, bats and rodents, and 3,590 in Carnivora.

Source: `docs/preregistration.md` §7, `reports/week6_report.md` §2.

**Trait coverage (full scale):** 4,191 species have an EltonTraits record, 4,529 a COMBINE record, 4,025 every trait, and 64 none. Source: `reports/tables/mammalia/coverage_by_order.csv`, `reports/week6_report.md` §2.

**Orders with ≥ 50 species:**

| order | species |
|---|---:|
| Rodentia | 1,722 |
| Chiroptera | 989 |
| Primates | 406 |
| Artiodactyla | 329 |
| Eulipotyphla | 325 |
| Carnivora | 284 |
| Diprotodontia | 139 |
| Lagomorpha | 85 |
| Didelphimorphia | 79 |
| Dasyuromorphia | 69 |

Source: `docs/preregistration_amendment_fullscale.md`, frozen species set.

**ADW (secondary A):** 279 of 284 pilot species have an account (238 under the current name, 26 under the tree name, 15 under a synonym). 246 have usable text (≥ 50 tokens; median 754 tokens); 505 requests. Source: `reports/week7_report.md` §4, `reports/tables/secondary/adw_coverage*.csv`.

## Methods

The full parameter table is in `reports/methods_parameters.md`.

**Distances:**
- **text:** 1 − cosine similarity.
- **phylogeny:** patristic distance in My on the MCC tree.
- **ecology:** Gower distance with the six trait groups weighted equally; pairwise-available rule; no imputation.
- **length:** |Δ log tokens|.
- **pageviews:** |Δ log(1 + median monthly views)|.

Source: `docs/preregistration_amendment_h2.md` §4.

**Tests:**
- Spearman on upper triangles with average ranks.
- Permutation of the rows and columns of the phylogenetic matrix (the text matrix for MRM).
- 9,999 permutations, seed 20261004, a fresh generator per test.
- p = (1 + count ≥ observed) / (1 + permutations).

Source: `docs/preregistration.md` §4, `docs/preregistration_amendment_h2.md` §5.

**Decision rules:**
- **H1** is supported if r > 0 with p < 0.05 on the primary spec and the MCC tree.
- **H2 rule (a):** partial r > 0 with p < 0.05.
- **H2 rule (b):** partial r below 50% of the raw r on the H2 species; undefined if the raw r ≤ 0.

Source: `docs/preregistration.md` §4, `docs/preregistration_amendment_h2.md` §6.

**Full-scale implementation:**
- Condensed float32 vectors, ranks computed once, chunked permutations. Tested to give the same statistic and null as the pilot code.
- Vectorized patristic distance, bit-identical to dendropy.
- Robustness variants and posterior trees use 999 permutations.

Source: `docs/preregistration_amendment_fullscale.md` §5 and §8, `tests/test_condensed.py`, `tests/test_phylo.py`.

**Compute:** Apple M3 Pro CPU, 8 threads, about 25 ms per permutation at 10.6 million pairs. H1 pass 224 s, H2 pass 293 s, MRM 247 s. Embeddings 84–206 s per specification on the Apple GPU. Source: `reports/week6_performance.md`.

**Post-hoc D:**
- 10 equal-size deciles by `n_tokens`, with ties broken by species order; within-decile pairs only.
- r_LM is the pair-weighted mean of the within-decile Spearman r.
- Null: species permuted within deciles; 9,999 permutations, seed 20261004, two-sided.

Source: `docs/preregistration_amendment_secondary.md` §D.

**Secondary A:**
- Strict masking with the pilot vocabulary; bge-large / chunk.
- H1 and the primary H2 model per source.
- Δ (ADW − Wikipedia) with a paired species bootstrap (2,000 draws, seed 20261004).

Source: `docs/preregistration_amendment_secondary.md` §A.

## Results

### Pilot: H1

- **Primary H1:** r = 0.033, p = 0.0285 on 284 species (40,186 pairs). **H1 supported.** Source: `reports/tables/h1_mantel.csv`.
- **The primary spec is the weakest of the 12 variants.** All 12 have r > 0 and p < 0.05.
  - Strict: 0.033–0.109.
  - Taxonomy: 0.215–0.378.
  - None: 0.204–0.373.

  Source: `reports/tables/h1_mantel.csv`, `reports/week4_notes.md`.
- **Pearson primary:** r = 0.098. Source: `reports/tables/h1_mantel.csv`.
- **Masking:** the nearest text neighbour is in the same family for 90% of species unmasked, 87% with taxonomy masking and 44% with strict masking; chance is 11%. Source: `reports/tables/final_masking_effect.csv`, `reports/week3_notes.md`.
- **Posterior trees (100):** median r = 0.080 (2.5th–97.5th percentile 0.004–0.161); 98 trees have r > 0 and 89 have p < 0.05. The MCC tree is near the low end. Source: `reports/tables/h1_posterior.csv`.
- **Within families (descriptive):** Canidae r = 0.200 (p = 0.003), Mustelidae −0.005, Felidae −0.078. Source: `reports/tables/h1_mantel_families.csv`.

### Pilot: H2

- **Raw r on the 283 H2 species:** 0.035 (p = 0.0241).
- **Partial r (ecology + length):** 0.020 (p = 0.1272), 57.5% of the raw r; the absolute change is 0.015.
- **Rule (a) is not met and rule (b) is not met.**

Source: `reports/tables/h2_partial_mantel.csv`.

**Secondary models:**

| model | partial r | p |
|---|---:|---:|
| S1, ecology | 0.024 | 0.0845 |
| S2, length | 0.031 | 0.0413 |
| S3, pageviews | 0.037 | 0.0185 |
| S4, ecology + length + pageviews | 0.027 | 0.0623 |

Source: `reports/tables/h2_partial_mantel.csv`.

**MRM:** R² = 0.0317.

| predictor | β | p (two-sided) |
|---|---:|---:|
| phylogeny | 0.027 | 0.13 |
| ecology | 0.135 | 0.0004 |
| length | 0.121 | 0.0002 |
| pageviews | −0.082 | 0.006 |

Source: `reports/tables/h2_mrm.csv`.

**Robustness and posterior trees:**
- Unweighted Gower gives partial r = 0.021 (p = 0.117); the complete-trait subset gives 0.021 (p = 0.122). Source: `reports/tables/h2_partial_mantel.csv`.
- On the 100 posterior trees, rule (a) holds on 85, the median partial r is 0.062, and the median raw − partial change is 0.018. Source: `reports/tables/h2_posterior.csv`, `reports/week5_notes.md`.

**Confound diagnostics (marginal correlations with text distance):**
- ecology 0.138;
- length 0.088;
- pageviews −0.023 (p = 0.40);
- phylogeny with pageviews 0.105.

Source: `reports/tables/confound_mantel.csv`.

### Full scale

**Replication:**
- 6a reproduces the pilot exactly: r = 0.032706, p = 0.0285, difference 7 × 10⁻¹⁷.
- 6b gives r = 0.028857 (p = 0.048), a difference of −0.0038, within the 0.01 threshold. 94 pilot texts were masked differently.

Source: `reports/tables/mammalia/replication.csv`.

**H1:** r = −0.028, p = 1.000 on 4,611 species (10.6 million pairs). **Not supported.** Source: `reports/tables/mammalia/h1_mantel.csv`.

**H1 variants (999 permutations):**

| variant | r | p |
|---|---:|---:|
| mpnet | 0.0006 | 0.444 |
| truncate | −0.025 | 0.999 |
| taxonomy mask | 0.186 | 0.001 |
| Pearson | −0.022 | 0.984 |

Source: `reports/tables/mammalia/h1_mantel.csv`.

**H2:**
- Raw r on the 4,519 H2 species: −0.029.
- Partial r: −0.034 (p = 1.000).
- **Rule (a) is not met; rule (b) is undefined** because the raw r ≤ 0.
- Every secondary model is negative (−0.033 to −0.015).
- Ecology without flight gives −0.031; the taxonomy-mask partial r is 0.162.

Source: `reports/tables/mammalia/h2_partial_mantel.csv`.

**MRM:** R² = 0.034.

| predictor | β | p (two-sided) |
|---|---:|---:|
| phylogeny | −0.023 | 0.006 |
| ecology | 0.058 | 0.0001 |
| length | 0.085 | 0.0001 |
| pageviews | −0.190 | 0.0001 |

Source: `reports/tables/mammalia/h2_mrm.csv`.

**Posterior trees (20):**
- H1: median r = −0.015 (range −0.060 to 0.011); 0/20 trees have p < 0.05.
- H2: median partial r = −0.022; rule (a) holds on none.

Source: `reports/tables/mammalia/posterior_summary.csv`.

**Per order (descriptive):**

| order | H1 r | p |
|---|---:|---:|
| Artiodactyla | 0.283 | ≤ 0.0011 |
| Chiroptera | 0.214 | ≤ 0.0011 |
| Lagomorpha | 0.143 | ≤ 0.0011 |
| Primates | 0.139 | ≤ 0.0011 |
| Carnivora | 0.029 | 0.048 |
| Rodentia | −0.006 | 0.689 |
| Eulipotyphla | −0.008 | 0.615 |

Source: `reports/tables/mammalia/per_order.csv`.

**Diagnostic, not registered:**
- Under strict masking, same-order pairs average a text distance of 0.412 against 0.403 for different-order pairs. Under taxonomy masking the figures are 0.401 against 0.441.
- The Week 6 report and Amendment 4 print the second taxonomy value as 0.442; this is a known erratum.

Source: recomputed by `scripts/34_verify_numbers.py`.

### Exploratory (both scales)

- **Mean text distance falls with popularity.** From the least- to the most-viewed quartile of pairs it goes from 0.356 to 0.239 in the pilot and from 0.439 to 0.357 at full scale. Source: `exploratory_well_known.csv`, `mammalia/exploratory_pageview_quartiles.csv`.
- **The text–phylogeny r is stronger among obscure pairs:**
  - Q1 vs Q4: pilot 0.133 vs 0.035; full scale 0.043 vs 0.024.
  - The least-viewed quartile has the larger r in 9 of 10 orders.

  Source: same tables, `reports/week6_report.md` §8.
- **30 My cap:** H1 r = 0.115 in the pilot and 0.074 at full scale. Only 0.83% of mammal pairs are within 30 My, against 13% in Carnivora. Source: `exploratory_nonlinear.csv` (both scales), `reports/week6_report.md` §8.

### Post-hoc D

- **r_LM = +0.0166, two-sided p = 0.0245.** Null mean 0.000, SD 0.0074; 1,060,761 within-decile pairs.
- With 5 bins r_LM is +0.0168; with 20 bins, +0.0167.
- Same-order pairs are closer in every decile (gap +0.002 to +0.022; weighted mean +0.008).

Source: `reports/tables/mammalia/posthoc_length_matched*.csv`.

### Secondary A

| test | Wikipedia | ADW | Δ (ADW − Wikipedia), 95% bootstrap interval | reading |
|---|---|---|---|---|
| H1 (246 species) | r = 0.019, p = 0.1753 | r = 0.075, p = 0.0002 | +0.056 [+0.019, +0.095] | ADW higher |
| H2 partial (245 species) | r = 0.012, p = 0.2719 | r = 0.064, p = 0.0008 | +0.052 [+0.013, +0.090] | ADW higher |

Source: `reports/tables/secondary/register_comparison.csv`, `register_bootstrap.csv`.

### Secondary B and C

- **B:** not run as of 2026-11-11. The hand labels are not yet committed and no model request has been sent; the results index rows are placeholders. Source: `reports/final_results_index.md`.
- **C:** built and tested on 3 fake rows; not run because no signage was collected. Moved to future work. Source: `reports/final_results_index.md`, `tests/test_zoo.py`.

## Discussion (facts the discussion can draw on)

- **The pilot and the full scale disagree in sign on H1.** Under Amendment 3 §2 neither overrides the other. Source: `reports/week6_report.md` §1.
- **Strict masking removes most of what separates orders:**
  - Full scale, taxonomy mask: r = 0.186 against −0.028 with strict masking.
  - Pilot, same comparison: 0.215 against 0.033.

  Source: the two `h1_mantel.csv` tables.
- **Within informative orders the signal is positive** (bats, primates, artiodactyls, lagomorphs). It is near zero in the orders with the most stubs and the least-viewed articles (rodents, shrews). Source: `per_order.csv`, `stubs_by_order.csv`.
- **Attention organizes the texts:**
  - Pageviews have the largest absolute MRM coefficient at full scale (−0.190).
  - Length beats phylogeny in both MRMs.
  - Popular pairs are described more alike at both scales.

  Source: the MRM and quartile tables.
- **Length composition:** after length matching, the full-scale relation is small and positive (post-hoc). Source: `posthoc_length_matched.csv`.
- **Register:** ADW text recovers the tree better than Wikipedia on the same Carnivora, but register and content are confounded. Source: `register_bootstrap.csv`, `docs/preregistration_amendment_secondary.md` §4.

## Limitations

See `reports/limitations_inventory.md` for the deduplicated list with sources.

## Future work

- Analysis B, if the hand labels are not completed: everything up to the model call is built and tested (`scripts/31_annotate.py`, `32_charisma_analysis.py`, `tests/test_charisma.py`). Source: `reports/week7_report.md` §5.
- Analysis C, zoo signage: the pipeline is built and tested on 3 fake rows (`scripts/28_zoo_signage.py`, `tests/test_zoo.py`); it needs a transcribed `data/raw/zoo_signage/signage.csv`. Source: `reports/final_results_index.md`.
- The full 12-variant grid and 100 posterior trees at full scale: not run (Amendment 3 limits these to 3 variants and 20 trees). Source: `docs/preregistration_amendment_fullscale.md` §8.
- A nonlinear phylogeny model: the cap result is exploratory at both scales and was chosen after seeing the data. Source: `docs/preregistration_amendment_h2.md` §9, `reports/week5_notes.md`.
