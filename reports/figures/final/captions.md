# Final figures: draft captions

Every figure is drawn by `scripts/35_final_figures.py` from the saved tables in `reports/tables/`. Figure 2 also reads the saved pilot distance matrices. Each figure is a 300 dpi PNG plus a PDF with embedded fonts.

The style is the same in every figure: 7.2 in wide and 9 pt Helvetica Neue. Blue is the pilot (284 Carnivora), orange is the full scale (4,611 mammals), and grey is a reference or comparison value. A filled marker means one-sided permutation p < 0.05 and a hollow marker means p ≥ 0.05, so no reading depends on colour alone.

Unless a caption says otherwise:
- **Primary specification:** BAAI/bge-large-en-v1.5 embeddings, chunked (the mean of overlapping passages), on strictly masked text (scientific names and common-name head nouns replaced).
- **Text distance:** cosine.
- **Phylogenetic distance:** patristic distance on the Upham et al. (2019) MCC tree.
- **Statistics:** Spearman correlations with permutation tests (seed 20261004).

The numbers quoted here are checked against the tables by `scripts/34_verify_numbers.py`.

---

**Figure 1. Study pipeline.** Data sources:
- MDD taxonomy and synonyms;
- the Upham et al. 2019 MCC and posterior trees;
- English Wikipedia articles;
- EltonTraits 1.0 and COMBINE traits;
- Wikipedia pageviews, 2023-01 to 2025-12.

Names are reconciled to MDD and to the tree tips. Articles are reduced to the lead plus description sections, and stubs (under 50 tokens) are dropped, leaving 284 Carnivora in the pilot and 4,611 mammals at full scale. Texts are masked at three levels and embedded with two models and two passage rules.

Five species-by-species distance matrices are built:
- text: cosine;
- phylogeny: patristic, in millions of years;
- ecology: Gower distance over trait groups;
- length: absolute difference in log token count;
- pageviews: absolute difference in log(1 + median monthly views).

The tests are:
- **H1:** Mantel r between text and phylogeny;
- **H2:** partial Mantel r given ecology and length;
- **MRM:** on all four predictors.

Every test is Spearman with 9,999 permutations. The full-scale robustness variants use 999. *File: `fig01_pipeline`.*

**Figure 2. Masking effect (pilot, descriptive).** n = 284 Carnivora.
- **Left:** the share of species whose nearest neighbour in text space is in the same family. Unmasked it is 90%, with taxonomy masking 87%, and with strict masking 44%, against 11% expected by chance (the mean share of other species in the same family).
- **Right:** H1 Mantel r at each mask level: 0.204, 0.215 and 0.033. One-sided p ≤ 0.0001, ≤ 0.0001 and 0.0285; 9,999 permutations.

Bars show the primary model and rule. Grey dots show the other three model/rule combinations (bge-large truncated; all-mpnet-base-v2 chunked and truncated). The nearest-neighbour shares repeat the Week 3 sanity check, recomputed from the saved matrices (`reports/tables/final_masking_effect.csv`). *File: `fig02_masking_effect`.*

**Figure 3. Pre-registered decisions, pilot vs full scale.** All values are on the primary specification and the MCC tree, with one-sided permutation p.
- **H1:** pilot r = 0.033, p = 0.0285, n = 284 (supported). Full scale r = −0.028, p = 1.0000, n = 4,611 (not supported); 9,999 permutations each.
- **H2 baseline:** the raw r on species with ecology data is 0.035 (p = 0.0241, n = 283) in the pilot and −0.029 (p = 0.9999, n = 4,519) at full scale.
- **H2 partial r (ecology + length):**
  - Pilot: 0.020, p = 0.1272. Rule (a), partial r > 0 with p < 0.05, is not met. Rule (b), partial r below half the raw r, is not met either: the partial r is 57.5% of the raw r.
  - Full scale: −0.034, p = 1.0000. Rule (a) is not met, and rule (b) is undefined because the raw r is ≤ 0.

*File: `fig03_pilot_vs_full`.*

**Figure 4. H2: the text–phylogeny correlation before and after each control.** Spearman partial Mantel r between text and phylogenetic distance on the H2 species (pilot n = 283, full scale n = 4,519), with one-sided p from 9,999 permutations.

The rows are:
- the raw r;
- three single controls: ecology (S1), article length (S2) and pageviews (S3);
- the primary H2 model, ecology + length (shaded);
- ecology + length + pageviews (S4).

In the pilot, every control leaves the r between 0.020 and 0.037. At full scale every value is negative (−0.034 to −0.015). *File: `fig04_h2_controls`.*

**Figure 5. Within-order correlations at full scale (descriptive).** For each order with at least 50 species, the left panel shows the H1 Mantel r (circles) and the H2 partial r given ecology and length (squares), computed within the order. There are 9,999 permutations per order, and hollow markers have one-sided p ≥ 0.05. The dashed line is the all-mammal H1 r (−0.028).

The right panel shows the share of each order's matched Wikipedia articles dropped as stubs (under 50 tokens). That share ranges from under 1% in Carnivora to 31% in Eulipotyphla.

H1 r is positive with p ≤ 0.0011 in:
- Artiodactyla (0.283);
- Chiroptera (0.214);
- Lagomorpha (0.143);
- Primates (0.139).

It is near zero in Rodentia (−0.006) and Eulipotyphla (−0.008). These per-order breakdowns are descriptive and do not test H1. *File: `fig05_per_order`.*

**Figure 6. Robustness to tree uncertainty.** The primary H1 r and H2 partial r recomputed on trees drawn from the Upham et al. posterior, with the seeded draws fixed in the amendments.

| sample | H1 r | H2 partial r |
|---|---|---|
| pilot, 100 trees, 9,999 permutations per tree | median 0.080; p < 0.05 on 89/100 | median 0.062; p < 0.05 on 85/100 |
| full scale, 20 trees, 999 permutations per tree | median −0.015 (range −0.060 to 0.011); p < 0.05 on 0/20 | median −0.022; p < 0.05 on 0/20 |

Each dot is one tree; hollow dots have p ≥ 0.05. Bars mark medians, and diamonds mark the MCC tree, which decides the registered tests. *File: `fig06_posterior_trees`.*

**Figure 7. Pageview quartiles (exploratory).** Species pairs are split into quartiles by the mean of the two species' log pageviews. The pilot has 10,046–10,047 pairs per quartile and the full scale about 2.66 million.

**Left:** the mean text distance falls from the least- to the most-viewed quartile, in the pilot from 0.356 to 0.239 and at full scale from 0.439 to 0.357.

**Right:** the Spearman r between text and phylogenetic distance within the least- and most-viewed quartiles.

| sample | Q1 (least viewed) | Q4 (most viewed) | permutations |
|---|---|---|---|
| pilot | 0.133, p = 0.0001 | 0.035, p = 0.0455 | 9,999 |
| full scale | 0.043, p = 0.001 | 0.024, p = 0.006 | 999 |

The analysis is exploratory, as registered. *File: `fig07_pageview_quartiles`.*

**Figure 8. Length-matched analysis (post-hoc, amendment 4, section D).** The 4,611 mammals are split into 10 equal-size deciles of token count, ranging from 50–64 tokens to 555 tokens or more, and only pairs within the same decile are compared (1,060,761 pairs).

- **Left:** the Spearman r between text and phylogenetic distance within each decile. The dashed line is the pair-weighted mean r_LM = +0.0166, with two-sided p = 0.0245 from 9,999 permutations of species within deciles. The null distribution has mean 0.000 and SD 0.0074. The dotted line is the unmatched H1 r (−0.028).
- **Right:** the mean text distance of different-order pairs minus that of same-order pairs. It is positive in every decile (+0.002 to +0.022), so same-order pairs are described more alike at every length.

This analysis was designed after the full-scale results were seen and does not change H1. *File: `fig08_length_matched_posthoc`.*

**Figure 9. Register comparison, Wikipedia vs Animal Diversity Web (secondary analysis A).** The pilot species with usable ADW text: n = 246 for H1, and n = 245 for H2, which needs ecology data.

**Left:** H1 Mantel r and the primary H2 partial r for each source on the same species, with one-sided p from 9,999 permutations.

| test | Wikipedia | ADW |
|---|---|---|
| H1 | 0.019, p = 0.1753 | 0.075, p = 0.0002 |
| H2 partial r | 0.012, p = 0.2719 | 0.064, p = 0.0008 |

**Right:** the difference ADW − Wikipedia, with 95% intervals from a paired species bootstrap (2,000 draws).
- H1: Δ = +0.056 [+0.019, +0.095].
- H2: Δ = +0.052 [+0.013, +0.090].

The two sources differ in content (the ADW sections used) as well as register. *File: `fig09_adw_vs_wikipedia`.*

**Figure 10. Content charisma (secondary analysis B). [Placeholder: Analysis B has not run.]** This figure is drawn automatically once `reports/tables/secondary/charisma_agreement.csv` exists.

- **Left:** for each of the four rubric dimensions, Spearman agreement between the author and the model on the 50 validation texts, and between model runs 1 and 2 on all pilot species. Dimensions with either agreement below 0.5 are dropped (registered rule).
- **Right:**
  - B1: partial Mantel r of text distance on charisma-content distance, given length and pageviews.
  - B2: the H2 partial r with charisma-content distance added as a control, next to the primary H2 partial r from the same run.

Model: claude-sonnet-5-5. The extended MRM (B3) is listed in `final_results_index.md`. *File: `fig10_charisma` (not yet produced).*
