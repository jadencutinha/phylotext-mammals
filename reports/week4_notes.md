## 6. Observations

These are hand-written notes on the run of 2026-10-04. `scripts/13_week4_report.py` appends this file, `reports/week4_notes.md`, to the report unchanged, so check the numbers against the tables above after any rebuild. Nothing here is an H2 conclusion: every correlation in this report is marginal.

**H1 holds, and the effect is very small.** The primary specification gives r = 0.033 with p = 0.0285. The permutation error on that p-value is about ±0.003, so it is below 0.05 whatever the seed, but it is not far below. An r of 0.033 means phylogenetic distance orders the text distances barely better than chance once names are masked.

**The primary specification is the weakest of the 12.** Every other variant has a larger r, and all 12 are positive and significant. bge-large / chunk is the lowest combination at every mask level, and strict is the lowest mask level for every model and rule. The choice of primary was made before these numbers existed (see the pre-registration, section 1), so this is where the result landed and not a selection.

**Most of the unmasked signal was names.** Taxonomy masking changes almost nothing (r moves by about 0.01 at most, in either direction). Strict masking leaves 16% to 29% of the unmasked r. This matches the Week 3 family check: common-name words such as "mongoose", "seal" and "fox" carried the family structure.

**Pearson is higher than Spearman in every variant** (0.098 against 0.033 for the primary). The scatter suggests why: median text distance rises over roughly the first 25 to 30 million years and is flat after that, and most of the 40,186 pairs sit in the flat part. This is a reading of the figure and has not been tested.

**The MCC tree gives a lower r than most posterior trees.** Across the 100 pre-registered posterior trees the median r is 0.080 (2.5th to 97.5th percentile 0.004 to 0.161), 98 are positive, and 89 have p < 0.05. Only 14 give a smaller r than the MCC tree's 0.033. So the direction of H1 does not depend on the MCC tree, but its size is uncertain by a factor of several, and the headline r is towards the low end of what the posterior supports. A likely reason for the spread is that a rank correlation depends on the order of the deep splits, where most pairs sit and where the trees disagree on dates (the largest distance ranges from 67 to 96 My across trees); this has not been tested.

**Inside families the picture is mixed.** Canidae shows a clear positive correlation (r = 0.200, p = 0.003). Mustelidae and Felidae show none (r = −0.005 and −0.078). These are descriptive, three tests with no correction.

**Things that look like the charisma effect.**

- Text distance is more strongly related to the difference in article length (r = 0.088) than to phylogeny (r = 0.033) in the primary specification. Species with articles of similar length have more similar embeddings. Week 3 saw the same thing for bge-large / chunk as a pile-up of long articles among each other's nearest neighbours.
- The difference in pageviews is not related to text distance (r = −0.023, p = 0.40). Length and pageviews are strongly correlated across species (0.76), so length is the part of attention that reaches the embeddings. This is the reason to keep the two attention matrices apart in Week 5 and not rely on the combined one (r = 0.049, p = 0.08).
- Phylogeny is itself related to the pageview difference (r = 0.105) and less to the length difference (r = 0.049): related species tend to get similar amounts of attention. That overlap is small, but it is the same size as or larger than the H1 effect.
- The `attention_mean_log_pageviews.npy` matrix for the "both well-known" check is built and has not been analysed.

**Ecology is the strongest marginal correlate of text distance** (r = 0.138 with groups weighted equally, 0.144 with columns weighted equally), about four times the phylogeny correlation. Ecology and phylogeny overlap only weakly (r = 0.081), and ecology is unrelated to any attention matrix (all |r| below 0.03, none significant). The two weightings of the ecology matrix agree in every row.

**What Week 5 has to deal with.** The H1 effect is smaller than the text–ecology and text–length correlations and about the size of the overlaps among the predictors. Whether anything of it is left after controlling for them is the H2 question and is not answered here.

**Data notes.**

- *Felis catus* has no COMBINE record and the eight species with no EltonTraits record have only COMBINE traits, so the domestic cat shares no trait with those eight. It is left out of the ecology matrices (283 species). Under amendment A1.2 the Week 5 raw text–phylogeny test is recomputed on those 283.
- The eight species missing from EltonTraits were split or described after the MSW3 taxonomy it follows (or are extinct, *Cryptoprocta spelea*). Their ecological distances rest on habitat breadth, land or water, and body mass where reported.
- Including redirects in the pageview counts matters for a few species. Redirects carry more than half the views for five articles; *Arctocephalus pusillus* has 90% of its views under earlier titles ("Brown fur seal", "Cape fur seal"). Some redirects are not names of the species ("Nandiniidae" redirects to African palm civet), and their views are counted too.
- The domestic cat's article ("Cat") has the most pageviews by a factor of more than two.
- The diagnostics in section 5 report two-sided p-values; the pre-registration did not say which side. H1 tests are one-sided as registered.
- The trend line in the scatter figure is the median text distance per 5 My bin, as in Week 3.
- Posterior tree tips are labelled `Genus_species`, without the `_FAMILY_ORDER` suffix the MCC tree uses. The names match one to one, so each tree was pruned with the Week 2 reconciliation after dropping the suffix.
