## 7. Observations

These are hand-written notes on the run of 2026-10-09. `scripts/17_week5_report.py` appends this file, `reports/week5_notes.md`, to the report unchanged, so check the numbers against the tables above after any rebuild.

**Neither decision rule applies.** On the primary specification the raw r on 283 species is 0.035 (p = 0.024) and the partial r given ecology and length is 0.020 (p = 0.127). Rule (a) fails: there is no evidence, on the registered test, that text distance tracks phylogeny beyond the measured ecology and article length. Rule (b) also fails: the partial r is 57% of the raw r, above the 50% line. So the registered outcome is "not shown to survive, and not largely confounded". This is not evidence that phylogeny has no effect.

**The controls removed a little of a small effect.** The absolute change is 0.015. The raw p-value was already close to 0.05, so a small drop was enough to lose significance. Ecology does most of it (r falls to 0.024, 69% of raw); length does less (0.031, 88%).

**Controlling for pageviews raises r.** Given pageviews alone the partial r is 0.037, above the raw 0.035, and adding pageviews to ecology and length lifts the partial r from 0.020 to 0.027 (p = 0.062). In the MRM the pageview coefficient is negative (−0.082, p = 0.006) while length is positive (0.121). Length and pageviews are correlated at 0.76 across species, so this looks like a suppressor effect between two overlapping measures of attention and should not be read as "popular pairs are described less alike". The choice to leave pageviews out of the primary controls was made after Week 4 and before these numbers; with pageviews in, rule (a) still fails.

**MRM gives the same ranking as the Week 4 marginals.** Ecology (0.135) and length (0.121) are each about five times the phylogeny coefficient (0.027, p = 0.13). The four predictors together explain 3% of the variance in ranked text distance. The pseudo-t for phylogeny is 5.3, which would look decisive under an ordinary regression; the permutation p of 0.13 is the reminder that 39,903 pairs from 283 species are not independent.

**The primary specification is again the weakest, and the only one of the 12 where rule (a) fails.** The other three strict-mask variants keep 72% to 78% of their raw r with p ≤ 0.002. Without strict masking the controls barely matter (93% to 98% kept), because most of that signal is shared names, which ecology and length do not measure. Rule (b) holds in none of the 12.

**Gower weighting and the trait-complete subset change nothing.** Partial r is 0.021 in both (p = 0.12).

**On the posterior trees the picture is more favourable to lineage than on the MCC tree.** Rule (a) holds on 85 of 100 trees, the median partial r is 0.062, and the partial r is below the raw r on every tree by a similar amount (median 0.018, 2.5th to 97.5th percentile 0.005 to 0.028). The ratio is under 50% on 8 trees and undefined on 1; its median is 0.79. The ratio is as noisy as expected: where the raw r is near zero it runs from −3.9 to 0.9, which is why the absolute change is the better summary. As in Week 4, the MCC tree sits at the low end of the posterior for both raw and partial r. The registered decision uses the MCC tree, so it stands; but the honest reading is that the MCC-tree result is fragile in both directions and the size of the controls' effect (about 0.015 to 0.02 off r) is the stable finding.

**Exploratory: pairs of well-known species are described more alike.** Mean text distance falls steadily from 0.356 in the least-viewed quartile of pairs to 0.239 in the most-viewed, while mean phylogenetic distance is flat (61 to 63 My). The text–phylogeny correlation is 0.13 among the least-viewed pairs and 0.035 among the most-viewed. This fits the "charisma" story (popular species get long, similarly structured articles that the embedding puts close together, whatever their lineage), but it is not a test of it: popular species have longer articles, Week 3 showed long articles crowd together for this model and rule, and the quartiles contain different mixes of species. Exploratory.

**Exploratory: the weak primary result depends on treating phylogenetic distance as linear in rank.** Capping distance at 30 My triples the Spearman H1 r (0.115) and gives a partial r of 0.078 (p = 0.0001). Pearson on untransformed, log and capped distances gives partial r of 0.070, 0.081 and 0.090, all at the smallest possible p, each 67% to 71% of its raw r. Text distance appears to respond to recent divergence and not to how deep an old split is, and 87% of pairs are beyond 30 My where it is flat. These transforms were chosen after seeing the Week 4 scatter, so they are hypotheses for the all-mammal data, not results.

**Limits.** As stated in the amendment: partial Mantel p-values can be too small when matrices are autocorrelated (Guillot and Rousset 2013), which cuts against the significant results here (the variants and posterior trees) and not against the non-significant primary one; the control is linear in ranks; and ecology is 19 coarse trait columns, so phylogeny can stand in for unmeasured ecology.

## 8. Go/no-go for scaling to all mammals (Week 6)

**Recommendation: go, with the data build, embeddings and MCC-tree primary analyses on the laptop, and the posterior trees and the variant grid on Adroit.** Four pieces of code have to change first, and one bottleneck is not computational.

Estimates use the measurements in section 6 (M3 Pro, 12 cores, 36 GB) and assume all 5,911 tree tips, which is an upper bound: species without a usable Wikipedia article will drop out.

| step | estimate at ~5,900 species | where |
|---|---|---|
| Wikipedia text and pageviews | about 12 pageview requests per species at the rate seen in Week 4 (3 per second): roughly 7 hours, plus a similar order for article text. Not measured for text. | laptop (needs internet; cached afterwards) |
| embeddings, 2 models × 2 rules × 3 masks | about 1.5 hours if articles average 1.7 passages as in Carnivora | laptop |
| one distance matrix | 0.28 GB as float64, 0.14 GB as float32; about 19 matrices, so 5 GB or 2.7 GB on disk | laptop |
| one Mantel test, 9,999 permutations | 0.7 hours (0.25 s per permutation) | either |
| one partial Mantel test (2 controls) | about 1 hour | either |
| H1 and H2 on the primary specification (raw, primary, four secondary, MRM) | about 8 hours | laptop, overnight |
| 12 variants, H1 both methods and H2 with MRM | about 45 hours serial | Adroit, or laptop over a few nights |
| 100 posterior trees, raw and partial | about 170 hours serial | Adroit job array |

The permutation tests are independent single-core jobs of a few GB each, so on Adroit the whole grid is an array of roughly 300 one-hour jobs. On the laptop alone the full Week 5 design at 9,999 permutations is about nine days of serial compute; running several at once would shorten that, but the speed-up has not been measured.

**What has to change in the code:**

1. **Permutation batch size.** `mantel`, `partial_mantel` and `mrm` gather 128 permuted matrices at a time. At 5,911 species that is 36 GB in one allocation. The batch size must be set from the matrix size (2 to 4 at this scale). The benchmark used 2.
2. **Tree distances.** `cophenetic_matrix` fills the matrix with a Python loop over pairs, 17.5 million calls per tree at this scale. It has not been timed at full size and needs a vectorized replacement before any all-mammal tree is processed, let alone 100 posterior trees.
3. **Storage.** float32 matrices, memory-mapped on load, and one cache of ranked and standardized matrices so that each is ranked once (4 seconds each) and not once per test.
4. **Shared permutations.** Every test uses the same seed and so the same permutations. Tests against the same permuted matrix can share one pass (the code already computes several correlations per permuted matrix), which would cut the H2 run on one text matrix from about six passes to two.

**What could change in the design (needs an amendment before Week 6 results):** 999 permutations for robustness variants and posterior trees (about 4 minutes per test, so the posterior run drops to under a day on the laptop), with 9,999 kept for the primary tests; and fewer posterior trees.

**The bottleneck that is not compute: manual review.** The trait stage stops when more than 5 species are unmatched and Carnivora had 9 of 284. At that rate all mammals would have about 190, plus name-reconciliation residuals and Wikipedia stubs, which will be far more common among rodents and bats than among carnivores. Week 6 needs written rules for these cases in place of case-by-case review, and the final species count is unknown until the articles are fetched.

**Scientific case for going.** The Carnivora result is small, sensitive to the tree, and decided by a p-value near the threshold. More species is the direct remedy. Two things to settle in a Week 6 pre-registration before any all-mammal number is seen: whether the primary correlation stays Spearman on untransformed distance, given the exploratory result above, and whether the Carnivora analysis is the confirmatory one with all mammals as replication, or the reverse.
