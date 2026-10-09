# Amendment 2 to the pre-registration: H2 (Week 5)

Written 2026-10-09 18:30 UTC, before any partial Mantel test or regression on distance matrices was run on the project's data. The commit that adds this file is the timestamp of record. It amends `docs/preregistration.md`, whose sections 1 to 8 and Amendment 1 stand unless an item below says otherwise. Later changes to this plan go in the "Changes" section at the end of this file, with a date and a reason.

## 1. What had been seen before this was written

Everything in `reports/week4_report.md`:

- the H1 result on the MCC tree (r = 0.033, one-sided p = 0.0285, 284 species) and all 12 variants;
- H1 on the 100 posterior trees (median r = 0.080);
- the marginal Mantel correlations of text distance with ecology (0.138), length distance (0.088), pageview distance (−0.023) and combined attention (0.049);
- the marginal correlations among the predictors (phylogeny with pageviews 0.105, ecology 0.081, length 0.049; ecology with attention all below 0.03);
- the scatter of text distance against phylogenetic distance.

No partial correlation, regression on distance matrices, or other analysis that conditions one matrix on another has been computed on these data. `partial_mantel` and `mrm` exist in `src/lineage_charisma/stats.py` and have been run on synthetic matrices only.

**Two choices below were made knowing the Week 4 marginals:**

- The primary controls are ecology and length distance, without pageview distance. Section 7 of the pre-registration planned controls for "ecology and human attention" and did not say which attention matrix enters the test. Week 4 showed that length distance is related to text distance and pageview distance is not. Pageview distance is kept as a secondary control and in the MRM.
- The nonlinear check in section 8 was prompted by the Week 4 scatter.

## 2. H2

**H2: text distance remains positively associated with phylogenetic distance after ecological distance and article-length distance are controlled for.**

## 3. Species

The 283 species in `data/processed/matrices/ecology_species.txt`: all 284 less *Felis catus*, which shares no trait with eight other species and so has no ecological distance to them. Every H2 analysis uses these 283 unless it says otherwise.

**Baseline.** The raw text–phylogeny Mantel test is recomputed on the same 283 species (Spearman, 9,999 permutations, seed 20261004, one-sided). This raw r is the baseline for every comparison in H2. Where a robustness analysis uses a different species set or a different tree, its baseline is the raw r on that same set and tree. H1 is not re-decided by this baseline.

## 4. Matrices

| name | file | definition |
|---|---|---|
| text | primary specification, bge-large / chunk / strict | 1 − cosine similarity |
| phylogeny | `phylo.npy` | patristic distance on the MCC tree, My |
| ecology | `ecology.npy` | Gower distance, six trait groups weighted equally |
| length | `attention_length.npy` | absolute difference in log article length in tokens |
| pageviews | `attention_pageviews.npy` | absolute difference in log(1 + median monthly pageviews) |

The combined attention matrix is not used in H2.

## 5. Primary H2 model

Partial Mantel test of text distance on phylogenetic distance, controlling for ecology and length jointly.

| setting | value |
|---|---|
| response | text distance |
| focal predictor | phylogenetic distance |
| controls | ecology (group-weighted) and length, together |
| correlation | Spearman: every matrix's upper triangle is ranked (average ranks for ties), then the partial correlation of the ranks is taken |
| statistic | the correlation between the residuals of text and the residuals of phylogeny, each regressed linearly on the two controls |
| null distribution | 9,999 permutations; each reorders the rows and columns of the phylogenetic matrix together, and the partial correlation is recomputed with text and controls left in place |
| seed | 20261004 (`stats.seed`), a fresh generator per test |
| p-value | one-sided, (1 + number of permuted partial r ≥ observed) / (1 + 9,999) |
| alpha | 0.05 |

## 6. Decision rules

Both rules are applied to the primary model in section 5 and nothing else. They are separate; both, one, or neither can hold, and the report states which.

- **(a) "Lineage signal survives"**: the primary partial r is greater than 0 with one-sided p < 0.05.
- **(b) "Signal is largely confounded"**: the primary partial r is less than 50% of the raw r on the same 283 species. A partial r at or below zero satisfies (b). If the raw r on the 283 species is itself at or below zero, (b) is not defined and the report says so.

The raw r, the partial r, and the partial r as a share of the raw r are reported whatever the outcome. If (a) does not hold, the wording is "no evidence that text distance tracks phylogeny beyond the measured ecology and article length", not "phylogeny has no effect".

## 7. Secondary models (descriptive)

Same settings as the primary model, differing only in the controls. They show which control moves r and do not decide H2.

| model | controls |
|---|---|
| S1 | ecology |
| S2 | length |
| S3 | pageviews |
| S4 | ecology + length + pageviews |

Each is reported with its partial r, one-sided p, and partial r as a share of the raw r. No multiple-testing correction is applied.

### MRM

Multiple regression on distance matrices (Lichstein 2007):

text ~ phylogeny + ecology + length + pageviews

- All five upper triangles are ranked and then standardized, so the coefficients are standardized coefficients on rank-transformed distances.
- Null distribution: 9,999 permutations of the rows and columns of the text matrix together, seed 20261004.
- Reported: each standardized coefficient, its two-sided permutation p-value (on the pseudo-t statistic), and R² with its permutation p-value.
- MRM puts the four predictors on one scale and is the main descriptive answer to "lineage or charisma". It does not decide H2.

## 8. Robustness (reported in full, not decisive)

For each item the raw r on the same species and tree, the primary partial model and the MRM are run, except where noted. No correction for multiple testing, and none of these changes the outcome of rules (a) and (b).

1. **All 12 model × rule × mask variants** of the text matrix (the primary specification is one of them).
2. **Unweighted Gower** (`ecology_unweighted.npy`) in place of the group-weighted ecology matrix, per A1.1.
3. **The 272 complete-trait species**: those with no missing value in any of the 19 trait columns of `data/processed/traits.csv`. All matrices are cut to these species.
4. **The 100 posterior trees** selected in A1.3 (same indices). For each tree: the raw r and the primary partial r on the 283 species, 9,999 permutations. Reported as two distributions (median, 2.5th and 97.5th percentiles, minimum, maximum) and the share of trees meeting rule (a) and rule (b). MRM is not run per tree.

## 9. Exploratory analyses

Everything in this section is labelled exploratory in tables, figures and the report. None of it bears on H1 or H2.

### "Both well-known"

Uses `attention_mean_log_pageviews.npy` (the mean of the two species' log pageviews) on all 284 species and the primary text matrix.

- The 40,186 pairs are split into quartiles of mean log pageviews. Reported per quartile: number of pairs, mean and median text distance, mean phylogenetic distance.
- Within the top quartile and within the bottom quartile: the Spearman correlation between text distance and phylogenetic distance over the pairs in that quartile. A permutation p-value is given by reordering the phylogenetic matrix (9,999 permutations, seed 20261004) with the set of pairs held fixed.
- Questions: are pairs of popular species described more alike, and does the text–phylogeny correlation differ between well-known and obscure pairs?

### Nonlinear phylogeny check

The Week 4 scatter suggested that text distance rises over roughly the first 25 to 30 My of patristic distance and is flat after that. Two transformed phylogenetic matrices are built:

- **log**: log(1 + patristic distance);
- **capped**: min(patristic distance, 30 My).

For each, H1 is rerun on the 284 species and the primary H2 model on the 283 species.

**A log transform cannot change a Spearman correlation.** It keeps the order of the distances, so the ranks, and therefore every Spearman result, are identical to the untransformed ones. The log variant is therefore run with Pearson correlation and compared with the Pearson result on untransformed distances. The capped variant changes the ranks (all pairs beyond 30 My tie) and is run with both Spearman and Pearson. About 13% of pairs are at or below 30 My.

### Unchanged from before

The 277-species `matrices_min100/` set stays exploratory and is not planned for Week 5.

## 10. Known limitations, stated in advance

- **Partial Mantel tests can have inflated type I error.** When the matrices are autocorrelated (here, through the tree), partial Mantel p-values can be too small (Guillot and Rousset 2013, *Methods in Ecology and Evolution* 4: 336–344). A significant partial r is weaker evidence than its p-value suggests. This is why MRM is reported alongside it, and why effect sizes are reported next to every p-value. MRM uses the same kind of permutation and is not immune to the problem.
- **The control is linear in ranks.** If text distance depends on ecology or length in a way that is not monotone, some of that dependence stays in the residuals. In synthetic tests, a partial correlation fell below zero when two matrices shared a driver but were not linearly related.
- **The controls are measured with error.** Ecology is 19 coarse trait columns and attention is two variables. Phylogeny can stand in for unmeasured ecology, so a surviving partial r means "not accounted for by the measured controls", not "caused by lineage".
- **The raw effect is small.** The H1 p-value was 0.0285 on 284 species; on 283 the raw result may fall either side of 0.05. A ratio of two small correlations, as in rule (b), is noisy; the posterior-tree distribution is reported to show how noisy.

## 11. Changes

### Change 1 (2026-10-09, before any H2 results)

Made after review of the first version of this file and before any partial Mantel test or MRM was run on the project's data. Sections 1 to 10 are left as first committed.

**Why.** The raw r is small, so the ratio in rule (b) is noisy. Rule (b) stays as written in section 6 and is still decided on the MCC tree alone. Two further quantities are reported next to it. Neither is decisive.

1. **Absolute change in r for the primary model**: raw r minus partial r, on the 283 species and the MCC tree.
2. **Distribution of the partial-to-raw ratio across the 100 posterior trees**: for each tree, the primary partial r divided by the raw r on the same tree and species. Reported as median, 2.5th and 97.5th percentiles, minimum and maximum, and the share of trees where the ratio is under 50%. The ratio is not defined for a tree whose raw r is at or below zero; such trees are counted and reported separately, the summary statistics are taken over the remaining trees, and the share under 50% is given as a share of all 100 trees with the undefined ones listed as their own category.
