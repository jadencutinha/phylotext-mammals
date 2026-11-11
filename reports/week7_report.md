# Week 7 report: secondary analyses, post-hoc check, and locking results

Week 7 of 8. All analysis code is in place, and every result that can be produced without new data has been produced. Two items are still open: the charisma annotation, which waits for the author's 50 hand labels, and zoo signage, which waits for the zoo visit. Each has a pipeline that runs with one command. The index of every result is in `reports/final_results_index.md`, and the figure plan for Week 8 is in `reports/week8_figure_plan.md`.

## 1. Week 6 completed

Week 6 is complete; see `reports/week6_report.md`, `reports/week6_performance.md` and the frozen species list (4,611 species, SHA-256 `4730f407…763f0e`). In brief:

- **Replication.** The exact replication with the pilot vocabulary reproduced the pilot H1 (r = 0.032706, p = 0.0285) bit for bit. The full-vocabulary replication gave r = 0.0289, within the 0.01 threshold.
- **H1 at full scale is not supported** (r = −0.028, p = 1.0, 4,611 species).
- **H2:** rule (a) fails; rule (b) is undefined because the raw r is not positive.
- **Taxonomy mask:** r = 0.19.
- **Within orders:** r is positive in bats (0.21), artiodactyls (0.28), primates (0.14) and lagomorphs (0.14), and near zero in rodents and shrews.
- **Attention (exploratory, replicated):** popular species pairs are described more alike, and the tree signal is stronger among obscure species.

## 2. Pre-registration Amendment 4

`docs/preregistration_amendment_secondary.md` was committed before any secondary result existed. It covers:

- three **secondary** analyses, which cannot change the H1/H2 conclusions: A, register comparison; B, LLM-annotated content charisma; C, zoo signage, descriptive only;
- section D, one **disclosed post-hoc** analysis, designed after the Week 6 results and kept separate from the registered outcomes.

**Change 1** was committed before any ADW request, LLM request or hand label. It sets:

- **ADW retrieval:** an 8-second crawl delay (from ADW's `robots.txt`), the project's contact User-Agent, and full caching. Text is used under CC BY-NC-SA 3.0, and full texts are kept out of the repository.
- **LLM settings:** Anthropic `claude-sonnet-5-5` with default sampling (the model rejects temperature 0), thinking off, schema-constrained JSON output, and no fallback model.
- **Signage:** implementation details fixed before any data exist.

## 3. Post-hoc: the length-matched analysis (section D)

Species were split into 10 equal-size deciles of text length, and only pairs within the same decile were compared. The statistic is the pair-weighted mean of the within-decile Spearman correlations. The null distribution permutes species within deciles (9,999 permutations, seed 20261004, two-sided). Table: `reports/tables/mammalia/posthoc_length_matched*.csv`.

- **r_LM = +0.0166, two-sided p = 0.0245.** The null distribution has mean 0.000 and SD 0.0074. With 5 bins r_LM is +0.0168, and with 20 bins +0.0167.
- In **every** decile, same-order pairs are more alike than different-order pairs (gap +0.002 to +0.022; pair-weighted mean +0.008). This reverses the pooled pattern (0.412 vs 0.403).
- **Reading, as fixed in advance:** the full-scale negative r is consistent with a length-composition effect. Among species with texts of similar length, closer relatives are described slightly more alike.

The effect is small, about half the pilot r. It is post-hoc and does not change H1.

## 4. Secondary A: Wikipedia vs Animal Diversity Web

**Retrieval and coverage.** ADW's terms were checked first: the text is CC BY-NC-SA 3.0, `robots.txt` sets `Crawl-delay: 8`, and no page addresses automated access. The author chose to fetch at the crawl delay.

- The fetch made 505 requests over about 3 hours, roughly 22 s per request including server time.
- 279 of 284 pilot species have an account: 238 under the current name, 26 under the tree name, 15 under a synonym.
- **246 have usable text** (50 tokens or more; median 754 tokens), against the registered minimum of 100.
- Coverage is lowest in Herpestidae (19 of 34 usable) and Viverridae (20 of 33), whose ADW accounts are often almost empty (`reports/tables/secondary/adw_coverage*.csv`).

**Results on the identical subset** (`register_comparison.csv`, `register_bootstrap.csv`, `reports/figures/week7/register_comparison.png`):

| test | Wikipedia | ADW | Δ (ADW − Wikipedia), 95% paired species bootstrap |
|---|---|---|---|
| H1 r (246 species) | 0.019, p = 0.175 | 0.075, p = 0.0002 | **+0.056 [+0.019, +0.095]**: ADW higher |
| primary H2 partial r (245 species) | 0.012, p = 0.272 | 0.064, p = 0.0008 | **+0.052 [+0.013, +0.090]**: ADW higher |

On the same Carnivora species, ADW text recovers phylogeny better than Wikipedia, before and after ecology and each source's own length are controlled for.

As stated in advance, the two sources differ in content as well as register: the ADW rule takes range, habitat, physical description, behaviour and diet sections. So the difference cannot be attributed to register alone. Wikipedia's r on this subset (0.019) is below the pilot value (0.033) because the subset drops the mongooses and civets that lack usable ADW text. The comparison is fair because both sources are measured on exactly the same species.

## 5. Secondary B: content charisma (in progress)

Everything up to the model call is built and tested.

- **Rubric and prompt:** fixed in the amendment. `src/lineage_charisma/charisma.py` holds the same text, and a test checks that it matches the amendment exactly.
- **Validation sheet:** 50 seeded species with masked texts only. The item-to-species key is in a separate file. The files are `data/interim/annotation_validation.csv` and `annotation_instructions.md`.
- **Runner** (`scripts/31_annotate.py`):
  - It refuses to call the model until the hand labels are complete and committed.
  - It reads the key from `ANTHROPIC_API_KEY`, caches every response, validates every response and retries malformed ones up to 5 times.
  - It stops if run 1 projects above US$10.
- **Analysis** (`scripts/32_charisma_analysis.py`): agreement, the 0.5 drop rule, the charisma-content distance matrix, and B1–B3.
- **Estimated cost:** about US$1.60 for both runs; US$2.60 under pessimistic assumptions.

**Next:** the author labels the 50 texts and commits the sheet. Then run 1, run 2 and the analysis follow (about 20 minutes), and the index is regenerated.

## 6. Secondary C: zoo signage (ready for data)

`data/raw/zoo_signage/template.csv` has the six registered columns. `scripts/28_zoo_signage.py` reconciles printed names to MDD, masks with the full-scale vocabulary, embeds, compares each sign with Wikipedia (and with ADW where present), and computes H1 on the overlap once at least 10 signed species fall in the frozen set.

On 3 fake rows it ran end to end: all three names resolved, including the printed trinomial *Gorilla gorilla gorilla*; similarity was computed; and H1 was correctly skipped. Sign transcriptions stay in `data/raw/` and never reach report tables.

**Next:** transcribe the signs into `data/raw/zoo_signage/signage.csv` after the visit and run the script.

## 7. Locked results and Week 8 preparation

- **`reports/final_results_index.md`:** 103 rows, generated by `scripts/33_results_index.py` from the result tables. It covers pilot and full-scale H1/H2, posterior trees, per-order breakdowns, robustness, exploratory analyses, the replication checks, the post-hoc analysis and secondary analyses A–C. Each row gives the role (4 primary, 24 secondary, 49 robustness, 23 exploratory, 2 checks, 1 post-hoc), spec, statistic, p, decision outcome and source file. The charisma and signage rows are marked pending and are filled in by rerunning the script.
- **`reports/week8_figure_plan.md`:** 11 main figures and the supporting set, with paths.
  - Every existing figure is at 200 dpi and needs a 300 dpi export.
  - Fonts are already consistent.
  - In the simulated palette check, blue/orange/grey is safe for all three colour-vision deficiencies. Aqua is not (ΔE 11 against grey under deuteranopia), so the variant plots, whose primary-spec colour is aqua, and the Week 3 mask plots need a palette change or marker shapes.
  - Two figures need to be made: the study-design diagram and the length-matched figure.

## 8. Tests

`tests/test_condensed.py`, `test_fullscale_alignment.py`, `test_adw.py`, `test_charisma.py` and `test_zoo.py` cover:

- the condensed permutation code against the pilot code;
- chunk- and thread-independence;
- alignment of every full-scale matrix to the frozen list;
- ADW section extraction on saved sample pages (structure with short attributed excerpts);
- LLM response validation (malformed JSON, missing or extra keys, out-of-range, non-integer and boolean scores) and the registered prompt text;
- the charisma distance (symmetry, zero diagonal, standardized Euclidean, alignment to `species_order.txt`);
- signage ingestion and name resolution on the fake rows.

All tests pass.

## 9. What remains for Week 8

1. Charisma: hand labels, model runs, analysis, and a rerun of the results index.
2. Zoo signage: collect, transcribe, run.
3. Figures at publication quality (the plan above), the final write-up, and the presentation.
