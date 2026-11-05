# Amendment 4 to the pre-registration: secondary analyses (Week 7)

Written 2026-11-05, before any Animal Diversity Web page was fetched, before any LLM annotation or hand label existed, and before any zoo signage was collected. The commit that adds this file is the timestamp of record. It covers the three secondary analyses promised in the proposal: a register comparison (A), content-based charisma (B) and zoo signage (C). Later changes go in the "Changes" section at the end, with a date and a reason, and always before the results they affect exist.

## 1. Status of these analyses

All three are **secondary**. Whatever they show, they do not change the H1 or H2 conclusions of the pilot (`reports/week4_report.md`, `reports/week5_report.md`) or of the full-scale analysis (`reports/week6_report.md`). They are reported in their own section, labelled secondary, and do not use the words "supports H1" or "supports H2". Analysis C is descriptive only.

Section D adds a fourth analysis of a different kind: a **disclosed post-hoc analysis**, designed after the Week 6 results were seen, to examine the full-scale negative r. It is not part of the registered H1/H2 outcomes and is reported apart from them, labelled post-hoc in every table, figure and sentence. Its specification is fixed here, before it is run.

## 2. What had been seen before this was written

- Every pilot result. On the primary specification (bge-large-en-v1.5 / chunk / strict, Spearman, MCC tree): H1 r = 0.033, p = 0.0285, 284 species. H2: raw r = 0.035 and partial r (ecology + length) = 0.020, p = 0.127, on 283 species. In the MRM, ecology and length were significant and phylogeny was not. Popular pairs are described more alike, and the tree signal is stronger among obscure species.
- Every full-scale result in `reports/week6_report.md`. On 4,611 species, H1 r = −0.028, H2 rule (a) fails and rule (b) is undefined. Taxonomy masking gives r = 0.19. Within-order r is positive in bats, primates and artiodactyls. The pageview pattern replicates.
- The unregistered diagnostics in section 1 of the Week 6 report, which motivate section D. Under strict masking, same-order pairs have mean text distance 0.412 and different-order pairs 0.403; under taxonomy masking, 0.401 and 0.442. Splitting species at the median token count, the strict-mask gap is 0.441 vs 0.444 among pairs of short texts and 0.348 vs 0.350 among pairs of long texts. No within-bin correlation, stratified test or other length-matched statistic has been computed.
- Nothing from Animal Diversity Web (ADW): no page, no coverage count, no text. No LLM output, no hand label, no signage.

## 3. Shared settings

Unless an item says otherwise, every test below uses the pilot primary settings: text distance = 1 − cosine similarity of `BAAI/bge-large-en-v1.5` embeddings under the chunk rule (512-token passages, 32-token overlap, mean-pooled); Spearman correlation on upper triangles with average ranks; 9,999 permutations; seed 20261004 with a fresh generator per test; one-sided p for correlations predicted to be positive; alpha 0.05. Phylogenetic, ecological and pageview distances are the pilot matrices (`data/processed/matrices/`), cut to the species analysed. The pilot code (`src/lineage_charisma/stats.py`) runs the tests.

## A. Register comparison: Wikipedia vs Animal Diversity Web

**Question.** On the same species, does text written in a scientific register (ADW species accounts) recover phylogeny better than Wikipedia?

**Retrieval** (before any analysis, subject to ADW's terms of use; if bulk retrieval is not clearly permitted, the analysis is not run and that is reported):

1. For each of the 284 pilot species (`data/processed/matrices/species_order.txt`), request `https://animaldiversity.org/accounts/<Genus>_<species>/` for the current MDD name. If there is no species account under that name, try in order: the Upham tip binomial, then each MDD synonym binomial in the order MDD lists them. The first page that is a species account (not a genus or higher-taxon page) and whose page title binomial equals the queried name is accepted. No other names and no manual matching are used.
2. Requests are cached on disk, rate-limited to at most one per second, sent with a descriptive User-Agent that includes contact details, and resumable.

**Text rule (fixed).** The ADW text is the concatenation, in this order, of the body text of these sections when present: *Geographic Range*, *Habitat*, *Physical Description*, *Behavior*, *Food Habits*. The section heading, figure captions, tables, reference lists, contributor lines and links are removed. These five sections are chosen to match what the Wikipedia text covers (the lead's range, habitat and diet, plus the description section). Reproduction, development, lifespan, communication, predation, ecosystem roles, economic importance and conservation status are not used. The cleaned text goes through the same normalization as the Wikipedia text (`wiki.clean_text`).

**Species set and minimum coverage.** These are the pilot species whose ADW text has at least 50 tokens (bge-large tokenizer, the pilot stub rule) and that are in the pilot species order. **At least 100 species are required.** Below 100, the coverage is reported by family and the analysis stops. For H2, the subset is further cut to the pilot H2 species (`ecology_species.txt`).

**Preprocessing.** Both texts are masked at the **strict** level with the **pilot** mask vocabulary (default `masking` settings in `config.yaml`, the same procedure that built `data/processed/corpus.parquet`). Both are then embedded with the primary model and chunk rule. Wikipedia texts are the existing pilot corpus texts, re-embedded only if needed; they are not re-fetched.

**Tests on the identical subset, for each source separately:**

- H1: Mantel r of text distance against phylogenetic distance, permuting phylogeny.
- Primary H2 model: partial Mantel r of text distance on phylogeny, controlling for ecology and length jointly. The length distance for each source is computed from that source's own text, by the same definition as `attention_length.npy`: the absolute difference in the log of the unmasked text's token count.

**Comparison.** The main statistic is the difference Δ_H1 = r_ADW − r_Wikipedia for H1. A second statistic is the same difference for the primary H2 partial r, Δ_H2. Each gets a **paired species bootstrap interval**:

- Draw species with replacement, the same draw for both sources (2,000 resamples; seed 20261004; `numpy.random.default_rng(seed)`).
- For each resample, compute both sources' statistics on the resampled species. Pairs made of two copies of the same species are dropped.
- Take the difference, and report the 2.5th and 97.5th percentiles.

The point estimate is the difference on the original subset. Both sources' r and permutation p-values are reported next to it.

**Reading.** If the 95% interval for Δ_H1 lies entirely above 0, the ADW text recovers phylogeny better on these species. If it lies entirely below 0, Wikipedia does. Otherwise there is no detectable difference. The same reading applies to Δ_H2. Nothing here is corrected for multiple tests. With about 100 to 280 species the interval is expected to be wide, and the report says so.

**Reported:** coverage (accounts found, by match route; texts at or above 50 tokens), coverage by family, and the median token count per source. Also r and p for each source and test, Δ with its interval, and the number of resamples in which a statistic was undefined (for example, a control fully explaining the text), which are dropped and counted.

## B. Content-based charisma (LLM annotation)

**Question.** Does the charisma *content* of a description explain text similarity beyond article length and pageviews?

**Unit and input.** Each pilot species' Wikipedia text at the **strict mask level** (`text_strict` in `data/processed/corpus.parquet`), the same text the primary embedding sees. The masked text keeps the rater from scoring the species' fame rather than what the text says. All 284 pilot species are annotated; H2 analyses use the 283 H2 species.

**Rubric** (fixed; the same wording goes to the model and into the hand-labelling instructions). Each description is scored 0 to 3 on four dimensions, using only what the text says:

| dimension | 0 | 1 | 2 | 3 |
|---|---|---|---|---|
| **warmth** (affective tone) | neutral, technical; no evaluative or emotional words | one or two mildly positive or emotive words | several affective or evaluative phrases; a noticeably warm or dramatic tone in places | affect-laden throughout; the description is written to make the reader feel something about the animal |
| **anthropomorphism** (mental-state language) | no mental states, intentions or personality attributed | one passing attribution (e.g. "curious", "shy") | several attributions of feelings, intentions or personality | the animal is described mainly through human-like mental life, motives, emotions or social roles |
| **aesthetic** (beauty, cuteness, impressiveness) | appearance described only in measurements and neutral colour terms | one aesthetic or impressive term ("striking", "large and powerful") | several; appearance is framed as beautiful, cute or impressive | aesthetic or impressive framing is central to the description |
| **human_relevance** (danger, usefulness, cultural role) | no mention of humans, human uses, danger or culture | one brief mention (e.g. hunted for fur, a single cultural reference) | several mentions, or one developed passage, on danger, economic use, conflict or cultural role | human relations (danger, use, symbolism, conflict) are a major theme of the description |

**Exact prompt** (system message, then user message; `{TEXT}` is the masked description):

System:
```
You are annotating encyclopedia descriptions of mammal species for a research study. Score only the language of the text you are given, not what you know about the animal. Some names in the text have been replaced by the token [TAXON]; do not try to guess them. Respond with JSON only.
```

User:
```
Score the description below from 0 to 3 on each of four dimensions, using these definitions.

warmth (affective tone): 0 = neutral, technical, no evaluative or emotional words; 1 = one or two mildly positive or emotive words; 2 = several affective or evaluative phrases, a noticeably warm or dramatic tone in places; 3 = affect-laden throughout, written to make the reader feel something about the animal.

anthropomorphism (mental-state language): 0 = no mental states, intentions or personality attributed; 1 = one passing attribution (e.g. "curious", "shy"); 2 = several attributions of feelings, intentions or personality; 3 = the animal is described mainly through human-like mental life, motives, emotions or social roles.

aesthetic (beauty, cuteness, impressiveness): 0 = appearance described only in measurements and neutral colour terms; 1 = one aesthetic or impressive term (e.g. "striking", "large and powerful"); 2 = several, with appearance framed as beautiful, cute or impressive; 3 = aesthetic or impressive framing is central to the description.

human_relevance (danger, usefulness, cultural role): 0 = no mention of humans, human uses, danger or culture; 1 = one brief mention (e.g. hunted for fur, a single cultural reference); 2 = several mentions, or one developed passage, on danger, economic use, conflict or cultural role; 3 = human relations (danger, use, symbolism, conflict) are a major theme of the description.

Return exactly this JSON object and nothing else, with an integer from 0 to 3 for each key:
{"warmth": <0-3>, "anthropomorphism": <0-3>, "aesthetic": <0-3>, "human_relevance": <0-3>}

Description:
"""
{TEXT}
"""
```

**Model and runs.** The provider and model are chosen by the author before any annotation. They are recorded, with the model version string, in a dated Change entry before the first validation request. Each species is scored in **two independent runs** at temperature 0, or the provider's most deterministic setting where 0 is not available, using structured (JSON) output where the provider offers it. Every request and response is cached. A response that is not valid JSON with exactly these four keys, each an integer from 0 to 3, is logged and the request is retried, up to 5 times per run. A species still without a valid response after 5 retries is reported and dropped from analysis B. The score used for analysis is the mean of the two runs.

**Validation (before the full analysis).**

1. **Sample.** 50 pilot species drawn with `numpy.random.default_rng(20261004).choice(284, size=50, replace=False)`, indexing the pilot species order. They are written to `data/interim/annotation_validation.csv` with blank score columns. The author labels these 50 masked texts by hand with the rubric. The model's scores are not shown before the labels are complete, and the labels are committed before the model is run on the sample.
2. **Agreement**, per dimension, as Spearman correlations: (a) author vs model (mean of the two runs) on the 50 species; (b) model run 1 vs model run 2 on all annotated species.
3. **Drop rule.** A dimension is **dropped** if either agreement (a) or (b) is **below 0.5**, or is undefined because a rater gave the same score to every species. Retained dimensions are used; dropped ones are reported with their agreement and not analysed further. If all four are dropped, analysis B stops and the agreement table is the result.

**Charisma-content distance.** Each retained dimension's mean score is standardized to mean 0 and SD 1 across the 284 species. The distance between two species is the Euclidean distance between their vectors of retained standardized scores; with one retained dimension, this is the absolute difference. The matrix is aligned to the pilot `species_order.txt` (zero diagonal, symmetric, checked by tests).

**Tests** (283 H2 species, Spearman, 9,999 permutations):

- **B1, the question:** partial Mantel r of text distance on charisma-content distance, controlling for length and pageviews jointly. The charisma matrix is permuted, and the test is one-sided (greater). Content charisma "explains text similarity beyond article length and pageviews" if the partial r > 0 with p < 0.05.
- **B2, extended H2 model:** partial Mantel r of text distance on phylogeny, controlling for ecology + length + charisma-content, with phylogeny permuted. Reported next to the primary H2 partial r (0.020), as the change in partial r when charisma content is added.
- **B3, extended MRM:** text ~ phylogeny + ecology + length + pageviews + charisma-content, ranked and standardized, with the text matrix permuted and two-sided p on pseudo-t (descriptive).

**Reported** for each dimension: the score distribution, both agreements, and whether it was retained. Also the marginal Mantel correlations of charisma-content distance with text, phylogeny, length and pageviews, plus B1–B3.

## C. Bronx Zoo signage: a third register (ingestion only this week)

**Data.** The author photographs and transcribes the species signs at the Bronx Zoo. Transcription format: `data/raw/zoo_signage/template.csv`, one row per sign, with columns:

| column | content |
|---|---|
| `species_common` | common name as printed |
| `species_scientific` | scientific name as printed (may be outdated or abbreviated) |
| `exhibit` | exhibit or building name |
| `photo_filename` | file name of the photo in `data/raw/zoo_signage/photos/` |
| `transcription` | the sign's body text, typed verbatim; line breaks become spaces; headings, the common and scientific name lines, and map or range labels are left out |
| `notes` | anything unusual (sign partly covered, several species on one sign, text in a graphic) |

A sign covering several species is entered once per species, with the same transcription. Only mammals are analysed.

**Processing.** Printed scientific names are reconciled to MDD with the existing engine, using the current name, then MDD synonyms, then a fuzzy match at the existing thresholds. A sign whose name does not resolve to one MDD species is listed and excluded. Signs are masked at the strict level with the **full-scale** vocabulary (`--clade mammalia`), because zoo species span many orders. They are embedded with the primary model and chunk rule. The Wikipedia texts used for comparison are the full-scale corpus texts (same vocabulary). The ADW texts are those from analysis A, re-masked with the full-scale vocabulary when a Carnivora species also has a sign. No stub rule is applied to signs, which are expected to be short; their token counts are reported.

**Analysis** (descriptive; no decision rule):

1. **Within-species cross-register similarity.** For each signed species, compute the cosine similarity between its sign and its Wikipedia text, and, where an ADW text exists, sign–ADW and Wikipedia–ADW. Compare these with the mean cosine similarity between that species' sign and the *other* signed species' Wikipedia texts, the between-species baseline. Report the median within-species similarity, the median baseline and the share of species whose within-species similarity exceeds their baseline.
2. **H1 on the overlap.** For the signed species that are in the full-scale frozen set, compute the Spearman Mantel r of text distance against phylogenetic distance separately for sign text and for Wikipedia text on those same species. Do the same for ADW if at least 10 signed species have ADW texts. Report the r values, with permutation p-values (999 permutations, same seed) given for completeness only.
3. If fewer than 10 signed species resolve, item 2 is not computed and only item 1 is reported.

Everything in C is reported descriptively because the sample is small and is selected by what the zoo exhibits.

## D. Length-matched analysis of the full-scale negative r (disclosed post-hoc)

**Status.** Post-hoc: designed on 2026-11-05 after the Week 6 results and the diagnostics in section 2 were seen. It does not change H1 or H2 at either scale. It is reported in its own section of the Week 7 report under the heading "post-hoc", and never as a registered outcome.

**Question.** Is the negative full-scale text–phylogeny correlation on the primary specification (r = −0.028) produced by how text length is distributed across orders? Rodents, bats and shrews have short, sparse texts and contribute most same-order pairs. If so, comparing only species with similar text lengths should remove the negative relation.

**Data.** The full-scale frozen set (4,611 species, `data/processed/mammalia/matrices/species_order.txt`), the primary text matrix (bge-large-en-v1.5 / chunk / strict, full-scale masks) and the MCC phylogenetic matrix, exactly as used for H1 in Week 6. The length of a species' text is `n_tokens` in `data/processed/mammalia/corpus.parquet`: the token count of the unmasked text (bge-large tokenizer, special tokens not counted), the quantity behind the stub rule and the length distance.

**Length bins.** The species are split into **10 bins of equal size** (deciles) by `n_tokens`. Ties are broken by position in `species_order.txt`, so the bins are deterministic:

```python
ranks = pd.Series(n_tokens, index=species_order).rank(method="first")
bin_of_species = pd.qcut(ranks, 10, labels=False)
```

With 4,611 species the bins hold 461 or 462 species each.

**Pairs compared.** Only **within-bin pairs**: both species in the same length decile. Pairs that cross bins are not used. In bin b with n_b species there are n_b(n_b − 1)/2 pairs, about 1.06 million in total across the 10 bins.

**Test statistic (D1, the main post-hoc statistic).** The pair-weighted mean of the within-bin Spearman correlations between text distance and phylogenetic distance:

  r_LM = Σ_b w_b · r_b,  with w_b = (pairs in bin b) / (pairs in all bins)

Here r_b is the Spearman correlation (average ranks) between the two distances over the pairs in bin b. Each bin is ranked on its own, so differences in text-distance level between short-text and long-text bins cannot enter the statistic.

**Null distribution.** The phylogenetic matrix's species labels are **permuted within bins**: each species can be swapped only with species in its own length decile, and rows and columns move together. This keeps each bin's set of phylogenetic distances and each species' text length fixed, and breaks the pairing between text and phylogeny within a bin. **9,999 permutations, seed 20261004**: one generator, `numpy.random.default_rng(20261004)`, and for each permutation the bins are permuted in bin order 0 to 9 with `rng.permutation`. The p-value is **two-sided**, because the post-hoc question is about sign: p = (1 + number of permutations with |r_LM*| ≥ |r_LM|) / 10,000.

**Also reported (descriptive, no further tests).**

- r_b, n_b, the token range and the pair count for each bin.
- For each bin, the mean text distance of same-order and of different-order pairs (MDD `order`), and their difference. Also the pair-weighted mean of that difference across bins, next to the unmatched values 0.412 and 0.403.
- r_LM with 5 and with 20 bins (same rules), as a check on the choice of 10. These are reported without p-values.
- The unmatched H1 r (−0.028) for comparison.

**Reading.** All readings are post-hoc and descriptive of the primary specification. None of them is a claim about H1.

- **r_LM > 0 with p < 0.05:** the full-scale negative r is consistent with a length-composition effect. Among species with similar text lengths, closer relatives are described more alike.
- **r_LM not distinguishable from 0:** once length is matched, the strict-mask texts carry no detectable phylogenetic signal across mammals, consistent with the mask removing the cues that separate orders.
- **r_LM < 0 with p < 0.05:** text length does not explain the negative relation.

**Limitations.** Equal-size decile bins leave some length variation within each bin, most in the top bin (the longest articles). Length is correlated with pageviews, so matching on length partly matches on attention as well. This analysis was chosen after seeing the data, and its p-value does not have the meaning of a pre-registered test.

## 4. Known limitations, stated in advance

- **ADW coverage and authorship.** ADW accounts are uneven in length and are often written by students. Accounts are more likely to exist for well-known species, so the ADW subset may be biased toward popular carnivores, the group in which the pilot found the weaker tree signal.
- **Different sections, different content.** The ADW text rule and the Wikipedia rule do not cover identical content. A difference between sources reflects both register and content, and cannot be attributed to register alone.
- **LLM ratings are a measurement model.** Agreement with one human labeller on 50 texts is a minimum check. A model can be consistent with itself and still be wrong in the same way each time. The 0.5 threshold is a floor, not evidence of validity.
- **Masked texts lose some cues.** Rating strict-masked texts removes common names, so some aesthetic or cultural content tied to a name is not visible to either rater. This is deliberate, matching what the embedding sees.
- **Zoo signage is a convenience sample** selected by the zoo's collection, is short, and cannot support inference.

## 5. Changes

None yet.
