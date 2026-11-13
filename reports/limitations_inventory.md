# Limitations inventory

This file lists every limitation stated in the pre-registration, its amendments and the weekly reports, merged where two sources say the same thing. Each entry gives its source or sources. Entries marked *(week 8)* record the status at the results lock on 2026-11-11.

## Statistical inference

1. **Partial Mantel tests can have inflated type I error.** When the matrices are autocorrelated through the tree, partial Mantel p-values can be too small (Guillot and Rousset 2013). A significant partial r is therefore weaker evidence than its p-value suggests. MRM uses the same kind of permutation and is not immune. This cuts against the significant results (the pilot variants and posterior trees), not against the non-significant primary H2 result.
   - Sources: Amendment 2 §10; Amendment 3 §10; `reports/week5_notes.md` ("Limits"); `reports/week6_report.md` §9.
2. **Pairs are not independent.** The Mantel permutation accounts for this, but the number of pairs (40,186 in the pilot, about 10.6 million at full scale) overstates the information, which is closer to the number of species. Very small correlations reach significance, so effect sizes are interpreted on their own terms. In the pilot MRM, a pseudo-t of 5.3 for phylogeny had a permutation p of 0.13.
   - Sources: Amendment 3 §10; `reports/week5_notes.md` (MRM); `reports/week6_report.md` §9.
3. **The raw effect is small, and ratios of small correlations are noisy.** The pilot H1 p of 0.0285 has a permutation error of about ±0.003. Rule (b)'s partial-to-raw ratio ran from −3.9 to 0.9 across posterior trees. The absolute change in r (about 0.015 to 0.02) is the stable summary.
   - Sources: Amendment 2 §10 and Change 1; `reports/week4_notes.md`; `reports/week5_notes.md`.
4. **No correction for multiple testing.** Robustness variants, secondary models, per-order tests, within-family tests and secondary analyses are reported without correction, and none decides H1 or H2.
   - Sources: prereg §5–6; Amendment 2 §7–8; Amendment 3 §7–8; Amendment 4 §A.

## The tree

5. **The MCC tree sits at the low end of the posterior** for the pilot raw and partial r: only 14 of 100 trees give a smaller H1 r. The size of the pilot effect is therefore uncertain by a factor of several. The registered decision uses the MCC tree alone. A likely reason, not tested, is disagreement among trees on deep split dates.
   - Sources: `reports/week4_notes.md`; `reports/week5_notes.md`.
6. **Only 20 posterior trees at full scale.** The percentiles are coarse (interpolated between the two most extreme trees), and 999 permutations per tree set a floor of p = 0.001.
   - Sources: Amendment 3 §8; `reports/tables/mammalia/posterior_summary.csv.meta.json`; `reports/week6_report.md` §5.
7. **Order-level structure dominates the full-scale distances.** Most variance in patristic distance lies between orders, so the full-scale r largely reflects whether texts separate orders. Only 0.83% of mammal pairs are within 30 My of each other, against 13% in Carnivora.
   - Sources: Amendment 3 §10; `reports/week6_report.md` §8.

## Controls and confounds

8. **The control is linear in ranks.** Non-monotone dependence of text distance on ecology or length stays in the residuals. In synthetic tests, a partial correlation fell below zero when two matrices shared a driver but were not linearly related.
   - Source: Amendment 2 §10.
9. **The controls are measured with error.** Ecology is 19 (pilot) or 20 (full scale) coarse trait columns, and attention is two variables. Phylogeny can stand in for unmeasured ecology, so a surviving partial r means "not accounted for by the measured controls", not "caused by lineage".
   - Sources: Amendment 2 §10; `reports/week5_notes.md`.
10. **The ecology control absorbs lineage for bats.** `terrestrial_volant` coincides exactly with Chiroptera. Controlling for ecology therefore partly controls for being a bat, which lowers the full-scale partial r for phylogenetic reasons. The variant without flight is reported (−0.031 against −0.034).
    - Sources: Amendment 3 §4, §10 and Change 3; `reports/tables/mammalia/h2_partial_mantel.csv`.
11. **Length and pageviews overlap.** They correlate at 0.76 across pilot species. In the pilot MRM the pageview coefficient is negative while length is positive, which looks like suppression between two measures of attention, not "popular pairs are described less alike". Matching on length (post-hoc D) also partly matches on attention.
    - Sources: `reports/week4_notes.md`; `reports/week5_notes.md`; Amendment 4 §D.
12. **Pageview counts include redirects.** Some redirects are not names of the species (for example "Nandiniidae" redirects to African palm civet). For a few species, redirects carry most of the views (90% for *Arctocephalus pusillus*).
    - Source: `reports/week4_notes.md`.
13. **Missing traits.** Nothing is imputed: the pairwise-available Gower rule is used.
    - Pilot: *Felis catus* shares no trait with 8 species and is dropped from H2 (283 species). The 8 species missing from EltonTraits rely on COMBINE traits alone.
    - Full scale: 64 species have no traits and 92 are dropped from H2 (4,519 species).
    - Sources: `reports/week4_notes.md`; Amendment 2 §3; `reports/week6_report.md` §2.

## Text, masking and embeddings

14. **The primary specification is the weakest of the variants and measures something narrow.** Strict masking removes the common-name words that carry most family and order structure, together with some real content tied to names. In the pilot, the primary spec has the lowest r of all 12 variants. At full scale, the taxonomy-mask variant gives r = 0.186 against −0.028. The primary was chosen before any H1 statistic existed.
    - Sources: prereg §1; `reports/week4_notes.md`; `reports/week5_notes.md`; `reports/week6_report.md` §9.
15. **Masking has side effects.**
    - Every masked text shares the `[TAXON]` token, which lowers all distances (mean 0.305 against 0.350 for bge-large / chunk). Texts with similar mask density are slightly closer.
    - Names of taxa outside the clade are left in the pilot (prey, for example), while the full-scale vocabulary masks them; this explains the 6b replication gap.
    - Common-name modifiers stay in ("sea [TAXON]").
    - Sources: `reports/week3_notes.md`; `reports/week6_report.md` §4.
16. **Description length drives proximity.** Long articles crowd together as nearest neighbours under bge-large / chunk. Length is tied to family and order, so it can look like phylogenetic signal or hide it. At full scale, length composition by order produces the negative r; post-hoc D addresses this.
    - Sources: `reports/week3_notes.md`; `reports/week4_notes.md`; `reports/week6_report.md` §1.
17. **The two embedding models agree only moderately.** At the taxonomy level, bge-large and mpnet distances have a rank correlation of 0.49–0.62.
    - Source: `reports/week3_notes.md`.
18. **The text is the lead plus one description section only.** Other parts of the article are not used.
    - Sources: prereg §2; Amendment 3 §3.

## Species coverage and matching

19. **Stubs are not random.** The 50-token rule removes 31% of Eulipotyphla, 23% of Rodentia and 19% of Chiroptera articles, against 0.4% of Carnivora. The frozen set therefore over-represents well-documented species, non-randomly with respect to both phylogeny and attention.
    - Sources: Amendment 3 §10; `reports/tables/mammalia/stubs_by_order.csv`; `reports/week6_report.md` §9.
20. **Name matching at scale uses written rules, not case-by-case review.** 88 residual species were excluded, and some could have been matched by hand. Known edge cases accepted as MDD synonymy:
    - *Cervus canadensis* is matched to the tip `Cervus_wallichi`.
    - *Pseudomys mimulus* is matched to `Pseudomys_pilligaensis`.
    - A few accepted titles are genus articles or redirects.
    - Sources: Amendment 3 §10 and Change 1; `reports/week6_report.md` §2.
21. **The pilot is one order.** Carnivora has unusually long, well-viewed articles (median of median monthly views 3,590, against 141 across mammals), so the pilot does not represent mammals in general. The pilot and full scale disagree in sign on H1, and neither overrides the other.
    - Sources: Amendment 3 §2; `reports/week6_report.md` §1–2.

## Post-hoc analysis D

22. **D was chosen after seeing the data.** Its p-value does not have the meaning of a pre-registered test, and it does not change H1. Equal-size deciles leave some length variation within each bin, most in the top bin.
    - Source: Amendment 4 §D.

## Secondary analyses

23. **Register comparison (A): register and content are confounded.** The ADW sections and the Wikipedia sections cover different content, so a difference cannot be attributed to register alone.
    - ADW accounts are uneven, often written by students, and more likely to exist for well-known species.
    - Coverage is lowest in Herpestidae (19 of 34 usable) and Viverridae (20 of 33).
    - Wikipedia's r on the ADW subset (0.019) is below the pilot r (0.033) because the subset drops those families.
    - The interval reflects 245–246 species.
    - Sources: Amendment 4 §4; `reports/week7_report.md` §4.
24. **LLM ratings (B) are a measurement model.** Agreement with one human labeller on 50 texts is a minimum check, and the 0.5 threshold is a floor, not evidence of validity. A model can be consistently wrong.
    - Masked texts hide name-linked content from both raters.
    - The model rejects temperature 0, so default sampling is used, and run-to-run agreement measures that variability.
    - Sources: Amendment 4 §4 and Change 1.
25. *(week 8)* **Analysis B is not finished at the lock.** The 50 hand labels are not committed, and no model request has been sent. Its index rows and figure are placeholders, and it is reported as incomplete if it cannot run this week.
    - Source: `reports/final_results_index.md`.
26. **Zoo signage (C) is a convenience sample.** It is selected by the zoo's collection, the signs are short, and it cannot support inference. *(week 8)* No signage was collected. The pipeline is built and tested on 3 fake rows only, and C is moved to future work.
    - Sources: Amendment 4 §4; `reports/final_results_index.md`.

## Reporting

27. *(week 8)* **One known numeric erratum.** The Week 6 report §1 and Amendment 4 §2 give the taxonomy-mask different-order mean text distance as 0.442. Recomputed from the matrices it is 0.44147, which rounds to 0.441. The committed text is left unchanged, and the correction is recorded in `scripts/34_verify_numbers.py` and `reports/number_check.md`. It does not affect any test.
28. **Some numbers in the weekly notes were computed by hand and are not in result tables.** The Week 3 mask-density check is one example, and so are the Week 6 diagnostics in §1. The Week 6 diagnostics are now recomputed by `scripts/34_verify_numbers.py`; the Week 3 hand checks are not.
    - Sources: `reports/week3_notes.md`; `reports/number_check.md`.
