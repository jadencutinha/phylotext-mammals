# Candidate figures for the final report and poster (Week 8)

Checked 2026-11-05 against four publication criteria: consistent fonts, a colourblind-safe palette, labelled axes, and 300 dpi.

**Applies to every figure:**

- **Resolution:** every PNG in `reports/figures/` is saved at 200 dpi (`plotstyle.save`). Changing that one setting to 300 and rerunning the figure scripts fixes resolution everywhere; for a poster, also export PDF or SVG.
- **Fonts:** consistent. Every figure uses `plotstyle.apply_style()` (Helvetica Neue, found on this machine; 9 pt body, 10 pt titles).
- **Palette:** the project palette was checked by simulating protanopia, deuteranopia and tritanopia (Machado et al. 2009, full severity; CIELAB ΔE between colours).
  - Blue `#2a78d6` against orange `#eb6834` stays far apart under every type (ΔE ≥ 97), and so does either colour against grey.
  - **Aqua `#1baf7a` is the problem:** against grey it collapses under deuteranopia (ΔE 11), against orange it is weak under protanopia (ΔE 25), and against blue it is weak under tritanopia (ΔE 22).
  - Figures that use aqua need it replaced (for example with a dark purple) or a second cue added (marker shape or direct labels).
- **Titles:** many figures carry long titles that name the internal specification ("bge-large-en-v1.5 / chunk / strict"). For publication, move them to captions and use plain labels such as "primary specification".

## Main narrative (final report)

| # | figure | path | needs regenerating because |
|---|---|---|---|
| 1 | Study design: pipeline from MDD, tree and Wikipedia to the three matrices and the tests | *not yet made* | new figure (diagram) |
| 2 | Species coverage: stub rate by order (who gets dropped) | `reports/figures/week6/stubs_by_order.png` | 300 dpi; the "not shown" footnote is small |
| 3 | Pilot H1: text vs phylogenetic distance, primary spec | `reports/figures/week4/h1_scatter_primary.png` | 300 dpi; plain-language title |
| 4 | Pilot H1 across the 12 specifications | `reports/figures/week4/h1_variants.png` | 300 dpi; **aqua (strict mask) is the colour of the primary spec**: change it or add marker shapes |
| 5 | Pilot H2: raw r and r after each control | `reports/figures/week5/h2_controls_primary.png` | 300 dpi |
| 6 | Full-scale H2: raw r and r after each control | `reports/figures/week6/h2_controls.png` | 300 dpi; best shown side by side with #5 on a shared x-axis |
| 7 | Within-order correlations, all mammals | `reports/figures/week6/per_order.png` | 300 dpi; Dasyuromorphia's two dots overlap (draw the H1 dot on top or offset it) |
| 8 | Length-matched analysis (post-hoc): r by length decile, same-order vs different-order distance | *not yet made* (`reports/tables/mammalia/posthoc_length_matched_bins.csv`) | new figure; label it post-hoc |
| 9 | Pageview quartiles: popular pairs more alike, obscure pairs track the tree | `reports/figures/week6/pageview_quartiles.png` (all mammals and by order); `reports/figures/week5/exploratory_well_known.png` (pilot) | 300 dpi; grey Q1 dots against blue Q4 are safe; label it exploratory |
| 10 | Register comparison, Wikipedia vs ADW | `reports/figures/week7/register_comparison.png` | 300 dpi; the x-axis upper limit leaves extra white space |
| 11 | Charisma content: agreement table and B1–B3 | *pending* (charisma results) | new figure once results exist |

## Supporting figures (report appendix or a poster side panel)

| figure | path | needs regenerating because |
|---|---|---|
| Pilot MRM coefficients | `reports/figures/week5/h2_mrm_primary.png` | 300 dpi |
| Full-scale MRM coefficients | `reports/figures/week6/h2_mrm.png` | 300 dpi |
| Pilot H1 on 100 posterior trees | `reports/figures/week4/h1_posterior.png` | 300 dpi |
| Pilot H2 on 100 posterior trees | `reports/figures/week5/h2_posterior.png` | 300 dpi |
| Full-scale posterior trees (20) | `reports/figures/week6/posterior.png` | 300 dpi; the range label "trees −0.060–0.011" reads badly with a minus sign: write "from −0.060 to 0.011" |
| Pilot H2 across 12 specifications | `reports/figures/week5/h2_variants.png` | 300 dpi; **aqua**, as #4 |
| Pilot transformed phylogeny (log, 30 My cap) | `reports/figures/week5/exploratory_nonlinear.png` | 300 dpi |
| Pilot H1 permutation null | `reports/figures/week4/h1_null_primary.png` | 300 dpi |
| Week 3 sanity: family ratio by mask level | `reports/figures/family_ratio_by_mask__*.png` | 300 dpi; **aqua**, as #4; check the axis labels |
| Week 3 sanity: distance heatmaps (taxonomy mask) | `reports/figures/heatmaps__*__taxonomy.png` | 300 dpi; check that the colour bars are labelled; sequential blue ramp is safe |
| Week 3 sanity: text vs phylogeny, length vs neighbours, within vs between family | `reports/figures/text_vs_phylo__*.png`, `length_vs_neighbors__*.png`, `family_within_between__*.png` | 300 dpi; probably appendix only |

## Poster shortlist

#1 (design), #3 or #6 (the controls), #7 (within orders), #9 (attention), #10 (register), and #8 if space allows. Each should be redrawn at poster size with larger fonts (at least 14 pt for labels at a typical 36 × 48 in print).
