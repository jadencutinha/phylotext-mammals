# Final results index

Every pre-registered result, the replication checks and the one disclosed post-hoc analysis, generated from the result tables by `scripts/33_results_index.py` (rerun it after any table changes). Roles: **primary** (the registered H1 and H2 decisions), **secondary** (registered tests that cannot change the H1/H2 conclusions), **robustness** (registered one-at-a-time variants), **exploratory** (registered, labelled exploratory or descriptive), **check** (registered replication checks), and **post-hoc** (designed after seeing the data, disclosed in amendment 4, section D). p-values are one-sided unless the decision column says otherwise; MRM p-values are two-sided.

Counts: primary 4, secondary 24, robustness 49, exploratory 23, check 2, post-hoc 1.

## Pilot H1

| result | role | spec | n | statistic | p | perms | decision | source |
|---|---|---|---:|---|---:|---:|---|---|
| H1 Mantel r | primary | bge-large-en-v1.5 / chunk / strict | 284 | r = 0.033 | 0.0285 | 9999 | supported | `reports/tables/h1_mantel.csv` |
| H1 Mantel r (Pearson) | robustness | bge-large-en-v1.5 / chunk / strict | 284 | r = 0.098 | 0.0001 | 9999 |  | `reports/tables/h1_mantel.csv` |
| H1 Mantel r | robustness | bge-large-en-v1.5 / truncate / strict | 284 | r = 0.064 | 0.0001 | 9999 |  | `reports/tables/h1_mantel.csv` |
| H1 Mantel r (Pearson) | robustness | bge-large-en-v1.5 / truncate / strict | 284 | r = 0.131 | 0.0001 | 9999 |  | `reports/tables/h1_mantel.csv` |
| H1 Mantel r | robustness | all-mpnet-base-v2 / chunk / strict | 284 | r = 0.083 | 0.0001 | 9999 |  | `reports/tables/h1_mantel.csv` |
| H1 Mantel r (Pearson) | robustness | all-mpnet-base-v2 / chunk / strict | 284 | r = 0.142 | 0.0001 | 9999 |  | `reports/tables/h1_mantel.csv` |
| H1 Mantel r | robustness | all-mpnet-base-v2 / truncate / strict | 284 | r = 0.109 | 0.0001 | 9999 |  | `reports/tables/h1_mantel.csv` |
| H1 Mantel r (Pearson) | robustness | all-mpnet-base-v2 / truncate / strict | 284 | r = 0.176 | 0.0001 | 9999 |  | `reports/tables/h1_mantel.csv` |
| H1 Mantel r | robustness | bge-large-en-v1.5 / chunk / taxonomy | 284 | r = 0.215 | 0.0001 | 9999 |  | `reports/tables/h1_mantel.csv` |
| H1 Mantel r (Pearson) | robustness | bge-large-en-v1.5 / chunk / taxonomy | 284 | r = 0.338 | 0.0001 | 9999 |  | `reports/tables/h1_mantel.csv` |
| H1 Mantel r | robustness | bge-large-en-v1.5 / truncate / taxonomy | 284 | r = 0.258 | 0.0001 | 9999 |  | `reports/tables/h1_mantel.csv` |
| H1 Mantel r (Pearson) | robustness | bge-large-en-v1.5 / truncate / taxonomy | 284 | r = 0.398 | 0.0001 | 9999 |  | `reports/tables/h1_mantel.csv` |
| H1 Mantel r | robustness | all-mpnet-base-v2 / chunk / taxonomy | 284 | r = 0.355 | 0.0001 | 9999 |  | `reports/tables/h1_mantel.csv` |
| H1 Mantel r (Pearson) | robustness | all-mpnet-base-v2 / chunk / taxonomy | 284 | r = 0.482 | 0.0001 | 9999 |  | `reports/tables/h1_mantel.csv` |
| H1 Mantel r | robustness | all-mpnet-base-v2 / truncate / taxonomy | 284 | r = 0.378 | 0.0001 | 9999 |  | `reports/tables/h1_mantel.csv` |
| H1 Mantel r (Pearson) | robustness | all-mpnet-base-v2 / truncate / taxonomy | 284 | r = 0.515 | 0.0001 | 9999 |  | `reports/tables/h1_mantel.csv` |
| H1 Mantel r | robustness | bge-large-en-v1.5 / chunk / none | 284 | r = 0.204 | 0.0001 | 9999 |  | `reports/tables/h1_mantel.csv` |
| H1 Mantel r (Pearson) | robustness | bge-large-en-v1.5 / chunk / none | 284 | r = 0.327 | 0.0001 | 9999 |  | `reports/tables/h1_mantel.csv` |
| H1 Mantel r | robustness | bge-large-en-v1.5 / truncate / none | 284 | r = 0.252 | 0.0001 | 9999 |  | `reports/tables/h1_mantel.csv` |
| H1 Mantel r (Pearson) | robustness | bge-large-en-v1.5 / truncate / none | 284 | r = 0.399 | 0.0001 | 9999 |  | `reports/tables/h1_mantel.csv` |
| H1 Mantel r | robustness | all-mpnet-base-v2 / chunk / none | 284 | r = 0.358 | 0.0001 | 9999 |  | `reports/tables/h1_mantel.csv` |
| H1 Mantel r (Pearson) | robustness | all-mpnet-base-v2 / chunk / none | 284 | r = 0.489 | 0.0001 | 9999 |  | `reports/tables/h1_mantel.csv` |
| H1 Mantel r | robustness | all-mpnet-base-v2 / truncate / none | 284 | r = 0.373 | 0.0001 | 9999 |  | `reports/tables/h1_mantel.csv` |
| H1 Mantel r (Pearson) | robustness | all-mpnet-base-v2 / truncate / none | 284 | r = 0.513 | 0.0001 | 9999 |  | `reports/tables/h1_mantel.csv` |
| H1 on 100 posterior trees | robustness | primary spec | 284 | median r = 0.080 [0.004, 0.161] |  | 9999 | p < 0.05 on 89/100 trees | `reports/tables/h1_posterior.csv` |
| H1 within 6 families | exploratory | primary spec |  | r from -0.078 to 0.200 |  | 9999 | descriptive | `reports/tables/h1_mantel_families.csv` |

## Pilot H2

| result | role | spec | n | statistic | p | perms | decision | source |
|---|---|---|---:|---|---:|---:|---|---|
| H2 partial r (ecology + length) | primary | bge-large-en-v1.5 / chunk / strict | 283 | partial r = 0.020 (raw 0.035, share 0.58) | 0.1272 | 9999 | rule (a) fails; rule (b) fails | `reports/tables/h2_partial_mantel.csv` |
| H2 S1 partial r (ecology) | secondary | bge-large-en-v1.5 / chunk / strict | 283 | partial r = 0.024 (raw 0.035, share 0.69) | 0.0845 | 9999 |  | `reports/tables/h2_partial_mantel.csv` |
| H2 S2 partial r (length) | secondary | bge-large-en-v1.5 / chunk / strict | 283 | partial r = 0.031 (raw 0.035, share 0.88) | 0.0413 | 9999 |  | `reports/tables/h2_partial_mantel.csv` |
| H2 S3 partial r (pageviews) | secondary | bge-large-en-v1.5 / chunk / strict | 283 | partial r = 0.037 (raw 0.035, share 1.06) | 0.0185 | 9999 |  | `reports/tables/h2_partial_mantel.csv` |
| H2 S4 partial r (ecology + length + pageviews) | secondary | bge-large-en-v1.5 / chunk / strict | 283 | partial r = 0.027 (raw 0.035, share 0.77) | 0.0623 | 9999 |  | `reports/tables/h2_partial_mantel.csv` |
| H2 partial r, variant | robustness | bge-large-en-v1.5 / truncate / strict | 283 | partial r = 0.047 (raw 0.066, share 0.72) | 0.0020 | 9999 |  | `reports/tables/h2_partial_mantel.csv` |
| H2 partial r, variant | robustness | all-mpnet-base-v2 / chunk / strict | 283 | partial r = 0.061 (raw 0.084, share 0.72) | 0.0003 | 9999 |  | `reports/tables/h2_partial_mantel.csv` |
| H2 partial r, variant | robustness | all-mpnet-base-v2 / truncate / strict | 283 | partial r = 0.086 (raw 0.110, share 0.78) | 0.0001 | 9999 |  | `reports/tables/h2_partial_mantel.csv` |
| H2 partial r, variant | robustness | bge-large-en-v1.5 / chunk / taxonomy | 283 | partial r = 0.201 (raw 0.217, share 0.93) | 0.0001 | 9999 |  | `reports/tables/h2_partial_mantel.csv` |
| H2 partial r, variant | robustness | bge-large-en-v1.5 / truncate / taxonomy | 283 | partial r = 0.244 (raw 0.259, share 0.94) | 0.0001 | 9999 |  | `reports/tables/h2_partial_mantel.csv` |
| H2 partial r, variant | robustness | all-mpnet-base-v2 / chunk / taxonomy | 283 | partial r = 0.346 (raw 0.355, share 0.97) | 0.0001 | 9999 |  | `reports/tables/h2_partial_mantel.csv` |
| H2 partial r, variant | robustness | all-mpnet-base-v2 / truncate / taxonomy | 283 | partial r = 0.371 (raw 0.378, share 0.98) | 0.0001 | 9999 |  | `reports/tables/h2_partial_mantel.csv` |
| H2 partial r, variant | robustness | bge-large-en-v1.5 / chunk / none | 283 | partial r = 0.191 (raw 0.206, share 0.93) | 0.0001 | 9999 |  | `reports/tables/h2_partial_mantel.csv` |
| H2 partial r, variant | robustness | bge-large-en-v1.5 / truncate / none | 283 | partial r = 0.240 (raw 0.255, share 0.94) | 0.0001 | 9999 |  | `reports/tables/h2_partial_mantel.csv` |
| H2 partial r, variant | robustness | all-mpnet-base-v2 / chunk / none | 283 | partial r = 0.349 (raw 0.358, share 0.97) | 0.0001 | 9999 |  | `reports/tables/h2_partial_mantel.csv` |
| H2 partial r, variant | robustness | all-mpnet-base-v2 / truncate / none | 283 | partial r = 0.365 (raw 0.374, share 0.98) | 0.0001 | 9999 |  | `reports/tables/h2_partial_mantel.csv` |
| H2 partial r, unweighted gower | robustness | bge-large-en-v1.5 / chunk / strict, ecology_unweighted | 283 | partial r = 0.021 (raw 0.035, share 0.60) | 0.1174 | 9999 |  | `reports/tables/h2_partial_mantel.csv` |
| H2 partial r, complete traits | robustness | bge-large-en-v1.5 / chunk / strict | 272 | partial r = 0.021 (raw 0.038, share 0.55) | 0.1216 | 9999 |  | `reports/tables/h2_partial_mantel.csv` |
| MRM coefficient: phylogeny | secondary | primary spec, descriptive | 283 | β = 0.027 (R² = 0.0317) | 0.1318 | 9999 | descriptive (two-sided) | `reports/tables/h2_mrm.csv` |
| MRM coefficient: ecology | secondary | primary spec, descriptive | 283 | β = 0.135 (R² = 0.0317) | 0.0004 | 9999 | descriptive (two-sided) | `reports/tables/h2_mrm.csv` |
| MRM coefficient: length | secondary | primary spec, descriptive | 283 | β = 0.121 (R² = 0.0317) | 0.0002 | 9999 | descriptive (two-sided) | `reports/tables/h2_mrm.csv` |
| MRM coefficient: pageviews | secondary | primary spec, descriptive | 283 | β = -0.082 (R² = 0.0317) | 0.0059 | 9999 | descriptive (two-sided) | `reports/tables/h2_mrm.csv` |
| H2 on 100 posterior trees | robustness | primary spec | 283 | median partial r = 0.062 |  | 9999 | rule (a) on 85/100; rule (b) on 8/100 | `reports/tables/h2_posterior.csv` |

## Pilot exploratory

| result | role | spec | n | statistic | p | perms | decision | source |
|---|---|---|---:|---|---:|---:|---|---|
| text–phylogeny r, Q1 (least viewed) pageview quartile | exploratory | primary spec | 10047 | r = 0.133 | 0.0001 | 9999 |  | `reports/tables/exploratory_well_known.csv` |
| text–phylogeny r, Q4 (most viewed) pageview quartile | exploratory | primary spec | 10047 | r = 0.035 | 0.0455 | 9999 |  | `reports/tables/exploratory_well_known.csv` |
| H1 with phylogeny capped (spearman) | exploratory | primary spec | 284 | r = 0.115; partial r = 0.078 | 0.0001 | 9999 |  | `reports/tables/exploratory_nonlinear.csv` |
| H1 with phylogeny none (pearson) | exploratory | primary spec | 284 | r = 0.098; partial r = 0.070 | 0.0001 | 9999 |  | `reports/tables/exploratory_nonlinear.csv` |
| H1 with phylogeny log (pearson) | exploratory | primary spec | 284 | r = 0.120; partial r = 0.081 | 0.0001 | 9999 |  | `reports/tables/exploratory_nonlinear.csv` |
| H1 with phylogeny capped (pearson) | exploratory | primary spec | 284 | r = 0.136; partial r = 0.090 | 0.0001 | 9999 |  | `reports/tables/exploratory_nonlinear.csv` |

## Full-scale checks

| result | role | spec | n | statistic | p | perms | decision | source |
|---|---|---|---:|---|---:|---:|---|---|
| 6a exact (pilot vocabulary) | check | primary spec, 284 pilot species | 284 | r = 0.032706 (pilot 0.032706) | 0.0285 | 9999 | passed | `reports/tables/mammalia/replication.csv` |
| 6b full-scale vocabulary | check | primary spec, 284 pilot species | 284 | r = 0.028857 (pilot 0.032706) | 0.0480 | 9999 | passed | `reports/tables/mammalia/replication.csv` |

## Full-scale H1

| result | role | spec | n | statistic | p | perms | decision | source |
|---|---|---|---:|---|---:|---:|---|---|
| H1 Mantel r | primary | bge-large-en-v1.5 / chunk / strict | 4611 | r = -0.028 | 1.0000 | 9999 | not supported | `reports/tables/mammalia/h1_mantel.csv` |
| H1 Mantel r | robustness | all-mpnet-base-v2 / chunk / strict | 4611 | r = 0.001 | 0.4440 | 999 |  | `reports/tables/mammalia/h1_mantel.csv` |
| H1 Mantel r | robustness | bge-large-en-v1.5 / truncate / strict | 4611 | r = -0.025 | 0.9990 | 999 |  | `reports/tables/mammalia/h1_mantel.csv` |
| H1 Mantel r | robustness | bge-large-en-v1.5 / chunk / taxonomy | 4611 | r = 0.186 | 0.0010 | 999 |  | `reports/tables/mammalia/h1_mantel.csv` |
| H1 Mantel r (Pearson) | robustness | bge-large-en-v1.5 / chunk / strict | 4611 | r = -0.022 | 0.9840 | 999 |  | `reports/tables/mammalia/h1_mantel.csv` |
| H1 on 20 posterior trees | robustness | primary spec | 4611 | median r = -0.015 [-0.060, 0.011] |  | 999 | p < 0.05 on 0% of trees | `reports/tables/mammalia/posterior_summary.csv` |

## Full-scale H2

| result | role | spec | n | statistic | p | perms | decision | source |
|---|---|---|---:|---|---:|---:|---|---|
| H2 partial r (ecology + length) | primary | bge-large-en-v1.5 / chunk / strict | 4519 | partial r = -0.034 (raw -0.029) | 1.0000 | 9999 | rule (a) fails; rule (b) undefined (raw r ≤ 0) | `reports/tables/mammalia/h2_partial_mantel.csv` |
| H2 S1 partial r (ecology) | secondary | bge-large-en-v1.5 / chunk / strict | 4519 | partial r = -0.033 (raw -0.029) | 1.0000 | 9999 |  | `reports/tables/mammalia/h2_partial_mantel.csv` |
| H2 S2 partial r (length) | secondary | bge-large-en-v1.5 / chunk / strict | 4519 | partial r = -0.030 (raw -0.029) | 0.9999 | 9999 |  | `reports/tables/mammalia/h2_partial_mantel.csv` |
| H2 S3 partial r (pageviews) | secondary | bge-large-en-v1.5 / chunk / strict | 4519 | partial r = -0.015 (raw -0.029) | 0.9721 | 9999 |  | `reports/tables/mammalia/h2_partial_mantel.csv` |
| H2 S4 partial r (ecology + length + pageviews) | secondary | bge-large-en-v1.5 / chunk / strict | 4519 | partial r = -0.023 (raw -0.029) | 0.9984 | 9999 |  | `reports/tables/mammalia/h2_partial_mantel.csv` |
| H2 partial r, variant | robustness | all-mpnet-base-v2 / chunk / strict | 4519 | partial r = -0.012 (raw 0.000) | 0.9430 | 999 |  | `reports/tables/mammalia/h2_partial_mantel.csv` |
| H2 partial r, variant | robustness | bge-large-en-v1.5 / truncate / strict | 4519 | partial r = -0.035 (raw -0.026) | 1.0000 | 999 |  | `reports/tables/mammalia/h2_partial_mantel.csv` |
| H2 partial r, variant | robustness | bge-large-en-v1.5 / chunk / taxonomy | 4519 | partial r = 0.162 (raw 0.186) | 0.0010 | 999 |  | `reports/tables/mammalia/h2_partial_mantel.csv` |
| H2 partial r, unweighted gower | robustness | bge-large-en-v1.5 / chunk / strict, ecology_unweighted | 4519 | partial r = -0.038 (raw -0.029) | 1.0000 | 999 |  | `reports/tables/mammalia/h2_partial_mantel.csv` |
| H2 partial r, ecology without flight | robustness | bge-large-en-v1.5 / chunk / strict, ecology_no_volant | 4519 | partial r = -0.031 (raw -0.029) | 1.0000 | 999 |  | `reports/tables/mammalia/h2_partial_mantel.csv` |
| MRM coefficient: phylogeny | secondary | primary spec, descriptive | 4519 | β = -0.023 (R² = 0.0338) | 0.0058 | 9999 | descriptive (two-sided) | `reports/tables/mammalia/h2_mrm.csv` |
| MRM coefficient: ecology | secondary | primary spec, descriptive | 4519 | β = 0.058 (R² = 0.0338) | 0.0001 | 9999 | descriptive (two-sided) | `reports/tables/mammalia/h2_mrm.csv` |
| MRM coefficient: length | secondary | primary spec, descriptive | 4519 | β = 0.085 (R² = 0.0338) | 0.0001 | 9999 | descriptive (two-sided) | `reports/tables/mammalia/h2_mrm.csv` |
| MRM coefficient: pageviews | secondary | primary spec, descriptive | 4519 | β = -0.190 (R² = 0.0338) | 0.0001 | 9999 | descriptive (two-sided) | `reports/tables/mammalia/h2_mrm.csv` |
| H2 on 20 posterior trees | robustness | primary spec | 4519 | median partial r = -0.022 |  | 999 | rule (a) on 0%; rule (b) on 20% (raw r near 0) | `reports/tables/mammalia/posterior_summary.csv` |

## Full-scale per order

| result | role | spec | n | statistic | p | perms | decision | source |
|---|---|---|---:|---|---:|---:|---|---|
| Rodentia: H1 r / H2 partial r | exploratory | primary spec, descriptive | 1722 | H1 r = -0.006; partial r = -0.020 | 0.6887 | 9999 | descriptive (H1 p shown) | `reports/tables/mammalia/per_order.csv` |
| Chiroptera: H1 r / H2 partial r | exploratory | primary spec, descriptive | 989 | H1 r = 0.214; partial r = 0.204 | 0.0001 | 9999 | descriptive (H1 p shown) | `reports/tables/mammalia/per_order.csv` |
| Primates: H1 r / H2 partial r | exploratory | primary spec, descriptive | 406 | H1 r = 0.139; partial r = 0.159 | 0.0001 | 9999 | descriptive (H1 p shown) | `reports/tables/mammalia/per_order.csv` |
| Artiodactyla: H1 r / H2 partial r | exploratory | primary spec, descriptive | 329 | H1 r = 0.283; partial r = 0.175 | 0.0001 | 9999 | descriptive (H1 p shown) | `reports/tables/mammalia/per_order.csv` |
| Eulipotyphla: H1 r / H2 partial r | exploratory | primary spec, descriptive | 325 | H1 r = -0.008; partial r = -0.027 | 0.6151 | 9999 | descriptive (H1 p shown) | `reports/tables/mammalia/per_order.csv` |
| Carnivora: H1 r / H2 partial r | exploratory | primary spec, descriptive | 284 | H1 r = 0.029; partial r = 0.017 | 0.0480 | 9999 | descriptive (H1 p shown) | `reports/tables/mammalia/per_order.csv` |
| Diprotodontia: H1 r / H2 partial r | exploratory | primary spec, descriptive | 139 | H1 r = 0.063; partial r = 0.027 | 0.1024 | 9999 | descriptive (H1 p shown) | `reports/tables/mammalia/per_order.csv` |
| Lagomorpha: H1 r / H2 partial r | exploratory | primary spec, descriptive | 85 | H1 r = 0.143; partial r = 0.063 | 0.0011 | 9999 | descriptive (H1 p shown) | `reports/tables/mammalia/per_order.csv` |
| Didelphimorphia: H1 r / H2 partial r | exploratory | primary spec, descriptive | 79 | H1 r = 0.021; partial r = 0.023 | 0.3510 | 9999 | descriptive (H1 p shown) | `reports/tables/mammalia/per_order.csv` |
| Dasyuromorphia: H1 r / H2 partial r | exploratory | primary spec, descriptive | 69 | H1 r = -0.033; partial r = -0.034 | 0.7350 | 9999 | descriptive (H1 p shown) | `reports/tables/mammalia/per_order.csv` |

## Full-scale exploratory

| result | role | spec | n | statistic | p | perms | decision | source |
|---|---|---|---:|---|---:|---:|---|---|
| text–phylogeny r, Q1 (least viewed) pageview quartile | exploratory | primary spec | 2657211 | r = 0.043 | 0.0010 | 999 |  | `reports/tables/mammalia/exploratory_pageview_quartiles.csv` |
| text–phylogeny r, Q4 (most viewed) pageview quartile | exploratory | primary spec | 2657085 | r = 0.024 | 0.0060 | 999 |  | `reports/tables/mammalia/exploratory_pageview_quartiles.csv` |
| H1 with phylogeny capped (spearman) | exploratory | primary spec | 4611 | r = 0.074; partial r = 0.071 | 0.0010 | 999 |  | `reports/tables/mammalia/exploratory_nonlinear.csv` |
| H1 with phylogeny none (pearson) | exploratory | primary spec | 4611 | r = -0.022; partial r = -0.027 | 0.9840 | 999 |  | `reports/tables/mammalia/exploratory_nonlinear.csv` |
| H1 with phylogeny log (pearson) | exploratory | primary spec | 4611 | r = 0.017; partial r = 0.010 | 0.0190 | 999 |  | `reports/tables/mammalia/exploratory_nonlinear.csv` |
| H1 with phylogeny capped (pearson) | exploratory | primary spec | 4611 | r = 0.086; partial r = 0.083 | 0.0010 | 999 |  | `reports/tables/mammalia/exploratory_nonlinear.csv` |

## Full-scale post-hoc

| result | role | spec | n | statistic | p | perms | decision | source |
|---|---|---|---:|---|---:|---:|---|---|
| Length-matched within-decile r (amendment 4, D) | post-hoc | primary spec, 10 length deciles | 4611 | r_LM = 0.0166 | 0.0245 | 9999 | positive, p < 0.05: consistent with length composition (post-hoc) | `reports/tables/mammalia/posthoc_length_matched.csv` |

## Secondary A: register

| result | role | spec | n | statistic | p | perms | decision | source |
|---|---|---|---:|---|---:|---:|---|---|
| H1 r, Wikipedia | secondary | bge-large / chunk / strict, pilot masks | 246 | r = 0.019 | 0.1753 | 9999 |  | `reports/tables/secondary/register_comparison.csv` |
| H2 partial r, Wikipedia | secondary | bge-large / chunk / strict, pilot masks | 245 | r = 0.012 | 0.2719 | 9999 |  | `reports/tables/secondary/register_comparison.csv` |
| H1 r, ADW | secondary | bge-large / chunk / strict, pilot masks | 246 | r = 0.075 | 0.0002 | 9999 |  | `reports/tables/secondary/register_comparison.csv` |
| H2 partial r, ADW | secondary | bge-large / chunk / strict, pilot masks | 245 | r = 0.064 | 0.0008 | 9999 |  | `reports/tables/secondary/register_comparison.csv` |
| Δ H1 (ADW − Wikipedia) | secondary | paired species bootstrap, 2,000 draws | 246 | Δ = +0.056 [+0.019, +0.095] |  |  | ADW higher | `reports/tables/secondary/register_bootstrap.csv` |
| Δ H2 partial (ADW − Wikipedia) | secondary | paired species bootstrap, 2,000 draws | 245 | Δ = +0.052 [+0.013, +0.090] |  |  | ADW higher | `reports/tables/secondary/register_bootstrap.csv` |

## Secondary B: charisma

| result | role | spec | n | statistic | p | perms | decision | source |
|---|---|---|---:|---|---:|---:|---|---|
| Agreement, drop rule, B1–B3 | secondary | claude-sonnet-5-5 |  | pending |  |  | pending: waiting for hand labels | `docs/preregistration_amendment_secondary.md` |

## Secondary C: zoo signage

| result | role | spec | n | statistic | p | perms | decision | source |
|---|---|---|---:|---|---:|---:|---|---|
| Cross-register similarity; H1 on overlap | secondary | descriptive |  | pending |  |  | pending: no signage collected yet | `scripts/28_zoo_signage.py` |
