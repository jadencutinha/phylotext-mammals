# Methods parameters

Every parameter of the analysis, with the place where it is set and the place where it was fixed. "Pilot" means Carnivora; "full" means all mammals (`--clade mammalia`). Settings in `config.yaml` are given by key.

## Sources and versions

| parameter | value | set in | fixed by |
|---|---|---|---|
| Taxonomy | MDD v2.5 (`MDD_v*_*species.csv`, `Species_Syn_Current_v*.csv`) | `sources.mdd` | prereg §2, Amendment 3 §3 |
| Tree (decisive) | Upham et al. 2019 MCC tree, `MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre`, 5,911 tips | `sources.tree.mcc_url` | prereg §2 |
| Tree posterior | `Completed_5911sp_topoCons_NDexp.zip`, 10,000 trees (VertLife) | `sources.tree.posterior_url` | Amendment 1, A1.3 |
| Posterior draws | `numpy.random.default_rng(20261004).choice(N, size=100, replace=False)`; full scale uses the first 20 of these | `scripts/14_posterior_h1.py`, `23_posterior_fullscale.py` | A1.3; Amendment 3 §8 |
| Text | English Wikipedia: the cleaned lead plus the first section whose heading is in `wikipedia.description_headings` | `wikipedia.*` | prereg §2; Amendment 3 §3 |
| Wikipedia client | User-Agent `lineage-or-charisma/0.1 (…; Princeton COS398 research)`; at least 0.5 s between requests; 5 retries; 2 s backoff | `wikipedia.*` | implementation |
| Traits | EltonTraits 1.0 `MamFuncDat.txt` (md5 `59c3eee2…`); COMBINE `trait_data_reported.csv` (md5 `d894586c…`), reported values only | `sources.elton_traits`, `sources.combine` | prereg §7 |
| Pageviews | Wikimedia Pageviews API, `en.wikipedia`, all-access, agent = user, 2023-01 to 2025-12, redirects included | `attention.*` | prereg §7 |
| ADW (secondary A) | `https://animaldiversity.org/accounts/<Genus>_<species>/`; route: MDD name, then tree binomial, then MDD synonyms; one request per 8 s; contact User-Agent; responses cached under `data/raw/adw/`; CC BY-NC-SA 3.0 | `src/lineage_charisma/adw.py` | Amendment 4 §A, Change 1 |

## Name matching

| parameter | value | set in | fixed by |
|---|---|---|---|
| Fuzzy thresholds | accept ≥ 92; review ≥ 85; minimum margin over the runner-up 3 | `matching.*` | Week 2; Amendment 3 §3 |
| Full-scale residual rule | Automatic acceptance only for exact, MSW3/CMW legacy-name or MDD-synonym tree matches. Everything else goes to review and is excluded unless the pilot hand-reviewed the same match. Synonym-only Wikipedia titles also go to review. | `scripts/03b_residual_rules.py` | Amendment 3, Change 1 |
| Trait-name matching | Exact, legacy-name and synonym matches only (full scale). Other matches count as missing. | `clades.mammalia.traits.auto_accept_methods` | Amendment 3, Change 1 |
| Pilot reviewed trait gaps | 9 species: 8 not in EltonTraits and *Felis catus* not in COMBINE | `traits.reviewed_missing` | Week 4 |

## Corpus and masking

| parameter | value | set in | fixed by |
|---|---|---|---|
| Stub threshold | Fewer than 50 tokens: bge-large tokenizer, unmasked text, special tokens not counted. Stubs are dropped from every matrix. | `corpus.min_tokens`, `matrices.exclude_stubs` | prereg §2; Amendment 3 §3 |
| Token counter | `BAAI/bge-large-en-v1.5` tokenizer | `corpus.tokenizer` | Week 3 |
| Mask levels | none, taxonomy, strict (primary: strict) | `masking.levels` | prereg §1 and §3 |
| Mask token | `[TAXON]` | `masking.token` | Week 3 |
| Pilot taxonomy mask | MDD order, suborder, family and subfamily names with derived forms; genera; epithets anywhere; synonyms; abbreviations ("P. leo"); a hand list (carnivoran, pinniped, Arctoidea, …) | `masking.*` | Week 3, prereg §3 |
| Pilot strict mask | Taxonomy, plus the head noun of every MDD common name and its plural. The head noun is the last word, plus the last part of a hyphenated word. | `masking.head_noun_rule: pilot` | Week 3 |
| Pilot keep rules | "canines" and "Iceland" never masked; "canine" kept before tooth, teeth or tip | `masking.keep_terms`, `keep_when_followed_by` | Week 3 |
| Full-scale head nouns | The whole final word, with hyphenated words kept whole. Names with a preposition give no head noun. | `clades.mammalia.masking.head_noun_rule: final_noun` | Amendment 3, Change 2 |
| Full-scale epithets | Masked only inside a binomial or an abbreviated binomial | `epithets_anywhere: false` | Amendment 3, Change 2 |
| Full-scale higher ranks | subclass through tribe, with derived forms | `clades.mammalia.masking.name_ranks` | Amendment 3, Change 2 |
| Full-scale extra terms | The pilot list plus order and clade vernaculars (primate, rodent, marsupial, …) | `clades.mammalia.masking.extra_taxonomy_terms` | Amendment 3, Change 2 |
| Full-scale keep terms | tail, head, sides, human, bull, cow, buck, stag, ounce, jack, digger, runner, canary, pinyin, fiber, boss, fossorial (with plurals) | `clades.mammalia.masking.keep_terms` | Amendment 3, Change 2 |
| Masking guard | The run stops if a masked text still contains its own genus name | `src/lineage_charisma/masking.py` | Amendment 3 §4 |

## Embeddings

| parameter | value | set in | fixed by |
|---|---|---|---|
| Models | `BAAI/bge-large-en-v1.5` (primary); `sentence-transformers/all-mpnet-base-v2` | `embedding.models` | prereg §3 and §5 |
| Passage limit | The model's own maximum, including the two special tokens: 512 for bge-large, 384 for mpnet | `embedding.max_tokens: null` | Week 3 |
| Long-text rules | chunk (primary): evenly sized overlapping passages, embeddings mean-pooled, then L2-normalized. truncate: the first passage only. | `embedding.rules` | prereg §3 |
| Chunk overlap | 32 tokens | `embedding.chunk_overlap` | Week 3; Amendment 4 §3 |
| Batch size, device | 16; auto (cuda, then mps, then cpu) | `embedding.*` | implementation |
| Pilot combinations | 2 models × 2 rules × 3 masks = 12 | `embedding` | prereg §5 |
| Full-scale combinations | bge/chunk/strict (primary), mpnet/chunk/strict, bge/truncate/strict, bge/chunk/taxonomy | `clades.mammalia.embedding.combos` | Amendment 3 §8 |

## Distance matrices

| matrix | definition | fixed by |
|---|---|---|
| text | 1 − cosine similarity of the embeddings | prereg §2 |
| phylogeny | Patristic distance in My on the pruned MCC tree (vectorized, bit-identical to dendropy) | prereg §2; Amendment 3 §4 |
| ecology (primary) | Gower distance with the six groups weighted equally and columns in a group sharing its weight. Continuous traits use \|Δ\| / range over the species analysed. Pairwise-available rule; no imputation. | prereg §7 |
| ecology_unweighted | Gower with every column weighted equally (robustness) | A1.1 |
| ecology_no_volant | Full scale: the primary matrix without `terrestrial_volant` (robustness) | Amendment 3 §4 and §8 |
| length | \|Δ log(n_tokens)\| | prereg §7 |
| pageviews | \|Δ log(1 + median monthly pageviews)\| | prereg §7 |
| combined attention | Euclidean distance on standardized log pageviews and log length (Week 4 diagnostics only; not used in H2) | prereg §7; Amendment 2 §4 |
| mean log pageviews | The mean of the two species' log pageviews (similarity covariate; exploratory quartiles) | prereg §7 |
| charisma content (B) | Euclidean distance on the standardized mean scores of retained dimensions (pending) | Amendment 4 §B |

## Trait list

| group | columns | source | type |
|---|---|---|---|
| diet | Diet-Inv, -Vend, -Vect, -Vfish, -Vunk, -Scav, -Fruit, -Nect, -Seed, -PlantO (10) | EltonTraits | continuous, 0–100 |
| foraging stratum | ForStrat-Value | EltonTraits | categorical |
| activity time | nocturnal, crepuscular, diurnal | EltonTraits | binary, not exclusive |
| habitat breadth | habitat breadth | COMBINE | continuous |
| terrestrial or aquatic | terrestrial, marine, freshwater; plus `terrestrial_volant` at full scale | COMBINE | binary |
| body mass | log10 adult body mass | COMBINE | continuous |

The pilot has 19 columns and the full scale 20. `terrestrial_volant` coincides exactly with Chiroptera (Amendment 3, Change 3). Source: `src/lineage_charisma/traits.py`, `docs/ecology_traits*.md`.

## Statistics

| parameter | value | set in | fixed by |
|---|---|---|---|
| Correlation | Spearman on upper triangles, average ranks for ties; Pearson reported as a variant | `stats.method` | prereg §4 |
| Permutations | 9,999 for the primary, secondary, MRM, per-order and replication tests, the pilot posterior trees, and post-hoc D. 999 for full-scale robustness variants, the 20 posterior trees and full-scale exploratory tests. | `stats.permutations`; script arguments | prereg §4; Amendment 3 §8 |
| Seed | 20261004, a fresh `numpy.random.default_rng` per test | `stats.seed` | prereg §4 |
| Permutation scheme | Rows and columns of the phylogenetic matrix permuted together. MRM permutes the text matrix; B1 permutes the charisma matrix; D permutes within length deciles. | `src/lineage_charisma/stats.py`, `condensed.py` | prereg §4; Amendment 2 §5 and §7; Amendment 4 |
| p-value | (1 + count ≥ observed) / (1 + permutations). One-sided for H1, H2 and the secondary tests. Two-sided for MRM (on the pseudo-t), Week 4 diagnostics and D. | `stats.py` | prereg §4; Amendment 2 §7; Amendment 4 §D |
| Alpha | 0.05; no multiple-testing correction | `stats.alpha` | prereg §4 |
| H1 decision | r > 0 with p < 0.05, primary specification, MCC tree | | prereg §4 |
| H2 primary model | Partial Mantel of text on phylogeny controlling for ecology and length jointly, on ranked upper triangles | | Amendment 2 §5 |
| H2 rule (a) | Partial r > 0 with p < 0.05 | | Amendment 2 §6 |
| H2 rule (b) | Partial r < 50% of the raw r on the H2 species; undefined if the raw r ≤ 0 | | Amendment 2 §6 |
| H2 secondary | S1 ecology; S2 length; S3 pageviews; S4 ecology + length + pageviews | | Amendment 2 §7 |
| MRM | text ~ phylogeny + ecology + length + pageviews; ranked and standardized | | Amendment 2 §7 |
| Full-scale engine | Condensed float32 vectors, ranks computed once, int32 index gathers, chunked permutations (identical results for any chunk size or thread count) | `src/lineage_charisma/condensed.py` | Amendment 3 §5 |
| Replication threshold | \|r_6b − r_pilot\| ≤ 0.01; 6a must match exactly | | Amendment 3 §6 |
| Per-order minimum | ≥ 50 species in the frozen set | | Amendment 3 §7 |
| Exploratory transforms | log(1 + d) with Pearson; min(d, 30 My) with Spearman and Pearson | | Amendment 2 §9; Amendment 3 §9 |
| Post-hoc D | 10 equal deciles of `n_tokens` (`rank(method="first")`, `qcut`); within-decile pairs; pair-weighted mean r; 9,999 within-decile permutations; two-sided. 5 and 20 bins are reported as a check. | `scripts/26_length_matched.py` | Amendment 4 §D |
| ADW comparison | Pilot strict masks, bge/chunk; minimum 100 usable species; Δ = r_ADW − r_Wiki; paired species bootstrap with 2,000 draws, seed 20261004, duplicate-copy pairs dropped; 2.5th and 97.5th percentiles | `scripts/30_register_comparison.py` | Amendment 4 §A |
| LLM annotation (B) | `claude-sonnet-5-5` via `anthropic==2.6.0`; default sampling; thinking off (`between_tools`); JSON-schema output; max_tokens 256; 5 retries per run; 2 runs; no fallback model; stops at a projected US$10 | `scripts/31_annotate.py`, `src/lineage_charisma/charisma.py` | Amendment 4 §B, Change 1 |
| Validation sample (B) | `default_rng(20261004).choice(284, size=50, replace=False)`; drop rule: Spearman agreement < 0.5 (author vs model, or run 1 vs run 2) | `scripts/27_annotation_sheet.py` | Amendment 4 §B |
| Zoo signage (C) | MDD reconciliation (current, legacy, synonym, then fuzzy 92/85/3); full-scale strict masks; H1 on the overlap only if ≥ 10 species; 999 permutations | `scripts/28_zoo_signage.py` | Amendment 4 §C, Change 1 |

## Species counts at each stage

| stage | pilot (Carnivora) | full scale (Mammalia) | source |
|---|---:|---:|---|
| MDD species | 320 | 6,904 | Amendment 3 §1; `species_by_order.csv` |
| excluded: no tree tip | 32 | 1,190 | `data/interim/excluded_species.csv`; week 6 report §2 |
| excluded: extinct or domestic without a tip | 3 | 43 | same |
| residual, not resolved by a written rule | 0 | 88 | same |
| matched to a tip and a Wikipedia article | 285 | 5,583 | `matched.csv`; Amendment 3 frozen set |
| stubs (< 50 tokens) | 1 | 972 | `stub_species.csv`; `stubs_by_order.csv` |
| **analysed (H1)** | **284** | **4,611** | `species_order.txt` (full-scale SHA-256 `4730f407…763f0e`) |
| H2 species (ecological distance to every other species) | 283 | 4,519 | `ecology_species.txt` |
| complete-trait subset (robustness) | 272 | — | Amendment 2 §8 |
| orders analysed separately (≥ 50 species) | — | 10 | Amendment 3, frozen set |
| ADW accounts / usable (≥ 50 tokens) / H2 | 279 / 246 / 245 | — | `adw_coverage.csv`, `register_comparison.csv` |
| LLM validation sample | 50 | — | Amendment 4 §B |

## Software

The Python version and pinned packages are in `requirements.lock`; the `README` gives the run order. The pilot `stats.py` code is the reference for the full-scale engine (`tests/test_condensed.py`).
