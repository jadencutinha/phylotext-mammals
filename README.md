# lineage-or-charisma

Does text-embedding similarity between mammal species descriptions recover the evolutionary tree, or does it mainly reflect human attention and charisma? This repo compares a text-distance matrix (from Wikipedia description embeddings) with a phylogenetic-distance matrix (from the Upham et al. 2019 mammal tree) using Mantel and partial Mantel tests.

**Current state: Week 3, aligned text-distance and phylogenetic-distance matrices for the pilot clade at three name-masking levels, with a descriptive sanity report in [`reports/week3_sanity.md`](reports/week3_sanity.md).** No Mantel tests, trait or attention matrices yet. The target clade is set in `config.yaml` (currently order Carnivora). No code names a clade directly, so moving to all mammals only means changing `target_clade`.

## Setup

```bash
make setup          # creates .venv and installs the package + pytest
make test           # runs the unit tests (no network needed)
make all            # runs every stage; each stage skips work whose output is cached
make week2          # stages 1-5 only (data foundation)
make week3          # stages 5b-8 only (corpus, embeddings, matrices, sanity report)
make all FORCE=1    # rebuilds every stage, including the downloads and API pulls
```

Python 3.11+ is required. Dependencies: requests, dendropy, rapidfuzz, pandas, pyyaml, numpy, pyarrow, sentence-transformers (which brings PyTorch), matplotlib (plus pytest for dev).

The first `make embed` downloads the two embedding models from the Hugging Face hub (about 1.3 GB and 0.4 GB) and the first `make corpus` downloads one tokenizer file. Both are free and need no account. Everything after that runs offline from the caches.

## Pipeline

| Stage | Command | Reads | Writes |
|---|---|---|---|
| 1 | `make mdd` → `scripts/01_fetch_mdd.py` | MDD zip from the MDD GitHub site | `data/raw/mdd/` (full release), `data/raw/mdd_species_<clade>.csv`, `data/raw/mdd_synonyms_<clade>.csv` |
| 2 | `make tree` → `scripts/02_fetch_tree.py` | Upham 2019 MCC tree | `data/raw/tree/*.tre`, `data/interim/tree_tips_all.csv`, `data/interim/tree_tips_<clade>.csv` |
| 3 | `make reconcile` → `scripts/03_reconcile_names.py` | stage 1 + 2 outputs, Wikipedia + Wikidata APIs | `data/interim/matched.csv`, `residual_review.csv`, `excluded_species.csv`, `unused_tree_tips.csv` |
| 4 | `make prune` → `scripts/04_prune_tree.py` | `matched.csv`, tree | `data/interim/tree_pruned.nwk`, `tree_pruned_mdd_names.nwk`, `cophenetic_distance.csv`, `sanity_check.txt` |
| 5 | `make wiki` → `scripts/05_fetch_wikipedia.py` | `matched.csv`, Wikipedia API | `data/interim/wiki/raw/*.json`, `data/interim/wiki/clean/*.txt`, `data/interim/wiki_texts.csv` |
| 5b | `make corpus` → `scripts/05b_build_corpus.py` | `matched.csv`, `wiki/clean/*.txt`, MDD species and synonym tables | `data/processed/corpus.parquet`, `data/interim/stub_species.csv`, `data/interim/mask_terms.csv` |
| 6 | `make embed` → `scripts/06_embed.py` | `corpus.parquet`, embedding models | `data/processed/embeddings/{model_slug}__{rule}__{mask_level}.npy`, `.species.txt`, `.meta.json` |
| 7 | `make matrices` → `scripts/07_build_matrices.py` | embeddings, `tree_pruned_mdd_names.nwk`, `corpus.parquet` | `data/processed/matrices/species_order.txt`, `phylo.npy`, `text__{model_slug}__{rule}__{mask_level}.npy` |
| 8 | `make sanity` → `scripts/08_sanity_report.py` | matrices, `corpus.parquet`, `matched.csv` | `reports/week3_sanity.md`, `reports/figures/*.png` |

Every script takes `--force` and `--config PATH`. Stages 3 and 5 also take `--offline`, which answers every API request from the cache and fails if a response is missing. Every CSV written gets a `.meta.json` sidecar with its sha256, row count, creation time, and source (for example, the MDD version).

Optional: `make posterior` downloads the full 10,000-tree posterior from VertLife (`Completed_5911sp_topoCons_NDexp.zip`, about 1.2 GB) for later weeks. The tip labels are the same across the posterior, so this week's reconciliation applies to every tree in it.

### Data sources (checked 2026-09-24)

- **MDD v2.5** (released 2026-07-28): `https://raw.githubusercontent.com/mammaldiversity/mammaldiversity.github.io/master/assets/data/MDD.zip`. Stage 1 picks the highest-versioned `MDD_v*_*species.csv` and `Species_Syn_Current_v*.csv` inside the zip, so a new MDD release is picked up without code changes.
- **Upham et al. 2019 tree**: the MCC tree of the *completed, topology-constrained, node-dated (NDexp)* posterior, from the authors' repo `n8upham/MamPhy_v1` (`MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre`). The matching 10k posterior is on `data.vertlife.org/mammaltree/`.

## Week 3: corpus, embeddings, and distance matrices

The species ID everywhere is the MDD binomial (`Panthera leo`), the `species` column of `matched.csv`. Tree tips use the same name with an underscore.

### Corpus (stage 5b)

`corpus.parquet` has one row per matched species with the columns `species_id`, `family`, `text`, `n_tokens`, `is_stub`, plus `text_<level>` and `n_masked_<level>` for each masked level. `text`, `n_tokens` and `is_stub` always describe the unmasked text, so the stub set is the same at every mask level. `n_tokens` is counted with the tokenizer named in `corpus.tokenizer`, without special tokens. A species whose text is shorter than `corpus.min_tokens` (50) gets `is_stub = True` and a row in `stub_species.csv`.

Stubs are flagged and still embedded. `matrices.exclude_stubs: true` (the current setting) leaves them out of every matrix and of `species_order.txt`. For a sensitivity check, `make matrices-sensitivity` (`scripts/07_build_matrices.py --sensitivity`) drops every species under `matrices.sensitivity_min_tokens` (100) and writes a complete, separately ordered set of matrices to `data/processed/matrices_min100/`. It reuses the cached embeddings and does not touch the primary matrices.

### Name masking (stage 5b)

`masking.levels` in `config.yaml` lists the mask levels to build; the first is the primary one (currently `taxonomy`). Every level is embedded and gets its own matrices.

- **`none`**: the cleaned text as fetched.
- **`taxonomy`**: scientific names are replaced by `masking.token` (`[TAXON]`). The vocabulary is built from the MDD tables for the whole target clade, not just the species being described:
  - order, suborder, family, subfamily and genus names, matched case-insensitively with plurals, wherever they occur;
  - terms derived from family-group names (`Felidae` → `felid`, `feline`; `Lutrinae` → `lutrine`; `Caniformia` → `caniform`), plus `masking.extra_taxonomy_terms` for names MDD's columns do not supply (`pinniped`, `carnivoran`, tribes and superfamilies);
  - current species epithets wherever they occur;
  - a genus followed by one or two epithets as a single mask, where the epithets may be any current or MDD synonym name, and abbreviated forms such as `P. leo` and `P. l. persica`;
  - genus names from MDD synonyms (`Felis` for the lion's original combination);
  - any word containing a genus name of six or more letters (`Propoecilogale`).
- **`strict`**: taxonomy masking plus the head noun of every MDD common name (main and other common names), with plurals: `cat`, `bear`, `fox`, `wolf`, `otter`, `mongoose`, `seal`, and so on.

Rules that keep ordinary English out of the mask:

- Synonym epithets are masked only inside a scientific name, because many are English words (`major`, `minor`, `vision`).
- A historical genus spelled like a current common name (`Hyena`, `Coati`, `Serval`) is left to the strict level. Current genera and epithets that double as common names (`lynx`, `caracal`, `puma`, `serval`, `binturong`) are masked at the taxonomy level.
- `masking.keep_terms` are never masked (`canines`), and `masking.keep_when_followed_by` spares a term before listed words (`canine teeth` stays, `a canine native to` is masked).

Adjacent masks collapse into one, so a binomial becomes a single `[TAXON]`. Names of taxa outside the clade (prey species, for example) are not masked. `mask_terms.csv` lists every replaced surface form with its count, per level, for auditing. The corpus build fails if any masked text still contains its own genus name, and a test checks the same on the built corpus.

### Embeddings (stage 6)

`embedding.models`, `embedding.rules` and `masking.levels` in `config.yaml` are lists, and stage 6 embeds the corpus once per model, rule and mask level. Adding a model is one more line under `models`. The first entry of each list is the primary one used for the main tables in the sanity report.

Two long-text rules are applied identically to every species:

- **`chunk`**: the text is tokenized with the model's own tokenizer and split into passages that fit the model's maximum sequence length (512 tokens for `bge-large-en-v1.5`, 384 for `all-mpnet-base-v2`, both including the two special tokens). Consecutive passages share `embedding.chunk_overlap` tokens. Passages are sized evenly, so a text just over the limit becomes two similar halves, not one full passage and a few leftover tokens. Each passage is embedded, the passage embeddings are averaged, and the average is L2-normalized.
- **`truncate`**: only the first passage-worth of tokens is embedded, then L2-normalized.

`embedding.max_tokens` lowers the per-passage limit for every model; `null` uses each model's own limit. Passages are passed to the model as token ids, so no passage can be silently truncated by the library. The device is picked automatically (`cuda`, then `mps`, then `cpu`) unless `embedding.device` names one.

Each output is cached as `{model_slug}__{rule}__{mask_level}.npy` (float32, one row per species, stubs included) with a matching `.species.txt` giving the row order and a `.meta.json` recording the model, rule, mask level, limits, per-species token and passage counts, and the corpus hash. A cached file is reused unless `--force` is given or its recorded settings no longer match the config or corpus. `--model`, `--rule` and `--mask` restrict a run to one combination.

### Distance matrices (stage 7)

- **Text distance**: 1 − cosine similarity between the normalized embeddings.
- **Phylogenetic distance**: patristic distance from the pruned tree with dendropy, in millions of years of branch length. On this ultrametric tree that is twice the divergence time of the two species.

`species_order.txt` is the single canonical order (species IDs sorted alphabetically, one per line). Every matrix is reindexed to it before saving. `distances.save_matrix` refuses a matrix that is not square, symmetric, zero on the diagonal, finite, and non-negative, and it records a hash of the species order in the matrix's `.meta.json`. `distances.load_matrix` raises `SpeciesOrderError` if a matrix's shape or recorded hash does not match `species_order.txt`, so a stale or misaligned matrix cannot be loaded quietly. Load matrices through `load_matrix`, not `np.load`.

### Sanity report (stage 8)

`reports/week3_sanity.md` covers corpus statistics and the stub list, what each mask level removes, nearest neighbours in text and phylogenetic space for the species in `week3_sanity.showcase_species`, within-family against between-family text distance, a comparison of the three mask levels on that family check, heatmaps in tree tip order, text distance against phylogenetic distance, and description length against near-neighbour count. The nearest-neighbour and family checks are repeated for every model and rule. The report is descriptive only.

The report is regenerated on every run. Hand-written observations live in `reports/week3_notes.md`, which the script appends to the report unchanged; edit that file, not the report.

## Outputs of stage 3

- **`matched.csv`**: one row per species that has both a tree tip and a validated English Wikipedia article. The key columns are `species | tree_tip | wiki_title | match_method | confidence`. The rest are provenance: `wiki_validation`, `wiki_source`, `wiki_qid`, `matched_on`, `tree_binomial`, `mdd_id`, `family`, `genus`, `common_name`, `extinct`, `domestic`.
- **`residual_review.csv`**: anything the pipeline could not decide with confidence. It is meant to be small, and it is currently empty.
- **`excluded_species.csv`**: MDD species with no counterpart in the 2019 tree. Most are species split or described after 2019 (e.g. *Leopardus garleppi*, *Mustela richardsonii*). The rest are extinct (*Lutra nippon*) or domestic (*Canis familiaris*, *Mustela furo*). Each row records the nearest tip and its fuzzy score as evidence that no counterpart was missed.
- **`unused_tree_tips.csv`**: tree tips not used, each with a reason. `junior_synonym_of_mdd_species` means MDD now lumps that tip into a species that already matched another tip. `no_mdd_counterpart` covers fossil taxa such as *Smilodon populator*.

### How names are reconciled

Names are first normalized: case and whitespace fixed, underscores turned into spaces, subgenus parentheses dropped, `cf./aff./sp.` qualifiers removed, trinomials cut to binomials, and `incertae sedis` rejected. Matching to tree tips then runs in passes. Each pass only sees tips not yet claimed, and a tip claimed by more than one species goes to review instead of being assigned.

1. **exact**: the MDD name equals a tip binomial.
2. **legacy_name**: MDD's own `MSW3_sciName` / `CMW_sciName` crosswalk (e.g. *Urva javanica* → `Herpestes_javanicus`).
3. **synonym**: MDD synonym records. Each record's original combination is used, with subspecies raised to species rank (*Canis aureus algirensis* → *Canis algirensis*, never *Canis aureus*), plus the current genus with each synonym's root name.
4. **epithet_in_family**: the same epithet after normalizing Latin gender endings (*griseus/grisea*, *edwardsi/edwardsii*), within the same family or genus.
5. **fuzzy**: rapidfuzz ratio within the same family. At or above `fuzzy_accept_threshold`, with a clear gap to the runner-up, the match is accepted. Between `fuzzy_review_threshold` and the accept threshold it goes to review. This catches the tree's misspelled `Otaria_bryonia`.
6. **wiki_bridge**: a leftover MDD species and a leftover tree tip of the same family whose names resolve to the same Wikipedia article (*Canis lupaster* ↔ tree tip `Canis_anthus`, both → "African wolf").

Wikipedia titles are found through the MediaWiki API with redirects followed. The candidates tried, in order, are: the scientific name, the tree name, the legacy names, the common names, then the synonyms. A final fallback searches Wikidata for items whose taxon name (P225) matches. Every candidate article is **validated** against its Wikidata item. The grades, best first:

- `exact`: P225 equals the MDD name.
- `subspecies_article`: P225 is a trinomial whose last epithet gives the MDD name.
- `synonym`: P225 is one of the species' known synonyms.
- `genus_transfer`: a species-rank item with the same epithet in another genus.
- `name_redirect`: the scientific name redirects straight to a species-rank article.

The last two grades are only allowed when the query was the scientific, tree, or legacy name. When two species land on the same article, the higher grade keeps it and the other species goes on to its next candidate. This stops newly split species from taking their parent species' article.

Leftover review items are then settled automatically where the evidence decides them. If a review species' candidate tip resolves to a different article than the species itself, the candidate is rejected.

## Wikipedia preprocessing rule (stage 5)

The same rule is applied to every species:

1. Fetch the page's plain text with the TextExtracts API (`explaintext`, `exsectionformat=wiki`). This removes infoboxes, tables, references, and markup on the server side. The revision id and timestamp are stored for reproducibility.
2. Keep the **lead** (text before the first heading) plus the **first section whose heading matches** `wikipedia.description_headings` in config (Description, Physical description, Characteristics, ...), together with its subsections.
3. Clean both parts with `wiki.clean_text`: Unicode NFC, any leftover citation markers (`[1]`, `[citation needed]`), templates, HTML, and wiki-link markup removed, empty parentheses dropped, whitespace normalized.

`wiki_texts.csv` records, for each species, whether a description section was found and which heading it came from. Species without one keep only the lead, and the `has_description` flag makes that visible in later analysis.

## API politeness

There is one requests session with a descriptive User-Agent (`wikipedia.user_agent` in config). Requests are throttled to at most one per `min_interval_seconds`. Titles and Wikidata entities are batched 50 per request, and `maxlag=5` is sent. HTTP 429/5xx responses and `maxlag` errors are retried with exponential backoff, and `Retry-After` is honoured. Every response is cached in `data/interim/api_cache/`, so a live pull happens only once.

## Layout

```
config.yaml                 target clade, paths, sources, thresholds, wiki rule
src/lineage_charisma/
  config.py                 config loading and path resolution
  io_utils.py               caching, atomic writes, provenance sidecars, downloads
  taxonomy.py               MDD loading, name normalization, reconciliation engine
  phylo.py                  tree loading, tip parsing, pruning, cophenetic distances
  wiki.py                   API client, title resolution/validation, text preprocessing
  corpus.py                 corpus table, token counts, stub flagging
  masking.py                name-masking vocabulary from MDD, the masker, leak checks
  embed.py                  chunking, pooling, normalization, model wrapper, embedding cache
  distances.py              text and phylogenetic matrices, species-order alignment, validation
  sanity.py                 nearest neighbours, family and length summaries for the report
scripts/01..08_*.py         one script per pipeline stage
tests/                      unit tests (no network or model needed); tests/fixtures holds saved API responses
reports/                    week3_sanity.md (generated), week3_notes.md (hand-written), figures/
data/                       gitignored: raw/, interim/, processed/
```
