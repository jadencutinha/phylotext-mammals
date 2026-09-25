# lineage-or-charisma

Does text-embedding similarity between mammal species descriptions recover the evolutionary tree, or does it mainly reflect human attention and charisma? This repo compares a text-distance matrix (from Wikipedia description embeddings) with a phylogenetic-distance matrix (from the Upham et al. 2019 mammal tree) using Mantel and partial Mantel tests.

**Current state: Week 2, the data foundation and taxonomy reconciliation.** The target clade is set in `config.yaml` (currently order Carnivora). No code names a clade directly, so moving to all mammals only means changing `target_clade`.

## Setup

```bash
make setup          # creates .venv and installs the package + pytest
make test           # runs the unit tests (no network needed)
make all            # runs stages 1-5; each stage skips work whose output is cached
make all FORCE=1    # rebuilds every stage
```

Python 3.11+ is required. Dependencies: requests, dendropy, rapidfuzz, pandas, pyyaml (plus pytest for dev).

## Pipeline

| Stage | Command | Reads | Writes |
|---|---|---|---|
| 1 | `make mdd` → `scripts/01_fetch_mdd.py` | MDD zip from the MDD GitHub site | `data/raw/mdd/` (full release), `data/raw/mdd_species_<clade>.csv`, `data/raw/mdd_synonyms_<clade>.csv` |
| 2 | `make tree` → `scripts/02_fetch_tree.py` | Upham 2019 MCC tree | `data/raw/tree/*.tre`, `data/interim/tree_tips_all.csv`, `data/interim/tree_tips_<clade>.csv` |
| 3 | `make reconcile` → `scripts/03_reconcile_names.py` | stage 1 + 2 outputs, Wikipedia + Wikidata APIs | `data/interim/matched.csv`, `residual_review.csv`, `excluded_species.csv`, `unused_tree_tips.csv` |
| 4 | `make prune` → `scripts/04_prune_tree.py` | `matched.csv`, tree | `data/interim/tree_pruned.nwk`, `tree_pruned_mdd_names.nwk`, `cophenetic_distance.csv`, `sanity_check.txt` |
| 5 | `make wiki` → `scripts/05_fetch_wikipedia.py` | `matched.csv`, Wikipedia API | `data/interim/wiki/raw/*.json`, `data/interim/wiki/clean/*.txt`, `data/interim/wiki_texts.csv` |

Every script takes `--force` and `--config PATH`. Stages 3 and 5 also take `--offline`, which answers every API request from the cache and fails if a response is missing. Every CSV written gets a `.meta.json` sidecar with its sha256, row count, creation time, and source (for example, the MDD version).

Optional: `make posterior` downloads the full 10,000-tree posterior from VertLife (`Completed_5911sp_topoCons_NDexp.zip`, about 1.2 GB) for later weeks. The tip labels are the same across the posterior, so this week's reconciliation applies to every tree in it.

### Data sources (checked 2026-09-24)

- **MDD v2.5** (released 2026-07-28): `https://raw.githubusercontent.com/mammaldiversity/mammaldiversity.github.io/master/assets/data/MDD.zip`. Stage 1 picks the highest-versioned `MDD_v*_*species.csv` and `Species_Syn_Current_v*.csv` inside the zip, so a new MDD release is picked up without code changes.
- **Upham et al. 2019 tree**: the MCC tree of the *completed, topology-constrained, node-dated (NDexp)* posterior, from the authors' repo `n8upham/MamPhy_v1` (`MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre`). The matching 10k posterior is on `data.vertlife.org/mammaltree/`.

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
scripts/01..05_*.py         one script per pipeline stage
tests/                      unit tests; tests/fixtures holds saved API responses
data/                       gitignored: raw/, interim/, processed/
```
