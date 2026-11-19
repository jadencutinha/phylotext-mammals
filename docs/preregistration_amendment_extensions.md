# Amendment 5 to the pre-registration: extension analyses (Week 8)

Written 2026-11-19. The main results were locked on 2026-11-11 (`reports/final_results_index.md`), and these analyses were designed **after** that, with every pilot and full-scale H1/H2 result already known. The specification below is fixed before any of the analyses has been run: no commercial embedding request has been sent, no TF-IDF vector has been built, no section-level text has been embedded or tested, and no zoo signage has been collected. The commit that adds this file is the timestamp of record. Later changes go in the "Changes" section at the end, with a date and a reason, and always before the results they affect exist.

## 1. Status of these analyses

Analyses A–C are **secondary, designed after the main results were locked**. Whatever they show, they **cannot change any H1 or H2 conclusion** of the pilot (`reports/week4_report.md`, `reports/week5_report.md`) or of the full-scale analysis (`reports/week6_report.md`). Those conclusions rest on the registered primary specification alone. These analyses are reported in their own section, labelled "secondary (Amendment 5, designed after the lock)" in every table, figure and sentence, and never use the words "supports H1" or "supports H2". Nothing here is corrected for multiple testing, and no result here is read as confirmatory.

- **A. Commercial embeddings** asks whether the registered outcomes depend on the open-source encoder.
- **B. TF-IDF baseline** asks whether a text representation with no neural model gives the same outcomes.
- **C. Section-level analysis** is descriptive. It asks which parts of an article carry the text–phylogeny relation.

Section **D** is not a new analysis. It amends one name-matching rule of Amendment 4's zoo signage analysis (secondary C there, descriptive only) before any signage exists; the visit is on 2026-11-28.

## 2. What had been seen before this was written

- Every registered result. On the primary specification (`BAAI/bge-large-en-v1.5` / chunk / strict, Spearman, MCC tree):
  - **Pilot** (284 Carnivora): H1 r = 0.033, p = 0.0285. H2 (283 species): raw r = 0.035, partial r (ecology + length) = 0.020, p = 0.127. Rule (a) fails; rule (b) does not hold (partial/raw = 0.575).
  - **Full scale** (4,611 mammals): H1 r = −0.028, p = 1.000. H2 (4,519 species): raw r = −0.029, partial r = −0.034, p = 1.000. Rule (a) fails; rule (b) is undefined.
  - The registered variants, including `sentence-transformers/all-mpnet-base-v2` (full-scale H1 r = 0.0006) and the taxonomy mask (full-scale H1 r = 0.186).
- Every exploratory, post-hoc and secondary result in `reports/final_results_index.md`, including post-hoc D (r_LM = 0.0166, two-sided p = 0.0245) and secondary A of Amendment 4 (ADW against Wikipedia).
- Analysis B of Amendment 4 (LLM charisma) had not run: no hand label was committed and no model request had been sent.
- **For section C only:** the frequency of every section heading in the cached Wikipedia extracts of the 4,611 frozen species. These are heading strings and counts, computed on 2026-11-19 to write the mapping in section C. No text distance, embedding or statistic was computed from section text.
- **For section D only:** MDD v2.5's `common_name` and `other_common_names` columns. After the normalization in section D, the 6,904 species carry 15,027 distinct names, 302 of which belong to more than one species. Nothing about the zoo was seen: no sign had been photographed or transcribed.
- **For the cost estimate in section A only:** the `cl100k_base` token count of the strict-masked texts. Pilot: 176,300 tokens over 284 texts, longest 2,768. Full scale: 1,348,009 tokens over 4,611 texts, longest 3,446.

## 3. Shared settings

Unless an item says otherwise, every test uses the registered settings:

- Spearman correlation on upper triangles, average ranks.
- **9,999 permutations**, **seed 20261004**, with a fresh generator per test.
- One-sided p (greater) for correlations predicted to be positive, and alpha 0.05.
- Phylogenetic distance: the MCC patristic matrix.
- H2 controls: the existing ecology and length matrices. The length distance is the absolute difference in log unmasked token count. It measures the text, not the model, so it is the same matrix for every representation.
- Statistics are computed with the condensed implementation (`src/lineage_charisma/condensed.py`), which gives the same results as `stats.py` with the same seed (`tests/test_condensed.py`).

**Species sets, fixed:**

| scale | H1 species | H2 species | file |
|---|---|---|---|
| pilot | 284 Carnivora | 283 | `data/processed/matrices/species_order.txt`, `ecology_species.txt` |
| full scale | 4,611 mammals (frozen set, SHA-256 `4730f407…763f0e`) | 4,519 | `data/processed/mammalia/matrices/species_order.txt`, `ecology_species.txt` |

**Texts.** Analyses A and B use the strict-masked texts already in the corpora: `text_strict` in `data/processed/corpus.parquet` for the pilot (pilot mask vocabulary) and in `data/processed/mammalia/corpus.parquet` for the full scale (full-scale vocabulary). They are not re-fetched or re-masked. Analysis C re-reads the cached Wikipedia extracts (`data/interim/wiki/raw/`, the same revisions the corpora were built from) and masks each section group with the same masker and vocabulary as its scale.

**Registered outcomes compared against.** For each representation and scale:

- **H1 outcome:** "positive" if r > 0 and p < 0.05, otherwise "not positive".
- **H2 rule (a):** survives if partial r > 0 and p < 0.05.
- **H2 rule (b):** confounded if partial/raw < 0.5; undefined if raw r ≤ 0.

A representation **agrees** with the primary specification at a scale when it has the same H1 outcome and the same rule (a) outcome. The primary outcomes are:

| scale | H1 | rule (a) | rule (b) |
|---|---|---|---|
| pilot | positive | fails | not confounded (0.575) |
| full scale | not positive | fails | undefined |

## A. Commercial embeddings

**Question.** Are the registered H1 and H2 outcomes specific to `bge-large-en-v1.5`, or do two commercial embedding models give the same outcomes on the same texts?

**Models** (fixed; any change goes in a dated Change entry before the first request):

| provider | model | input limit | output | other settings |
|---|---|---|---|---|
| OpenAI | `text-embedding-3-large` | 8,191 tokens (`cl100k_base`) | 3,072 dimensions (default; the `dimensions` parameter is not set), `encoding_format: "float"` | none |
| Voyage AI | `voyage-3-large` | 32,000 tokens | model default dimension (1,024), `output_dtype: "float"` | `input_type: null` (no retrieval prompt, since documents are compared with documents); `truncation: false`, so an over-limit input fails instead of being cut silently |

**Chunking rule, adapted to each model's limit.** The registered chunk rule is kept, with the model's own limit in place of 512 tokens:

- Each text is tokenized with the model's tokenizer. For OpenAI this is `tiktoken` `cl100k_base`. For Voyage it is the tokenizer Voyage publishes for the model; if that cannot be loaded, `cl100k_base` is used with a limit of 16,000, half of 32,000, as a safety margin.
- A text longer than the limit is split into evenly sized passages with a 32-token overlap, the same `chunk_spans` function as the primary specification.
- Each passage is embedded. The species vector is the mean of its passage vectors, L2-normalized.
- At the counts in section 2 (longest text 3,446 `cl100k_base` tokens), no text is expected to need more than one passage for either model. The number of multi-passage texts is reported.

**Requests.**

- Every request and response is cached on disk under `data/interim/api_embeddings/<provider>/`, keyed by provider, model and the SHA-256 of the passage. A rerun reads the cache and never re-requests a cached passage.
- Batches stay under each provider's per-request limits.
- Failed requests (HTTP 429 or 5xx) are retried with exponential backoff, up to 5 times. A passage still without an embedding stops the run; no species is silently dropped.
- API keys are read from `OPENAI_API_KEY` and `VOYAGE_API_KEY`, are never printed, and are never written to a file.
- The response's `model` field and token usage are recorded with each cached response.

**Tests** (both models, both scales; 9,999 permutations, seed 20261004):

- **A-H1:** Mantel r of text distance (1 − cosine) against phylogenetic distance, with phylogeny permuted, on the H1 species.
- **A-H2:** raw r and partial Mantel r of text distance on phylogeny, controlling for ecology + length jointly, with phylogeny permuted, on the H2 species. Rules (a) and (b) are computed as for the primary analysis.
- **Descriptive:** the Spearman correlation between each pair of text-distance matrices (bge-large, mpnet, OpenAI, Voyage) on the upper triangle, without p-values.

**What agreement or disagreement would mean** (stated before any result):

1. **Both models agree with the primary outcomes at both scales.** The registered outcomes are not specific to the open-source encoder. They are more likely to reflect the texts than a property of one model.
2. **A model gives a positive full-scale H1 (r > 0, p < 0.05).** The full-scale negative r depends on the encoder. That would fit the post-hoc D account, in which length composition drives the sign under bge-large. The registered full-scale result stands for the registered specification. Any wording that the strict-masked texts carry no across-mammal phylogenetic signal in general would have to be qualified as encoder-dependent.
3. **A model gives a pilot H1 that is not positive.** The pilot H1 result (p = 0.0285, close to alpha) is encoder-dependent and should be described as fragile. The registered pilot decision is not changed.
4. **A model gives rule (a) "survives" at either scale.** A stronger encoder detects lineage signal beyond the measured controls. The registered H2 conclusion ("not supported") is unchanged, and the result is reported as a secondary finding that motivates future work.
5. **The two commercial models disagree with each other on any outcome.** The choice of encoder is a main source of variation in this design. This is reported as a limitation, and no model is preferred after the fact.
6. **Mixed pattern** (for example, agreement at one scale only). Each scale is reported as it falls, and no overall verdict is drawn.

Differences in r between models are reported as numbers, without inference.

**Cost and stop rule.**

- About 1.52 million input tokens per model across both scales (section 2).
- Estimated cost at list prices (OpenAI US$0.13 and Voyage US$0.18 per million tokens): about US$0.20 for OpenAI and US$0.28 for Voyage, about US$0.47 in all. Under pessimistic assumptions (Voyage's tokenizer counting 25% more tokens, and every request sent twice) it is about US$1.08.
- Prices are confirmed on the providers' pricing pages before the first request.
- Spending stops and is reported if the measured cost after the pilot projects above US$10 for the full run.

**Limitations.**

- Commercial models are versioned by the provider and can change or be withdrawn. The cached vectors are the record, and a rerun from cache is exact while a rerun against the API may not be.
- The training data of these models is undisclosed and very probably includes Wikipedia. That is true of the open-source models as well.
- Only the texts are sent to the providers: public Wikipedia text with names masked. No personal data is sent.

## B. Non-embedding baseline: TF-IDF

**Question.** Does a bag-of-words representation, with no neural model, give the same H1 and H2 outcomes as the embeddings?

**Vectors** (scikit-learn 1.9.1 `TfidfVectorizer`, as pinned in `requirements.lock`), fitted separately at each scale on that scale's H1 species and texts only:

| setting | value |
|---|---|
| input | `text_strict` (the same masked text the primary embedding sees) |
| mask token | every `[TAXON]` is replaced by a space before vectorizing, so the mask token is not a term |
| lowercase | yes |
| token pattern | `(?u)\b[^\W\d_]{2,}\b`: runs of two or more letters; digits, punctuation and single letters are not tokens |
| stopwords | scikit-learn's built-in English list (`ENGLISH_STOP_WORDS`, 318 words) |
| n-grams | unigrams only |
| stemming / lemmatizing | none |
| minimum document frequency | 5 documents (`min_df=5`) |
| maximum document frequency | none (`max_df=1.0`) |
| term weighting | `sublinear_tf=True` (1 + log tf), `smooth_idf=True`, `use_idf=True` |
| normalization | L2 |
| distance | 1 − cosine similarity of the L2-normalized vectors |

A species whose vector is all zeros after these rules is listed and dropped from that scale's tests. None is expected, because every analysed text has at least 50 tokens. Vocabulary size, the number of non-zero entries per document (median) and any dropped species are reported.

**Tests.** Pilot and full scale, 9,999 permutations, seed 20261004:

- **B-H1:** H1 Mantel.
- **B-H2:** raw and partial Mantel controlling for ecology + length, with rules (a) and (b).
- **Descriptive:** the Spearman correlation between the TF-IDF and bge-large text-distance matrices.

**Reading** (stated before any result):

- **TF-IDF agrees with the primary outcomes.** The outcomes do not need a neural representation. Whatever relation there is, it is visible in word use alone.
- **TF-IDF gives a positive H1 where the embeddings do not** (for example at full scale). Lexical overlap in non-name vocabulary carries phylogenetic signal that the embedding does not express. One reason could be that TF-IDF ignores the shared `[TAXON]` token and responds differently to text length.
- **TF-IDF gives no positive H1 where the embeddings do** (the pilot). The embedding signal is not just shared vocabulary.

No reading changes H1 or H2.

**Limitations.**

- TF-IDF distances on short texts are noisy, because few terms are shared. Stub-adjacent texts (50–100 tokens) are common at full scale in Rodentia, Chiroptera and Eulipotyphla.
- The document-frequency floor and stopword list are one fixed choice. No grid over them is run.

## C. Section-level analysis (descriptive)

**Question.** Which parts of a Wikipedia article carry the text–phylogeny relation? This is descriptive only. Every p-value is reported for completeness, with **no multiple-testing claim** and no reading of any group as significant.

**Source.** For each species, the cached extract for the same revision the corpus used, split into sections by `wiki.split_sections` (headings `== … ==` at levels 2 to 6). The **lead** (the text before the first heading) is not one of the five groups, because it summarizes every topic at once. It is reported as a **sixth, reference row** (`lead`), under the same text, threshold, embedding and test rules as the groups. Unlike the primary text, which is the lead plus one description section, the five groups together use every section of the article.

**Heading normalization.** Lowercase; strip surrounding whitespace; replace `&` with `and`; British spellings `behaviour` → `behavior` and `colouration` → `coloration`; collapse internal whitespace.

**Groups.** Each section is assigned to exactly one of: `appearance`, `behavior`, `distribution_habitat`, `diet`, `other`, or `excluded` (not used). The assignment follows these rules, in order:

1. **Exact table.** If the normalized heading is in this table, it takes that group.

   | group | headings |
   |---|---|
   | excluded | references, external links, see also, further reading, notes, literature cited, sources, bibliography, footnotes, citations, gallery, notes and references, references and notes, cited literature, works cited |
   | appearance | description, physical description, characteristics, physical characteristics, appearance, morphology, anatomy, anatomy and morphology, anatomy and physiology, size, description and anatomy, identification, coloration, pelage, dentition, skull |
   | behavior | behavior, ecology, biology, behavior and ecology, ecology and behavior, biology and ecology, biology and behavior, behavior and biology, life history, natural history, habits, reproduction, breeding, mating, reproduction and development, reproduction and life cycle, mating and reproduction, social behavior, social structure, communication, vocalizations, echolocation, hibernation, locomotion, senses, predators, predation, lifespan |
   | distribution_habitat | distribution and habitat, distribution, habitat, range and habitat, habitat and distribution, range, population and distribution, habitat and ecology, distribution and ecology, range and status, distribution and status |
   | diet | diet, feeding, food habits, hunting, foraging, diet and foraging, food and foraging, hunting and diet, diet and hunting, feeding ecology, diet and feeding |
   | other | taxonomy, etymology, taxonomy and etymology, taxonomy and evolution, taxonomy and phylogeny, taxonomy and naming, evolution, phylogeny, fossil record, classification, names, subspecies, hybrids, genetics, history, discovery, extinction, conservation, conservation status, status, threats, threats and conservation, status and conservation, population status, relationship with humans, relationships with humans, interaction with humans, interactions with humans, human interaction, human interactions, in culture, in popular culture, cultural significance, in captivity, parasites, physiology, attacks on humans, cultural references |

2. **Keyword rule.** Otherwise the heading is split into words (letters only). If any word is in a keyword list below, the heading takes the group of the keyword that occurs **first** in the heading. The lists are checked in the order shown, so `excluded` wins a tie at the same word.
   - excluded: references, notes, sources, bibliography, footnotes, citations, gallery, literature, links
   - appearance: description, appearance, characteristics, morphology, anatomy, size, coloration, pelage, fur, dentition, skull, measurements, physical, identification
   - behavior: behavior, ecology, biology, habits, social, communication, reproduction, breeding, mating, echolocation, vocalizations, vocalization, hibernation, locomotion, activity, life, lifecycle, lifespan, predators, predation, senses, development
   - distribution_habitat: distribution, habitat, habitats, range, occurrence, geographic
   - diet: diet, feeding, food, foraging, hunting, prey, nutrition
3. **Inheritance.** If neither rule applies, a section at level 3 or deeper takes the group of its nearest enclosing section, and that inherited group may be `excluded`. A level-2 section takes `other`.

The complete heading → group table that these rules produce over the 4,611 frozen species is written to `reports/tables/extensions/section_heading_map.csv` (heading, count, group, and which rule assigned it).

**Group text.**

- Each section's body is cleaned with `wiki.clean_text`. A species' bodies in one group are joined in article order with a blank line between them. The `lead` text is the cleaned lead alone, the same text as the lead part of the corpus text.
- The unmasked token count uses the bge-large tokenizer, special tokens not counted, as in the corpus stub rule.
- A species is in a group's analysis only if that group's text has **at least 50 tokens**, the registered stub threshold.
- Group texts are masked at the **strict** level with the scale's vocabulary (pilot or full scale), the same masker used to build the corpora.
- The pipeline stops if a masked group text still contains its own genus name.

**Embedding.** Primary model and rule only: `BAAI/bge-large-en-v1.5`, chunk rule (512-token passages, 32-token overlap, mean-pooled, L2-normalized). Distance is 1 − cosine.

**Tests** (each of the five groups and the `lead` row, each scale, on that row's species):

- **C-H1:** Mantel r of the group's text distance against phylogenetic distance, with phylogeny permuted, 9,999 permutations, seed 20261004, one-sided p for completeness.
- **Like-for-like reference:** the primary H1 r (bge-large / chunk / strict corpus text) recomputed on the same species, so each group's r has a comparison on identical species. Its p uses the same settings.
- A group is not tested at a scale if it has fewer than **50 species** there. Its count is reported.
- **Reported per group:** number of species, median tokens, r, p, and the reference r.
- `other` is reported too, with a note that it mixes taxonomy (dense in `[TAXON]` tokens), conservation and human relations.
- `lead` is reported in the same table, marked as a reference row. Its texts are part of the corpus text, so its r is not independent of the like-for-like reference.

**Reading.** Descriptive only. A group whose r is above the reference r on the same species carries relatively more phylogenetic signal than the corpus text. Groups, including the `lead` row, are not compared with each other by test, and their species sets differ.

**Limitations.**

- Section structure varies with article quality. Species with many sections are better documented and more viewed, so each group's species set is selected by attention, and differently from group to group.
- Headings are an imperfect guide to content. The keyword rule is a fixed approximation, and compound headings take the first topic named.
- Group texts are shorter than the full text, so their distances are noisier, and short groups lose more species to the 50-token threshold.
- Taxonomy sections are dominated by names, which the mask replaces. What remains is largely the mask token, which could inflate similarity in `other`.

## D. Zoo signage: matching signs with no printed scientific name (amends Amendment 4, section C)

Dated **2026-11-19**, before the Bronx Zoo visit on 2026-11-28 and before any sign was photographed or transcribed.

**What changes.** Under Amendment 4, section C and its Change 1, a sign is reconciled to MDD through its printed scientific name only, and a row with a blank `species_scientific` is skipped. From this amendment on, a row whose `species_scientific` is blank (the sign prints no scientific name) is matched through its printed common name (`species_common`) instead. Rows with a printed scientific name are resolved exactly as before. A printed scientific name that does not resolve is still excluded, with no common-name fallback.

**Common-name rule.**

1. **Names searched:** each MDD species' `common_name` and every entry of its `other_common_names` (entries separated by `|`), from the full-scale MDD table (`data/raw/mdd_species_mammalia.csv`).
2. **Normalization** of the printed name and every MDD name:
   - Lowercase.
   - Replace every hyphen or dash (U+002D, U+2010–U+2015, U+2212) with a space.
   - Delete every other character that is not a letter, a digit or whitespace (apostrophes, periods, commas, parentheses and the like).
   - Collapse runs of whitespace to one space and trim.
   - Accents are letters and are kept, so "Père" and "Pere" differ.
3. **Match:** the normalized printed name must **equal** a normalized MDD name exactly; there is no fuzzy or partial matching.
   - If the names it equals all belong to one species, the row resolves to that species.
   - If they belong to **more than one species**, whether through main or other names, the row is **ambiguous**. It is excluded and counted, and no precedence between main and other names is applied.
   - If it equals no MDD name, the row is excluded and counted.
4. A row with neither a scientific nor a common name is skipped and counted, as before.

**Reported:**
- The number of signed species matched **by scientific name** and **by common name**. A species with rows matched both ways counts as matched by scientific name.
- The number of common-name rows that were ambiguous or had no match, each listed by its printed name.

Everything downstream (one text per species, masking with the full-scale vocabulary, embedding, similarities and the H1 overlap) is unchanged.

**Limitations.**
- Common names are less specific than scientific names. A sign's common name may follow a different authority from MDD (for example, a name MDD has since split between two species), and such rows are lost as no-match or ambiguous rather than guessed.
- Not folding accents is deliberate, to keep the rule exact, and may lose a few rows.

## 4. Outputs

- Tables go to `reports/tables/extensions/` (pilot) and `reports/tables/extensions/mammalia/` (full scale). Each has a `.meta.json` recording "secondary (Amendment 5, designed after the 2026-11-11 lock)", the settings and the permutations.
- Embeddings and TF-IDF matrices go to `data/processed/…/extensions/` (not committed).
- The zoo pipeline (`scripts/28_zoo_signage.py`) adds a `route` column (scientific name or common name) to `zoo_resolution.csv`, adds the section D counts to `zoo_summary.json`, and refuses to use the common-name route unless this file is committed.
- The scripts (`36_api_embeddings.py`, `37_tfidf_baseline.py`, `38_section_groups.py`) refuse to compute any statistic, and script 36 refuses to send any request, unless this file is committed with no uncommitted changes.

## 5. Changes

(none yet)
