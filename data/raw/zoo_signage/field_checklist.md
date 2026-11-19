# Bronx Zoo signage: field checklist (visit 2026-11-28)

This is for secondary analysis C (Amendment 4, section C, and its Change 1, with the common-name rule of Amendment 5, section D). The rules there are fixed, and this checklist only puts them into practice. Do not change a rule here; if something on the day does not fit, write it in `notes`.

After the visit, the whole pipeline is one command: `make zoo`, then `make final` to fill the index rows. See "After the visit" at the end.

## Before you go

- [ ] Phone charged; a power bank; enough free storage for about 300 photos.
- [ ] Camera settings: no filters, HDR on, and location tagging on or off (either is fine; photos are never committed).
- [ ] Have this checklist, and a notes app or paper, to record each building's name as posted at its entrance.
- [ ] Target: **every sign for a mammal species** you can reach. The H1 comparison on the overlap needs at least **10 signed species in the frozen set**, so aim for at least 20 mammal species to leave a margin for names that do not resolve.

## What to photograph

For every **mammal species sign** (skip birds, reptiles, amphibians, fish and invertebrates; only mammals are analysed):

- [ ] **The whole sign, straight on**, filling the frame, with no glare. Step sideways or shade the sign with your hand if needed.
- [ ] **A close-up of the body text** if it is small; take several close-ups for a long sign.
- [ ] **The printed scientific name**, legible. If the sign prints none, make sure the **printed common name** is legible instead: it is then the only name the pipeline can match (see "Special cases").
- [ ] Signs whose text sits inside a graphic, map or illustration: photograph them and note "text in graphic".
- [ ] Digital or rotating screens: photograph every page that names the species.
- [ ] **The building or exhibit entrance sign**, once per building, so the `exhibit` column matches what is posted.

Do not photograph people. Retake any photo with a face in it, or crop it before saving.

## How to name the photos

Rename them the same evening, while you still remember which is which:

```
20261128_<exhibit>_s<NNN>[<part>].jpg
```

| part | rule | example |
|---|---|---|
| `20261128` | the visit date | |
| `<exhibit>` | the exhibit name as posted, lowercase, letters and digits only, no spaces | `jungleworld` for "JungleWorld" |
| `s<NNN>` | the sign number, counting up through the day from `s001`; one number per physical sign | `s014` |
| `<part>` | for several photos of one sign: `a`, `b`, `c`, … (`a` is the whole sign) | `s014a`, `s014b` |

Examples: `20261128_jungleworld_s014a.jpg`, `20261128_jungleworld_s014b.jpg`.

Put every photo in `data/raw/zoo_signage/photos/`. That folder is under `data/`, which git ignores, so photos never go into the repository.

## How to fill in the template

1. Copy `data/raw/zoo_signage/template.csv` to `data/raw/zoo_signage/signage.csv`, keeping the header row exactly.
2. Add **one row per species per sign**. The columns:

| column | what to type | example |
|---|---|---|
| `species_common` | the common name, **as printed**, with the sign's own spelling, hyphens and accents. It is matched against MDD when no scientific name is printed. | `Malayan Tapir` |
| `species_scientific` | the scientific name, **as printed**, including outdated names, trinomials or abbreviations; do not correct it from memory or a book | `Tapirus indicus` |
| `exhibit` | the exhibit or building name as posted at its entrance | `JungleWorld` |
| `photo_filename` | the file name of photo `a`; list several with `;` | `20261128_jungleworld_s014a.jpg;20261128_jungleworld_s014b.jpg` |
| `transcription` | the sign's **body text, verbatim**: same words, spelling and punctuation. Line breaks become single spaces. **Leave out** headings, the common-name line, the scientific-name line, and map or range labels. Keep every other sentence, including captions that are full sentences. | `Tapirs use their flexible snouts to grab leaves …` |
| `notes` | anything unusual, using the fixed phrases below where they fit | `multi-species sign` |

Fixed phrases for `notes` (several can be joined with `; `):

- `partly covered`
- `multi-species sign`
- `text in graphic`
- `digital screen`
- `no scientific name printed`
- `illegible in photo`
- `name printed differently from the common usage`

3. Save as **UTF-8 CSV**. Excel's "CSV UTF-8" and Google Sheets' "Download → CSV" both work; any field containing a comma is quoted automatically. Do not add, remove or rename columns: the pipeline stops if a template column is missing.

## Special cases (all decided in Amendment 4; the checklist only restates them)

- **One sign for several species.** One row per species, each with the **same** `transcription` and the same photo names. Note `multi-species sign`.
- **Several signs for one species.** One row per sign. The pipeline joins the species' distinct transcriptions in row order.
- **No scientific name printed.** Type the common name exactly as printed, leave `species_scientific` **blank**, and note `no scientific name printed`. Do not fill in a scientific name from elsewhere.
  - The pipeline matches the common name to MDD's main and other common names. The match is exact after lowercasing, turning hyphens into spaces and dropping other punctuation.
  - A name shared by several species (for example "Red Panda", which MDD splits into two species) is excluded and counted, and so is a name MDD does not list. Accents are not ignored, so type them as printed.
- **A trinomial** (for example *Panthera tigris jacksoni*): type it as printed. The pipeline reads it as the binomial.
- **Text you cannot read in any photo.** Type what you can, and mark the gap as `[illegible]`. Note `illegible in photo`.
- **Non-mammals.** Do not enter them.

## After the visit

1. Transcribe within a few days, checking each row against its photo.
2. Run the pipeline, from the repository root:

   ```
   make zoo
   ```

   This reconciles the printed names to MDD, masks with the full-scale vocabulary, embeds with the primary model, and writes the following to `reports/tables/secondary/`:

   | file | contents |
   |---|---|
   | `zoo_resolution.csv` | names and match outcomes only; no sign text |
   | `zoo_similarity.csv` | cross-register similarities |
   | `zoo_summary.json` | summary |
   | `zoo_h1.csv` | written only if at least 10 signed species are in the frozen set |

   It also reports how many species were matched by scientific name and how many by common name (in `zoo_summary.json`). The common-name route runs only once Amendment 5 is committed. Everything it needs is already on this machine: the full-scale matrices and embeddings, the ADW texts and the cached embedding model. It sends no paid request.
3. Read the `[zoo] not resolved:` lines it prints, which give each row's route. Fix a row only if the **typed** name differs from the **photo**, which is a transcription error. Names that do not resolve as printed stay excluded, including ambiguous common names, as the amendments require.
4. Run `make final` to refresh the results index and number check.
5. Run `make test`, then commit the tables under `reports/tables/secondary/`; `signage.csv` and the photos stay under `data/`, which is not committed.
