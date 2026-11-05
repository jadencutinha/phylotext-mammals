# Ecological traits

Written by `scripts/10_build_ecology.py`; do not edit by hand. The trait list is fixed in `docs/preregistration.md` (section 7 and amendment A1.1) and in `TRAITS` in `src/lineage_charisma/traits.py`.

## Sources

- **EltonTraits 1.0** (Wilman et al. 2014, *Ecology* 95:2027), mammal file `MamFuncDat.txt`: diet, foraging stratum, activity time.
- **COMBINE** (Soria et al. 2021, *Ecology* 102:e03344), `trait_data_reported.csv`: habitat breadth, terrestrial or aquatic, body mass. This is the file of reported values. COMBINE's imputed file is not used.

## Traits

20 columns in 6 groups, for 4611 species. "Weight" is each column's share of the distance when a pair has every trait: equal by group in the primary matrix, equal by column in the robustness variant.

| column | group | type | source | source column | meaning | group weight | column weight | missing | missing rate |
|---|---|---|---|---|---|---:|---:|---:|---:|
| `diet_inv` | diet | continuous | elton | `Diet-Inv` | percentage of diet: invertebrates | 0.017 | 0.050 | 420 | 9.1% |
| `diet_vend` | diet | continuous | elton | `Diet-Vend` | percentage of diet: mammals and birds | 0.017 | 0.050 | 420 | 9.1% |
| `diet_vect` | diet | continuous | elton | `Diet-Vect` | percentage of diet: reptiles and amphibians | 0.017 | 0.050 | 420 | 9.1% |
| `diet_vfish` | diet | continuous | elton | `Diet-Vfish` | percentage of diet: fish | 0.017 | 0.050 | 420 | 9.1% |
| `diet_vunk` | diet | continuous | elton | `Diet-Vunk` | percentage of diet: vertebrates, type unknown | 0.017 | 0.050 | 420 | 9.1% |
| `diet_scav` | diet | continuous | elton | `Diet-Scav` | percentage of diet: carrion | 0.017 | 0.050 | 420 | 9.1% |
| `diet_fruit` | diet | continuous | elton | `Diet-Fruit` | percentage of diet: fruit | 0.017 | 0.050 | 420 | 9.1% |
| `diet_nect` | diet | continuous | elton | `Diet-Nect` | percentage of diet: nectar and pollen | 0.017 | 0.050 | 420 | 9.1% |
| `diet_seed` | diet | continuous | elton | `Diet-Seed` | percentage of diet: seeds | 0.017 | 0.050 | 420 | 9.1% |
| `diet_planto` | diet | continuous | elton | `Diet-PlantO` | percentage of diet: other plant material | 0.017 | 0.050 | 420 | 9.1% |
| `foraging_stratum` | foraging stratum | categorical | elton | `ForStrat-Value` | M marine, G ground, S scansorial, Ar arboreal, A aerial | 0.167 | 0.050 | 420 | 9.1% |
| `activity_nocturnal` | activity time | binary | elton | `Activity-Nocturnal` | active at night | 0.056 | 0.050 | 420 | 9.1% |
| `activity_crepuscular` | activity time | binary | elton | `Activity-Crepuscular` | active at dawn and dusk | 0.056 | 0.050 | 420 | 9.1% |
| `activity_diurnal` | activity time | binary | elton | `Activity-Diurnal` | active by day | 0.056 | 0.050 | 420 | 9.1% |
| `habitat_breadth` | habitat breadth | continuous | combine | `habitat_breadth_n` | number of IUCN habitat types used | 0.167 | 0.050 | 241 | 5.2% |
| `terrestrial` | terrestrial or aquatic | binary | combine | `terrestrial_non-volant` | lives on land | 0.042 | 0.050 | 92 | 2.0% |
| `marine` | terrestrial or aquatic | binary | combine | `marine` | lives in the sea | 0.042 | 0.050 | 92 | 2.0% |
| `freshwater` | terrestrial or aquatic | binary | combine | `freshwater` | lives in fresh water | 0.042 | 0.050 | 92 | 2.0% |
| `volant` | terrestrial or aquatic | binary | combine | `terrestrial_volant` | flies | 0.042 | 0.050 | 92 | 2.0% |
| `log_body_mass` | body mass | continuous | combine | `adult_mass_g` | log10 of adult body mass in grams | 0.167 | 0.050 | 164 | 3.6% |

## How the distance is computed

- Continuous traits contribute |difference| / range, where the range is taken over the species with a value. Categorical and binary traits contribute 0 for a match and 1 for a mismatch.
- A pair's distance is the weighted mean of those contributions over the traits both species have (Gower's pairwise-available rule). Nothing is imputed.
- `ecology.npy` weights the six groups equally (primary). `ecology_unweighted.npy` weights every column equally (robustness variant).
- Body mass is log10 of `adult_mass_g`. The three activity flags are not exclusive: a species can be both nocturnal and crepuscular.
- `volant` is COMBINE's `terrestrial_volant` (full scale only, docs/preregistration_amendment_fullscale.md section 4). `ecology_no_volant.npy` is the group-weighted matrix without it, a declared robustness variant.

## Coverage

4519 of 4611 species have an ecological distance to each other and are in the matrices (`ecology_species.txt`). 4025 species have every trait.
- `ecology.npy`: over all 10,628,355 pairs of the 4611 species, a pair shares on average 87.8% of the total trait weight; 23.8% of pairs are missing at least one trait and 300433 pairs share none.
- `ecology_unweighted.npy`: over all 10,628,355 pairs of the 4611 species, a pair shares on average 86.2% of the total trait weight; 23.8% of pairs are missing at least one trait and 300433 pairs share none.
