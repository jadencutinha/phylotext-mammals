## 8. Observations

These are hand-written notes on the run of 2026-09-30 (284 Carnivora species after excluding one stub; primary combination bge-large / chunk / taxonomy). `scripts/08_sanity_report.py` appends this file, `reports/week3_notes.md`, to the report unchanged, so check the numbers against the tables above after any rebuild.

**Masking scientific names changes almost nothing.** With taxonomy masking no text names its own genus, species epithet or family (283, 277 and 59 texts did before). The family check is unchanged: the pooled within-to-between ratio is 0.794 with and without taxonomy masking for bge-large / chunk, and a species' nearest text neighbour is in its own family for 87% of species against 90% unmasked. The masked and unmasked matrices have a rank correlation of 0.91–0.95 and share 3.7–4.0 of each species' top five. The scientific names were not what held families together.

**Masking common names removes most of the family structure.** Under strict masking the ratio rises to 0.89–0.92 in every model and rule, and the nearest text neighbour is in the same family for only 43–46% of species. A random other species shares the family 11% of the time, so some structure survives without any names, but most of what the unmasked matrices showed was carried by words like "mongoose", "seal", "skunk" and "fox".

- The families that lose the most are the ones whose members share one common word: Herpestidae goes from 0.71 to 1.01 (no tighter than its surroundings), Mephitidae from 0.75 to 0.95, Phocidae from 0.65 to 0.85.
- Ursidae, Hyaenidae and Felidae stay the tightest families under strict masking (0.72, 0.72, 0.81). Those are also the families with the longest articles, so this is not clean evidence of content-based signal.
- Eupleridae is the one family that gets slightly tighter (0.85 to 0.81). Its members have different common names (fossa, falanouc, vontsira), so there was no shared word to lose.
- Prionodontidae has two species and moves from 0.29 to 0.61 with taxonomy masking alone, because both texts leaned on the genus name. Ignore it as a family-level result.

**Under strict masking the neighbour lists are dominated by other long articles.** The wolf's top five become tiger, fisher, stoat, caracal and wolverine; the meerkat's become serval, caracal, ocelot, cheetah and fennec fox. Some ecology survives (the sea otter's neighbours are walrus, elephant seal and leopard seal; the polar bear's are Arctic fox, wolf and walrus). The giant panda's nearest neighbour is the red panda with no masking and with taxonomy masking, and drops to fifth under strict.

**Strict masking has side effects of its own.** It replaces a median of 13 terms per text (4% of words) and lowers every distance (mean 0.305 against 0.350 for bge-large / chunk), because every text now shares the same token. Texts with similar mask density are slightly closer (rank correlation 0.17 for bge-large, 0.09 for mpnet, between the difference in density and text distance), which is small next to the changes above. This check was run by hand and is not in the tables.

**Text distance still separates close relatives from everything else, then flattens.** For the primary combination the median text distance rises from 0.22 for pairs under 5 My apart to about 0.33 at 25–30 My, and stays between 0.30 and 0.37 out to 81 My.

**The description-length effect is still there under taxonomy masking.** For bge-large / chunk, pairs where both texts exceed one passage average 0.269 against 0.380 where neither does, and the longest third of texts has a median of 26 near neighbours against about 5 for the rest. Truncation shrinks the gap (0.333 against 0.380) and mpnet does not show it (0.546 against 0.564 for chunk; 0.588 against 0.564 for truncate). Length is tied to family (median 194 tokens for Herpestidae, 1061 for Ursidae), so it can look like phylogenetic signal or hide it.

**The two models agree only moderately.** At the taxonomy level, bge-large and mpnet distances have a rank correlation of 0.49–0.62 and share 2.5–2.7 of each species' top five. Chunk and truncate within one model agree at 0.88–0.90.

**What the masking does not cover.**

- Names of taxa outside the clade (prey species, for example) are left in.
- Common names that are also current genus or epithet names (lynx, caracal, puma, fossa, serval, binturong) are masked at the taxonomy level, because the instruction was to match scientific names case-insensitively. Historical genus names spelled like common names (Hyena, Coati) are left to the strict level.
- "canine" is masked unless followed by "teeth", "tooth" or "tip"; "canines" is never masked.
- Common-name modifiers stay in under strict ("red", "giant", "Arctic", "sea"), so "sea [TAXON]" still marks an otariid.

**Open questions for Week 4.**

- Taxonomy is the primary mask level as you asked, but it behaves like no masking. The strict level is the one that tests whether description content, as opposed to names, tracks the tree.
- Description length belongs in the Week 4 design as a covariate alongside the attention measures, most of all under strict masking.
- The 100-token sensitivity matrices exist in `data/processed/matrices_min100/` (277 species) and have not been analysed.
