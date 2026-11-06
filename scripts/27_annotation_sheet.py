from __future__ import annotations

import pandas as pd

from lineage_charisma.charisma import DIMENSIONS, RUBRIC, VALIDATION_SIZE, validation_sample
from lineage_charisma.config import base_parser, load_config
from lineage_charisma.corpus import read_corpus
from lineage_charisma.embed import read_species_order
from lineage_charisma.io_utils import atomic_write_text, write_csv

# docs/preregistration_amendment_secondary.md, section B, validation item 1


def main() -> None:
    args = base_parser("Write the hand-labelling sheet for the charisma rubric (50 seeded pilot species, strict-masked texts).").parse_args()
    cfg = load_config(args.config, None)   # analysis B is on the pilot species
    seed = int(cfg["stats"]["seed"])
    out_csv, out_md = cfg.interim / "annotation_validation.csv", cfg.interim / "annotation_instructions.md"
    out_key = cfg.interim / "annotation_validation_key.csv"   # item -> species; kept apart so the labeller sees only masked text
    if out_csv.exists() and not args.force:
        filled = pd.read_csv(out_csv)[list(DIMENSIONS)].notna().any().any()
        raise SystemExit(f"{out_csv.relative_to(cfg.root)} exists{' and has scores in it' if filled else ''}; not overwritten (use --force only to start over)")

    order = read_species_order(cfg.species_order_path())
    sample = validation_sample(order, seed)
    corpus = read_corpus(cfg.corpus_path()).set_index("species_id")
    sheet = pd.DataFrame({"item": range(1, VALIDATION_SIZE + 1), "species_id": sample})
    for d in DIMENSIONS:
        sheet[d] = pd.array([pd.NA] * len(sheet), dtype="Int64")
    sheet["notes"] = ""
    write_csv(sheet[["item", "species_id"]], out_key, seed=seed, note="do not open while labelling")
    sheet = sheet.drop(columns="species_id")
    write_csv(sheet, out_csv, seed=seed, sample="numpy.random.default_rng(seed).choice(284, size=50, replace=False) over the pilot species_order.txt",
              text_level="strict", plan="docs/preregistration_amendment_secondary.md, section B")

    lines = [
        "# Charisma rubric: hand-labelling instructions",
        "",
        "Score each of the 50 descriptions below from 0 to 3 on the four dimensions, using **only the wording of the text**, "
        "not what you know about the animal. Names have been replaced by `[TAXON]`; do not try to guess them. "
        "Enter whole numbers 0 to 3 in the matching columns of `annotation_validation.csv` (row `item` = the number below). "
        "Leave `annotation_validation_key.csv` closed until you have finished. "
        "Do this before looking at any model output; model scores are not produced until these labels are committed.",
        "",
        "## Rubric (the same definitions the model receives)",
        "",
        *RUBRIC.split("\n"),
        "",
        "## Descriptions",
        "",
    ]
    for item, species in zip(range(1, VALIDATION_SIZE + 1), sample):
        lines += [f"### {item}", "", corpus.at[species, "text_strict"], ""]
    atomic_write_text(out_md, "\n".join(lines))
    print(f"[annotation] wrote {out_csv.relative_to(cfg.root)} ({len(sheet)} items, blank score columns), {out_md.relative_to(cfg.root)} "
          f"and the item-to-species key {out_key.relative_to(cfg.root)}")


if __name__ == "__main__":
    main()
