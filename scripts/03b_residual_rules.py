from __future__ import annotations

import pandas as pd

from lineage_charisma.config import base_parser, load_config
from lineage_charisma.io_utils import read_csv, read_meta, write_csv

# docs/preregistration_amendment_fullscale.md, section 3 and Change 1: at full scale only these tree matches are
# accepted without review; the others go to residual_review.csv unless the pilot already hand-reviewed the same match
AUTO_ACCEPT = {"exact", "legacy_name", "synonym"}


def main() -> None:
    args = base_parser("Full scale: move tree matches made by fuzzy, epithet or Wikipedia-bridge rules to the residual pile.").parse_args()
    cfg = load_config(args.config, args.clade)
    if not cfg.output_subdir:
        raise SystemExit("this step applies to full-scale runs only (e.g. --clade mammalia); the pilot matches were hand-reviewed")
    matched, residual = read_csv(cfg.matched_path()), read_csv(cfg.residual_path())
    pilot = read_csv(load_config(args.config, clade=None).matched_path()).set_index("species")["tree_tip"]

    reviewed = matched["species"].map(pilot).eq(matched["tree_tip"])
    move = ~matched["match_method"].isin(AUTO_ACCEPT) & ~reviewed
    kept_by_pilot = matched[~matched["match_method"].isin(AUTO_ACCEPT) & reviewed]
    moved = matched[move].assign(
        issue="unreviewed_" + matched.loc[move, "match_method"],
        candidate_tip=matched.loc[move, "tree_tip"],
        candidate_score=matched.loc[move, "confidence"],
        note="auto-accepted by the engine; full-scale rule sends non-exact, non-synonym matches to review",
    )[residual.columns]

    keep = {k: v for k, v in read_meta(cfg.matched_path()).items() if k in ("clade", "matching")}
    meta = {**keep, "residual_rule": f"auto-accept {sorted(AUTO_ACCEPT)}; others to residual unless the pilot reviewed the same tip"}
    write_csv(matched[~move], cfg.matched_path(), **meta)
    write_csv(pd.concat([residual, moved], ignore_index=True).sort_values("species"), cfg.residual_path(), **meta)
    print(f"[residual-rules] moved {len(moved)} matches to {cfg.residual_path().relative_to(cfg.root)}: {moved['issue'].value_counts().to_dict()}")
    print(f"[residual-rules] kept {len(kept_by_pilot)} non-exact matches reviewed in the pilot: {kept_by_pilot['species'].tolist()}")
    print(f"[residual-rules] matched species: {int((~move).sum())}")


if __name__ == "__main__":
    main()
