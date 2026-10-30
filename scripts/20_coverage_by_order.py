from __future__ import annotations

import pandas as pd

from lineage_charisma.config import base_parser, load_config
from lineage_charisma.embed import read_species_order
from lineage_charisma.io_utils import read_csv, write_csv
from lineage_charisma.traits import trait_list


def main() -> None:
    parser = base_parser("Full scale: trait and pageview coverage of the frozen species set, by order.")
    args = parser.parse_args()
    cfg = load_config(args.config, args.clade)
    order = read_species_order(cfg.species_order_path())
    mdd_order = read_csv(cfg.mdd_species_path()).set_index("species")["order"]
    traits = trait_list(bool(cfg["traits"].get("include_volant", False)))
    columns = [t.column for t in traits]

    table = read_csv(cfg.processed / "traits.csv").set_index("species_id").loc[order, columns]
    eco = set(read_species_order(cfg.matrices_dir() / "ecology_species.txt"))
    elton = [t.column for t in traits if t.source == "elton"]
    combine = [t.column for t in traits if t.source == "combine"]
    att = read_csv(cfg.processed / "attention.csv").set_index("species_id").loc[order]

    frame = pd.DataFrame({
        "order": mdd_order.loc[order].to_numpy(),
        "any_elton": table[elton].notna().any(axis=1).to_numpy(),
        "any_combine": table[combine].notna().any(axis=1).to_numpy(),
        "all_traits": table.notna().all(axis=1).to_numpy(),
        "no_traits": table.isna().all(axis=1).to_numpy(),
        "in_ecology_species": [s in eco for s in order],
        "has_pageviews": (att["total_pageviews"] > 0).to_numpy(),
        "all_months_viewed": (att["months_with_views"] == att["months_with_views"].max()).to_numpy(),
        "median_monthly_pageviews": att["median_monthly_pageviews"].to_numpy(),
    }, index=order)
    agg = dict(
        n_species=("order", "size"),
        with_elton_record=("any_elton", "sum"), with_combine_record=("any_combine", "sum"), all_traits=("all_traits", "sum"),
        no_traits=("no_traits", "sum"), in_ecology_species=("in_ecology_species", "sum"),
        with_pageviews=("has_pageviews", "sum"), viewed_every_month=("all_months_viewed", "sum"),
        median_of_median_monthly_pageviews=("median_monthly_pageviews", "median"),
    )
    counts = frame.groupby("order").agg(**agg).sort_values("n_species", ascending=False)
    counts.loc["total"] = frame.assign(order="total").groupby("order").agg(**agg).loc["total"]
    counts = counts.astype({c: int for c in counts.columns if c != "median_of_median_monthly_pageviews"})
    missing = pd.DataFrame({"order": frame["order"], **{c: table[c].isna().to_numpy() for c in columns}}).groupby("order").mean()

    tables = cfg.tables_dir()
    write_csv(counts.rename_axis("order").reset_index(), tables / "coverage_by_order.csv", n_species=len(order), traits=columns,
              pageview_window=[cfg["attention"]["start"], cfg["attention"]["end"]])
    write_csv(missing.rename_axis("order").reset_index(), tables / "trait_missing_rate_by_order.csv", traits=columns)
    print(counts.to_string())
    print(f"[coverage] wrote {(tables / 'coverage_by_order.csv').relative_to(cfg.root)} and trait_missing_rate_by_order.csv")


if __name__ == "__main__":
    main()
