from __future__ import annotations

import pandas as pd

from lineage_charisma.adw import ADWClient, extract_text, parse_account, query_names
from lineage_charisma.config import base_parser, load_config
from lineage_charisma.corpus import load_token_counter
from lineage_charisma.embed import read_species_order
from lineage_charisma.io_utils import read_csv, write_csv
from lineage_charisma.taxonomy import normalize_name

# docs/preregistration_amendment_secondary.md, section A and Change 1: one request per 8 s, every response cached,
# full texts kept in data/ (not committed); only the coverage table goes to reports/
MIN_SPECIES = 100


def main() -> None:
    parser = base_parser("Register comparison: fetch ADW species accounts for the pilot species and apply the fixed text rule.")
    parser.add_argument("--offline", action="store_true", help="use cached pages only")
    args = parser.parse_args()
    cfg = load_config(args.config, None)
    order = read_species_order(cfg.species_order_path())
    mdd = read_csv(cfg.mdd_species_path()).set_index("species")
    synonyms = read_csv(cfg.mdd_synonyms_path())
    matched = read_csv(cfg.matched_path()).set_index("species")
    client = ADWClient(cfg.root / "data" / "raw" / "adw" / "pages", cfg["wikipedia"]["user_agent"], offline=args.offline)
    count_tokens = load_token_counter(cfg["corpus"]["tokenizer"])
    min_tokens = int(cfg["corpus"]["min_tokens"])

    rows = []
    for k, species in enumerate(order, 1):
        mdd_id = int(mdd.at[species, "mdd_id"])
        combos = synonyms.loc[synonyms["mdd_id"] == mdd_id, "normalized_original_combination"].dropna().tolist()
        names = query_names(species, matched.at[species, "tree_binomial"], combos)
        row = {"species_id": species, "family": mdd.at[species, "family"], "adw_name": None, "route": "not_found", "n_names_tried": 0, "text": ""}
        for i, name in enumerate(names):
            html, meta = client.fetch(name)
            row["n_names_tried"] = i + 1
            account = parse_account(html)
            if account.is_species_account and normalize_name(account.title) == name:
                tree_name = normalize_name(matched.at[species, "tree_binomial"])
                route = "current_name" if name == species else "tree_name" if name == tree_name else "synonym"
                row.update(adw_name=name, route=route, text=extract_text(account), url=meta["url"])
                break
        rows.append(row)
        if k % 20 == 0:
            found = sum(r["adw_name"] is not None for r in rows)
            print(f"[adw] {k}/{len(order)} species, {found} accounts found, {client.requests_made} requests this run", flush=True)

    table = pd.DataFrame(rows)
    table["n_tokens"] = [count_tokens([t])[0] if t else 0 for t in table["text"]]
    table["n_words"] = table["text"].str.split().str.len().fillna(0).astype(int)
    table["usable"] = table["n_tokens"] >= min_tokens
    texts_path = cfg.root / "data" / "interim" / "adw" / "adw_texts.csv"
    write_csv(table.loc[table["usable"], ["species_id", "adw_name", "route", "url", "n_tokens", "text"]], texts_path,
              licence="Animal Diversity Web text, CC BY-NC-SA 3.0; not for redistribution in this repository",
              rule="geographic_range, habitat, physical_description, behavior, food_habits paragraphs, in that order")
    coverage = table.drop(columns=["text"])
    out = cfg.root / "reports" / "tables" / "secondary"
    write_csv(coverage, out / "adw_coverage.csv", min_tokens=min_tokens, source="Animal Diversity Web (animaldiversity.org), CC BY-NC-SA 3.0")
    by_family = coverage.groupby("family").agg(species=("species_id", "size"), accounts=("adw_name", lambda s: int(s.notna().sum())),
                                                usable=("usable", "sum")).sort_values("species", ascending=False)
    write_csv(by_family.reset_index(), out / "adw_coverage_by_family.csv")
    n_usable = int(table["usable"].sum())
    print(f"[adw] accounts found {int(table['adw_name'].notna().sum())}/{len(order)} {table['route'].value_counts().to_dict()}; "
          f"usable (>= {min_tokens} tokens) {n_usable}; median tokens {int(table.loc[table['usable'], 'n_tokens'].median()) if n_usable else 0}; "
          f"{client.requests_made} requests this run")
    print(by_family.to_string())
    if n_usable < MIN_SPECIES:
        print(f"[adw] coverage {n_usable} is below the registered minimum of {MIN_SPECIES}: analysis A stops here (amendment 4, section A)")


if __name__ == "__main__":
    main()
