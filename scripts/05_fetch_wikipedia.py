from __future__ import annotations

import pandas as pd

from lineage_charisma.config import base_parser, load_config
from lineage_charisma.io_utils import atomic_write_text, is_cached, read_csv, read_json, skip_if_cached, slugify_species, write_csv, write_json
from lineage_charisma.wiki import WikiClient, page_from_response, preprocess

SUMMARY_COLUMNS = [
    "species", "wiki_title", "pageid", "revid", "rev_timestamp", "wikibase_item",
    "has_description", "description_heading", "n_chars", "n_words", "raw_path", "clean_path",
]


def main() -> None:
    parser = base_parser("Fetch and clean the Wikipedia lead + Description section for each matched species.")
    parser.add_argument("--offline", action="store_true", help="use only cached API responses")
    args = parser.parse_args()
    cfg = load_config(args.config)
    wiki_dir = cfg.wiki_dir()
    raw_dir, clean_dir = wiki_dir / "raw", wiki_dir / "clean"
    raw_dir.mkdir(exist_ok=True)
    clean_dir.mkdir(exist_ok=True)
    summary_path = cfg.interim / "wiki_texts.csv"
    if skip_if_cached(summary_path, args.force):
        return

    matched = read_csv(cfg.matched_path())
    headings = cfg["wikipedia"]["description_headings"]
    client = WikiClient.from_config(cfg, offline=args.offline)
    rows, failures = [], []
    for i, rec in enumerate(matched.itertuples(index=False), 1):
        slug = slugify_species(rec.species)
        raw_path, clean_path = raw_dir / f"{slug}.json", clean_dir / f"{slug}.txt"
        if args.force or not is_cached(raw_path):
            data = client.fetch_page(rec.wiki_title)
            write_json(raw_path, data)
        else:
            data = read_json(raw_path)
        try:
            page = page_from_response(data)
        except ValueError as exc:
            failures.append((rec.species, rec.wiki_title, str(exc)))
            continue
        cleaned = preprocess(page["extract"], headings)
        atomic_write_text(clean_path, cleaned["text"] + "\n")
        rows.append({
            "species": rec.species,
            "wiki_title": page["title"],
            "pageid": page["pageid"],
            "revid": page["revid"],
            "rev_timestamp": page["rev_timestamp"],
            "wikibase_item": page["wikibase_item"],
            "has_description": cleaned["has_description"],
            "description_heading": cleaned["description_heading"],
            "n_chars": cleaned["n_chars"],
            "n_words": cleaned["n_words"],
            "raw_path": str(raw_path.relative_to(cfg.root)),
            "clean_path": str(clean_path.relative_to(cfg.root)),
        })
        if i % 50 == 0:
            print(f"[wiki] {i}/{len(matched)}")

    summary = pd.DataFrame(rows, columns=SUMMARY_COLUMNS)
    if failures:
        for f in failures:
            print(f"[fail] {f}")
        raise SystemExit(f"{len(failures)} pages failed; summary not written")
    write_csv(summary, summary_path, headings=headings, rule="lead + first matching description section (with subsections), TextExtracts plaintext, clean_text()")
    mismatched = summary[summary["wiki_title"] != matched.set_index("species").loc[summary["species"], "wiki_title"].values]
    print(f"[wiki] {len(summary)} pages cached; with description section: {int(summary['has_description'].sum())}; "
          f"median words: {int(summary['n_words'].median())}; min words: {int(summary['n_words'].min())}")
    print(f"[wiki] description headings used: {summary['description_heading'].value_counts(dropna=False).to_dict()}")
    if len(mismatched):
        print(f"[warn] {len(mismatched)} titles redirected since reconciliation: {mismatched['species'].tolist()}")


if __name__ == "__main__":
    main()
