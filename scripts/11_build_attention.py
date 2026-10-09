from __future__ import annotations

import numpy as np
import pandas as pd

from lineage_charisma.attention import PageviewsClient, attention_matrices, attention_variables, month_range, pairwise_mean_matrix, resolve_articles, summed_monthly_views
from lineage_charisma.config import base_parser, load_config
from lineage_charisma.corpus import read_corpus
from lineage_charisma.distances import order_hash, save_matrix
from lineage_charisma.embed import read_species_order
from lineage_charisma.io_utils import read_csv, write_csv, write_meta
from lineage_charisma.wiki import WikiClient

MEAN_MATRIX = "attention_mean_log_pageviews.npy"


def main() -> None:
    parser = base_parser("Fetch Wikipedia pageviews and build the attention distance matrices.")
    parser.add_argument("--offline", action="store_true", help="use only cached API responses")
    args = parser.parse_args()
    cfg = load_config(args.config)
    acfg = cfg["attention"]
    start, end = acfg["start"], acfg["end"]
    months = month_range(start, end)

    order = read_species_order(cfg.species_order_path())
    matched = read_csv(cfg.matched_path()).set_index("species").loc[order]
    corpus = read_corpus(cfg.corpus_path()).set_index("species_id").loc[order]

    wiki = WikiClient.from_config(cfg, offline=args.offline)
    views = PageviewsClient.from_config(cfg, offline=args.offline)
    articles = resolve_articles(wiki, matched["wiki_title"])
    n_titles = sum(1 + len(redirects) for _, redirects in articles.values())
    print(f"[attention] {len(articles)} articles and {n_titles - len(articles)} redirects to them; window {start} to {end} ({len(months)} months), "
          f"{acfg['access']}, agent={acfg['agent']}")

    rows, monthly = [], {}
    for i, species in enumerate(order, 1):
        title = matched.loc[species, "wiki_title"]
        current, redirects = articles[title]
        own = summed_monthly_views(views, [current], start, end)
        total = summed_monthly_views(views, [current, *redirects], start, end)
        monthly[species] = total
        rows.append({
            "species_id": species, "wiki_title": title, "current_title": current, "renamed": current != title, "n_redirects": len(redirects),
            "median_monthly_pageviews": float(total.median()), "total_pageviews": int(total.sum()), "months_with_views": int((total > 0).sum()),
            "redirect_share": float(1 - own.sum() / total.sum()) if total.sum() else 0.0, "n_tokens": int(corpus.loc[species, "n_tokens"]),
        })
        if i % 50 == 0 or i == len(order):
            print(f"[attention] {i}/{len(order)} species ({views.network_calls} pageview requests so far)")

    table = pd.DataFrame(rows)
    variables = attention_variables(table["median_monthly_pageviews"], table["n_tokens"])
    table = pd.concat([table, variables], axis=1)
    meta = dict(window=[start, end], project=acfg["project"], access=acfg["access"], agent=acfg["agent"], redirects_included=True)
    write_csv(table, cfg.processed / "attention.csv", **meta)
    write_csv(pd.DataFrame(monthly).T.rename_axis("species_id").reset_index(), cfg.processed / "pageviews_monthly.csv", **meta)

    out_dir = cfg.matrices_dir()
    for name, matrix in attention_matrices(variables).items():
        save_matrix(out_dir / f"{name}.npy", matrix, order, kind=name, **meta)
        off = matrix[np.triu_indices(len(order), k=1)]
        print(f"[attention] {name}.npy: mean {off.mean():.3f}, max {off.max():.3f}")
    # not a distance (non-zero diagonal), so it is saved without the distance-matrix validation in save_matrix
    mean_path = out_dir / MEAN_MATRIX
    np.save(mean_path, pairwise_mean_matrix(variables["log_pageviews"]))
    write_meta(mean_path, shape=[len(order)] * 2, species_order_sha256=order_hash(order), kind="pairwise mean of log pageviews; similarity-type covariate, not a distance", **meta)

    shared = table[table.duplicated("current_title", keep=False)]
    print(f"[attention] median monthly pageviews: min {table['median_monthly_pageviews'].min():.0f}, median {table['median_monthly_pageviews'].median():.0f}, "
          f"max {table['median_monthly_pageviews'].max():.0f}; {int((table['months_with_views'] < len(months)).sum())} species have a month with no views; "
          f"{int(table['renamed'].sum())} titles now redirect elsewhere; {len(shared)} species share an article with another species")


if __name__ == "__main__":
    main()
