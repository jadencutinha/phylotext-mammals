from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lineage_charisma.config import base_parser, load_config
from lineage_charisma.corpus import read_corpus
from lineage_charisma.distances import canonical_order
from lineage_charisma.embed import read_species_order, write_species_order
from lineage_charisma.io_utils import read_csv, write_csv
from lineage_charisma.plotstyle import AXIS, BLUE, INK_2, MUTED, apply_style, save

MIN_ORDER_SIZE = 50   # docs/preregistration_amendment_fullscale.md, section 7
FIGURE_MIN_ARTICLES = 10


def by_order(series: pd.Series, name: str) -> pd.Series:
    return series.value_counts().rename(name)


def plot_stubs(stubs: pd.DataFrame, min_tokens: int, path) -> None:
    shown = stubs[stubs["n_articles"] >= FIGURE_MIN_ARTICLES].sort_values("stub_rate")
    fig, ax = plt.subplots(figsize=(7.2, 0.28 * len(shown) + 1.4), constrained_layout=True)
    ys = np.arange(len(shown))
    ax.barh(ys, shown["stub_rate"] * 100, color=BLUE, height=0.62)
    for y, (_, r) in zip(ys, shown.iterrows()):
        ax.text(r["stub_rate"] * 100 + 0.6, y, f"{int(r['n_stubs']):,} of {int(r['n_articles']):,}", va="center", fontsize=7.5, color=INK_2)
    ax.set_yticks(ys, shown.index)
    ax.set_xlabel(f"share of matched species whose Wikipedia text has fewer than {min_tokens} tokens (%)")
    ax.set_xlim(0, max(5.0, shown["stub_rate"].max() * 100 * 1.25))
    ax.grid(axis="x")
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_color(AXIS)
    ax.set_title(f"Stub articles by mammalian order (orders with at least {FIGURE_MIN_ARTICLES} matched articles)")
    hidden = stubs[stubs["n_articles"] < FIGURE_MIN_ARTICLES]
    if len(hidden):
        fig.text(0.01, -0.01, f"Not shown: {len(hidden)} orders with fewer than {FIGURE_MIN_ARTICLES} articles "
                 f"({int(hidden['n_stubs'].sum())} stubs among {int(hidden['n_articles'].sum())} articles).", fontsize=7.5, color=MUTED)
    save(fig, path)


def main() -> None:
    parser = base_parser("Full scale: freeze species_order.txt and write the species and stub coverage tables by order.")
    parser.add_argument("--candidate", action="store_true",
                        help="write the list to data/interim/<clade>/species_order_candidate.txt instead of freezing it, so it can be reviewed and committed first")
    args = parser.parse_args()
    cfg = load_config(args.config, args.clade)
    if not cfg.output_subdir:
        raise SystemExit("this step is for full-scale runs (e.g. --clade mammalia); the pilot order is written by scripts/07_build_matrices.py")
    apply_style()
    min_tokens = int(cfg["corpus"]["min_tokens"])

    mdd = read_csv(cfg.mdd_species_path()).set_index("species")
    matched, residual, excluded = read_csv(cfg.matched_path()), read_csv(cfg.residual_path()), read_csv(cfg.excluded_species_path())
    corpus = read_corpus(cfg.corpus_path())
    corpus["order"] = mdd.loc[corpus["species_id"], "order"].to_numpy()
    wiki = read_csv(cfg.interim / "wiki_texts.csv").set_index("species")
    corpus["has_description"] = wiki.loc[corpus["species_id"], "has_description"].astype(bool).to_numpy()
    if set(corpus["species_id"]) != set(matched["species"]):
        raise SystemExit("[freeze] the corpus species are not exactly the matched species; rebuild the corpus (scripts/05b_build_corpus.py)")

    order = canonical_order(corpus.loc[~corpus["is_stub"], "species_id"])
    path = cfg.interim / "species_order_candidate.txt" if args.candidate else cfg.species_order_path()
    if not args.candidate and path.exists() and read_species_order(path) != order and not args.force:
        raise SystemExit(f"[freeze] {path.relative_to(cfg.root)} already holds a different frozen list; it is not overwritten without --force")
    write_species_order(path, order)

    table = pd.concat([
        by_order(mdd["order"], "mdd_species"),
        by_order(excluded.loc[excluded["reason"] != "no_tree_tip", "order"], "excluded_extinct_or_domestic_no_tip"),
        by_order(excluded.loc[excluded["reason"] == "no_tree_tip", "order"], "excluded_no_tree_tip"),
        *[by_order(residual.loc[residual["issue"] == issue, "order"], f"residual_{issue}") for issue in sorted(residual["issue"].unique())],
        by_order(matched["order"], "matched_with_article"),
        by_order(corpus.loc[corpus["is_stub"], "order"], "stubs"),
        by_order(corpus.loc[~corpus["is_stub"], "order"], "frozen_species"),
    ], axis=1).fillna(0).astype(int).sort_values("mdd_species", ascending=False).rename_axis("order")
    table["frozen_share_of_mdd"] = table["frozen_species"] / table["mdd_species"]
    table["in_per_order_analyses"] = table["frozen_species"] >= MIN_ORDER_SIZE
    tables = cfg.tables_dir()
    write_csv(table.reset_index(), tables / "species_by_order.csv", n_frozen=len(order), min_order_size=MIN_ORDER_SIZE)

    tokens = corpus.groupby("order")["n_tokens"]
    stubs = pd.DataFrame({"n_articles": corpus["order"].value_counts(), "n_stubs": corpus.loc[corpus["is_stub"], "order"].value_counts(),
                          "median_tokens": tokens.median(), "share_with_description_section": corpus.groupby("order")["has_description"].mean()}).fillna({"n_stubs": 0}).astype({"n_stubs": int})
    stubs["stub_rate"] = stubs["n_stubs"] / stubs["n_articles"]
    stubs = stubs.sort_values("stub_rate", ascending=False).rename_axis("order")
    write_csv(stubs.reset_index(), tables / "stubs_by_order.csv", min_tokens=min_tokens, tokenizer=cfg["corpus"]["tokenizer"])
    fig_path = cfg.week6_figures_dir() / "stubs_by_order.png"
    plot_stubs(stubs, min_tokens, fig_path)

    big = table[table["in_per_order_analyses"]]
    print(f"[freeze] {len(order)} species {'in the candidate list' if args.candidate else 'frozen'} -> {path.relative_to(cfg.root)}")
    print(f"[freeze] {len(matched)} matched with an article, {int(corpus['is_stub'].sum())} stubs under {min_tokens} tokens")
    print(f"[freeze] orders with at least {MIN_ORDER_SIZE} species: " + ", ".join(f"{o} ({n})" for o, n in big["frozen_species"].items()))
    print(f"[freeze] wrote {(tables / 'species_by_order.csv').relative_to(cfg.root)}, stubs_by_order.csv and {fig_path.relative_to(cfg.root)}")


if __name__ == "__main__":
    main()
