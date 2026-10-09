from __future__ import annotations

import pandas as pd

from lineage_charisma.config import base_parser, load_config
from lineage_charisma.distances import load_matrix, phylo_matrix_path, text_matrix_path
from lineage_charisma.io_utils import write_csv
from lineage_charisma.plotstyle import combo_label
from lineage_charisma.stats import mantel, submatrix

ECOLOGY = {"ecology": "ecology.npy", "ecology_unweighted": "ecology_unweighted.npy"}
ATTENTION = {"attention_pageviews": "attention_pageviews.npy", "attention_length": "attention_length.npy", "attention_combined": "attention_combined.npy"}


def main() -> None:
    args = base_parser("Marginal Mantel correlations among phylogeny, ecology, attention and the primary text distance.").parse_args()
    cfg = load_config(args.config)
    scfg = cfg["stats"]
    perms, seed, method = int(scfg["permutations"]), int(scfg["seed"]), scfg["method"]
    out_dir, order_path = cfg.matrices_dir(), cfg.species_order_path()
    primary = cfg.primary_combo

    needed = [out_dir / f for f in (*ECOLOGY.values(), *ATTENTION.values())]
    absent = [p.name for p in needed if not p.exists()]
    if absent:
        raise SystemExit(f"missing matrices {absent}; run make ecology and make attention first")

    matrices, orders = {}, {}
    matrices["phylogeny"], order = load_matrix(phylo_matrix_path(out_dir), order_path)
    matrices["text"], _ = load_matrix(text_matrix_path(out_dir, *primary), order_path)
    for name, filename in ATTENTION.items():
        matrices[name], _ = load_matrix(out_dir / filename, order_path)
    for name, filename in ECOLOGY.items():
        matrices[name], orders[name] = load_matrix(out_dir / filename, out_dir / "ecology_species.txt")

    pairs = (
        [("phylogeny", e, "predictor overlap") for e in ECOLOGY] + [("phylogeny", a, "predictor overlap") for a in ATTENTION]
        + [(e, a, "predictor overlap") for e in ECOLOGY for a in ATTENTION]
        + [("text", e, "text against predictor") for e in ECOLOGY] + [("text", a, "text against predictor") for a in ATTENTION]
    )
    rows = []
    for x, y, role in pairs:
        # ecology matrices cover only the species with trait data; every other matrix is cut down to that set for the pair
        species = next((orders[m] for m in (x, y) if m in orders), order)
        mx = matrices[x] if x in orders else submatrix(matrices[x], order, species)
        my = matrices[y] if y in orders else submatrix(matrices[y], order, species)
        res = mantel(mx, my, method=method, permutations=perms, seed=seed, alternative="two-sided")
        rows.append({"x": x, "y": y, "role": role, "method": method, "r": res.r, "p": res.p, "n_species": res.n, "n_pairs": res.n_pairs})
        print(f"[diagnostics] {x:20s} vs {y:20s} r = {res.r:+.4f}  p = {res.p:.4f}  (n = {res.n})")
    write_csv(pd.DataFrame(rows), cfg.tables_dir() / "confound_mantel.csv", permutations=perms, seed=seed, alternative="two-sided",
              text_matrix=combo_label(*primary))


if __name__ == "__main__":
    main()
