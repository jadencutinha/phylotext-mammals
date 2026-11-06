from __future__ import annotations

import hashlib

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from lineage_charisma.config import base_parser, load_config
from lineage_charisma.corpus import load_token_counter, read_corpus
from lineage_charisma.distances import cosine_distance_matrix, load_matrix, phylo_matrix_path, text_matrix_path
from lineage_charisma.embed import SentenceTransformerEncoder, embed_texts, read_species_order
from lineage_charisma.io_utils import read_csv, write_csv, write_json
from lineage_charisma.masking import PILOT_NAME_RANKS, Masker, build_vocabulary, term_leaks
from lineage_charisma.plotstyle import BLUE, INK, MUTED, ORANGE, SURFACE, apply_style, plt, save
from lineage_charisma.stats import mantel, partial_mantel, submatrix

# docs/preregistration_amendment_secondary.md, section A (secondary)
MIN_SPECIES = 100
N_BOOT = 2000


def pilot_masker(cfg) -> Masker:
    m = cfg.get("masking") or {}
    vocab = build_vocabulary(read_csv(cfg.mdd_species_path()), read_csv(cfg.mdd_synonyms_path()), extra_terms=m.get("extra_taxonomy_terms", []),
                             keep_terms=m.get("keep_terms", []), keep_when_followed_by=m.get("keep_when_followed_by", {}),
                             head_noun_rule=m.get("head_noun_rule", "pilot"), epithets_anywhere=m.get("epithets_anywhere", True),
                             name_ranks=m.get("name_ranks", PILOT_NAME_RANKS))
    return Masker(vocab, m.get("token", "[TAXON]"))


def length_distance(n_tokens: np.ndarray) -> np.ndarray:
    """As attention_length.npy: absolute difference in log token count of the unmasked text."""
    v = np.log(np.asarray(n_tokens, dtype=float))
    return np.abs(v[:, None] - v[None, :])


def _std_ranks(v: np.ndarray) -> np.ndarray:
    r = rankdata(v)
    r = r - r.mean()
    return r / np.linalg.norm(r)


def pair_statistics(idx: np.ndarray, text: np.ndarray, phylo: np.ndarray, controls: list[np.ndarray] | None) -> float:
    """Spearman r (no controls) or Spearman partial r (controls) over the pairs of a species resample, dropping pairs of two copies."""
    a, b = np.triu_indices(len(idx), k=1)
    keep = idx[a] != idx[b]
    i, j = idx[a][keep], idx[b][keep]
    x, y = _std_ranks(phylo[i, j]), _std_ranks(text[i, j])
    if not controls:
        return float(x @ y)
    zs = [_std_ranks(c[i, j]) for c in controls]
    k = len(zs)
    zz = np.array([[1.0 if p == q else zs[p] @ zs[q] for q in range(k)] for p in range(k)])
    inv = np.linalg.inv(zz)
    rxz, ryz = np.array([x @ z for z in zs]), np.array([y @ z for z in zs])
    den = (1 - rxz @ inv @ rxz) * (1 - ryz @ inv @ ryz)
    if den <= 1e-12:
        return float("nan")
    return float((x @ y - rxz @ inv @ ryz) / np.sqrt(den))


def bootstrap(stat, n: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    draws = [rng.choice(n, size=n, replace=True) for _ in range(N_BOOT)]
    return np.array([stat(d, "adw") for d in draws]), np.array([stat(d, "wikipedia") for d in draws])


def plot(results: pd.DataFrame, deltas: pd.DataFrame, path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2), constrained_layout=True)
    for ax, test in zip(axes, ("H1", "H2 partial")):
        sub = results[results["test"] == test].set_index("source")
        ys = [1, 0]
        ax.barh(ys, [sub.loc["wikipedia", "r"], sub.loc["adw", "r"]], color=[BLUE, ORANGE], height=0.55, edgecolor=SURFACE)
        for y, src in zip(ys, ("wikipedia", "adw")):
            ax.annotate(f"r = {sub.loc[src, 'r']:.3f}, p = {sub.loc[src, 'p']:.4f}", xy=(max(sub.loc[src, 'r'], 0), y), xytext=(4, 0),
                        textcoords="offset points", va="center", color=INK, fontsize=8)
        d = deltas.set_index("test").loc[test]
        ax.set_title(f"{test}, {int(sub['n_species'].iloc[0])} species\nΔ (ADW − Wikipedia) = {d['delta']:+.3f} [{d['ci_low']:+.3f}, {d['ci_high']:+.3f}]", fontsize=9)
        ax.set_yticks(ys, ["Wikipedia", "Animal Diversity Web"])
        ax.axvline(0, color=MUTED, linewidth=0.8)
        ax.set_xlabel("Spearman r, text vs phylogenetic distance")
        ax.set_xlim(min(0, results["r"].min()) - 0.02, results["r"].max() * 1.9 + 0.02)
        ax.grid(True, axis="x")
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)
    fig.suptitle("Secondary: register comparison on the same Carnivora species (95% paired species-bootstrap intervals)", x=0.01, ha="left",
                 fontweight="bold", fontsize=10)
    save(fig, path)


def main() -> None:
    args = base_parser("Register comparison (amendment 4, section A): Wikipedia vs ADW on the same pilot species.").parse_args()
    cfg = load_config(args.config, None)
    apply_style()
    scfg = cfg["stats"]
    perms, seed, method = int(scfg["permutations"]), int(scfg["seed"]), scfg["method"]
    model, rule, level = cfg.primary_combo
    order = read_species_order(cfg.species_order_path())
    adw = read_csv(cfg.root / "data" / "interim" / "adw" / "adw_texts.csv").set_index("species_id")
    subset = [s for s in order if s in adw.index]
    out, figs = cfg.root / "reports" / "tables" / "secondary", cfg.root / "reports" / "figures" / "week7"
    figs.mkdir(parents=True, exist_ok=True)
    if len(subset) < MIN_SPECIES:
        write_json(out / "register_comparison_stopped.json", {"usable_species": len(subset), "minimum": MIN_SPECIES})
        raise SystemExit(f"[register] {len(subset)} usable species is below the registered minimum {MIN_SPECIES}; analysis A stops")

    masker = pilot_masker(cfg)
    masked = [masker.mask(adw.at[s, "text"], level) for s in subset]
    leaks = [s for s, t in zip(subset, masked) if term_leaks(t, s.split()[0])]
    if leaks:
        raise SystemExit(f"[register] masked ADW text still names its own genus for {leaks}")
    cache = cfg.root / "data" / "interim" / "adw" / f"embeddings_{hashlib.sha256(('|'.join([model, rule, *masked])).encode()).hexdigest()[:16]}.npy"
    if cache.exists():
        emb = np.load(cache)
    else:
        ecfg = cfg["embedding"]
        emb = embed_texts(masked, SentenceTransformerEncoder(model, device=ecfg.get("device", "auto"), max_tokens=ecfg.get("max_tokens")), rule,
                          overlap=ecfg["chunk_overlap"], batch_size=ecfg["batch_size"]).embeddings
        np.save(cache, emb)

    phylo_full = load_matrix(phylo_matrix_path(cfg.matrices_dir()), cfg.species_order_path())[0]
    wiki_full = load_matrix(text_matrix_path(cfg.matrices_dir(), *cfg.primary_combo), cfg.species_order_path())[0]
    corpus = read_corpus(cfg.corpus_path()).set_index("species_id")
    text = {"wikipedia": submatrix(wiki_full, order, subset), "adw": cosine_distance_matrix(emb)}
    tokens = {"wikipedia": corpus.loc[subset, "n_tokens"].to_numpy(), "adw": adw.loc[subset, "n_tokens"].to_numpy()}
    phylo = submatrix(phylo_full, order, subset)

    eco_order = read_species_order(cfg.matrices_dir() / "ecology_species.txt")
    eco_full = load_matrix(cfg.matrices_dir() / "ecology.npy", cfg.matrices_dir() / "ecology_species.txt")[0]
    h2 = [s for s in subset if s in set(eco_order)]
    pos = [subset.index(s) for s in h2]
    eco = submatrix(eco_full, eco_order, h2)
    text_h2 = {k: v[np.ix_(pos, pos)] for k, v in text.items()}
    length_h2 = {k: length_distance(tokens[k][pos]) for k in text}
    phylo_h2 = phylo[np.ix_(pos, pos)]

    rows = []
    for src in ("wikipedia", "adw"):
        r1 = mantel(text[src], phylo, method=method, permutations=perms, seed=seed)
        rows.append({"status": "secondary", "test": "H1", "source": src, "r": r1.r, "p": r1.p, "n_species": r1.n, "n_pairs": r1.n_pairs,
                     "median_tokens": float(np.median(tokens[src]))})
        r2 = partial_mantel(phylo_h2, text_h2[src], [eco, length_h2[src]], method=method, permutations=perms, seed=seed)
        rows.append({"status": "secondary", "test": "H2 partial", "source": src, "r": r2.r, "p": r2.p, "n_species": r2.n, "n_pairs": r2.n_pairs,
                     "raw_r": r2.r_marginal, "median_tokens": float(np.median(tokens[src][pos]))})
    results = pd.DataFrame(rows)

    stat_h1 = lambda d, src: pair_statistics(d, text[src], phylo, None)
    stat_h2 = lambda d, src: pair_statistics(d, text_h2[src], phylo_h2, [eco, length_h2[src]])
    deltas = []
    for test, stat, n in (("H1", stat_h1, len(subset)), ("H2 partial", stat_h2, len(h2))):
        a, w = bootstrap(stat, n, seed)
        d = a - w
        ok = ~np.isnan(d)
        point = results[(results["test"] == test) & (results["source"] == "adw")]["r"].iloc[0] - results[(results["test"] == test) & (results["source"] == "wikipedia")]["r"].iloc[0]
        lo, hi = np.percentile(d[ok], [2.5, 97.5])
        verdict = "ADW higher" if lo > 0 else "Wikipedia higher" if hi < 0 else "no detectable difference"
        deltas.append({"status": "secondary", "test": test, "n_species": n, "delta": point, "ci_low": lo, "ci_high": hi, "n_boot": N_BOOT,
                       "n_undefined": int((~ok).sum()), "reading": verdict})
    deltas = pd.DataFrame(deltas)
    meta = dict(status="secondary (amendment 4, section A)", spec=f"{model} / {rule} / {level}, pilot masks", permutations=perms, seed=seed,
                source="Animal Diversity Web text, CC BY-NC-SA 3.0")
    write_csv(results, out / "register_comparison.csv", **meta)
    write_csv(deltas, out / "register_bootstrap.csv", **meta, bootstrap="paired species bootstrap, same draw for both sources, pairs of two copies dropped")
    plot(results, deltas, figs / "register_comparison.png")
    print(results.to_string(index=False))
    print(deltas.to_string(index=False))


if __name__ == "__main__":
    main()
