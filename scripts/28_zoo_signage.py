from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd

from lineage_charisma import zoo
from lineage_charisma.config import base_parser, load_config
from lineage_charisma.distances import cosine_distance_matrix, load_matrix, phylo_matrix_path, text_matrix_path
from lineage_charisma.embed import SentenceTransformerEncoder, embed_texts, embedding_paths, load_embeddings, read_species_order
from lineage_charisma.io_utils import read_csv, write_csv, write_json
from lineage_charisma.masking import PILOT_NAME_RANKS, Masker, build_vocabulary, term_leaks
from lineage_charisma.stats import mantel, submatrix

# docs/preregistration_amendment_secondary.md, section C: descriptive only; runs end to end once signage.csv has rows
PERMS = 999


def full_scale_masker(cfg) -> Masker:
    m = cfg.get("masking") or {}
    vocab = build_vocabulary(read_csv(cfg.mdd_species_path()), read_csv(cfg.mdd_synonyms_path()), extra_terms=m.get("extra_taxonomy_terms", []),
                             keep_terms=m.get("keep_terms", []), keep_when_followed_by=m.get("keep_when_followed_by", {}),
                             head_noun_rule=m.get("head_noun_rule", "pilot"), epithets_anywhere=m.get("epithets_anywhere", True),
                             name_ranks=m.get("name_ranks", PILOT_NAME_RANKS))
    return Masker(vocab, m.get("token", "[TAXON]"))


def embed_cached(texts: dict[str, str], cfg, cache_dir: Path, encoder_factory) -> dict[str, np.ndarray]:
    """Embed with the primary model and rule; each text's vector is cached under the sha256 of (model, rule, text)."""
    model, rule, _ = cfg.primary_combo
    cache_dir.mkdir(parents=True, exist_ok=True)
    keys = {s: hashlib.sha256(f"{model}|{rule}|{t}".encode()).hexdigest() for s, t in texts.items()}
    todo = [s for s, k in keys.items() if not (cache_dir / f"{k}.npy").exists()]
    if todo:
        ecfg = cfg["embedding"]
        result = embed_texts([texts[s] for s in todo], encoder_factory(), rule, overlap=ecfg["chunk_overlap"], batch_size=ecfg["batch_size"])
        for s, v in zip(todo, result.embeddings):
            np.save(cache_dir / f"{keys[s]}.npy", v)
    return {s: np.load(cache_dir / f"{k}.npy") for s, k in keys.items()}


def main() -> None:
    parser = base_parser("Zoo signage (amendment 4, section C): reconcile, mask, embed, and run the descriptive comparison.")
    parser.add_argument("--input", type=Path, default=None, help="transcription sheet (default data/raw/zoo_signage/signage.csv)")
    parser.add_argument("--output-dir", type=Path, default=None, help="where tables go (default reports/tables/secondary)")
    args = parser.parse_args()
    cfg = load_config(args.config, "mammalia")    # zoo species span many orders: full-scale names, masks and matrices
    root = cfg.root
    source = args.input or root / "data" / "raw" / "zoo_signage" / "signage.csv"
    out_dir = args.output_dir or root / "reports" / "tables" / "secondary"
    work = cfg.interim / "zoo_signage"
    if not source.exists():
        raise SystemExit(f"{source} not found: copy data/raw/zoo_signage/template.csv to signage.csv and add one row per sign")

    signs = zoo.read_signage(source)
    usable = signs[(signs["transcription"] != "") & (signs["species_scientific"] != "")]
    print(f"[zoo] {len(signs)} rows in {source.name}; {len(signs) - len(usable)} without a transcription or scientific name are skipped")
    species_table, synonyms = read_csv(cfg.mdd_species_path()), read_csv(cfg.mdd_synonyms_path())
    resolved = zoo.resolve_signage(usable, species_table, synonyms, cfg["matching"])
    table = resolved.drop(columns=["transcription"])     # sign text stays in data/raw; only names and outcomes are written to reports
    write_csv(table, out_dir / "zoo_resolution.csv", source=source.name, rule="current name, then MDD legacy and synonym names, then fuzzy at the engine thresholds")
    ok = resolved[resolved["resolved"]]
    print(f"[zoo] resolved {len(ok)} of {len(resolved)} signs to MDD species: {ok['method'].value_counts().to_dict()}")
    for r in resolved[~resolved["resolved"]].itertuples():
        print(f"[zoo]   not resolved: row {r.row} {r.species_scientific!r} ({r.method}: {r.note})")
    if ok.empty:
        print("[zoo] nothing to analyse")
        return

    # one text per species: the distinct transcriptions of its signs, joined in row order (amendment 4, Change 1, part C)
    per_species = ok.groupby("species_id", sort=False)["transcription"].agg(lambda t: " ".join(dict.fromkeys(t)))
    masker = full_scale_masker(cfg)
    masked = {s: masker.mask(t, "strict") for s, t in per_species.items()}
    leaks = [s for s, t in masked.items() if term_leaks(t, s.split()[0])]
    if leaks:
        raise SystemExit(f"[zoo] masked sign text still contains its own genus name for {leaks}")
    model = cfg.primary_combo[0]
    factory = lambda: SentenceTransformerEncoder(model, device=cfg["embedding"].get("device", "auto"), max_tokens=cfg["embedding"].get("max_tokens"))
    sign_emb = embed_cached(masked, cfg, work / "embeddings", factory)

    wiki_emb_all, wiki_species = load_embeddings(embedding_paths(cfg.embeddings_dir(), *cfg.primary_combo))
    position = {s: i for i, s in enumerate(wiki_species)}
    wiki_emb = {s: wiki_emb_all[position[s]] for s in masked if s in position}
    adw_emb = {}
    adw_path = cfg.root / "data" / "interim" / "adw" / "adw_texts.csv"
    if adw_path.exists():   # written by the register comparison (section A); re-masked here with the full-scale vocabulary
        adw = read_csv(adw_path).set_index("species_id")["text"]
        adw_emb = embed_cached({s: masker.mask(adw[s], "strict") for s in masked if s in adw.index}, cfg, work / "adw_embeddings", factory)

    species = list(masked)
    sims = zoo.similarity_table(species, sign_emb, wiki_emb, adw_emb)
    sims["sign_tokens"] = [len(per_species[s].split()) for s in species]
    summary = zoo.similarity_summary(sims)

    frozen = read_species_order(cfg.species_order_path())
    overlap = [s for s in frozen if s in set(species)]
    h1 = []
    if len(overlap) >= zoo.MIN_H1_SPECIES:
        phylo = submatrix(load_matrix(phylo_matrix_path(cfg.matrices_dir()), cfg.species_order_path())[0], frozen, overlap)
        wiki_text = submatrix(load_matrix(text_matrix_path(cfg.matrices_dir(), *cfg.primary_combo), cfg.species_order_path())[0], frozen, overlap)
        sign_text = cosine_distance_matrix(np.stack([sign_emb[s] for s in overlap]))
        for name, matrix in (("sign", sign_text), ("wikipedia", wiki_text)):
            res = mantel(matrix, phylo, permutations=PERMS, seed=int(cfg["stats"]["seed"]))
            h1.append({"register": name, "n_species": res.n, "r": res.r, "p": res.p, "permutations": PERMS})
        with_adw = [s for s in overlap if s in adw_emb]
        if len(with_adw) >= zoo.MIN_H1_SPECIES:
            res = mantel(cosine_distance_matrix(np.stack([adw_emb[s] for s in with_adw])), submatrix(phylo, overlap, with_adw),
                         permutations=PERMS, seed=int(cfg["stats"]["seed"]))
            h1.append({"register": "adw", "n_species": res.n, "r": res.r, "p": res.p, "permutations": PERMS})
    summary["n_in_frozen_set"] = len(overlap)
    summary["h1_computed"] = bool(h1)

    meta = dict(status="secondary, descriptive (amendment 4, section C)", source=source.name, text_spec="bge-large-en-v1.5 / chunk / strict, full-scale masks")
    write_csv(sims, out_dir / "zoo_similarity.csv", **meta)
    if h1:
        write_csv(pd.DataFrame(h1), out_dir / "zoo_h1.csv", **meta, note="p-values for completeness only")
    write_json(out_dir / "zoo_summary.json", summary)
    print(sims.to_string(index=False))
    print(f"[zoo] summary: {summary}")
    if not h1:
        print(f"[zoo] H1 on the overlap not computed: {len(overlap)} signed species are in the frozen set (need {zoo.MIN_H1_SPECIES})")


if __name__ == "__main__":
    main()
