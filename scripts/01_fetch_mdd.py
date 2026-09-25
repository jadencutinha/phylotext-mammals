from __future__ import annotations

import zipfile

from lineage_charisma.config import base_parser, load_config
from lineage_charisma.io_utils import download, is_cached, latest_versioned, read_meta, skip_if_cached, version_key, write_csv
from lineage_charisma.taxonomy import filter_clade, load_mdd, load_mdd_synonyms, mdd_species_table


def main() -> None:
    args = base_parser("Download the Mammal Diversity Database and extract the target clade.").parse_args()
    cfg = load_config(args.config)
    src = cfg["sources"]["mdd"]
    out_species, out_syn = cfg.mdd_species_path(), cfg.mdd_synonyms_path()
    if skip_if_cached([out_species, out_syn], args.force):
        return

    mdd_dir = cfg.raw / "mdd"
    mdd_dir.mkdir(parents=True, exist_ok=True)
    zip_path = mdd_dir / "MDD.zip"
    if args.force or not is_cached(zip_path):
        print(f"[download] {src['zip_url']}")
        download(src["zip_url"], zip_path, user_agent=cfg["wikipedia"]["user_agent"])
    with zipfile.ZipFile(zip_path) as zf:
        members = [m for m in zf.namelist() if not m.startswith("__MACOSX") and not m.endswith(".DS_Store")]
        zf.extractall(mdd_dir, members=members)

    species_csv = latest_versioned(mdd_dir, src["species_glob"])
    synonym_csv = latest_versioned(mdd_dir, src["synonym_glob"])
    version = "v" + ".".join(map(str, version_key(species_csv)))
    print(f"[mdd] using {species_csv.name} and {synonym_csv.name} ({version})")

    raw = load_mdd(species_csv)
    clade = filter_clade(raw, cfg.clade_rank, cfg.clade_name)
    if clade.empty:
        raise SystemExit(f"no MDD species for {cfg.clade_rank}={cfg.clade_name}")
    species = mdd_species_table(clade)
    synonyms = load_mdd_synonyms(synonym_csv, species["mdd_id"])

    zip_meta = read_meta(zip_path)
    meta = dict(mdd_version=version, mdd_zip_sha256=zip_meta.get("sha256"), source_url=src["zip_url"], clade=f"{cfg.clade_rank}={cfg.clade_name}")
    write_csv(species, out_species, source_file=species_csv.name, **meta)
    write_csv(synonyms, out_syn, source_file=synonym_csv.name, **meta)
    print(f"[mdd] {len(species)} species, {len(synonyms)} synonym records -> {out_species.name}, {out_syn.name}")
    print(f"[mdd] extinct={int(species['extinct'].sum())} domestic={int(species['domestic'].sum())} families={species['family'].nunique()}")


if __name__ == "__main__":
    main()
