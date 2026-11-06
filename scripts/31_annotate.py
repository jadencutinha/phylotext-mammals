from __future__ import annotations

import json
import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor

import pandas as pd

from lineage_charisma.charisma import DIMENSIONS, MAX_RETRIES, OUTPUT_SCHEMA, SYSTEM_PROMPT, InvalidResponse, user_message, validate_scores
from lineage_charisma.config import base_parser, load_config
from lineage_charisma.corpus import read_corpus
from lineage_charisma.embed import read_species_order
from lineage_charisma.io_utils import read_json, write_csv, write_json

# docs/preregistration_amendment_secondary.md, section B and Change 1
MODEL = "claude-sonnet-5-5"
MAX_TOKENS = 256
PRICE_PER_MTOK = {"input": 2.00, "output": 10.00}   # US$ for claude-sonnet-5-5
COST_STOP_USD = 10.0
WORKERS = 4


def labels_are_committed(cfg) -> tuple[bool, str]:
    sheet = cfg.interim / "annotation_validation.csv"
    labels = pd.read_csv(sheet)
    if labels[list(DIMENSIONS)].isna().any().any():
        return False, f"{sheet.name} still has blank scores"
    rel = str(sheet.relative_to(cfg.root))
    tracked = subprocess.run(["git", "ls-files", "--error-unmatch", rel], cwd=cfg.root, capture_output=True).returncode == 0
    clean = subprocess.run(["git", "diff", "--quiet", "HEAD", "--", rel], cwd=cfg.root).returncode == 0
    if not (tracked and clean):
        return False, f"{rel} must be committed (git add -f) with no uncommitted changes before any model request"
    return True, ""


def request_scores(client, text: str) -> tuple[dict, dict]:
    response = client.messages.create(
        model=MODEL,
        max_tokens=MAX_TOKENS,
        system=SYSTEM_PROMPT,
        thinking={"type": "between_tools"},
        output_config={"format": {"type": "json_schema", "schema": OUTPUT_SCHEMA}},
        messages=[{"role": "user", "content": user_message(text)}],
    )
    raw = next((b.text for b in response.content if b.type == "text"), "")
    record = {"model": response.model, "stop_reason": response.stop_reason, "raw": raw, "request_id": response._request_id,
              "input_tokens": response.usage.input_tokens, "output_tokens": response.usage.output_tokens}
    try:
        if response.stop_reason in ("refusal", "max_tokens"):
            raise InvalidResponse(f"stop_reason {response.stop_reason}")
        return validate_scores(raw), record
    except InvalidResponse as e:
        e.record = record     # logged with the failed attempt
        raise


def annotate(client, species: str, text: str, run: int, cache_dir) -> dict:
    path = cache_dir / f"run{run}" / f"{species.replace(' ', '_')}.json"
    if path.exists():
        return read_json(path)
    attempts = []
    for attempt in range(MAX_RETRIES + 1):
        try:
            scores, record = request_scores(client, text)
            out = {"species_id": species, "run": run, "scores": scores, "attempts": attempts, **record}
            path.parent.mkdir(parents=True, exist_ok=True)
            write_json(path, out)
            return out
        except InvalidResponse as e:
            attempts.append({"attempt": attempt, "error": str(e), **getattr(e, "record", {})})
        time.sleep(1.0)
    out = {"species_id": species, "run": run, "scores": None, "attempts": attempts}
    path.parent.mkdir(parents=True, exist_ok=True)
    write_json(path, out)
    return out


def usage_cost(records) -> float:
    tokens_in = sum(r.get("input_tokens", 0) + sum(a.get("input_tokens", 0) for a in r.get("attempts", [])) for r in records)
    tokens_out = sum(r.get("output_tokens", 0) + sum(a.get("output_tokens", 0) for a in r.get("attempts", [])) for r in records)
    return tokens_in / 1e6 * PRICE_PER_MTOK["input"] + tokens_out / 1e6 * PRICE_PER_MTOK["output"]


def main() -> None:
    parser = base_parser("LLM charisma annotation (amendment 4, section B): two runs of claude-sonnet-5-5 over the pilot species.")
    parser.add_argument("--run", type=int, choices=(1, 2), required=True)
    args = parser.parse_args()
    cfg = load_config(args.config, None)
    ok, why = labels_are_committed(cfg)
    if not ok:
        raise SystemExit(f"[annotate] refusing to call the model: {why}")
    if not os.environ.get("ANTHROPIC_API_KEY"):
        raise SystemExit("[annotate] ANTHROPIC_API_KEY is not set")
    if args.run == 2 and len(list((cfg.interim / "annotation" / "run1").glob("*.json"))) < 284:
        raise SystemExit("[annotate] run 1 must be complete before run 2")
    import anthropic

    client = anthropic.Anthropic()    # reads ANTHROPIC_API_KEY from the environment; the key is never written anywhere
    order = read_species_order(cfg.species_order_path())
    corpus = read_corpus(cfg.corpus_path()).set_index("species_id")
    cache_dir = cfg.interim / "annotation"
    with ThreadPoolExecutor(WORKERS) as pool:
        records = list(pool.map(lambda s: annotate(client, s, corpus.at[s, "text_strict"], args.run, cache_dir), order))
    cost = usage_cost(records)
    failed = [r["species_id"] for r in records if r["scores"] is None]
    retried = sum(bool(r["attempts"]) for r in records)
    table = pd.DataFrame([{"species_id": r["species_id"], **(r["scores"] or {d: None for d in DIMENSIONS})} for r in records])
    write_csv(table, cache_dir / f"scores_run{args.run}.csv", model=MODEL, run=args.run, cost_usd=round(cost, 4), failed=failed)
    write_json(cache_dir / f"run{args.run}_summary.json", {"model": MODEL, "models_returned": sorted({r.get("model") for r in records if r.get("model")}),
                                                           "cost_usd": cost, "failed": failed, "species_with_retries": retried})
    print(f"[annotate] run {args.run}: {len(records) - len(failed)}/{len(records)} valid, {retried} needed retries, failed {failed}; cost ${cost:.2f}")
    if args.run == 1 and 2 * cost > COST_STOP_USD:
        raise SystemExit(f"[annotate] projected total ${2 * cost:.2f} is above the registered stop of ${COST_STOP_USD:.0f}; stopping before run 2")


if __name__ == "__main__":
    main()
