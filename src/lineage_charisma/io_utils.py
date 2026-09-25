from __future__ import annotations

import hashlib
import json
import os
import re
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import pandas as pd
import requests


def is_cached(paths: Path | Iterable[Path]) -> bool:
    items = [paths] if isinstance(paths, (str, Path)) else list(paths)
    return all(Path(p).exists() and Path(p).stat().st_size > 0 for p in items)


def skip_if_cached(paths: Path | Iterable[Path], force: bool) -> bool:
    if force:
        return False
    if is_cached(paths):
        items = [paths] if isinstance(paths, (str, Path)) else list(paths)
        for p in items:
            print(f"[cache] {p} exists, skipping (use --force to rebuild)")
        return True
    return False


def sha256(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def meta_path(path: Path) -> Path:
    path = Path(path)
    return path.with_name(path.name + ".meta.json")


def write_meta(path: Path, **meta: Any) -> Path:
    path = Path(path)
    record = {
        "file": path.name,
        "sha256": sha256(path),
        "bytes": path.stat().st_size,
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        **meta,
    }
    mp = meta_path(path)
    mp.write_text(json.dumps(record, indent=2, default=str))
    return mp


def read_meta(path: Path) -> dict[str, Any]:
    mp = meta_path(path)
    return json.loads(mp.read_text()) if mp.exists() else {}


def _atomic_target(path: Path) -> tuple[int, str]:
    path.parent.mkdir(parents=True, exist_ok=True)
    return tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.", suffix=".tmp")


def atomic_write_bytes(path: Path, data: bytes) -> Path:
    path = Path(path)
    fd, tmp = _atomic_target(path)
    with os.fdopen(fd, "wb") as fh:
        fh.write(data)
    os.chmod(tmp, 0o644)
    os.replace(tmp, path)
    return path


def atomic_write_text(path: Path, text: str) -> Path:
    return atomic_write_bytes(path, text.encode("utf-8"))


def write_json(path: Path, obj: Any) -> Path:
    return atomic_write_text(path, json.dumps(obj, indent=2, ensure_ascii=False, default=str))


def read_json(path: Path) -> Any:
    return json.loads(Path(path).read_text())


def write_csv(df: pd.DataFrame, path: Path, **meta: Any) -> Path:
    path = Path(path)
    atomic_write_text(path, df.to_csv(index=False))
    write_meta(path, rows=len(df), columns=list(df.columns), **meta)
    return path


def read_csv(path: Path, **kwargs: Any) -> pd.DataFrame:
    return pd.read_csv(path, keep_default_na=False, na_values=[""], **kwargs)


def download(url: str, dest: Path, user_agent: str | None = None, timeout: int = 120) -> Path:
    dest = Path(dest)
    headers = {"User-Agent": user_agent} if user_agent else {}
    fd, tmp = _atomic_target(dest)
    try:
        with requests.get(url, headers=headers, stream=True, timeout=timeout) as resp:
            resp.raise_for_status()
            with os.fdopen(fd, "wb") as fh:
                for block in resp.iter_content(chunk_size=1 << 20):
                    fh.write(block)
        os.chmod(tmp, 0o644)
        os.replace(tmp, dest)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)
    write_meta(dest, source_url=url)
    return dest


_VERSION_RE = re.compile(r"v(\d+(?:\.\d+)*)")


def version_key(path: Path) -> tuple[int, ...]:
    m = _VERSION_RE.search(Path(path).name)
    return tuple(int(x) for x in m.group(1).split(".")) if m else ()


def latest_versioned(directory: Path, pattern: str) -> Path:
    matches = sorted(Path(directory).rglob(pattern), key=version_key)
    matches = [m for m in matches if "__MACOSX" not in m.parts and not m.name.startswith("._")]
    if not matches:
        raise FileNotFoundError(f"no file matching {pattern} under {directory}")
    return matches[-1]


def slugify_species(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_")
