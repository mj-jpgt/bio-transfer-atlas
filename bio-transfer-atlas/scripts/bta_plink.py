"""Resolve plink2 binary across Windows (.exe) and Linux."""
from __future__ import annotations

import os
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def plink2_bin() -> str:
    env = os.environ.get("PLINK2")
    if env and Path(env).exists():
        return env
    candidates = [
        ROOT / "tools/plink2/plink2",
        ROOT / "tools/plink2/plink2.exe",
    ]
    for c in candidates:
        if c.exists():
            return str(c)
    which = shutil.which("plink2")
    if which:
        return which
    raise FileNotFoundError("plink2 not found (tools/plink2/plink2[.exe] or PATH)")
