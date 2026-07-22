"""Cross-platform curl binary resolution (Windows curl.exe / Linux curl)."""
from __future__ import annotations

import shutil


def curl_bin() -> str:
    for name in ("curl", "curl.exe"):
        path = shutil.which(name)
        if path:
            return path
    raise RuntimeError("curl not found on PATH")
