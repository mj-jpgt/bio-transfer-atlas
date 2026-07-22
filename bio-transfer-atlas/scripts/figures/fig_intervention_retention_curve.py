#!/usr/bin/env python3
"""Thin entrypoint for intervention retention–MAD curve figure (plan artifact path)."""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from eval_intervention_retention_controls import main  # noqa: E402

if __name__ == "__main__":
    main()
