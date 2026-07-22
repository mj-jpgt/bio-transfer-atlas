"""
Phase 18 gate: validate partial16 intervention outputs.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from intervention_common import SUPERPOPS

TAG = "partial16"
RESULTS = ROOT / f"results/tables/intervention_results.{TAG}.csv"
CONFIDENCE = ROOT / f"results/tables/per_patient_confidence.{TAG}.parquet"


def main() -> None:
    assert RESULTS.exists(), f"Missing {RESULTS}"
    res = pd.read_csv(RESULTS)
    assert len(res) > 0, f"intervention_results.{TAG}.csv is empty"

    mad = res[res["metric"] == "mean_abs_delta_EUR"].copy()
    assert not mad.empty, "No mean_abs_delta_EUR rows"

    reweight = mad[mad["mode"] == "reweight_linear"]
    reweight_wins = (reweight["reduction"] > 0).any()
    assert reweight_wins, "reweight_linear did not reduce |delta_EUR| for any PGS"

    risk = mad[(mad["mode"] == "filter_10") & (mad["pgs_id"].notna())]
    rand = mad[(mad["mode"] == "random") & (mad["pgs_id"].notna())]
    merged = risk.merge(rand, on="pgs_id", suffixes=("_risk", "_rand"), how="inner")
    beats_random = (merged["reduction_risk"] > merged["reduction_rand"]).any()
    assert beats_random, "filter_10 did not beat random on delta_EUR reduction for any PGS"

    assert CONFIDENCE.exists(), f"Missing {CONFIDENCE}"
    conf = pd.read_parquet(CONFIDENCE)
    covered = set(conf["super_pop"].dropna().unique())
    missing = set(SUPERPOPS) - covered
    assert not missing, f"per_patient_confidence missing superpops: {missing}"

    print("GATE PASS: partial16 intervention")


if __name__ == "__main__":
    main()
