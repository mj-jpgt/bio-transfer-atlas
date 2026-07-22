"""
Phase 18 gate: validate genome-wide intervention outputs.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from intervention_common import GW_TAG, SUPERPOPS

RESULTS = ROOT / f"results/tables/intervention_results.{GW_TAG}.csv"
CONFIDENCE = ROOT / f"results/tables/per_patient_confidence.{GW_TAG}.parquet"


def main() -> None:
    assert RESULTS.exists(), f"Missing {RESULTS}"
    res = pd.read_csv(RESULTS)
    assert len(res) > 0, f"intervention_results.{GW_TAG}.csv is empty"

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

    # Science-deepen: MAD rows should carry bootstrap CI columns
    assert "value_lo" in mad.columns and "value_hi" in mad.columns, "missing MAD bootstrap CI columns"
    assert mad["value_lo"].notna().any(), "MAD value_lo all missing"

    # At least one metric where reweight_linear beats random (MAD reduction or concordance gain)
    rw = mad[mad["mode"] == "reweight_linear"][["pgs_id", "reduction"]].rename(columns={"reduction": "rw"})
    rnd = mad[mad["mode"] == "random"][["pgs_id", "reduction"]].rename(columns={"reduction": "rnd"})
    m2 = rw.merge(rnd, on="pgs_id", how="inner")
    assert (m2["rw"] > m2["rnd"]).any(), "reweight_linear did not beat random on MAD reduction for any PGS"

    # Tail ancestry composition metrics present
    tail = res[res["metric"] == "tail5_ancestry_L1"]
    assert not tail.empty, "missing tail5_ancestry_L1 metrics"

    print("GATE PASS: genome-wide intervention")


if __name__ == "__main__":
    main()
