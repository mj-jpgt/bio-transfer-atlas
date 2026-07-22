"""Gate for pathway enrichment FDR / LOLO science-deepen artifacts."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TAG = "genomewide"
TRAITS = ["CAD", "T2D", "BMI", "LDL"]


def main() -> None:
    fdr_path = ROOT / f"results/tables/pathway_enrichment_fdr_{TAG}.csv"
    lolo_path = ROOT / f"results/tables/pathway_lolo_sensitivity_{TAG}.csv"
    assert fdr_path.exists(), f"missing {fdr_path}"
    assert lolo_path.exists(), f"missing {lolo_path}"

    fdr = pd.read_csv(fdr_path)
    assert {"trait", "rid", "bh_q", "ld_block_perm_q", "fdr_significant"}.issubset(fdr.columns)
    for trait in TRAITS:
        sub = fdr[fdr["trait"] == trait]
        assert len(sub) > 0, f"no enrichment rows for {trait}"
        n_sig = int(sub["fdr_significant"].astype(str).str.lower().isin(["true", "1"]).sum())
        # Either some survive, or explicit zero is fine (publishable null)
        print(f"  {trait}: {n_sig} FDR+LD significant pathways (of {len(sub)})")

    lolo = pd.read_csv(lolo_path)
    assert {"trait", "fragile", "fisher_p_lolo"}.issubset(lolo.columns)
    assert len(lolo) > 0

    print("GATE PASS: pathway enrichment stats present")


if __name__ == "__main__":
    main()
