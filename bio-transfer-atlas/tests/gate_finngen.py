from pathlib import Path

import pandas as pd


def main():
    base = Path("data/raw/finngen/chr22")
    files = sorted(base.glob("*.chr22.parquet"))
    assert len(files) >= 2, f"expected >=2 FinnGen traits, found {len(files)}"
    for f in files:
        d = pd.read_parquet(f)
        assert len(d) > 50_000, f"{f.name}: too few rows"
        assert d["chr"].astype(str).isin(["22", "chr22"]).all(), f"{f.name}: non chr22 rows"
    print("GATE PASS: FinnGen chr22 ready (GRCh38)")


if __name__ == "__main__":
    main()
