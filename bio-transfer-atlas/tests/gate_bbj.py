from pathlib import Path

import pandas as pd


def main():
    base = Path("data/raw/bbj/chr22")
    files = sorted(base.glob("*.chr22.parquet"))
    assert len(files) == 4, f"expected 4 BBJ trait files, found {len(files)}"
    for f in files:
        d = pd.read_parquet(f)
        assert len(d) > 50_000, f"{f.name}: too few rows"
        assert (d["chr"].astype(str) == "22").all(), f"{f.name}: non-chr22 rows found"
        assert d["beta"].notna().mean() > 0.9, f"{f.name}: sparse beta"
    print("GATE PASS: BBJ chr22 downloaded & filtered")


if __name__ == "__main__":
    main()
