"""Post-process SuSiE outputs: mark top-decile PIP as in_cs, rebuild tiers, dump status."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    d = ROOT / "data/labels/susie"
    for f in d.glob("susie_*.parquet"):
        df = pd.read_parquet(f)
        if "pip" not in df.columns:
            continue
        thr = float(df["pip"].quantile(0.9)) if len(df) else 0.0
        thr = max(thr, 1e-6)
        df["in_cs"] = df["pip"] >= thr
        df.to_parquet(f, index=False)
        print(f.name, "n", len(df), "in_cs", int(df.in_cs.sum()), "thr", thr, "maxpip", float(df.pip.max()))

    # rebuild tiers via import
    import subprocess
    import sys

    subprocess.check_call([sys.executable, str(ROOT / "scripts/build_finemap_tier_labels.py"), "--tag", "genomewide"])
    t = pd.read_parquet(ROOT / "data/labels/finemap_tiers_genomewide.parquet")
    print(t.groupby(["tier_method", "finemap_tier"]).size())

    rows = []
    for c in range(1, 8):
        p = ROOT / f"data/interim/1000g_grch38/chr{c}.score.pgen"
        rows.append(
            {
                "chrom": c,
                "status": "ok" if p.exists() and p.stat().st_size > 1000 else "pending",
                "bytes": p.stat().st_size if p.exists() else 0,
            }
        )
    out = ROOT / "results/tables/score_pgen_chr1_7_rebuild_status.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    print(pd.DataFrame(rows))


if __name__ == "__main__":
    main()
