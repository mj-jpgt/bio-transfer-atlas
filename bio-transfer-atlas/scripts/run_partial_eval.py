"""
Run concat + per-chromosome and pooled mechanism ablations in parallel.

Example:
  python scripts/run_partial_eval.py --chroms 1,2,3,4,5,22 --workers 3
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def parse_chrom_list(spec: str) -> list[str]:
    return [str(int(x.strip())) for x in spec.split(",") if x.strip()]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--chroms", default="1,2,3,4,5,22")
    p.add_argument("--workers", type=int, default=3)
    p.add_argument("--tag", default="partial6")
    p.add_argument("--python", default=sys.executable)
    p.add_argument("--skip-concat", action="store_true")
    p.add_argument("--skip-per-chr", action="store_true")
    p.add_argument("--skip-pooled", action="store_true")
    return p.parse_args()


def run(cmd: list[str], desc: str) -> None:
    print(f"\n[{desc}] {' '.join(cmd)}", flush=True)
    r = subprocess.run(cmd, cwd=ROOT, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"FAILED ({r.returncode}): {desc}")


def run_chr_baseline(py: str, chrom: str, tag: str) -> dict:
    master = ROOT / f"data/modeling/master_variant_table.chr{chrom}.parquet"
    groups = ROOT / f"data/modeling/feature_groups.chr{chrom}.json"
    out_suffix = f".chr{chrom}.{tag}"
    cmd = [
        py,
        "scripts/compute_baselines.py",
        "--master",
        str(master),
        "--groups",
        str(groups),
        "--subsets",
        "associated",
        "--associated-only",
        "--splits",
        "split_variant",
        "--skip-riskclass",
        "--out-suffix",
        out_suffix,
    ]
    r = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True)
    if r.returncode != 0:
        return {"chrom": chrom, "ok": False, "error": r.stderr[-500:]}
    cls = ROOT / f"results/tables/ablation_classification{out_suffix}.csv"
    row = None
    if cls.exists():
        import pandas as pd

        df = pd.read_csv(cls)
        row = df[
            (df.subset == "associated")
            & (df.split == "split_variant")
            & (df.model == "hgb")
            & (df.feature_group == "AF_LD_SEL")
        ]
        if len(row):
            row = row.iloc[0].to_dict()
    return {"chrom": chrom, "ok": True, "metrics": row}


def main() -> None:
    args = parse_args()
    chroms = parse_chrom_list(args.chroms)
    py = args.python
    tag = args.tag

    if not args.skip_concat:
        run(
            [py, "scripts/concat_genomewide.py", "--chroms", ",".join(chroms), "--tag", tag],
            "Concat partial genome-wide tables",
        )

    per_chr_results = []
    if not args.skip_per_chr:
        print("\nPer-chromosome baselines (sequential) ...")
        for chrom in chroms:
            res = run_chr_baseline(py, chrom, tag)
            per_chr_results.append(res)
            status = "ok" if res["ok"] else "FAIL"
            print(f"  chr{res['chrom']}: {status}")
            if not res["ok"] and res.get("error"):
                print(f"    {res['error']}")

        import pandas as pd

        rows = [r["metrics"] for r in per_chr_results if r.get("metrics")]
        if rows:
            summary = pd.DataFrame(rows)
            summary["chrom"] = [r["chrom"] for r in per_chr_results if r.get("metrics")]
            out = ROOT / f"results/tables/ablation_per_chromosome_{tag}.csv"
            summary.to_csv(out, index=False, float_format="%.4f")
            print(f"Saved {out}")
            print(
                f"Per-chr AF_LD_SEL AUROC (associated, variant split): "
                f"mean={summary['AUROC'].mean():.3f}  "
                f"min={summary['AUROC'].min():.3f}  max={summary['AUROC'].max():.3f}"
            )

    if not args.skip_pooled:
        master = ROOT / f"data/modeling/master_variant_table_genomewide_{tag}.parquet"
        groups = ROOT / f"data/modeling/feature_groups_genomewide_{tag}.json"
        suffix = f".genomewide_{tag}"
        common = [
            "--master", str(master),
            "--groups", str(groups),
            "--subsets", "associated",
            "--associated-only",
            "--splits", "split_variant",
        ]
        run(
            [py, "scripts/compute_baselines.py", *common, "--skip-riskclass", "--out-suffix", suffix],
            "Pooled baselines (associated only)",
        )
        run(
            [py, "scripts/evaluate_suite.py", *common, "--n-boot", "50", "--out-suffix", suffix],
            "Pooled evaluation suite (bootstrap CIs)",
        )

    manifest = ROOT / f"data/modeling/genomewide_manifest_{tag}.json"
    if manifest.exists():
        print(f"\nManifest: {manifest.read_text(encoding='utf-8')}")


if __name__ == "__main__":
    main()
