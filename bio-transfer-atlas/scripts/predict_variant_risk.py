"""
Phase 18.2: Predict per-variant portability risk (chunked for genome-wide master).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from intervention_common import GW_MASTER, GW_MODEL, GW_PREDS, ROOT, SEED

MDIR = ROOT / "data/modeling"
BATCH_SIZE = 500_000


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Score variant×trait portability risk.")
    p.add_argument("--master", default=str(GW_MASTER))
    p.add_argument("--model", default=str(GW_MODEL))
    p.add_argument("--meta", default=str(GW_MODEL.with_suffix(".meta.json")))
    p.add_argument("--out", default=str(GW_PREDS))
    p.add_argument("--batch-size", type=int, default=BATCH_SIZE)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    meta = json.loads(Path(args.meta).read_text())
    feats = meta["features"]
    model = joblib.load(args.model)

    master_path = Path(args.master)
    if not master_path.exists():
        raise FileNotFoundError(master_path)

    pf = pq.ParquetFile(master_path)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.exists():
        out_path.unlink()

    writer: pq.ParquetWriter | None = None
    total = 0
    cols = ["variant_id", "trait"] + feats

    for batch in pf.iter_batches(batch_size=args.batch_size, columns=cols):
        df = batch.to_pandas()
        X = df[feats].to_numpy(dtype=np.float32)
        risk = model.predict_proba(X)[:, 1].astype(np.float32)
        out = df[["variant_id", "trait"]].copy()
        out["predicted_risk"] = risk
        table = pa.Table.from_pandas(out, preserve_index=False)
        if writer is None:
            writer = pq.ParquetWriter(out_path, table.schema)
        writer.write_table(table)
        total += len(out)
        if total % (args.batch_size * 4) == 0:
            print(f"  scored {total:,} rows ...")

    if writer is not None:
        writer.close()
    print(f"Saved {total:,} predictions -> {out_path}")


if __name__ == "__main__":
    main()
