"""
Concat per-chromosome parquets into genome-wide tables for a chosen chromosome set.

Memory-efficient: streams one chromosome at a time; reassigns global split_variant
across all unique variants in the merged table.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[1]
SEED = 719


def parse_chrom_list(spec: str) -> list[str]:
    return [str(int(x.strip())) for x in spec.split(",") if x.strip()]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--chroms", required=True, help="Comma list, e.g. 1,2,3,4,5,22")
    p.add_argument("--tag", default="", help="Optional suffix for output names, e.g. partial6")
    p.add_argument(
        "--skip-master",
        action="store_true",
        help="Skip master concat when output already exists (labels/pathways/meta only)",
    )
    return p.parse_args()


def chrom_sort_key(path: Path) -> int:
    stem = path.stem
    for part in stem.split("."):
        if part.startswith("chr") and part[3:].isdigit():
            return int(part[3:])
    return 0


def pick_files(parent: Path, pattern: str, chroms: list[str]) -> list[Path]:
    wanted = {f"chr{c}" for c in chroms}
    files = []
    for p in parent.glob(pattern):
        if any(f".{w}." in p.name or p.name.endswith(f".{w}.parquet") for w in wanted):
            files.append(p)
        elif any(p.name.endswith(f".chr{c}.parquet") for c in chroms):
            files.append(p)
    return sorted(set(files), key=chrom_sort_key)


def assign_global_variant_split(variant_ids: np.ndarray) -> dict[str, str]:
    rng = np.random.default_rng(SEED)
    uniq = np.array(sorted(set(variant_ids)), dtype=object)
    perm = rng.permutation(len(uniq))
    n_tr = int(0.70 * len(uniq))
    n_va = int(0.15 * len(uniq))
    folds = np.full(len(uniq), "test", dtype=object)
    folds[perm[:n_tr]] = "train"
    folds[perm[n_tr : n_tr + n_va]] = "val"
    return dict(zip(uniq, folds))


def master_column_order(files: list[Path]) -> list[str]:
    cols: list[str] = []
    seen: set[str] = set()
    for f in files:
        for name in pq.ParquetFile(f).schema.names:
            if name not in seen:
                cols.append(name)
                seen.add(name)
    return cols


def normalize_master_df(df: pd.DataFrame, column_order: list[str]) -> pd.DataFrame:
    for col in column_order:
        if col not in df.columns:
            df[col] = np.nan
    df = df[column_order].copy()
    for col in df.columns:
        if pd.api.types.is_float_dtype(df[col]):
            df[col] = df[col].astype(np.float32)
        elif isinstance(df[col].dtype, pd.CategoricalDtype):
            df[col] = df[col].astype(str)
    return df


def stream_concat_master(
    files: list[Path], out_path: Path, split_map: dict[str, str], column_order: list[str]
) -> int:
    if out_path.exists():
        out_path.unlink()
    writer = None
    n_rows = 0
    batch_size = 250_000
    for f in files:
        pf = pq.ParquetFile(f)
        chr_rows = 0
        for batch in pf.iter_batches(batch_size=batch_size):
            df = normalize_master_df(batch.to_pandas(), column_order)
            df["split_variant"] = df["variant_id"].map(split_map).astype(str)
            table = pa.Table.from_pandas(df, preserve_index=False)
            if writer is None:
                writer = pq.ParquetWriter(str(out_path), table.schema)
            else:
                table = table.cast(writer.schema, safe=False)
            writer.write_table(table)
            chr_rows += len(df)
            n_rows += len(df)
            del df, table
        print(f"  merged {f.name}: {chr_rows:,} rows")
    if writer:
        writer.close()
    return n_rows


def stream_concat_plain(files: list[Path], out_path: Path) -> int:
    if out_path.exists():
        out_path.unlink()
    writer = None
    n_rows = 0
    batch_size = 250_000
    for f in files:
        pf = pq.ParquetFile(f)
        chr_rows = 0
        for batch in pf.iter_batches(batch_size=batch_size):
            table = pa.Table.from_batches([batch])
            if writer is None:
                writer = pq.ParquetWriter(str(out_path), table.schema)
            else:
                table = table.cast(writer.schema, safe=False)
            writer.write_table(table)
            chr_rows += table.num_rows
            n_rows += table.num_rows
        print(f"  merged {f.name}: {chr_rows:,} rows")
    if writer:
        writer.close()
    return n_rows


def collect_variant_ids(files: list[Path]) -> np.ndarray:
    uniq: set[str] = set()
    for f in files:
        pf = pq.ParquetFile(f)
        for batch in pf.iter_batches(batch_size=500_000, columns=["variant_id"]):
            uniq.update(batch.column("variant_id").to_pylist())
    return np.array(sorted(uniq), dtype=object)


def copy_feature_groups(chroms: list[str], out_path: Path) -> None:
    src = ROOT / "data/modeling" / f"feature_groups.chr{chroms[0]}.json"
    if not src.exists():
        src = ROOT / "data/modeling/feature_groups.json"
    out_path.write_text(src.read_text(encoding="utf-8"), encoding="utf-8")


def main() -> None:
    args = parse_args()
    chroms = parse_chrom_list(args.chroms)
    suffix = f"_{args.tag}" if args.tag else ""

    label_dir = ROOT / "data/labels"
    model_dir = ROOT / "data/modeling"
    result_dir = ROOT / "results/tables"

    master_files = [model_dir / f"master_variant_table.chr{c}.parquet" for c in chroms]
    missing = [str(p) for p in master_files if not p.exists()]
    if missing:
        raise FileNotFoundError(f"Missing master tables: {missing}")

    print(f"Collecting variant IDs from {len(master_files)} chromosomes ...")
    column_order = master_column_order(master_files)
    variant_ids = collect_variant_ids(master_files)
    split_map = assign_global_variant_split(variant_ids)
    print(f"  {len(split_map):,} unique variants -> global split_variant")

    master_out = model_dir / f"master_variant_table_genomewide{suffix}.parquet"
    if args.skip_master and master_out.exists():
        n_master = pq.read_metadata(master_out).num_rows
        print(f"Skipping master concat; reusing {master_out.name} ({n_master:,} rows)")
    else:
        print(f"Writing {master_out.name} ({len(column_order)} columns) ...")
        n_master = stream_concat_master(master_files, master_out, split_map, column_order)
        print(f"Saved {master_out} ({n_master:,} rows)")

    labels_files = [label_dir / f"gwas_concordance_labels_multisource.chr{c}.parquet" for c in chroms]
    labels_out = label_dir / f"gwas_concordance_labels_multisource_genomewide{suffix}.parquet"
    print(f"Writing {labels_out.name} ...")
    n_labels = stream_concat_plain([f for f in labels_files if f.exists()], labels_out)
    print(f"Saved {labels_out} ({n_labels:,} rows)")

    pathway_files = pick_files(result_dir, "pathway_risk_table.chr*.parquet", chroms)
    if pathway_files:
        pathway_out = result_dir / f"pathway_risk_table_genomewide{suffix}.parquet"
        print(f"Writing {pathway_out.name} ...")
        n_path = stream_concat_plain(pathway_files, pathway_out)
        print(f"Saved {pathway_out} ({n_path:,} rows)")

    fg_out = model_dir / f"feature_groups_genomewide{suffix}.json"
    copy_feature_groups(chroms, fg_out)
    print(f"Saved {fg_out}")

    meta = {
        "chromosomes": chroms,
        "n_master_rows": n_master,
        "n_unique_variants": len(split_map),
        "n_label_rows": n_labels,
        "master_path": str(master_out.relative_to(ROOT)),
        "labels_path": str(labels_out.relative_to(ROOT)),
        "feature_groups_path": str(fg_out.relative_to(ROOT)),
    }
    meta_path = model_dir / f"genomewide_manifest{suffix}.json"
    meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(f"Saved {meta_path}")
    print("Done.")


if __name__ == "__main__":
    main()
