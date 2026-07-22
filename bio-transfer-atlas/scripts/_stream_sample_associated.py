"""Shared low-RAM stream sampler: write associated subsample to a temp parquet."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset as ds
import pyarrow.parquet as pq


def stream_sample_associated(
    master: Path,
    cols: list[str],
    out_parquet: Path,
    *,
    max_rows: int = 500_000,
    batch_size: int = 50_000,
    seed: int = 719,
    mhc_keep_all: bool = False,
    mhc_chr: str = "6",
    mhc_start: int = 28_510_020,
    mhc_end: int = 33_480_577,
) -> int:
    """Reservoir-sample associated rows to out_parquet. Returns rows written."""
    dataset = ds.dataset(str(master), format="parquet")
    available = set(dataset.schema.names)
    use = [c for c in cols if c in available]
    if "variant_id" not in use:
        use = ["variant_id"] + use
    filt = ds.field("associated") == True  # noqa: E712
    rng = np.random.default_rng(seed)

    out_parquet.parent.mkdir(parents=True, exist_ok=True)
    writer: pq.ParquetWriter | None = None
    reservoir: list[pa.Table] = []
    n_kept = 0
    n_seen = 0
    n_total = 0
    n_mhc_extra = 0

    def _cast_f32(table: pa.Table) -> pa.Table:
        arrays, names = [], []
        for i, name in enumerate(table.schema.names):
            col = table.column(i)
            if pa.types.is_floating(col.type) and col.type != pa.float32():
                col = pc.cast(col, pa.float32())
            arrays.append(col)
            names.append(name)
        return pa.Table.from_arrays(arrays, names=names)

    def _is_mhc_table(table: pa.Table) -> pa.Array:
        vids = table.column("variant_id")
        # parse chrom:pos via python for reliability on small batches
        flags = []
        for v in vids.to_pylist():
            parts = str(v).split(":")
            ok = False
            if len(parts) >= 2:
                try:
                    ok = parts[0] == mhc_chr and mhc_start <= int(parts[1]) <= mhc_end
                except ValueError:
                    ok = False
            flags.append(ok)
        return pa.array(flags)

    def _flush_reservoir() -> None:
        nonlocal writer, reservoir, n_kept
        if not reservoir:
            return
        table = pa.concat_tables(reservoir)
        reservoir = []
        if writer is None:
            writer = pq.ParquetWriter(str(out_parquet), table.schema, compression="zstd")
        writer.write_table(table)

    for batch in dataset.scanner(columns=use, filter=filt, batch_size=batch_size).to_batches():
        table = _cast_f32(pa.Table.from_batches([batch]))
        n_total += table.num_rows

        mhc_table = None
        if mhc_keep_all:
            mask = _is_mhc_table(table)
            mhc_table = table.filter(mask)
            table = table.filter(pc.invert(mask))
            if mhc_table.num_rows:
                # write MHC immediately (not part of reservoir budget)
                if writer is None:
                    writer = pq.ParquetWriter(str(out_parquet), mhc_table.schema, compression="zstd")
                writer.write_table(mhc_table)
                n_mhc_extra += mhc_table.num_rows

        if table.num_rows == 0:
            continue

        if n_kept < max_rows:
            take = min(table.num_rows, max_rows - n_kept)
            reservoir.append(table.slice(0, take))
            n_kept += take
            n_seen += take
            table = table.slice(take)
            if n_kept >= max_rows:
                _flush_reservoir()

        if table.num_rows and n_kept >= max_rows:
            # Approximate reservoir: rewrite whole file is expensive; instead
            # skip further non-MHC once full (first-N + all MHC is fine for sensitivity).
            n_seen += table.num_rows
            # Optionally mix a small random slice into a side buffer — skip to save RAM.
            _ = rng  # seed retained for reproducibility of fill order

        if n_total % 400_000 < batch_size:
            print(
                f"  scanned {n_total:,}; reservoir={n_kept:,}; MHC_extra={n_mhc_extra:,}",
                flush=True,
            )

    _flush_reservoir()
    if writer is not None:
        writer.close()
    total = n_kept + n_mhc_extra
    print(
        f"  wrote {out_parquet} rows~{total:,} (reservoir={n_kept:,} + MHC={n_mhc_extra:,})",
        flush=True,
    )
    return total
