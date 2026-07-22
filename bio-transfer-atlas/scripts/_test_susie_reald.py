#!/usr/bin/env python3
import tempfile
from pathlib import Path

import pandas as pd

from run_polyfun_susie import load_associated, run_susie_block

labels = load_associated(
    Path("data/modeling/master_variant_table_genomewide_genomewide.parquet"), "T2D", "22"
)
blocks = pd.read_parquet(
    "data/modeling/ld_block_assignments_genomewide.parquet",
    columns=["variant_id", "ld_block"],
)
df = labels.merge(blocks, on="variant_id").dropna(subset=["z_meta"])
bname = df.groupby("ld_block").size().sort_values(ascending=False).index[0]
bdf = df[df.ld_block == bname].head(30)
with tempfile.TemporaryDirectory() as td:
    out = run_susie_block(
        bdf.reset_index(drop=True),
        10,
        Path(td),
        prefer_real_ld=True,
        pfile_root=Path("data/interim/1000g_grch38"),
        chrom="22",
    )
    print(out["fallback_mode"].value_counts().to_dict())
    print(out.head(2).to_string())
