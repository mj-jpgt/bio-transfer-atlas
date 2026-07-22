#!/usr/bin/env python3
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from run_polyfun_susie import _plink2, extract_ld_matrix, load_associated

print("plink", _plink2())
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
print("vids", bdf.variant_id.head(3).tolist())
with tempfile.TemporaryDirectory() as td:
    R, vids, mode = extract_ld_matrix(
        bdf, Path(td), Path("data/interim/1000g_grch38"), "22"
    )
    print("mode", mode, "shape", R.shape, "max", float(np.nanmax(R)), "offdiag", float(R[0, 1]))
    print("files", list(Path(td).iterdir()))
