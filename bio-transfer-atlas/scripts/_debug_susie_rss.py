#!/usr/bin/env python3
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from run_polyfun_susie import extract_ld_matrix, load_associated
import subprocess

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
zmap = dict(zip(bdf.variant_id.astype(str), bdf.z_meta.astype(float)))
with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    R, vids, mode = extract_ld_matrix(bdf, td, Path("data/interim/1000g_grch38"), "22")
    z = np.array([zmap[v] for v in vids], dtype=float)
    print("mode", mode, "n", len(vids), "z finite", np.isfinite(z).sum())
    # make R PSD-ish
    R = 0.9 * R + 0.1 * np.eye(len(vids))
    np.fill_diagonal(R, 1.0)
    pd.DataFrame({"z": z}).to_csv(td / "z.csv", index=False)
    pd.DataFrame(R).to_csv(td / "R.csv", index=False, header=False)
    rscript = f"""
suppressPackageStartupMessages(library(susieR))
z <- as.numeric(read.csv('{td.as_posix()}/z.csv')$z)
R <- as.matrix(read.csv('{td.as_posix()}/R.csv', header=FALSE))
print(dim(R)); print(range(R)); print(summary(z))
res <- tryCatch(susie_rss(z=z, R=R, L=10, estimate_residual_variance=TRUE), error=function(e) {{print(e); NULL}})
print(is.null(res))
if (!is.null(res)) print(summary(res$pip))
"""
    (td / "run.R").write_text(rscript)
    r = subprocess.run(["Rscript", str(td / "run.R")], capture_output=True, text=True)
    print("rc", r.returncode)
    print("STDOUT", r.stdout[-1500:])
    print("STDERR", r.stderr[-1500:])
