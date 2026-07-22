#!/usr/bin/env python3
from pathlib import Path
import pyarrow.parquet as pq
import pandas as pd

ROOT = Path("/lambda/nfs/geeg/fairness")
pgens = list((ROOT / "data/interim/1000g_grch38").glob("chr*.score.pgen"))
print("score_pgens", sorted(int(p.name.split(".")[0].replace("chr", "")) for p in pgens if p.name[3].isdigit() or p.name[4:5].isdigit() or True))
chroms = []
for p in pgens:
    try:
        chroms.append(int(p.name.split("chr")[1].split(".")[0]))
    except Exception:
        pass
print("score_pgen_chroms", sorted(chroms), "n=", len(chroms))

sub = list((ROOT / "data/features/af").glob("subpop_af_features.chr*.parquet"))
sch = []
for p in sub:
    try:
        sch.append(int(p.name.split("chr")[1].split(".")[0]))
    except Exception:
        pass
print("subpop_chroms", sorted(sch), "n=", len(sch))

susie = list((ROOT / "data/labels/susie").glob("*.parquet"))
print("susie_files", len(susie), [p.name for p in sorted(susie)[:12]])

rg = list((ROOT / "data/features/baselines").glob("rg*"))
print("rg", [p.name for p in rg])

mp = ROOT / "data/modeling/master_variant_table_genomewide_genomewide.parquet"
schema = pq.read_schema(mp)
chr_cols = [n for n in schema.names if n.lower() in ("chrom", "chr", "chromosome")]
print("master_chr_cols", chr_cols)
if chr_cols:
    c = pq.read_table(mp, columns=chr_cols[:1]).to_pandas().iloc[:, 0]
    print("master_chrom_nunique", c.nunique())
    print(c.astype(str).value_counts().sort_index().head(30).to_string())
print("score_matrix", (ROOT / "data/processed/scores_grch38/score_matrix_grch38_genomewide_genomewide.parquet").exists())
print("int_modes", len(list((ROOT / "data/processed/scores_grch38_intervention_genomewide").glob("score_matrix_*_genomewide.parquet"))))

import json
rg_meta = ROOT / "data/features/baselines/rg_real_meta.json"
if rg_meta.exists():
    print("RG_META", json.dumps(json.loads(rg_meta.read_text()), indent=2)[:1000])

# Master chrom span from variant_id
t = pq.read_table(mp, columns=["variant_id"])
ids = t.column(0).to_pandas().astype(str)
chrom = ids.str.split(":", n=1).str[0]
print("master_nrows", len(ids))
print("master_chrom_nunique", chrom.nunique())
print(chrom.value_counts().sort_index().head(30).to_string())
