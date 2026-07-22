#!/usr/bin/env python3
import pandas as pd
from pathlib import Path
p = Path("data/processed/scores_grch38_intervention_genomewide/score_matrix_fst_genomewide.parquet")
df = pd.read_parquet(p)
print(df.shape)
print([c for c in df.columns if str(c).startswith("PGS")])
b = Path("data/processed/scores_grch38/score_matrix_grch38_genomewide_genomewide.parquet")
bb = pd.read_parquet(b)
print("baseline", bb.shape, [c for c in bb.columns if str(c).startswith("PGS")])
