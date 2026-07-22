import pandas as pd
from pathlib import Path

p = Path("data/processed/scores_grch38/score_matrix_grch38_genomewide_genomewide.parquet")
df = pd.read_parquet(p)
print("shape", df.shape)
print("pgs", [c for c in df.columns if c.startswith("PGS")])
