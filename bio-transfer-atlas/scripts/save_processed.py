"""Save Phase 2 outputs: sample_metadata.parquet and ancestry_pcs.parquet"""
import pandas as pd
from pathlib import Path

root = Path(__file__).resolve().parents[1]
processed = root / "data" / "processed" / "1000g"
processed.mkdir(parents=True, exist_ok=True)

# --- Sample metadata ---
panel = pd.read_csv(
    root / "data/raw/1000g/metadata/integrated_call_samples_v3.20130502.ALL.panel",
    sep="\t",
    header=0,
    usecols=["sample", "pop", "super_pop", "gender"],
)
panel = panel.rename(columns={"sample": "sample_id"})
panel.to_parquet(root / "data/processed/sample_metadata.parquet", index=False)
print(f"sample_metadata.parquet: {len(panel)} rows, columns={list(panel.columns)}")
print(panel["super_pop"].value_counts().to_string())

# --- Ancestry PCs ---
eigenvec = pd.read_csv(
    root / "data/processed/1000g/chr22_pca.eigenvec",
    sep="\t",
    header=0,
)
eigenvec = eigenvec.rename(columns={"#IID": "sample_id"})
if "FID" in eigenvec.columns:
    eigenvec = eigenvec.drop(columns=["FID"])

# Merge with population labels
pcs = eigenvec.merge(panel[["sample_id", "pop", "super_pop"]], on="sample_id", how="left")
pcs.to_parquet(root / "data/processed/ancestry_pcs.parquet", index=False)
print(f"\nancestry_pcs.parquet: {len(pcs)} rows, {len(pcs.columns)} columns")
print("PC columns:", [c for c in pcs.columns if c.startswith("PC")])

# --- Eigenvalues (variance explained) ---
eigenval = pd.read_csv(
    root / "data/processed/1000g/chr22_pca.eigenval",
    header=None,
    names=["eigenvalue"],
)
eigenval["variance_explained"] = eigenval["eigenvalue"] / eigenval["eigenvalue"].sum()
eigenval["PC"] = [f"PC{i+1}" for i in range(len(eigenval))]
eigenval.to_parquet(root / "data/processed/1000g/pca_eigenvalues.parquet", index=False)
print(f"\nVariance explained by PC1–5:")
print(eigenval.head(5).to_string(index=False))
