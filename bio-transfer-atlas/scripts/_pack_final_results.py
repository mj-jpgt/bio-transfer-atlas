import tarfile
from pathlib import Path

root = Path("/lambda/nfs/geeg/fairness")
out = Path("/tmp/bta_final_m1m3.tgz")
paths = [
    "results/tables/ablation_ldblock_and_baselines_genomewide.csv",
    "results/tables/ablation_xgboost_gpu_companion.csv",
    "results/tables/vep_af_interaction_eval.csv",
    "results/tables/popcorn_rg_summary.csv",
    "results/tables/ldsc_rg_companion.csv",
    "results/tables/subpop_af_features_status.csv",
    "results/tables/shap_mechanism_attribution_genomewide.csv",
    "results/tables/score_pgen_chr1_7_rebuild_status.csv",
    "data/annotations/alphamissense_grch38.parquet",
    "data/features/selection/vep_af_interaction_features.parquet",
    "data/features/af/subpop_af_features.chr22.parquet",
    "data/modeling/feature_groups_genomewide_genomewide.json",
    "data/labels/finemap_tiers_genomewide.parquet",
]
with tarfile.open(out, "w:gz") as t:
    for p in paths:
        fp = root / p
        if fp.exists():
            t.add(fp, arcname=p)
    for fp in (root / "data/features/baselines").glob("*"):
        if fp.is_file():
            t.add(fp, arcname=str(fp.relative_to(root)))
    for fp in (root / "data/labels/susie").glob("*"):
        if fp.is_file():
            t.add(fp, arcname=str(fp.relative_to(root)))
print("wrote", out, out.stat().st_size)
