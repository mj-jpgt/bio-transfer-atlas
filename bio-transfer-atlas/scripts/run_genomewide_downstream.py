"""
Orchestrate genome-wide / partial downstream steps after per-chr masters exist.

Usage:
  python scripts/run_genomewide_downstream.py --step all --tag partial16 \\
    --chroms 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,22 --score-chroms 8,9,10,11,12,13,14,15
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CHROMS = ",".join(str(i) for i in range(1, 23))
TAG = "genomewide"
PY = sys.executable


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run genome-wide / partial downstream pipeline.")
    p.add_argument("--step", default="all")
    p.add_argument("--chroms", default=CHROMS)
    p.add_argument("--score-chroms", default="", help="Chroms with score.pgen for PLINK scoring")
    p.add_argument("--tag", default=TAG)
    p.add_argument("--memory-mb", type=int, default=640)
    p.add_argument("--skip-per-chr", action="store_true", help="Skip per-chr baselines in eval")
    return p.parse_args()


def run(cmd: list[str], desc: str) -> None:
    print(f"\n[{desc}]", flush=True)
    print("  " + " ".join(cmd), flush=True)
    r = subprocess.run(cmd, cwd=ROOT)
    if r.returncode != 0:
        raise SystemExit(f"FAILED: {desc}")


def main() -> None:
    args = parse_args()
    tag = args.tag
    steps = {s.strip() for s in args.step.split(",") if s.strip()}
    if "all" in steps:
        steps = {"concat", "eval", "pathway", "score", "atlas", "intervention", "gates"}

    master = ROOT / f"data/modeling/master_variant_table_genomewide_{tag}.parquet"
    groups = ROOT / f"data/modeling/feature_groups_genomewide_{tag}.json"
    pathway_table = ROOT / f"results/tables/pathway_risk_table_genomewide_{tag}.parquet"
    score_matrix = ROOT / f"data/processed/scores_grch38/score_matrix_grch38_genomewide_{tag}.parquet"
    model = ROOT / f"data/modeling/portability_model.{tag}.joblib"
    preds = ROOT / f"data/modeling/variant_portability_predictions.{tag}.parquet"
    int_root = ROOT / f"data/processed/pgs_grch38_intervention_{tag}"
    int_scores = ROOT / f"data/processed/scores_grch38_intervention_{tag}"

    score_chroms = args.score_chroms.strip() or args.chroms

    if "concat" in steps:
        run(
            [PY, "scripts/concat_genomewide.py", "--chroms", args.chroms, "--tag", tag],
            "Concat tables",
        )

    if "eval" in steps:
        cmd = [PY, "scripts/run_partial_eval.py", "--chroms", args.chroms, "--tag", tag]
        if args.skip_per_chr:
            cmd.append("--skip-per-chr")
        run(cmd, "Pooled evaluation")

    if "pathway" in steps:
        run(
            [
                PY, "scripts/publish_pathway_risk_top.py",
                "--pathway-table", str(pathway_table),
                "--out-csv", str(ROOT / f"results/tables/pathway_risk_top_{tag}.csv"),
                "--out-summary", str(ROOT / f"results/tables/pathway_risk_summary_{tag}.txt"),
            ],
            "Pathway top summaries",
        )

    if "score" in steps:
        run(
            [
                PY, "scripts/score_genomewide.py",
                "--chroms", score_chroms,
                "--memory-mb", str(args.memory_mb),
                "--out", str(score_matrix),
            ],
            "PGS scoring",
        )

    if "atlas" in steps:
        run(
            [PY, "scripts/compute_score_shifts.py", "--score-matrix", str(score_matrix), "--tag", tag],
            "Score-shift atlas tables",
        )
        run(
            [
                PY, "scripts/compute_distance_sensitivity.py",
                "--score-matrix", str(score_matrix), "--tag", tag,
            ],
            "Distance sensitivity",
        )
        run([PY, "scripts/make_atlas_figures.py", "--tag", tag], "Atlas figures")

    if "intervention" in steps:
        run(
            [
                PY, "scripts/train_portability_model.py",
                "--master", str(master),
                "--groups", str(groups),
                "--out", str(model),
            ],
            "train_portability_model",
        )
        run(
            [
                PY, "scripts/predict_variant_risk.py",
                "--master", str(master),
                "--model", str(model),
                "--meta", str(model.with_suffix(".meta.json")),
                "--out", str(preds),
            ],
            "predict_variant_risk",
        )
        run(
            [
                PY, "scripts/apply_intervention.py",
                "--master", str(master),
                "--predictions", str(preds),
                "--out-root", str(int_root),
            ],
            "apply_intervention",
        )
        # PLINK --score is CPU/RAM bound; parallelize chrom×mode×PGS jobs.
        # Cap per-job memory so 8 concurrent jobs fit in ~216GB.
        score_jobs = 8
        score_mem = min(args.memory_mb, 8192)
        run(
            [
                PY, "scripts/score_intervention.py",
                "--chroms", score_chroms,
                "--intervention-root", str(int_root),
                "--out-dir", str(int_scores),
                "--memory-mb", str(score_mem),
                "--jobs", str(score_jobs),
                "--threads", "4",
                "--tag", tag,
            ],
            "score_intervention",
        )
        run(
            [
                PY, "scripts/evaluate_intervention.py",
                "--baseline", str(score_matrix),
                "--labels", str(master),
                "--intervention-root", str(int_root),
                "--intervention-scores", str(int_scores),
                "--tag", tag,
                "--out-csv", str(ROOT / f"results/tables/intervention_results.{tag}.csv"),
                "--out-summary", str(ROOT / f"results/tables/intervention_summary.{tag}.txt"),
            ],
            "evaluate_intervention",
        )
        run(
            [
                PY, "scripts/per_patient_confidence.py",
                "--chroms", score_chroms,
                "--predictions", str(preds),
                "--out", str(ROOT / f"results/tables/per_patient_confidence.{tag}.parquet"),
            ],
            "per_patient_confidence",
        )

    if "gates" in steps:
        if tag == "genomewide":
            run([PY, "tests/gate_genomewide.py"], "Gate genome-wide eval")
            run([PY, "tests/gate_intervention_genomewide.py"], "Gate genome-wide intervention")
        else:
            gate_eval = ROOT / f"tests/gate_{tag}.py"
            gate_int = ROOT / f"tests/gate_intervention_{tag}.py"
            if gate_eval.exists():
                run([PY, str(gate_eval.relative_to(ROOT))], f"Gate {tag} eval")
            if gate_int.exists():
                run([PY, str(gate_int.relative_to(ROOT))], f"Gate {tag} intervention")


if __name__ == "__main__":
    main()
