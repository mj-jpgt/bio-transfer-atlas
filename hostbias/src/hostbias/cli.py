"""Command-line entry point for the complete methods/statistics stage."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Annotated

import typer

from hostbias.config import (
    ValidationError,
    load_and_validate,
    load_thresholds,
)
from hostbias.controls import evaluate_controls
from hostbias.endpoints import calculate_endpoints
from hostbias.labeling import label_contigs
from hostbias.schemas import (
    AlignmentRow,
    BinQcRow,
    ContigBinRow,
    ControlTruthRow,
    SampleGroupRow,
    SensitivityRow,
    read_tsv,
)
from hostbias.statistics import analyze, attach_groups
from hostbias.provenance import build_provenance, write_json_atomic
from hostbias.runtime_manifest import prepare_runtime
from hostbias.sentinel_runner import LocalExecutor, run_sentinel_panel
from hostbias.verdict import make_verdict, write_verdict


app = typer.Typer(
    no_args_is_help=True,
    help="Produce auditable HostBias Gate A aggregate endpoints and verdicts.",
)


@app.callback()
def main() -> None:
    """HostBias workflow-facing analysis commands."""


@app.command("validate")
def validate_command(
    config: Annotated[Path, typer.Option(exists=True, readable=True)],
) -> None:
    """Validate workflow configuration and its sample manifest."""

    try:
        inputs = load_and_validate(config)
    except ValidationError as error:
        raise typer.BadParameter(str(error), param_hint="--config") from error
    typer.echo(
        f"valid: {len(inputs.samples)} samples; "
        f"experiment={inputs.config['experiment']['id']}"
    )


@app.command("provenance")
def provenance_command(
    config: Annotated[Path, typer.Option(exists=True, readable=True)],
    output: Annotated[Path, typer.Option()],
) -> None:
    """Write a privacy-safe reproducibility manifest."""

    try:
        inputs = load_and_validate(config)
    except ValidationError as error:
        raise typer.BadParameter(str(error), param_hint="--config") from error
    write_json_atomic(build_provenance(inputs), output)
    typer.echo(output)


@app.command("sentinel-run")
def sentinel_run_command(
    manifest: Annotated[Path, typer.Option(exists=True, readable=True)],
    thresholds: Annotated[Path, typer.Option(exists=True, readable=True)],
    grch38_index: Annotated[Path, typer.Option()],
    scratch_root: Annotated[Path, typer.Option()],
    output_dir: Annotated[Path, typer.Option()],
    threads: Annotated[int, typer.Option(min=1)] = 16,
    fastq_dump: Annotated[str, typer.Option()] = "fastq-dump",
    bowtie2: Annotated[str, typer.Option()] = "bowtie2",
) -> None:
    """Run or resume the six aggregate-only Stage 0 sentinel checks."""

    report = run_sentinel_panel(
        manifest_path=manifest,
        thresholds_path=thresholds,
        grch38_index=grch38_index,
        scratch_root=scratch_root,
        output_dir=output_dir,
        threads=threads,
        executor=LocalExecutor(fastq_dump=fastq_dump, bowtie2=bowtie2),
    )
    typer.echo(
        json.dumps(
            {
                "status": report["status"],
                "eligible": report["eligible"],
                "report": str(output_dir / "sentinel_eligibility.json"),
            }
        )
    )
    if report["status"] != "complete":
        raise typer.Exit(2)


def _snapshot_specs(values: list[str]) -> dict[str, Path]:
    result: dict[str, Path] = {}
    for value in values:
        arm, separator, path = value.partition("=")
        if not separator or not arm or not path:
            raise typer.BadParameter(
                "snapshot must be ARM=PATH", param_hint="--snapshot"
            )
        if arm in result:
            raise typer.BadParameter(
                f"duplicate snapshot arm {arm!r}", param_hint="--snapshot"
            )
        snapshot_path = Path(path)
        if not snapshot_path.is_file():
            raise typer.BadParameter(
                f"snapshot does not exist: {snapshot_path}", param_hint="--snapshot"
            )
        result[arm] = snapshot_path
    return result


@app.command("prepare-runtime")
def prepare_runtime_command(
    snapshot: Annotated[
        list[str],
        typer.Option(help="Canonical runtime ENA snapshot as ARM=PATH; repeat per arm."),
    ],
    scope: Annotated[str, typer.Option(help="Either sentinel or primary.")] = "sentinel",
    frozen_manifest: Annotated[
        Path, typer.Option(exists=True, readable=True)
    ] = Path("config/stage0_samples.tsv"),
    config_template: Annotated[
        Path, typer.Option(exists=True, readable=True)
    ] = Path("config/config.example.yaml"),
    runtime_dir: Annotated[Path, typer.Option()] = Path("runtime"),
    evidence: Annotated[Path | None, typer.Option()] = None,
) -> None:
    """Create checksum/size-pinned runtime inputs and aggregate-safe evidence."""

    project_root = Path.cwd().resolve()
    evidence_output = evidence or (
        Path("results")
        / "aggregate"
        / "checkpoints"
        / f"runtime_manifest_{scope}.json"
    )
    try:
        report = prepare_runtime(
            frozen_manifest=frozen_manifest,
            snapshot_paths=_snapshot_specs(snapshot),
            config_template=config_template,
            runtime_manifest=runtime_dir / f"stage0_samples.{scope}.tsv",
            runtime_config=runtime_dir / f"config.{scope}.yaml",
            evidence_output=evidence_output,
            project_root=project_root,
            scope=scope,
        )
    except (ValidationError, ValueError) as error:
        raise typer.BadParameter(str(error)) from error
    typer.echo(
        f"valid: {len(report['ordered_accessions'])} {scope} runs; "
        f"config={runtime_dir / f'config.{scope}.yaml'}"
    )


def _write_json(path: Path, value: object) -> None:
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def run_analysis(
    alignments_path: Path,
    contig_bins_path: Path,
    bin_qc_path: Path,
    sample_groups_path: Path,
    control_alignments_path: Path,
    control_truth_path: Path,
    sensitivities_path: Path,
    output_dir: Path,
    thresholds_path: Path | None = None,
) -> str:
    """Run and persist the complete deterministic methods-stage analysis."""

    label_thresholds, endpoint_thresholds, gate_thresholds, effective = (
        load_thresholds(thresholds_path)
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    alignments = read_tsv(alignments_path, AlignmentRow)
    contig_bins = read_tsv(contig_bins_path, ContigBinRow)
    bin_qc = read_tsv(bin_qc_path, BinQcRow)
    groups = read_tsv(sample_groups_path, SampleGroupRow)
    sensitivities = read_tsv(sensitivities_path, SensitivityRow)
    control_alignments = read_tsv(control_alignments_path, AlignmentRow)
    control_truth = read_tsv(control_truth_path, ControlTruthRow)

    calls = label_contigs(alignments, label_thresholds)
    endpoints, bin_fractions = calculate_endpoints(
        calls,
        contig_bins,
        bin_qc,
        endpoint_thresholds,
        expected_samples=[row.sample_id for row in groups],
    )
    control_calls = label_contigs(control_alignments, label_thresholds)
    controls = evaluate_controls(
        control_calls,
        control_truth,
        min_sensitivity=effective["controls"]["min_sensitivity"],
        max_false_positive_bp_rate=effective["controls"][
            "max_false_positive_bp_rate"
        ],
    )
    grouped = attach_groups(endpoints, groups)
    statistics = analyze(
        grouped,
        bootstrap_iterations=effective["statistics"]["bootstrap_iterations"],
        permutation_iterations=effective["statistics"]["permutation_iterations"],
        seed=effective["statistics"]["seed"],
    )
    verdict = make_verdict(statistics, controls, sensitivities, gate_thresholds)

    _write_json(output_dir / "effective_thresholds.json", effective)
    _write_json(output_dir / "contig_calls.json", [row.to_dict() for row in calls])
    _write_json(
        output_dir / "sample_endpoints.json", [row.to_dict() for row in endpoints]
    )
    _write_json(
        output_dir / "bin_human_fractions.json",
        [row.to_dict() for row in bin_fractions],
    )
    _write_json(output_dir / "control_results.json", controls.to_dict())
    _write_json(output_dir / "analysis_statistics.json", statistics.to_dict())
    write_verdict(verdict, output_dir)
    return verdict.status


@app.command("analyze")
def analyze_command(
    alignments: Annotated[Path, typer.Option(exists=True, readable=True)],
    contig_bins: Annotated[Path, typer.Option(exists=True, readable=True)],
    bin_qc: Annotated[Path, typer.Option(exists=True, readable=True)],
    sample_groups: Annotated[Path, typer.Option(exists=True, readable=True)],
    control_alignments: Annotated[Path, typer.Option(exists=True, readable=True)],
    control_truth: Annotated[Path, typer.Option(exists=True, readable=True)],
    sensitivities: Annotated[Path, typer.Option(exists=True, readable=True)],
    output_dir: Annotated[Path, typer.Option()],
    thresholds: Annotated[
        Path | None, typer.Option(exists=True, readable=True)
    ] = None,
) -> None:
    """Compute labels, endpoints, controls, statistics, and final verdict."""

    status = run_analysis(
        alignments,
        contig_bins,
        bin_qc,
        sample_groups,
        control_alignments,
        control_truth,
        sensitivities,
        output_dir,
        thresholds,
    )
    typer.echo(status)


if __name__ == "__main__":
    app()
