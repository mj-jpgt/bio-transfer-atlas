"""Command-line entry point for the complete methods/statistics stage."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Annotated

import typer

from hostbias.alignment_bridge import run_alignment_bridge
from hostbias.assembly_qc import assembly_qc
from hostbias.config import (
    ValidationError,
    load_and_validate,
    load_thresholds,
)
from hostbias.controls import evaluate_controls
from hostbias.endpoints import calculate_endpoints
from hostbias.endpoint_bridge import aggregate_sample_endpoint
from hostbias.fetch_audit import audit_fetch
from hostbias.labeling import label_contigs
from hostbias.mag_bridge import (
    bins_to_scaffolds2bin,
    build_mag_contracts,
    depth_to_maxbin_abundance,
)
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
from hostbias.operations import (
    OperationError,
    launch_production,
    prepare_production_overlay,
    production_status,
)
from hostbias.reference_acquisition import Minimap2IndexBuilder, build_reference_panel
from hostbias.runtime_manifest import prepare_runtime
from hostbias.sentinel_runner import LocalExecutor, run_sentinel_panel
from hostbias.stage_audit import audit_stage
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


@app.command("reference-build")
def reference_build_command(
    metadata_sources: Annotated[
        Path, typer.Option(exists=True, readable=True)
    ] = Path("config/reference_metadata_sources.tsv"),
    donors: Annotated[Path, typer.Option(exists=True, readable=True)] = Path(
        "config/hprc_balanced_donors.tsv"
    ),
    panel: Annotated[Path, typer.Option(exists=True, readable=True)] = Path(
        "config/competitive_human_panel.tsv"
    ),
    reference_root: Annotated[Path, typer.Option()] = Path("references"),
    checkpoint: Annotated[Path, typer.Option()] = Path(
        "results/aggregate/checkpoints/competitive_human_reference_panel.json"
    ),
    threads: Annotated[int, typer.Option(min=1)] = 32,
    index_batch: Annotated[str, typer.Option()] = "64G",
    minimap2: Annotated[str, typer.Option()] = "minimap2",
) -> None:
    """Acquire and index the verified ancestry-balanced human panel."""

    report = build_reference_panel(
        metadata_sources_path=metadata_sources,
        donors_path=donors,
        panel_path=panel,
        reference_root=reference_root,
        checkpoint_path=checkpoint,
        threads=threads,
        index_batch=index_batch,
        index_builder=Minimap2IndexBuilder(minimap2),
    )
    typer.echo(
        json.dumps(
            {
                "status": report["status"],
                "reference_id": report["reference_id"],
                "checkpoint": str(checkpoint),
            }
        )
    )


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


@app.command("stage-audit")
def stage_audit_command(
    sample_id: Annotated[str, typer.Option()],
    normalized_r1: Annotated[Path, typer.Option(exists=True, readable=True)],
    normalized_r2: Annotated[Path, typer.Option(exists=True, readable=True)],
    source_r1: Annotated[Path, typer.Option(exists=True, readable=True)],
    source_r2: Annotated[Path, typer.Option(exists=True, readable=True)],
    strict_r1: Annotated[Path, typer.Option(exists=True, readable=True)],
    strict_r2: Annotated[Path, typer.Option(exists=True, readable=True)],
    expected_r1_sha256: Annotated[str, typer.Option()],
    expected_r2_sha256: Annotated[str, typer.Option()],
    expected_r1_bytes: Annotated[int, typer.Option(min=1)],
    expected_r2_bytes: Annotated[int, typer.Option(min=1)],
    output: Annotated[Path, typer.Option()],
    expected_pairs: Annotated[int, typer.Option(min=1)] = 8_000_000,
    expected_length: Annotated[int, typer.Option(min=1)] = 100,
) -> None:
    """Audit normalized and GRCh38-filtered pairs without exposing reads."""

    report = audit_stage(
        sample_id=sample_id,
        normalized_r1=normalized_r1,
        normalized_r2=normalized_r2,
        source_r1=source_r1,
        source_r2=source_r2,
        strict_r1=strict_r1,
        strict_r2=strict_r2,
        expected_r1_sha256=expected_r1_sha256,
        expected_r2_sha256=expected_r2_sha256,
        expected_r1_bytes=expected_r1_bytes,
        expected_r2_bytes=expected_r2_bytes,
        expected_pairs=expected_pairs,
        expected_length=expected_length,
    )
    write_json_atomic(report, output)
    typer.echo(json.dumps({"status": report["status"], "sample_id": sample_id}))


@app.command("fetch-audit")
def fetch_audit_command(
    manifest: Annotated[Path, typer.Option(exists=True, readable=True)],
    raw_root: Annotated[Path, typer.Option(exists=True, file_okay=False)],
    output: Annotated[Path, typer.Option()],
    expected_pairs: Annotated[int, typer.Option(min=1)] = 40,
    threads: Annotated[int, typer.Option(min=1)] = 4,
) -> None:
    """Verify atomic FASTQ acquisition and write aggregate-only evidence."""

    report = audit_fetch(
        manifest=manifest,
        raw_root=raw_root,
        expected_pairs=expected_pairs,
        threads=threads,
    )
    write_json_atomic(report, output)
    typer.echo(
        json.dumps(
            {
                "status": report["status"],
                "complete_pairs": report["observed"]["complete_pairs"],
            }
        )
    )


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


@app.command("production-prepare")
def production_prepare_command(
    base_config: Annotated[Path, typer.Option(exists=True, readable=True)],
    nfs_root: Annotated[Path, typer.Option()],
    run_id: Annotated[str, typer.Option()],
    output_config: Annotated[Path | None, typer.Option()] = None,
    evidence: Annotated[Path | None, typer.Option()] = None,
    allow_non_shared_filesystem: Annotated[bool, typer.Option()] = False,
) -> None:
    """Create a private shared-NFS overlay for the 40-sample production run."""

    output = output_config or Path("runtime") / "operations" / f"{run_id}.yaml"
    evidence_output = evidence or (
        Path("results") / "aggregate" / "operations" / f"{run_id}.ready.json"
    )
    try:
        report = prepare_production_overlay(
            base_config=base_config,
            output_config=output,
            evidence_output=evidence_output,
            nfs_root=nfs_root,
            run_id=run_id,
            require_shared_fs=not allow_non_shared_filesystem,
        )
    except (OperationError, ValidationError) as error:
        raise typer.BadParameter(str(error)) from error
    typer.echo(
        json.dumps(
            {
                "status": report["status"],
                "sample_count": report["sample_count"],
                "config": str(output),
                "evidence": str(evidence_output),
            }
        )
    )


@app.command("production-launch")
def production_launch_command(
    config: Annotated[Path, typer.Option(exists=True, readable=True)],
    stage: Annotated[str, typer.Option()],
    evidence: Annotated[Path | None, typer.Option()] = None,
    cores: Annotated[int, typer.Option(min=1)] = 30,
    jobs: Annotated[int, typer.Option(min=1)] = 8,
    mem_mb: Annotated[int, typer.Option(min=1)] = 204_800,
    disk_mb: Annotated[int, typer.Option(min=1)] = 2_500_000,
    latency_wait_seconds: Annotated[int, typer.Option(min=1)] = 120,
    snakemake: Annotated[str, typer.Option()] = "snakemake",
    dry_run: Annotated[bool, typer.Option()] = False,
) -> None:
    """Run or resume one production stage in the foreground."""

    evidence_output = evidence or (
        Path("results")
        / "aggregate"
        / "operations"
        / f"production.{stage}.launch.json"
    )
    try:
        report = launch_production(
            config_path=config,
            stage=stage,
            evidence_output=evidence_output,
            cores=cores,
            jobs=jobs,
            mem_mb=mem_mb,
            disk_mb=disk_mb,
            latency_wait_seconds=latency_wait_seconds,
            snakemake_executable=snakemake,
            dry_run=dry_run,
        )
    except (OperationError, ValidationError) as error:
        raise typer.BadParameter(str(error)) from error
    typer.echo(json.dumps({"status": report["status"], "evidence": str(evidence_output)}))
    if report["exit_code"] != 0:
        raise typer.Exit(report["exit_code"])


@app.command("production-status")
def production_status_command(
    config: Annotated[Path, typer.Option(exists=True, readable=True)],
    output: Annotated[
        Path, typer.Option()
    ] = Path("results/aggregate/operations/production.status.json"),
) -> None:
    """Write aggregate-only completion counts for every production stage."""

    try:
        report = production_status(config)
    except (OperationError, ValidationError) as error:
        raise typer.BadParameter(str(error)) from error
    write_json_atomic(report, output)
    typer.echo(json.dumps({"stages": report["stages"], "output": str(output)}))


@app.command("assembly-qc")
def assembly_qc_command(
    assembly: Annotated[Path, typer.Option(exists=True, readable=True)],
    sample_id: Annotated[str, typer.Option()],
    filter_mode: Annotated[str, typer.Option()],
    output: Annotated[Path, typer.Option()],
) -> None:
    """Write identifier-free assembly QC aggregates."""

    write_json_atomic(assembly_qc(assembly, sample_id, filter_mode), output)
    typer.echo(output)


@app.command("build-alignment-table")
def build_alignment_table_command(
    spec: Annotated[Path, typer.Option(exists=True, readable=True)],
    output_tsv: Annotated[Path, typer.Option()],
    output_manifest: Annotated[Path, typer.Option()],
) -> None:
    """Convert private minimap2 PAFs into the strict competitive-label table."""

    run_alignment_bridge(spec, output_tsv, output_manifest)
    typer.echo(output_tsv)


@app.command("aggregate-endpoint")
def aggregate_endpoint_command(
    alignments: Annotated[Path, typer.Option(exists=True, readable=True)],
    contig_bins: Annotated[Path, typer.Option(exists=True, readable=True)],
    bin_qc: Annotated[Path, typer.Option(exists=True, readable=True)],
    sample_id: Annotated[str, typer.Option()],
    filter_mode: Annotated[str, typer.Option()],
    output: Annotated[Path, typer.Option()],
    thresholds: Annotated[
        Path | None, typer.Option(exists=True, readable=True)
    ] = None,
) -> None:
    """Aggregate bin propagation without publishing sequence-derived identifiers."""

    payload = aggregate_sample_endpoint(
        alignments,
        contig_bins,
        bin_qc,
        sample_id,
        filter_mode,
        thresholds,
    )
    write_json_atomic(payload, output)
    typer.echo(output)


@app.command("maxbin-abundance")
def maxbin_abundance_command(
    depth: Annotated[Path, typer.Option(exists=True, readable=True)],
    output: Annotated[Path, typer.Option()],
) -> None:
    """Translate MetaBAT depth output into MaxBin's abundance input."""

    count = depth_to_maxbin_abundance(depth, output)
    typer.echo(f"{count} contigs")


@app.command("bins-to-map")
def bins_to_map_command(
    bin_dir: Annotated[Path, typer.Option(exists=True, file_okay=False)],
    output: Annotated[Path, typer.Option()],
    bin_prefix: Annotated[str, typer.Option()],
) -> None:
    """Translate a FASTA bin directory into DAS Tool scaffolds-to-bin TSV."""

    count = bins_to_scaffolds2bin(bin_dir, output, bin_prefix)
    typer.echo(f"{count} assignments")


@app.command("mag-contract")
def mag_contract_command(
    sample_id: Annotated[str, typer.Option()],
    dastool_map: Annotated[Path, typer.Option(exists=True, readable=True)],
    checkm2_report: Annotated[Path, typer.Option(exists=True, readable=True)],
    gunc_report: Annotated[Path, typer.Option(exists=True, readable=True)],
    gtdb_bacterial_summary: Annotated[
        Path, typer.Option(exists=True, readable=True)
    ],
    gtdb_archaeal_summary: Annotated[
        Path, typer.Option(exists=True, readable=True)
    ],
    contig_bins_output: Annotated[Path, typer.Option()],
    bin_qc_output: Annotated[Path, typer.Option()],
) -> None:
    """Build exact private endpoint inputs from MAG tool reports."""

    count = build_mag_contracts(
        sample_id,
        dastool_map,
        checkm2_report,
        gunc_report,
        (gtdb_bacterial_summary, gtdb_archaeal_summary),
        contig_bins_output,
        bin_qc_output,
    )
    typer.echo(f"{count} selected bins")


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
