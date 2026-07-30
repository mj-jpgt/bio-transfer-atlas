from __future__ import annotations

import json
import os
import shutil
import stat
import csv
from pathlib import Path

import pytest
import yaml

from hostbias.config import load_and_validate
from hostbias.operations import (
    OperationError,
    build_snakemake_command,
    launch_production,
    prepare_production_overlay,
    production_status,
)


PROJECT = Path(__file__).parents[1]


def make_production_config(
    name: str, *, work_dir: Path, results_dir: Path
) -> tuple[Path, Path]:
    runtime = PROJECT / "runtime" / name
    runtime.mkdir(parents=True, exist_ok=True)
    manifest = runtime / "samples.tsv"
    fields = (
        "sample_id",
        "cohort",
        "accession",
        "layout",
        "fastq1_url",
        "fastq1_md5",
        "fastq1_bytes",
        "fastq2_url",
        "fastq2_md5",
        "fastq2_bytes",
        "status",
        "rank",
    )
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for index in range(40):
            accession = f"SRR{9000000 + index}"
            writer.writerow(
                {
                    "sample_id": accession,
                    "cohort": "tanzania" if index < 20 else "netherlands",
                    "accession": accession,
                    "layout": "PAIRED",
                    "fastq1_url": f"https://example.org/{accession}_1.fastq.gz",
                    "fastq1_md5": f"{index + 1:032x}",
                    "fastq1_bytes": "100",
                    "fastq2_url": f"https://example.org/{accession}_2.fastq.gz",
                    "fastq2_md5": f"{index + 101:032x}",
                    "fastq2_bytes": "100",
                    "status": "primary",
                    "rank": str(index % 20 + 1),
                }
            )
    config = yaml.safe_load(
        (PROJECT / "config" / "config.example.yaml").read_text(encoding="utf-8")
    )
    config["paths"] = {
        "sample_manifest": manifest.relative_to(PROJECT).as_posix(),
        "work_dir": str(work_dir),
        "results_dir": str(results_dir),
    }
    config_path = runtime / "production.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    return config_path, runtime


def test_overlay_keeps_work_on_private_nfs_and_results_in_project(
    tmp_path: Path,
) -> None:
    nfs = tmp_path / "nfs"
    nfs.mkdir()
    runtime = PROJECT / "runtime" / tmp_path.name
    evidence = tmp_path / "evidence.json"
    output = runtime / "production.yaml"
    try:
        report = prepare_production_overlay(
            base_config=PROJECT / "config" / "config.example.yaml",
            output_config=output,
            evidence_output=evidence,
            nfs_root=nfs,
            run_id="gate-a-test",
            expected_samples=2,
            require_shared_fs=False,
        )

        generated = load_and_validate(output)
        work = Path(generated.config["paths"]["work_dir"])
        assert work.is_relative_to(nfs.resolve())
        if os.name == "posix":
            assert stat.S_IMODE(work.stat().st_mode) == 0o700
        assert Path(generated.config["paths"]["results_dir"]) == PROJECT / "results"
        assert Path(
            generated.config["references"]["grch38"]["fasta"]
        ).is_relative_to(nfs.resolve())
        assert report["storage"]["scratch_layout"] == "scratch/gate-a-test/work"
        serialized = evidence.read_text(encoding="utf-8")
        assert str(nfs.resolve()) not in serialized
        assert json.loads(serialized)["privacy"]["contains_absolute_paths"] is False
    finally:
        shutil.rmtree(runtime, ignore_errors=True)


def test_overlay_rejects_relative_or_missing_nfs_root(tmp_path: Path) -> None:
    output = PROJECT / "runtime" / "missing-nfs-test" / "production.yaml"
    try:
        with pytest.raises(OperationError, match="absolute path"):
            prepare_production_overlay(
                base_config=PROJECT / "config" / "config.example.yaml",
                output_config=output,
                evidence_output=tmp_path / "evidence.json",
                nfs_root=Path("relative"),
                run_id="gate-a-test",
                expected_samples=2,
                require_shared_fs=False,
            )
        with pytest.raises(OperationError, match="already exist"):
            prepare_production_overlay(
                base_config=PROJECT / "config" / "config.example.yaml",
                output_config=output,
                evidence_output=tmp_path / "evidence.json",
                nfs_root=tmp_path / "missing",
                run_id="gate-a-test",
                expected_samples=2,
                require_shared_fs=False,
            )
    finally:
        shutil.rmtree(output.parent, ignore_errors=True)


def test_overlay_rejects_escape_and_oversubscribed_resources(tmp_path: Path) -> None:
    nfs = tmp_path / "nfs"
    nfs.mkdir()
    outside = tmp_path / "outside.yaml"
    with pytest.raises(OperationError, match="output config must remain"):
        prepare_production_overlay(
            base_config=PROJECT / "config" / "config.example.yaml",
            output_config=outside,
            evidence_output=tmp_path / "evidence.json",
            nfs_root=nfs,
            run_id="gate-a-test",
            expected_samples=2,
            require_shared_fs=False,
        )

    config = yaml.safe_load(
        (PROJECT / "config" / "config.example.yaml").read_text(encoding="utf-8")
    )
    config["resources"]["assembly"]["threads"] = 31
    base = PROJECT / "runtime" / tmp_path.name / "oversubscribed.yaml"
    output = base.parent / "production.yaml"
    base.parent.mkdir(parents=True)
    base.write_text(yaml.safe_dump(config), encoding="utf-8")
    try:
        with pytest.raises(OperationError, match="core cap"):
            prepare_production_overlay(
                base_config=base,
                output_config=output,
                evidence_output=tmp_path / "evidence.json",
                nfs_root=nfs,
                run_id="gate-a-test",
                expected_samples=2,
                require_shared_fs=False,
            )
    finally:
        shutil.rmtree(base.parent, ignore_errors=True)


def test_launch_command_has_exact_resume_and_resource_semantics(tmp_path: Path) -> None:
    work = tmp_path / "scratch" / "work"
    work.mkdir(parents=True)
    results = tmp_path / "results"
    config, runtime = make_production_config(
        tmp_path.name + "-launch", work_dir=work, results_dir=results
    )
    observed: list[object] = []

    def runner(command: list[str], cwd: Path) -> int:
        observed.extend((command, cwd))
        return 0

    try:
        evidence = tmp_path / "launch.json"
        report = launch_production(
            config_path=config,
            stage="filter",
            evidence_output=evidence,
            dry_run=True,
            require_clean_git=False,
            runner=runner,
        )

        command = observed[0]
        assert isinstance(command, list)
        assert "--rerun-incomplete" in command
        assert "--keep-going" in command
        assert "--dry-run" in command
        assert command[-1] == "filter_stage"
        assert ["--cores", "30"] == command[
            command.index("--cores") : command.index("--cores") + 2
        ]
        assert report["status"] == "DRY_RUN_COMPLETE"
        serialized = evidence.read_text(encoding="utf-8")
        assert str(tmp_path) not in serialized
        assert json.loads(serialized)["resume"]["same_command_resumes"] is True
    finally:
        shutil.rmtree(runtime, ignore_errors=True)


def test_status_reports_only_aggregate_unit_counts(tmp_path: Path) -> None:
    work = tmp_path / "scratch" / "work"
    work.mkdir(parents=True)
    work.parent.chmod(0o700)
    results = tmp_path / "results"
    config, runtime = make_production_config(
        tmp_path.name + "-status", work_dir=work, results_dir=results
    )
    sample = "SRR9000000"
    raw = work / "raw"
    raw.mkdir()
    (raw / f"{sample}_R1.fastq.gz").touch()
    (raw / f"{sample}_R2.fastq.gz").touch()
    try:
        report = production_status(config)

        assert report["stages"]["fetch"] == {
            "state": "partial",
            "completed_units": 1,
            "expected_units": 40,
            "percent": 2.5,
        }
        assert report["stages"]["assemble"]["expected_units"] == 80
        serialized = json.dumps(report)
        assert sample not in serialized
        assert str(tmp_path) not in serialized
    finally:
        shutil.rmtree(runtime, ignore_errors=True)


def test_scheduler_rejects_vm_oversubscription() -> None:
    with pytest.raises(OperationError, match="30-core"):
        build_snakemake_command(config_path=Path("config.yaml"), stage="fetch", cores=31)
    with pytest.raises(OperationError, match="cap of 8"):
        build_snakemake_command(config_path=Path("config.yaml"), stage="fetch", jobs=9)
    with pytest.raises(OperationError, match="stage must"):
        build_snakemake_command(config_path=Path("config.yaml"), stage="unknown")
