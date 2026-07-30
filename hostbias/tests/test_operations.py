from __future__ import annotations

import json
import os
import shutil
import stat
from pathlib import Path

import pytest
import yaml

from hostbias.config import load_and_validate
from hostbias.operations import OperationError, prepare_production_overlay


PROJECT = Path(__file__).parents[1]


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
