from __future__ import annotations

import io
import json
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml
from typer.testing import CliRunner

from hostbias.cli import app
from hostbias.data_manifest import read_tsv
from hostbias.sentinel_runner import (
    LocalExecutor,
    SentinelExecutionError,
    _selected_sentinels,
    consume_mapped_pair_sam,
    run_sentinel_panel,
    validate_fastq_pair,
)


def write_fastq_pair(
    r1: Path, r2: Path, *, pairs: int = 3, mismatch_at: int | None = None
) -> None:
    first = []
    second = []
    for index in range(1, pairs + 1):
        first.extend([f"@read{index}/1\n", "ACGT\n", "+\n", "IIII\n"])
        mate_index = index + 1 if mismatch_at == index else index
        second.extend([f"@read{mate_index}/2\n", "TGCA\n", "+\n", "IIII\n"])
    r1.write_text("".join(first), encoding="ascii", newline="\n")
    r2.write_text("".join(second), encoding="ascii", newline="\n")


def write_inputs(root: Path) -> tuple[Path, Path, Path]:
    manifest = root / "samples.tsv"
    rows = [
        "run_accession\tarm\tbioproject\tsample_accession\trole\trank\t"
        "library_layout\tinstrument_platform\tlibrary_strategy\tbase_count\tread_count"
    ]
    for arm in ("a", "b"):
        for rank in range(1, 31):
            role = "primary" if rank <= 20 else "reserve"
            rows.append(
                f"{arm.upper()}{rank}\t{arm}\tPRJ{arm.upper()}\tS{arm.upper()}{rank}\t"
                f"{role}\t{rank}\tPAIRED\tILLUMINA\tWGS\t{1000-rank}\t10"
            )
    manifest.write_text("\n".join(rows) + "\n", encoding="utf-8", newline="\n")
    thresholds = root / "thresholds.yaml"
    thresholds.write_text(
        yaml.safe_dump(
            {
                "sentinel": {
                    "runs_per_arm": 3,
                    "streamed_spots_per_run": 3,
                    "minimum_grch38_mapped_pairs_per_million": 300000,
                    "minimum_passing_runs_per_arm": 2,
                }
            }
        ),
        encoding="utf-8",
        newline="\n",
    )
    index = root / "grch38"
    for suffix in (".1.bt2", ".2.bt2"):
        (root / f"grch38{suffix}").write_bytes(suffix.encode())
    return manifest, thresholds, index


class FakeExecutor:
    def __init__(self, *, fail_once: str | None = None) -> None:
        self.fetch_calls: list[str] = []
        self.align_calls: list[str] = []
        self.fail_once = fail_once

    def tool_versions(self) -> dict[str, str]:
        return {"fastq_dump": "mock fastq 1", "bowtie2": "mock bowtie2 1"}

    def fetch(
        self, accession: str, scratch: Path, spots: int, threads: int
    ) -> tuple[Path, Path]:
        assert spots == 3
        assert threads == 2
        self.fetch_calls.append(accession)
        r1 = scratch / f"{accession}_1.fastq"
        r2 = scratch / f"{accession}_2.fastq"
        write_fastq_pair(r1, r2)
        return r1, r2

    def align(
        self,
        r1: Path,
        r2: Path,
        grch38_index: Path,
        threads: int,
        scratch: Path,
    ) -> int:
        accession = r1.name.split("_", 1)[0]
        self.align_calls.append(accession)
        if self.fail_once == accession:
            self.fail_once = None
            raise SentinelExecutionError(
                "alignment",
                "mock bowtie failure",
                return_code=12,
                stderr=b"tool error ACGTACGTACGTACGTACGTACGTACGT\n",
            )
        return 0 if accession.endswith("3") else 1


def test_panel_runs_exact_first_three_and_resumes_from_checkpoints(
    tmp_path: Path,
) -> None:
    manifest, thresholds, index = write_inputs(tmp_path)
    scratch = tmp_path / "scratch"
    output = tmp_path / "aggregate"
    executor = FakeExecutor()
    report = run_sentinel_panel(
        manifest_path=manifest,
        thresholds_path=thresholds,
        grch38_index=index,
        scratch_root=scratch,
        output_dir=output,
        threads=2,
        executor=executor,
    )
    assert report["status"] == "complete"
    assert report["eligible"] is True
    assert executor.fetch_calls == ["A1", "A2", "A3", "B1", "B2", "B3"]
    assert not any(scratch.iterdir())
    assert not list(output.rglob("*.fastq"))
    assert not list(output.rglob("*.sam"))
    assert not list(output.rglob("*.bam"))

    resumed = FakeExecutor()
    second = run_sentinel_panel(
        manifest_path=manifest,
        thresholds_path=thresholds,
        grch38_index=index,
        scratch_root=scratch,
        output_dir=output,
        threads=2,
        executor=resumed,
    )
    assert second["eligible"] is True
    assert resumed.fetch_calls == []
    assert resumed.align_calls == []


def test_failed_run_is_cleaned_diagnosed_and_restartable(tmp_path: Path) -> None:
    manifest, thresholds, index = write_inputs(tmp_path)
    scratch = tmp_path / "scratch"
    output = tmp_path / "aggregate"
    failed = run_sentinel_panel(
        manifest_path=manifest,
        thresholds_path=thresholds,
        grch38_index=index,
        scratch_root=scratch,
        output_dir=output,
        threads=2,
        executor=FakeExecutor(fail_once="A2"),
    )
    assert failed["status"] == "operational_failure"
    assert failed["completed_runs"] == 5
    assert not any(scratch.iterdir())
    diagnostic = json.loads(
        (output / "diagnostics" / "A2.json").read_text(encoding="utf-8")
    )
    assert diagnostic["stage"] == "alignment"
    assert diagnostic["return_code"] == 12
    assert "ACGTACGT" not in diagnostic["sanitized_stderr_tail"]
    assert "<redacted-sequence>" in diagnostic["sanitized_stderr_tail"]
    assert diagnostic["sequence_data_recorded"] is False

    retry = FakeExecutor()
    completed = run_sentinel_panel(
        manifest_path=manifest,
        thresholds_path=thresholds,
        grch38_index=index,
        scratch_root=scratch,
        output_dir=output,
        threads=2,
        executor=retry,
    )
    assert completed["status"] == "complete"
    assert retry.fetch_calls == ["A2"]
    assert not (output / "diagnostics" / "A2.json").exists()


def test_pair_validator_reports_only_record_position(tmp_path: Path) -> None:
    r1 = tmp_path / "r1.fastq"
    r2 = tmp_path / "r2.fastq"
    write_fastq_pair(r1, r2, mismatch_at=2)
    with pytest.raises(SentinelExecutionError, match="record 2") as captured:
        validate_fastq_pair(r1, r2)
    assert "read" not in str(captured.value)


def test_sam_counter_counts_each_mapped_primary_query_once() -> None:
    sam = io.BytesIO(
        b"@HD\tVN:1.6\n"
        b"r1\t99\tchr1\n"
        b"r1\t147\tchr1\n"
        b"r2\t77\t*\n"
        b"r2\t141\t*\n"
        b"r3\t99\tchr1\n"
        b"r3\t355\tchr2\n"
    )
    assert consume_mapped_pair_sam(sam) == 2


def test_local_fetch_limits_spots_and_requires_exact_pair_shape(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    observed_command: list[str] = []

    def fake_run(command: list[str], **_: object) -> SimpleNamespace:
        observed_command.extend(command)
        output = Path(command[command.index("-O") + 1])
        accession = command[-1]
        write_fastq_pair(
            output / f"{accession}_1.fastq", output / f"{accession}_2.fastq"
        )
        return SimpleNamespace(returncode=0, stderr=b"")

    monkeypatch.setattr("hostbias.sentinel_runner.subprocess.run", fake_run)
    scratch = tmp_path / "run"
    scratch.mkdir()
    LocalExecutor().fetch("SRR1", scratch, 1_000_000, 8)
    assert observed_command[observed_command.index("-N") + 1] == "1"
    assert observed_command[observed_command.index("-X") + 1] == "1000000"
    assert "--split-files" in observed_command


def test_local_alignment_streams_sam_without_writing_alignment(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    observed_command: list[str] = []

    class FakeProcess:
        def __init__(self, command: list[str], **kwargs: object) -> None:
            observed_command.extend(command)
            self.stdout = io.BytesIO(b"r1\t99\tchr1\nr1\t147\tchr1\n")
            kwargs["stderr"].write(b"mock bowtie summary\n")  # type: ignore[union-attr]

        def wait(self) -> int:
            return 0

        def kill(self) -> None:
            return None

    monkeypatch.setattr("hostbias.sentinel_runner.subprocess.Popen", FakeProcess)
    mapped = LocalExecutor().align(
        tmp_path / "r1.fastq",
        tmp_path / "r2.fastq",
        tmp_path / "grch38",
        8,
        tmp_path,
    )
    assert mapped == 1
    assert "--very-sensitive-local" in observed_command
    assert "--reorder" in observed_command
    assert "--no-unal" in observed_command
    assert not list(tmp_path.glob("*.sam"))
    assert not list(tmp_path.glob("*.bam"))


def test_cli_exposes_restartable_sentinel_command() -> None:
    result = CliRunner().invoke(app, ["sentinel-run", "--help"])
    assert result.exit_code == 0
    assert "--grch38-index" in result.stdout
    assert "--scratch-root" in result.stdout
    assert "--output-dir" in result.stdout


def test_frozen_manifest_selects_exact_preregistered_sentinels() -> None:
    manifest = read_tsv(Path(__file__).parents[1] / "config" / "stage0_samples.tsv")
    selected = _selected_sentinels(manifest, 3)
    assert [row["run_accession"] for row in selected] == [
        "SRR13348579",
        "SRR13348648",
        "SRR13348823",
        "SRR5127734",
        "SRR5127631",
        "SRR5127544",
    ]
