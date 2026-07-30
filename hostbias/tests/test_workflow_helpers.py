from __future__ import annotations

import gzip
import hashlib
import subprocess
import sys
from pathlib import Path


PROJECT = Path(__file__).parents[1]
FETCH = PROJECT / "workflow" / "scripts" / "fetch_pair.py"
VALIDATE = PROJECT / "workflow" / "scripts" / "validate_fastq_pair.py"
SUBSAMPLE = PROJECT / "workflow" / "scripts" / "subsample_fastq_pair.py"


def fastq_bytes(mate: int, names: tuple[str, ...] = ("read-a", "read-b")) -> bytes:
    lines = []
    for name in names:
        lines.extend((f"@{name}/{mate}", "ACGT", "+", "IIII"))
    return gzip.compress(("\n".join(lines) + "\n").encode("ascii"), mtime=0)


def many_fastq_bytes(mate: int, count: int = 20) -> bytes:
    return fastq_bytes(mate, tuple(f"read-{index:02d}" for index in range(count)))


def run_process(*arguments: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, *arguments],
        check=False,
        capture_output=True,
        text=True,
    )


def test_fetch_pair_publishes_only_checksum_verified_mates(tmp_path: Path) -> None:
    source1 = tmp_path / "source1.fastq.gz"
    source2 = tmp_path / "source2.fastq.gz"
    source1.write_bytes(fastq_bytes(1))
    source2.write_bytes(fastq_bytes(2))
    output1 = tmp_path / "out" / "R1.fastq.gz"
    output2 = tmp_path / "out" / "R2.fastq.gz"

    result = run_process(
        str(FETCH),
        "--url1",
        source1.as_uri(),
        "--md5-1",
        hashlib.md5(source1.read_bytes(), usedforsecurity=False).hexdigest(),
        "--output1",
        str(output1),
        "--bytes-1",
        str(source1.stat().st_size),
        "--url2",
        source2.as_uri(),
        "--md5-2",
        hashlib.md5(source2.read_bytes(), usedforsecurity=False).hexdigest(),
        "--output2",
        str(output2),
        "--bytes-2",
        str(source2.stat().st_size),
    )

    assert result.returncode == 0, result.stderr
    assert output1.read_bytes() == source1.read_bytes()
    assert output2.read_bytes() == source2.read_bytes()


def test_fetch_pair_checksum_failure_publishes_neither_mate(tmp_path: Path) -> None:
    source1 = tmp_path / "source1.fastq.gz"
    source2 = tmp_path / "source2.fastq.gz"
    source1.write_bytes(fastq_bytes(1))
    source2.write_bytes(fastq_bytes(2))
    output1 = tmp_path / "out" / "R1.fastq.gz"
    output2 = tmp_path / "out" / "R2.fastq.gz"

    result = run_process(
        str(FETCH),
        "--url1",
        source1.as_uri(),
        "--md5-1",
        hashlib.md5(source1.read_bytes(), usedforsecurity=False).hexdigest(),
        "--output1",
        str(output1),
        "--bytes-1",
        str(source1.stat().st_size),
        "--url2",
        source2.as_uri(),
        "--md5-2",
        "f" * 32,
        "--output2",
        str(output2),
        "--bytes-2",
        str(source2.stat().st_size),
    )

    assert result.returncode != 0
    assert not output1.exists()
    assert not output2.exists()


def test_fetch_pair_size_failure_publishes_neither_mate(tmp_path: Path) -> None:
    source1 = tmp_path / "source1.fastq.gz"
    source2 = tmp_path / "source2.fastq.gz"
    source1.write_bytes(fastq_bytes(1))
    source2.write_bytes(fastq_bytes(2))
    output1 = tmp_path / "out" / "R1.fastq.gz"
    output2 = tmp_path / "out" / "R2.fastq.gz"

    result = run_process(
        str(FETCH),
        "--url1",
        source1.as_uri(),
        "--md5-1",
        hashlib.md5(source1.read_bytes(), usedforsecurity=False).hexdigest(),
        "--bytes-1",
        str(source1.stat().st_size + 1),
        "--output1",
        str(output1),
        "--url2",
        source2.as_uri(),
        "--md5-2",
        hashlib.md5(source2.read_bytes(), usedforsecurity=False).hexdigest(),
        "--bytes-2",
        str(source2.stat().st_size),
        "--output2",
        str(output2),
    )

    assert result.returncode != 0
    assert "size mismatch" in result.stderr
    assert not output1.exists()
    assert not output2.exists()


def test_fastq_pair_validator_accepts_synchronized_exact_length(tmp_path: Path) -> None:
    r1 = tmp_path / "R1.fastq.gz"
    r2 = tmp_path / "R2.fastq.gz"
    r1.write_bytes(fastq_bytes(1))
    r2.write_bytes(fastq_bytes(2))

    result = run_process(
        str(VALIDATE),
        "--r1",
        str(r1),
        "--r2",
        str(r2),
        "--expected-pairs",
        "2",
        "--expected-length",
        "4",
    )

    assert result.returncode == 0, result.stderr
    assert "2 pairs" in result.stdout


def test_fastq_pair_validator_rejects_desynchronization(tmp_path: Path) -> None:
    r1 = tmp_path / "R1.fastq.gz"
    r2 = tmp_path / "R2.fastq.gz"
    r1.write_bytes(fastq_bytes(1))
    r2.write_bytes(fastq_bytes(2, names=("read-a", "read-c")))

    result = run_process(str(VALIDATE), "--r1", str(r1), "--r2", str(r2))

    assert result.returncode != 0
    assert "mate names differ" in result.stderr


def test_pair_subsampling_is_exact_synchronized_and_byte_deterministic(
    tmp_path: Path,
) -> None:
    r1 = tmp_path / "R1.fastq.gz"
    r2 = tmp_path / "R2.fastq.gz"
    r1.write_bytes(many_fastq_bytes(1))
    r2.write_bytes(many_fastq_bytes(2))
    output_sets = []
    for run_name in ("first", "second"):
        output1 = tmp_path / f"{run_name}.R1.fastq.gz"
        output2 = tmp_path / f"{run_name}.R2.fastq.gz"
        result = run_process(
            str(SUBSAMPLE),
            "--r1",
            str(r1),
            "--r2",
            str(r2),
            "--output1",
            str(output1),
            "--output2",
            str(output2),
            "--pairs",
            "7",
            "--seed",
            "20260729",
            "--expected-length",
            "4",
        )
        assert result.returncode == 0, result.stderr
        output_sets.append((output1.read_bytes(), output2.read_bytes()))
    assert output_sets[0] == output_sets[1]

    validation = run_process(
        str(VALIDATE),
        "--r1",
        str(tmp_path / "first.R1.fastq.gz"),
        "--r2",
        str(tmp_path / "first.R2.fastq.gz"),
        "--expected-pairs",
        "7",
        "--expected-length",
        "4",
    )
    assert validation.returncode == 0, validation.stderr


def test_workflow_has_explicit_privacy_and_resource_boundaries() -> None:
    snakefile = (PROJECT / "workflow" / "Snakefile").read_text(encoding="utf-8")

    for rule in (
        "rule fetch_fastq_pair:",
        "rule trim_and_normalize:",
        "rule filter_host:",
        "rule assemble:",
    ):
        assert rule in snakefile
    assert "resources:" in snakefile
    assert "samtools view -u -f 12 -F 2304" in snakefile
    assert "--un-conc-gz" in snakefile
    assert "--max_len1" in snakefile and "--max_len2" in snakefile
    assert "--length_limit" not in snakefile
    assert "subsample_fastq_pair.py" in snakefile
    assert "seqtk sample" not in snakefile
    assert "results/aggregate" not in snakefile
