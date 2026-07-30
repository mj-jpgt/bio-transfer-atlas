from __future__ import annotations

import gzip
import hashlib
import subprocess
import sys
from pathlib import Path


PROJECT = Path(__file__).parents[1]
FETCH = PROJECT / "workflow" / "scripts" / "fetch_pair.py"
VALIDATE = PROJECT / "workflow" / "scripts" / "validate_fastq_pair.py"


def fastq_bytes(mate: int, names: tuple[str, ...] = ("read-a", "read-b")) -> bytes:
    lines = []
    for name in names:
        lines.extend((f"@{name}/{mate}", "ACGT", "+", "IIII"))
    return gzip.compress(("\n".join(lines) + "\n").encode("ascii"), mtime=0)


def run(*arguments: str) -> subprocess.CompletedProcess[str]:
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

    result = run(
        str(FETCH),
        "--url1",
        source1.as_uri(),
        "--md5-1",
        hashlib.md5(source1.read_bytes(), usedforsecurity=False).hexdigest(),
        "--output1",
        str(output1),
        "--url2",
        source2.as_uri(),
        "--md5-2",
        hashlib.md5(source2.read_bytes(), usedforsecurity=False).hexdigest(),
        "--output2",
        str(output2),
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

    result = run(
        str(FETCH),
        "--url1",
        source1.as_uri(),
        "--md5-1",
        hashlib.md5(source1.read_bytes(), usedforsecurity=False).hexdigest(),
        "--output1",
        str(output1),
        "--url2",
        source2.as_uri(),
        "--md5-2",
        "f" * 32,
        "--output2",
        str(output2),
    )

    assert result.returncode != 0
    assert not output1.exists()
    assert not output2.exists()


def test_fastq_pair_validator_accepts_synchronized_exact_length(tmp_path: Path) -> None:
    r1 = tmp_path / "R1.fastq.gz"
    r2 = tmp_path / "R2.fastq.gz"
    r1.write_bytes(fastq_bytes(1))
    r2.write_bytes(fastq_bytes(2))

    result = run(
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

    result = run(str(VALIDATE), "--r1", str(r1), "--r2", str(r2))

    assert result.returncode != 0
    assert "mate names differ" in result.stderr


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
    assert "results/aggregate" not in snakefile
