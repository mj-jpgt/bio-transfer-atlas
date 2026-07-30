from __future__ import annotations

import gzip
import hashlib
import importlib.util
import subprocess
import sys
from pathlib import Path

import yaml


PROJECT = Path(__file__).parents[1]
FETCH = PROJECT / "workflow" / "scripts" / "fetch_pair.py"
VALIDATE = PROJECT / "workflow" / "scripts" / "validate_fastq_pair.py"
SUBSAMPLE = PROJECT / "workflow" / "scripts" / "subsample_fastq_pair.py"
MAG_RULES = PROJECT / "workflow" / "rules" / "mag.smk"


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


def load_fetch_module():
    spec = importlib.util.spec_from_file_location("fetch_pair", FETCH)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


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


def test_fetch_pair_converts_ena_ftp_to_https() -> None:
    module = load_fetch_module()

    assert module.ena_https_url(
        "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR000/001/example.fastq.gz"
    ) == "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR000/001/example.fastq.gz"
    assert module.ena_https_url("ftp://example.org/file") == "ftp://example.org/file"


def test_fetch_pair_resumes_a_persistent_partial_file(tmp_path: Path, monkeypatch) -> None:
    module = load_fetch_module()
    monkeypatch.setattr(module.shutil, "which", lambda executable: None)
    payload = b"0123456789"
    destination = tmp_path / ".R1.fastq.gz.partial"
    destination.write_bytes(payload[:4])
    requests = []

    class Response:
        status = 206

        def __init__(self):
            self.payload = payload[4:]

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return None

        def read(self, size):
            remaining, self.payload = self.payload, b""
            return remaining

    def urlopen(request, timeout):
        requests.append((request.full_url, request.headers, timeout))
        return Response()

    monkeypatch.setattr(module.urllib.request, "urlopen", urlopen)
    module.fetch(
        "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/example.fastq.gz",
        destination,
        len(payload),
    )

    assert destination.read_bytes() == payload
    assert requests == [
        (
            "https://ftp.sra.ebi.ac.uk/vol1/fastq/example.fastq.gz",
            {"User-agent": "hostbias/0.1", "Range": "bytes=4-"},
            120,
        )
    ]


def test_fetch_pair_uses_resumable_segmented_aria2(tmp_path: Path, monkeypatch) -> None:
    module = load_fetch_module()
    destination = tmp_path / ".R1.fastq.gz.partial"
    payload = b"0123456789"
    calls = []

    monkeypatch.setattr(
        module.shutil, "which", lambda executable: "/usr/bin/aria2c"
    )

    def run(command, check, capture_output, text):
        calls.append((command, check, capture_output, text))
        destination.write_bytes(payload)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(module.subprocess, "run", run)
    module.fetch(
        "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/example.fastq.gz",
        destination,
        len(payload),
    )

    command, check, capture_output, text = calls[0]
    assert command[0] == "aria2c"
    assert "--continue=true" in command
    assert "--file-allocation=none" in command
    assert "--max-connection-per-server=4" in command
    assert "--split=4" in command
    assert command[-1] == (
        "https://ftp.sra.ebi.ac.uk/vol1/fastq/example.fastq.gz"
    )
    assert (check, capture_output, text) == (False, True, True)


def test_segmented_fetch_retries_http_success_with_wrong_size(
    tmp_path: Path, monkeypatch
) -> None:
    module = load_fetch_module()
    destination = tmp_path / ".R1.fastq.gz.partial"
    payload = b"0123456789"
    calls = []
    sleeps = []

    def run(command, check, capture_output, text):
        calls.append(command)
        destination.write_bytes(b"service unavailable" if len(calls) == 1 else payload)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(module.subprocess, "run", run)
    monkeypatch.setattr(module.time, "sleep", sleeps.append)
    module.fetch_with_aria2(
        "https://ftp.sra.ebi.ac.uk/example.fastq.gz",
        destination,
        len(payload),
        attempts=2,
    )

    assert destination.read_bytes() == payload
    assert len(calls) == 2
    assert sleeps == [1]


def test_segmented_fetch_preserves_controlled_partial_across_retry(
    tmp_path: Path, monkeypatch
) -> None:
    module = load_fetch_module()
    destination = tmp_path / ".R1.fastq.gz.partial"
    control = module.aria2_control_path(destination)
    payload = b"0123456789"
    calls = []

    def run(command, check, capture_output, text):
        calls.append(command)
        if len(calls) == 1:
            destination.write_bytes(payload[:4])
            control.write_bytes(b"aria2 session")
            return subprocess.CompletedProcess(command, 1, "", "interrupted")
        assert destination.read_bytes() == payload[:4]
        assert control.exists()
        destination.write_bytes(payload)
        control.unlink()
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(module.subprocess, "run", run)
    monkeypatch.setattr(module.time, "sleep", lambda seconds: None)
    module.fetch_with_aria2(
        "https://ftp.sra.ebi.ac.uk/example.fastq.gz",
        destination,
        len(payload),
        attempts=2,
    )

    assert destination.read_bytes() == payload
    assert len(calls) == 2


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
        "rule assembly_qc:",
        "rule map_competitive:",
        "rule alignment_contract:",
        "rule aggregate_endpoint:",
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


def test_production_stage_targets_and_vm_caps_are_explicit() -> None:
    snakefile = (PROJECT / "workflow" / "Snakefile").read_text(encoding="utf-8")
    for target in (
        "rule fetch_stage:",
        "rule normalize_stage:",
        "rule filter_stage:",
        "rule assemble_stage:",
        "rule downstream_bridge:",
    ):
        assert target in snakefile
    assert "--memory {resources.mem_mb}000000" in snakefile
    assert "--memory 0.9" not in snakefile

    profile = yaml.safe_load(
        (PROJECT / "profiles" / "vm" / "config.yaml").read_text(encoding="utf-8")
    )
    assert profile["cores"] == 30
    assert profile["jobs"] == 8
    assert "mem_mb=204800" in profile["resources"]
    assert profile["rerun-incomplete"] is True
    assert "minimap2 -c -x {params.preset:q} --secondary=yes -N 5" in snakefile
    assert "work/analysis_inputs" not in snakefile
    assert "PYTHONPATH=src python -m hostbias.cli" in snakefile


def test_mag_consensus_uses_shared_coverage_and_three_binners() -> None:
    rules = MAG_RULES.read_text(encoding="utf-8")
    for command in (
        "jgi_summarize_bam_contig_depths",
        "metabat2 -i",
        "run_MaxBin.pl",
        "cut_up_fasta.py",
        "concoct_coverage_table.py",
        "merge_cutup_clustering.py",
        "DAS_Tool",
    ):
        assert command in rules
    assert rules.count("reads.sorted.bam") >= 3
    assert "results/" not in rules
    assert f'{chr(34)}rules/mag.smk{chr(34)}' in (
        PROJECT / "workflow" / "Snakefile"
    ).read_text(encoding="utf-8")
    assert (
        (PROJECT / "workflow" / "Snakefile")
        .read_text(encoding="utf-8")
        .rstrip()
        .endswith('include: "rules/mag.smk"')
    )


def test_mag_qc_taxonomy_and_contract_are_private_and_pinned() -> None:
    rules = MAG_RULES.read_text(encoding="utf-8")
    for command in (
        "checkm2 predict",
        "gunc run",
        "gtdbtk classify_wf",
        "mag_translate.py contract",
    ):
        assert command in rules
    assert "GTDBTK_DATA_PATH=" in rules
    assert 'os.environ.get("CHECKM2DB"' in rules
    assert 'os.environ.get("GUNC_DB"' in rules
    assert 'os.environ.get("GTDBTK_DATA_PATH"' in rules
    assert "results/" not in rules


def test_checkm2_environment_avoids_unsatisfiable_bioconda_build() -> None:
    environment = yaml.safe_load(
        (PROJECT / "envs" / "checkm2.yaml").read_text(encoding="utf-8")
    )
    dependencies = environment["dependencies"]
    assert "python=3.12" in dependencies
    assert "tensorflow=2.17=cpu*" in dependencies
    assert "diamond=2.1.11" in dependencies
    assert "scikit-learn=1.6.1" in dependencies
    assert "checkm2=1.1.0" not in dependencies
    assert {
        "pip": [
            "git+https://github.com/chklovski/CheckM2.git@"
            "777bd767d93ac65decaf42a46d43ecc1bf7c41a3"
        ]
    } in dependencies
