from pathlib import Path


def test_fetch_validation_sentinel_is_ancient_for_protected_fastqs() -> None:
    snakefile = Path(__file__).parents[1] / "workflow" / "Snakefile"
    source = snakefile.read_text(encoding="utf-8")
    assert 'validation=ancient(f"{WORK}/validated.ok")' in source
    assert 'r1=protected(f"{WORK}/raw/{{sample}}_R1.fastq.gz")' in source
