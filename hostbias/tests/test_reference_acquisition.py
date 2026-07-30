from __future__ import annotations

import csv
import gzip
import hashlib
from pathlib import Path
from typing import Mapping, Sequence

import pytest

from hostbias.data_manifest import canonical_tsv
from hostbias.reference_acquisition import (
    METADATA_FIELDS,
    AcquisitionError,
    _expected_md5_from_sidecar,
    build_reference_panel,
)
from hostbias.reference_catalog import (
    DONOR_FIELDS,
    PANEL_FIELDS,
    SUPERPOPULATIONS,
    build_balanced_donor_catalog,
)


def tsv_bytes(rows: Sequence[Mapping[str, object]], fields: Sequence[str]) -> bytes:
    return canonical_tsv(rows, fields)


def csv_bytes(rows: list[dict[str, str]], fields: list[str], delimiter: str = ",") -> bytes:
    import io

    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=fields, delimiter=delimiter, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue().encode()


class FakeHttpClient:
    def __init__(self, payloads: dict[str, bytes], sidecars: dict[str, str]) -> None:
        self.payloads = payloads
        self.sidecars = sidecars
        self.download_calls: list[str] = []

    def text(self, url: str) -> str:
        return self.sidecars[url]

    def content_length(self, url: str) -> int:
        return len(self.payloads[url])

    def download(self, url: str, destination: Path, expected_bytes: int) -> None:
        payload = self.payloads[url]
        assert len(payload) == expected_bytes
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(payload)
        self.download_calls.append(url)


class FakeIndexBuilder:
    def __init__(self) -> None:
        self.build_calls = 0

    def version(self) -> str:
        return "minimap2-mock-1"

    def build(
        self,
        union_fasta: Path,
        output_index: Path,
        *,
        threads: int,
        index_batch: str,
    ) -> dict[str, object]:
        self.build_calls += 1
        output_index.parent.mkdir(parents=True, exist_ok=True)
        output_index.write_bytes(
            b"mock-mmi:" + hashlib.sha256(union_fasta.read_bytes()).hexdigest().encode()
        )
        return {
            "preset": "asm5",
            "index_batch": index_batch,
            "threads": threads,
            "seconds": 0.01,
        }


def build_fixture(root: Path) -> tuple[Path, Path, Path, dict[str, bytes], dict[str, str]]:
    assemblies: list[dict[str, str]] = []
    samples: list[dict[str, str]] = []
    igsr: list[dict[str, str]] = []
    donors: list[dict[str, str]] = []
    panel: list[dict[str, str]] = []
    payloads: dict[str, bytes] = {}
    sidecars: dict[str, str] = {}

    chm_url = "https://example.test/chm13.fa.gz"
    chm_md5_url = "https://example.test/chm13.md5"
    chm_payload = gzip.compress(b">chr1\nACGT\n", mtime=0)
    payloads[chm_url] = chm_payload
    chm_md5 = hashlib.md5(chm_payload, usedforsecurity=False).hexdigest()
    sidecars[chm_md5_url] = f"{chm_md5}  ./chm13.fa.gz\n"
    panel.append(
        {
            "reference_id": "chm13v2.0",
            "kind": "chm13",
            "donor_id": "CHM13",
            "haplotype": "0",
            "population_code": "NA",
            "superpopulation": "NA",
            "panel_role": "reference",
            "genbank_accession": "GCA_TEST",
            "source_url": chm_url,
            "source_md5_url": chm_md5_url,
            "expected_bytes": str(len(chm_payload)),
            "expected_md5": chm_md5,
            "expected_sha256": hashlib.sha256(chm_payload).hexdigest(),
            "include_in_union": "true",
            "local_filename": "chm13.fa.gz",
        }
    )

    for group in SUPERPOPULATIONS:
        population = f"{group}P"
        for rank in (1, 2):
            donor = f"{group}{rank}"
            role = "primary" if rank == 1 else "holdout"
            samples.append(
                {
                    "sample_id": donor,
                    "biosample_id": f"SAMN{group}{rank}",
                    "population_abbreviation": population,
                }
            )
            igsr.append({"sample": donor, "pop": population, "super_pop": group})
            donors.append(
                {
                    "donor_id": donor,
                    "biosample_accession": f"SAMN{group}{rank}",
                    "population_code": population,
                    "superpopulation": group,
                    "panel_role": role,
                    "selection_rank": str(rank),
                    "haplotype_count": "2",
                    "metadata_basis": "HPRC+IGSR exact donor join",
                }
            )
            for haplotype in ("1", "2"):
                filename = f"{donor}_hap{haplotype}.fa.gz"
                s3_url = f"s3://human-pangenomics/test/{filename}"
                source_url = (
                    "https://s3-us-west-2.amazonaws.com/human-pangenomics/"
                    f"test/{filename}"
                )
                md5_url = source_url + ".md5"
                payload = gzip.compress(
                    f">contig\n{group}ACGT{haplotype}\n".encode(), mtime=0
                )
                digest = hashlib.md5(payload, usedforsecurity=False).hexdigest()
                payloads[source_url] = payload
                sidecars[md5_url] = f"{digest}  {filename}\n"
                accession = f"GCA_{group}_{rank}_{haplotype}"
                assemblies.append(
                    {
                        "sample_id": donor,
                        "haplotype": haplotype,
                        "genbank_accession": accession,
                        "assembly_md5": s3_url + ".md5",
                        "assembly": s3_url,
                    }
                )
                panel.append(
                    {
                        "reference_id": f"hprc_{donor}_h{haplotype}",
                        "kind": "hprc_assembly",
                        "donor_id": donor,
                        "haplotype": haplotype,
                        "population_code": population,
                        "superpopulation": group,
                        "panel_role": role,
                        "genbank_accession": accession,
                        "source_url": source_url,
                        "source_md5_url": md5_url,
                        "expected_bytes": str(len(payload)),
                        "expected_md5": digest,
                        "expected_sha256": "",
                        "include_in_union": str(role == "primary").lower(),
                        "local_filename": filename,
                    }
                )

    assembly_bytes = csv_bytes(
        assemblies,
        ["sample_id", "haplotype", "genbank_accession", "assembly_md5", "assembly"],
    )
    sample_bytes = csv_bytes(
        samples, ["sample_id", "biosample_id", "population_abbreviation"]
    )
    igsr_bytes = csv_bytes(igsr, ["sample", "pop", "super_pop"], delimiter="\t")
    metadata_payloads = {
        "https://example.test/assemblies.csv": assembly_bytes,
        "https://example.test/samples.csv": sample_bytes,
        "https://example.test/panel.tsv": igsr_bytes,
    }
    payloads.update(metadata_payloads)
    metadata = []
    for source_id, url in (
        ("hprc_release2_assemblies", "https://example.test/assemblies.csv"),
        ("hprc_release2_samples", "https://example.test/samples.csv"),
        ("igsr_1000g_20130502_panel", "https://example.test/panel.tsv"),
    ):
        payload = metadata_payloads[url]
        metadata.append(
            {
                "source_id": source_id,
                "authority": "fixture",
                "url": url,
                "expected_bytes": str(len(payload)),
                "sha256": hashlib.sha256(payload).hexdigest(),
                "retrieved_date_utc": "2026-07-30",
                "fields_used": "fixture",
            }
        )

    metadata_path = root / "metadata.tsv"
    donors_path = root / "donors.tsv"
    panel_path = root / "panel.tsv"
    metadata_path.write_bytes(tsv_bytes(metadata, METADATA_FIELDS))
    donors_path.write_bytes(tsv_bytes(donors, DONOR_FIELDS))
    panel_path.write_bytes(tsv_bytes(panel, PANEL_FIELDS))
    return metadata_path, donors_path, panel_path, payloads, sidecars


def test_catalog_selection_uses_exact_join_and_balances_groups(tmp_path: Path) -> None:
    metadata_path, _, _, payloads, _ = build_fixture(tmp_path)
    metadata = list(csv.DictReader(metadata_path.open(), delimiter="\t"))
    assemblies = list(
        csv.DictReader(payloads[metadata[0]["url"]].decode().splitlines())
    )
    samples = list(csv.DictReader(payloads[metadata[1]["url"]].decode().splitlines()))
    igsr = list(
        csv.DictReader(
            payloads[metadata[2]["url"]].decode().splitlines(), delimiter="\t"
        )
    )
    selected = build_balanced_donor_catalog(assemblies, samples, igsr)
    assert len(selected) == 10
    assert [(row["superpopulation"], row["panel_role"]) for row in selected] == [
        pair for group in SUPERPOPULATIONS for pair in ((group, "primary"), (group, "holdout"))
    ]


def test_complete_acquisition_is_verified_aggregate_only_and_restartable(
    tmp_path: Path,
) -> None:
    metadata, donors, panel, payloads, sidecars = build_fixture(tmp_path)
    client = FakeHttpClient(payloads, sidecars)
    indexer = FakeIndexBuilder()
    reference_root = tmp_path / "references"
    checkpoint = tmp_path / "aggregate" / "reference.json"
    report = build_reference_panel(
        metadata_sources_path=metadata,
        donors_path=donors,
        panel_path=panel,
        reference_root=reference_root,
        checkpoint_path=checkpoint,
        threads=4,
        client=client,
        index_builder=indexer,
    )
    assert report["status"] == "PASS"
    assert report["download_totals"]["references"] == 21
    assert report["union_fasta"]["included_reference_count"] == 11
    assert report["union_fasta"]["headers_namespaced"] is True
    assert indexer.build_calls == 1
    union = (reference_root / "panel" / "hprc-balanced-chm13-v1.fa").read_text()
    assert ">chm13v2.0|chr1" in union
    assert ">hprc_AFR1_h1|contig" in union
    assert ">hprc_AFR2_h1|contig" not in union
    aggregate_text = checkpoint.read_text(encoding="utf-8")
    assert str(tmp_path) not in aggregate_text
    assert "ACGT" not in aggregate_text

    resumed_client = FakeHttpClient(payloads, sidecars)
    resumed_indexer = FakeIndexBuilder()
    resumed = build_reference_panel(
        metadata_sources_path=metadata,
        donors_path=donors,
        panel_path=panel,
        reference_root=reference_root,
        checkpoint_path=checkpoint,
        threads=4,
        client=resumed_client,
        index_builder=resumed_indexer,
    )
    assert resumed["download_totals"]["reused"] == 21
    assert resumed_indexer.build_calls == 0
    assert resumed_client.download_calls == []


def test_corrupt_reference_is_redownloaded(tmp_path: Path) -> None:
    metadata, donors, panel, payloads, sidecars = build_fixture(tmp_path)
    reference_root = tmp_path / "references"
    build_reference_panel(
        metadata_sources_path=metadata,
        donors_path=donors,
        panel_path=panel,
        reference_root=reference_root,
        checkpoint_path=tmp_path / "first.json",
        threads=2,
        client=FakeHttpClient(payloads, sidecars),
        index_builder=FakeIndexBuilder(),
    )
    corrupt = reference_root / "downloads" / "AFR1_hap1.fa.gz"
    corrupt.write_bytes(b"corrupt")
    retry_client = FakeHttpClient(payloads, sidecars)
    build_reference_panel(
        metadata_sources_path=metadata,
        donors_path=donors,
        panel_path=panel,
        reference_root=reference_root,
        checkpoint_path=tmp_path / "second.json",
        threads=2,
        client=retry_client,
        index_builder=FakeIndexBuilder(),
    )
    assert retry_client.download_calls == [
        "https://s3-us-west-2.amazonaws.com/human-pangenomics/test/AFR1_hap1.fa.gz"
    ]


def test_md5_sidecar_requires_exact_filename() -> None:
    assert (
        _expected_md5_from_sidecar("a" * 32 + "  ./target.fa.gz\n", "target.fa.gz")
        == "a" * 32
    )
    with pytest.raises(AcquisitionError, match="found 0"):
        _expected_md5_from_sidecar("a" * 32 + "  other.fa.gz\n", "target.fa.gz")
