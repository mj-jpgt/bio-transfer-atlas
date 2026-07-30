#!/usr/bin/env python3
"""Uniformly subsample synchronized FASTQ pairs with one deterministic RNG."""

from __future__ import annotations

import argparse
import gzip
import io
import json
import random
from contextlib import ExitStack, contextmanager
from pathlib import Path
from typing import Iterator, TextIO

from validate_fastq_pair import canonical_name, open_fastq, validate_pair


@contextmanager
def deterministic_gzip_text(path: Path) -> Iterator[TextIO]:
    """Write gzip with a fixed timestamp so identical inputs are byte-identical."""

    raw = path.open("wb")
    compressed = gzip.GzipFile(
        filename="", mode="wb", fileobj=raw, compresslevel=1, mtime=0
    )
    text = io.TextIOWrapper(compressed, encoding="ascii", newline="\n")
    try:
        yield text
    finally:
        text.close()


def read_record(handle: TextIO) -> tuple[str, tuple[str, str, str, str]] | None:
    header = handle.readline()
    if not header:
        return None
    sequence = handle.readline()
    separator = handle.readline()
    quality = handle.readline()
    if not sequence or not separator.startswith("+") or not quality:
        raise ValueError(f"invalid FASTQ record at {header.rstrip()!r}")
    if len(sequence.rstrip("\r\n")) != len(quality.rstrip("\r\n")):
        raise ValueError(f"sequence/quality length mismatch at {header.rstrip()!r}")
    return canonical_name(header), (
        header.rstrip("\r\n"),
        sequence.rstrip("\r\n"),
        separator.rstrip("\r\n"),
        quality.rstrip("\r\n"),
    )


def fastp_pair_count(path: Path) -> int:
    report = json.loads(path.read_text(encoding="utf-8"))
    try:
        total_reads = report["summary"]["after_filtering"]["total_reads"]
    except (KeyError, TypeError) as error:
        raise ValueError("fastp JSON is missing summary.after_filtering.total_reads") from error
    if not isinstance(total_reads, int) or total_reads <= 0 or total_reads % 2:
        raise ValueError("fastp total_reads must be a positive even integer")
    return total_reads // 2


def subsample_pair(
    *,
    r1: Path,
    r2: Path,
    output1: Path,
    output2: Path,
    pairs: int,
    seed: int,
    expected_length: int,
    total_pairs: int | None = None,
) -> int:
    """Select exactly ``pairs`` uniformly without replacement, in input order."""

    if pairs <= 0:
        raise ValueError("pairs must be positive")
    total = (
        validate_pair(r1, r2, expected_length=expected_length)
        if total_pairs is None
        else total_pairs
    )
    if total < pairs:
        raise ValueError(f"requested {pairs} pairs but only {total} passed preprocessing")

    output1.parent.mkdir(parents=True, exist_ok=True)
    output2.parent.mkdir(parents=True, exist_ok=True)
    rng = random.Random(seed)
    needed = pairs
    remaining = total
    success = False
    try:
        with ExitStack() as stack:
            first = stack.enter_context(open_fastq(r1))
            second = stack.enter_context(open_fastq(r2))
            first_out = stack.enter_context(deterministic_gzip_text(output1))
            second_out = stack.enter_context(deterministic_gzip_text(output2))
            while remaining:
                record1 = read_record(first)
                record2 = read_record(second)
                if record1 is None or record2 is None:
                    raise ValueError("FASTQ ended before validated pair count")
                if record1[0] != record2[0]:
                    raise ValueError(
                        f"mate names differ: {record1[0]!r} versus {record2[0]!r}"
                    )
                if (
                    len(record1[1][1]) != expected_length
                    or len(record2[1][1]) != expected_length
                ):
                    raise ValueError(
                        f"pair {record1[0]!r} is not exactly {expected_length} bp"
                    )
                if rng.randrange(remaining) < needed:
                    first_out.write("\n".join(record1[1]) + "\n")
                    second_out.write("\n".join(record2[1]) + "\n")
                    needed -= 1
                remaining -= 1
            if read_record(first) is not None or read_record(second) is not None:
                raise ValueError(
                    f"FASTQ contains more than the reported {total} paired records"
                )
        if needed:
            raise AssertionError(f"subsampling ended with {needed} pairs unselected")
        validate_pair(
            output1,
            output2,
            expected_pairs=pairs,
            expected_length=expected_length,
        )
        success = True
        return pairs
    finally:
        if not success:
            output1.unlink(missing_ok=True)
            output2.unlink(missing_ok=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--r1", required=True, type=Path)
    parser.add_argument("--r2", required=True, type=Path)
    parser.add_argument("--output1", required=True, type=Path)
    parser.add_argument("--output2", required=True, type=Path)
    parser.add_argument("--pairs", required=True, type=int)
    parser.add_argument("--seed", required=True, type=int)
    parser.add_argument("--expected-length", required=True, type=int)
    parser.add_argument("--fastp-json", type=Path)
    return parser.parse_args()


if __name__ == "__main__":
    arguments = parse_args()
    observed = subsample_pair(
        r1=arguments.r1,
        r2=arguments.r2,
        output1=arguments.output1,
        output2=arguments.output2,
        pairs=arguments.pairs,
        seed=arguments.seed,
        expected_length=arguments.expected_length,
        total_pairs=(
            fastp_pair_count(arguments.fastp_json) if arguments.fastp_json else None
        ),
    )
    print(f"deterministically selected {observed} synchronized pairs")
