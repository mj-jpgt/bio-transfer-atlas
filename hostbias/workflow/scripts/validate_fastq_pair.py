#!/usr/bin/env python3
"""Stream-validation for synchronized paired FASTQ files."""

from __future__ import annotations

import argparse
import gzip
from contextlib import ExitStack
from pathlib import Path
from typing import TextIO


def open_fastq(path: Path) -> TextIO:
    if path.suffix == ".gz" or path.name.endswith(".gz.partial"):
        return gzip.open(path, "rt", encoding="ascii")
    return path.open("r", encoding="ascii")


def canonical_name(header: str) -> str:
    if not header.startswith("@"):
        raise ValueError(f"invalid FASTQ header: {header.rstrip()!r}")
    name = header[1:].split()[0]
    if name.endswith("/1") or name.endswith("/2"):
        name = name[:-2]
    return name


def read_record(handle: TextIO) -> tuple[str, str] | None:
    header = handle.readline()
    if not header:
        return None
    sequence = handle.readline().rstrip("\r\n")
    separator = handle.readline()
    quality = handle.readline().rstrip("\r\n")
    if not separator.startswith("+") or len(sequence) != len(quality):
        raise ValueError(f"invalid FASTQ record at {header.rstrip()!r}")
    return canonical_name(header), sequence


def validate_pair(
    r1: Path,
    r2: Path,
    expected_pairs: int | None = None,
    expected_length: int | None = None,
) -> int:
    count = 0
    with ExitStack() as stack:
        first = stack.enter_context(open_fastq(r1))
        second = stack.enter_context(open_fastq(r2))
        while True:
            record1 = read_record(first)
            record2 = read_record(second)
            if record1 is None and record2 is None:
                break
            if record1 is None or record2 is None:
                raise ValueError("mate files contain different numbers of records")
            if record1[0] != record2[0]:
                raise ValueError(f"mate names differ: {record1[0]!r} versus {record2[0]!r}")
            if expected_length is not None and (
                len(record1[1]) != expected_length or len(record2[1]) != expected_length
            ):
                raise ValueError(
                    f"{record1[0]} does not have expected read length {expected_length}"
                )
            count += 1
    if expected_pairs is not None and count != expected_pairs:
        raise ValueError(f"expected {expected_pairs} pairs, observed {count}")
    return count


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--r1", required=True, type=Path)
    parser.add_argument("--r2", required=True, type=Path)
    parser.add_argument("--expected-pairs", type=int)
    parser.add_argument("--expected-length", type=int)
    return parser.parse_args()


if __name__ == "__main__":
    arguments = parse_args()
    observed = validate_pair(
        arguments.r1,
        arguments.r2,
        arguments.expected_pairs,
        arguments.expected_length,
    )
    print(f"valid paired FASTQ: {observed} pairs")
