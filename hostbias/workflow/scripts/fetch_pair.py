#!/usr/bin/env python3
"""Atomically fetch and checksum a paired FASTQ input."""

from __future__ import annotations

import argparse
import hashlib
import os
import tempfile
import urllib.request
from pathlib import Path


def md5_file(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch(url: str, destination: Path) -> None:
    request = urllib.request.Request(url, headers={"User-Agent": "hostbias/0.1"})
    with urllib.request.urlopen(request, timeout=120) as response:
        with destination.open("wb") as output:
            while chunk := response.read(1024 * 1024):
                output.write(chunk)
            output.flush()
            os.fsync(output.fileno())


def fetch_pair(
    url1: str,
    md5_1: str,
    output1: Path,
    url2: str,
    md5_2: str,
    output2: Path,
) -> None:
    if output1.parent != output2.parent:
        raise ValueError("paired outputs must share a directory for atomic publication")
    output1.parent.mkdir(parents=True, exist_ok=True)
    temporary_paths: list[Path] = []
    try:
        for url, expected, output in (
            (url1, md5_1.lower(), output1),
            (url2, md5_2.lower(), output2),
        ):
            descriptor, name = tempfile.mkstemp(dir=output.parent, prefix=f".{output.name}.")
            os.close(descriptor)
            temporary = Path(name)
            temporary_paths.append(temporary)
            fetch(url, temporary)
            observed = md5_file(temporary)
            if observed != expected:
                raise ValueError(
                    f"checksum mismatch for {output.name}: expected {expected}, observed {observed}"
                )
        os.replace(temporary_paths[0], output1)
        os.replace(temporary_paths[1], output2)
    finally:
        for temporary in temporary_paths:
            temporary.unlink(missing_ok=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--url1", required=True)
    parser.add_argument("--md5-1", required=True)
    parser.add_argument("--output1", required=True, type=Path)
    parser.add_argument("--url2", required=True)
    parser.add_argument("--md5-2", required=True)
    parser.add_argument("--output2", required=True, type=Path)
    return parser.parse_args()


if __name__ == "__main__":
    arguments = parse_args()
    fetch_pair(
        arguments.url1,
        arguments.md5_1,
        arguments.output1,
        arguments.url2,
        arguments.md5_2,
        arguments.output2,
    )
