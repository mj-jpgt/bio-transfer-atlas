#!/usr/bin/env python3
"""Atomically fetch and checksum a paired FASTQ input."""

from __future__ import annotations

import argparse
import hashlib
import os
import shutil
import subprocess
import time
import urllib.request
from pathlib import Path
from urllib.parse import urlsplit, urlunsplit


def md5_file(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def ena_https_url(url: str) -> str:
    """Use ENA's HTTPS endpoint because its FTP endpoint truncates large transfers."""
    parsed = urlsplit(url)
    if parsed.scheme == "ftp" and parsed.hostname == "ftp.sra.ebi.ac.uk":
        return urlunsplit(("https", parsed.netloc, parsed.path, parsed.query, parsed.fragment))
    return url


def aria2_control_path(destination: Path) -> Path:
    return destination.with_name(f"{destination.name}.aria2")


def fetch_with_aria2(url: str, destination: Path, expected_size: int) -> None:
    """Download one object with resumable segmented HTTP ranges."""
    control = aria2_control_path(destination)
    observed = destination.stat().st_size if destination.exists() else 0
    if observed == expected_size and not control.exists():
        return
    if observed > expected_size:
        destination.unlink()
        control.unlink(missing_ok=True)

    command = [
        "aria2c",
        "--allow-overwrite=true",
        "--auto-file-renaming=false",
        "--continue=true",
        "--file-allocation=none",
        "--max-connection-per-server=8",
        "--split=8",
        "--min-split-size=16M",
        "--max-tries=12",
        "--retry-wait=3",
        "--connect-timeout=30",
        "--timeout=120",
        "--summary-interval=0",
        "--console-log-level=warn",
        f"--dir={destination.parent}",
        f"--out={destination.name}",
        url,
    ]
    result = subprocess.run(command, check=False, capture_output=True, text=True)
    observed = destination.stat().st_size if destination.exists() else 0
    if result.returncode != 0 or observed != expected_size or control.exists():
        detail = (result.stderr or result.stdout).strip()[-1000:]
        raise ValueError(
            f"size mismatch for {destination.name}: expected {expected_size}, "
            f"observed {observed}; aria2c exit code {result.returncode}: {detail}"
        )


def fetch_with_urllib(
    url: str, destination: Path, expected_size: int, attempts: int
) -> None:
    """Portable single-stream fallback for local files and minimal environments."""
    is_network = urlsplit(url).scheme in {"http", "https", "ftp"}
    attempts = attempts if is_network else 1
    last_error: Exception | None = None

    for attempt in range(attempts):
        observed = destination.stat().st_size if destination.exists() else 0
        if observed == expected_size:
            return
        if observed > expected_size:
            destination.unlink()
            observed = 0

        headers = {"User-Agent": "hostbias/0.1"}
        if observed:
            headers["Range"] = f"bytes={observed}-"
        request = urllib.request.Request(url, headers=headers)

        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                status = getattr(response, "status", None)
                resume = observed > 0 and status == 206
                with destination.open("ab" if resume else "wb") as output:
                    while chunk := response.read(1024 * 1024):
                        output.write(chunk)
                    output.flush()
                    os.fsync(output.fileno())
        except Exception as error:
            last_error = error

        observed = destination.stat().st_size if destination.exists() else 0
        if observed == expected_size:
            return
        if attempt + 1 < attempts:
            time.sleep(min(2**attempt, 30))

    detail = f"; last transfer error: {last_error}" if last_error else ""
    raise ValueError(
        f"size mismatch for {destination.name}: expected {expected_size}, "
        f"observed {observed}{detail}"
    )


def fetch(url: str, destination: Path, expected_size: int, attempts: int = 12) -> None:
    """Download to a persistent partial file and verify its exact byte count."""
    url = ena_https_url(url)
    is_network = urlsplit(url).scheme in {"http", "https", "ftp"}
    if is_network and shutil.which("aria2c"):
        fetch_with_aria2(url, destination, expected_size)
        return
    fetch_with_urllib(url, destination, expected_size, attempts)


def fetch_pair(
    url1: str,
    md5_1: str,
    bytes_1: int,
    output1: Path,
    url2: str,
    md5_2: str,
    bytes_2: int,
    output2: Path,
) -> None:
    if output1.parent != output2.parent:
        raise ValueError("paired outputs must share a directory for atomic publication")
    output1.parent.mkdir(parents=True, exist_ok=True)
    temporary_paths = [
        output1.with_name(f".{output1.name}.partial"),
        output2.with_name(f".{output2.name}.partial"),
    ]
    for url, expected, expected_size, output, temporary in (
        (url1, md5_1.lower(), bytes_1, output1, temporary_paths[0]),
        (url2, md5_2.lower(), bytes_2, output2, temporary_paths[1]),
    ):
        fetch(url, temporary, expected_size)
        observed = md5_file(temporary)
        if observed != expected:
            temporary.unlink(missing_ok=True)
            raise ValueError(
                f"checksum mismatch for {output.name}: expected {expected}, observed {observed}"
            )
    os.replace(temporary_paths[0], output1)
    os.replace(temporary_paths[1], output2)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--url1", required=True)
    parser.add_argument("--md5-1", required=True)
    parser.add_argument("--bytes-1", required=True, type=int)
    parser.add_argument("--output1", required=True, type=Path)
    parser.add_argument("--url2", required=True)
    parser.add_argument("--md5-2", required=True)
    parser.add_argument("--bytes-2", required=True, type=int)
    parser.add_argument("--output2", required=True, type=Path)
    return parser.parse_args()


if __name__ == "__main__":
    arguments = parse_args()
    fetch_pair(
        arguments.url1,
        arguments.md5_1,
        arguments.bytes_1,
        arguments.output1,
        arguments.url2,
        arguments.md5_2,
        arguments.bytes_2,
        arguments.output2,
    )
