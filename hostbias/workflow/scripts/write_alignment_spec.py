#!/usr/bin/env python3
"""Write a private, checksum-pinned alignment bridge run specification."""

from __future__ import annotations

import argparse
import json
import os
import tempfile
from pathlib import Path

from hostbias.provenance import sha256_file


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--filter-mode", choices=("source", "strict"), required=True)
    parser.add_argument("--assembly", type=Path, required=True)
    parser.add_argument("--mapper-version-file", type=Path, required=True)
    parser.add_argument("--mapper-preset", required=True)
    parser.add_argument("--human-paf", type=Path, required=True)
    parser.add_argument("--human-reference-id", required=True)
    parser.add_argument("--human-reference-sha256", required=True)
    parser.add_argument("--gtdb-paf", type=Path, required=True)
    parser.add_argument("--gtdb-reference-id", required=True)
    parser.add_argument("--gtdb-reference-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = arguments()
    payload = {
        "schema_version": "1.0",
        "sample_id": args.sample_id,
        "filter_mode": args.filter_mode,
        "assembly": {
            "path": str(args.assembly.resolve()),
            "sha256": sha256_file(args.assembly),
        },
        "mapper": {
            "version": args.mapper_version_file.read_text(encoding="utf-8").strip(),
            "preset": args.mapper_preset,
        },
        "references": [
            {
                "domain": "human",
                "reference_id": args.human_reference_id,
                "reference_sha256": args.human_reference_sha256,
                "paf_path": str(args.human_paf.resolve()),
            },
            {
                "domain": "gtdb",
                "reference_id": args.gtdb_reference_id,
                "reference_sha256": args.gtdb_reference_sha256,
                "paf_path": str(args.gtdb_paf.resolve()),
            },
        ],
    }
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=output.parent, prefix=f".{output.name}.", suffix=".tmp"
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, output)
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise


if __name__ == "__main__":
    main()
