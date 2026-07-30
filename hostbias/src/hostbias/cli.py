"""Command-line entry points for Hostbias."""

from __future__ import annotations

import argparse
from pathlib import Path

from hostbias.config import ValidationError, load_and_validate
from hostbias.provenance import build_provenance, write_json_atomic


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="hostbias")
    parser.add_argument("--version", action="version", version="%(prog)s 0.1.0")
    subparsers = parser.add_subparsers(dest="command", required=True)
    validate_parser = subparsers.add_parser(
        "validate", help="validate configuration and manifests"
    )
    validate_parser.add_argument("--config", required=True, type=Path)
    provenance_parser = subparsers.add_parser(
        "provenance", help="write a reproducibility manifest"
    )
    provenance_parser.add_argument("--config", required=True, type=Path)
    provenance_parser.add_argument("--output", required=True, type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        inputs = load_and_validate(args.config)
        if args.command == "validate":
            print(f"valid: {len(inputs.samples)} samples; experiment={inputs.config['experiment']['id']}")
            return 0
        if args.command == "provenance":
            write_json_atomic(build_provenance(inputs), args.output)
            print(args.output)
            return 0
    except ValidationError as error:
        parser.exit(2, f"hostbias: error: {error}\n")
    parser.error(f"unknown command: {args.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
