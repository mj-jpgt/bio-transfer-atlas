#!/usr/bin/env python3
"""Python 3.9+-compatible wrappers around strict MAG translators."""

from __future__ import annotations

import argparse
from pathlib import Path

from hostbias.mag_bridge import (
    bins_to_scaffolds2bin,
    build_mag_contracts,
    depth_to_maxbin_abundance,
)


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser()
    commands = root.add_subparsers(dest="command", required=True)

    abundance = commands.add_parser("abundance")
    abundance.add_argument("--depth", type=Path, required=True)
    abundance.add_argument("--output", type=Path, required=True)

    bin_map = commands.add_parser("bins-to-map")
    bin_map.add_argument("--bin-dir", type=Path, required=True)
    bin_map.add_argument("--output", type=Path, required=True)
    bin_map.add_argument("--bin-prefix", required=True)

    contract = commands.add_parser("contract")
    contract.add_argument("--sample-id", required=True)
    contract.add_argument("--dastool-map", type=Path, required=True)
    contract.add_argument("--checkm2-report", type=Path, required=True)
    contract.add_argument("--gunc-report", type=Path, required=True)
    contract.add_argument("--gtdb-bacterial-summary", type=Path, required=True)
    contract.add_argument("--gtdb-archaeal-summary", type=Path, required=True)
    contract.add_argument("--contig-bins-output", type=Path, required=True)
    contract.add_argument("--bin-qc-output", type=Path, required=True)
    return root


def main() -> None:
    args = parser().parse_args()
    if args.command == "abundance":
        depth_to_maxbin_abundance(args.depth, args.output)
    elif args.command == "bins-to-map":
        bins_to_scaffolds2bin(
            args.bin_dir,
            args.output,
            args.bin_prefix,
        )
    else:
        build_mag_contracts(
            args.sample_id,
            args.dastool_map,
            args.checkm2_report,
            args.gunc_report,
            (args.gtdb_bacterial_summary, args.gtdb_archaeal_summary),
            args.contig_bins_output,
            args.bin_qc_output,
        )


if __name__ == "__main__":
    main()
