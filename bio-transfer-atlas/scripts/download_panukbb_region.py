"""
Lean Pan-UKB region fetch via tabix linear index (16 kb bins).

Downloads ONLY a genomic window (default: Duffy/ACKR1 neighborhood), never a
full chromosome. Safe for ~16 GB RAM laptops.
"""
from __future__ import annotations

import argparse
import gzip
import struct
import subprocess
import zlib
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
S3 = "https://pan-ukb-us-east-1.s3.amazonaws.com"
TBI_DIR = ROOT / "data/raw/panukbb/tbi"
TBI_DIR.mkdir(parents=True, exist_ok=True)

TRAITS = {
    "WBC": "sumstats_flat_files/continuous-30000-both_sexes-irnt.tsv.bgz",
    "T2D": "sumstats_flat_files/phecode-250.2-both_sexes.tsv.bgz",
    "CAD": "sumstats_flat_files/phecode-411.4-both_sexes.tsv.bgz",
    "BMI": "sumstats_flat_files/continuous-21001-both_sexes-irnt.tsv.bgz",
    "LDL": "sumstats_flat_files/continuous-LDLC-both_sexes-medadj_irnt.tsv.bgz",
}
ANCS = ["AFR", "AMR", "CSA", "EAS", "EUR", "MID"]

# Duffy-null rs2814778 neighborhood (GRCh38)
DEFAULT_CHR = "1"
DEFAULT_START = 159_000_000
DEFAULT_END = 159_400_000


def normalize_chr(value: str) -> str:
    s = str(value).strip()
    return s[3:] if s.lower().startswith("chr") else s


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Lean Pan-UKB genomic-window download.")
    p.add_argument("--trait", default="WBC", choices=sorted(TRAITS))
    p.add_argument("--chrom", default=DEFAULT_CHR)
    p.add_argument("--start", type=int, default=DEFAULT_START)
    p.add_argument("--end", type=int, default=DEFAULT_END)
    p.add_argument(
        "--max-mb",
        type=float,
        default=80.0,
        help="Abort if estimated byte range exceeds this many MB (safety).",
    )
    return p.parse_args()


def curl_range(url: str, start: int, end: int | None, out_path: Path) -> Path:
    from bta_curl import curl_bin

    rng = f"{start}-{end if end is not None else ''}"
    cmd = [curl_bin(), "-sL", "-r", rng, "-o", str(out_path), url]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"curl range {rng} failed: {r.stderr[-300:]}")
    return out_path


def curl_full(url: str, out_path: Path) -> Path:
    from bta_curl import curl_bin

    cmd = [curl_bin(), "-sL", "-o", str(out_path), url]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"curl failed: {r.stderr[-300:]}")
    return out_path


def region_coffsets(tbi_path: Path, chrom: str, start: int, end: int) -> tuple[int, int]:
    """Return (coffset_start, coffset_end_exclusive) from TBI linear index."""
    raw = gzip.open(tbi_path, "rb").read()
    p = 0
    magic = raw[p : p + 4]
    p += 4
    if magic != b"TBI\x01":
        raise RuntimeError(f"bad tbi magic {magic!r}")
    (n_ref, _fmt, _col_seq, _col_beg, _col_end, _meta, _skip, l_nm) = struct.unpack(
        "<8i", raw[p : p + 32]
    )
    p += 32
    names_blob = raw[p : p + l_nm]
    p += l_nm
    names = [n.decode() for n in names_blob.split(b"\x00") if n]
    target_idx = next((i for i, nm in enumerate(names) if normalize_chr(nm) == chrom), None)
    if target_idx is None:
        raise RuntimeError(f"chr{chrom} not in tbi: {names}")

    bin_size = 16_384  # TBI linear index window
    i0 = max(0, start // bin_size)
    i1 = max(i0 + 1, (end // bin_size) + 2)

    for ref in range(n_ref):
        (n_bin,) = struct.unpack("<i", raw[p : p + 4])
        p += 4
        for _ in range(n_bin):
            (_bin_id, n_chunk) = struct.unpack("<Ii", raw[p : p + 8])
            p += 8
            p += 16 * n_chunk
        (n_intv,) = struct.unpack("<i", raw[p : p + 4])
        p += 4
        ioffs = struct.unpack(f"<{n_intv}Q", raw[p : p + 8 * n_intv])
        p += 8 * n_intv
        if ref != target_idx:
            continue
        # virtual offsets -> compressed file offsets
        def coff(v: int) -> int:
            return int(v >> 16)

        # find first nonzero at/before i0
        c_start = None
        for i in range(min(i0, n_intv - 1), -1, -1):
            if ioffs[i] > 0:
                c_start = coff(ioffs[i])
                break
        if c_start is None:
            nz = [coff(o) for o in ioffs if o > 0]
            c_start = min(nz) if nz else 0
        c_end = None
        for i in range(min(i1, n_intv - 1), n_intv):
            if ioffs[i] > 0 and coff(ioffs[i]) > c_start:
                c_end = coff(ioffs[i])
                break
        if c_end is None:
            # pad ~8 MB past start as a soft end (region scripts should stay small)
            c_end = c_start + 8 * 1024 * 1024
        return c_start, c_end
    raise RuntimeError("region coffsets not found")


def bgzf_decompress(data: bytes) -> bytes:
    out = bytearray()
    n = len(data)
    p = 0
    while p + 18 <= n:
        if data[p] != 0x1F or data[p + 1] != 0x8B:
            break
        xlen = struct.unpack("<H", data[p + 10 : p + 12])[0]
        extra = data[p + 12 : p + 12 + xlen]
        bsize = None
        q = 0
        while q + 4 <= len(extra):
            si1, si2, slen = extra[q], extra[q + 1], struct.unpack("<H", extra[q + 2 : q + 4])[0]
            if si1 == 66 and si2 == 67:
                bsize = struct.unpack("<H", extra[q + 4 : q + 6])[0]
            q += 4 + slen
        if bsize is None:
            break
        block_len = bsize + 1
        if p + block_len > n:
            break
        cdata = data[p + 12 + xlen : p + block_len - 8]
        try:
            out.extend(zlib.decompress(cdata, -15))
        except Exception:
            break
        p += block_len
    return bytes(out)


def read_header(url: str, out_dir: Path) -> list[str]:
    tmp = out_dir / "_hdr_region.bgz"
    curl_range(url, 0, 200_000, tmp)
    txt = bgzf_decompress(tmp.read_bytes()).decode("utf-8", errors="replace")
    tmp.unlink(missing_ok=True)
    return txt.split("\n", 1)[0].rstrip("\r").split("\t")


def main() -> None:
    args = parse_args()
    chrom = normalize_chr(args.chrom)
    if args.end <= args.start:
        raise SystemExit("--end must be > --start")
    if (args.end - args.start) > 5_000_000:
        raise SystemExit("Refuse windows >5 Mb on this lean path; split the region.")

    key = TRAITS[args.trait]
    url = f"{S3}/{key}"
    tbi_url = f"{S3}/sumstats_flat_files_tabix/{Path(key).name}.tbi"
    out_dir = ROOT / "data/raw/panukbb" / f"regions"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_pq = out_dir / f"{args.trait}.chr{chrom}_{args.start}_{args.end}.parquet"
    if out_pq.exists() and out_pq.stat().st_size > 10_000:
        print(f"cached -> {out_pq}")
        return

    tbi_path = TBI_DIR / f"{args.trait}.region.tbi"
    print(f"Downloading tabix index for {args.trait} ...", flush=True)
    curl_full(tbi_url, tbi_path)
    c0, c1 = region_coffsets(tbi_path, chrom, args.start, args.end)
    est_mb = (c1 - c0) / 1e6
    print(f"Byte range {c0:,}-{c1:,} (~{est_mb:.1f} MB compressed)", flush=True)
    if est_mb > args.max_mb:
        raise SystemExit(
            f"Estimated download {est_mb:.1f} MB > --max-mb {args.max_mb}. "
            "Narrow the window or raise --max-mb deliberately."
        )

    header = read_header(url, out_dir)
    slice_path = out_dir / f"_{args.trait}.chr{chrom}_{args.start}_{args.end}.bgz"
    print("Fetching region slice ...", flush=True)
    curl_range(url, c0, c1, slice_path)
    print(f"  got {slice_path.stat().st_size / 1e6:.1f} MB", flush=True)

    keep = ["chr", "pos", "ref", "alt"]
    for a in ANCS:
        for pre in ("beta", "se", "pval", "low_confidence"):
            c = f"{pre}_{a}"
            if c in header:
                keep.append(c)
    for c in ("beta_meta", "se_meta", "pval_meta", "af_meta"):
        if c in header:
            keep.append(c)
    keep = list(dict.fromkeys([c for c in keep if c in header]))

    # Decompress small slice fully (must stay under max-mb)
    text = bgzf_decompress(slice_path.read_bytes()).decode("utf-8", errors="replace")
    lines = text.splitlines()
    # drop partial first line
    if lines and not lines[0].startswith(header[0]) and "\t" in lines[0]:
        lines = lines[1:]
    from io import StringIO

    buf = StringIO("\n".join(lines))
    df = pd.read_csv(buf, sep="\t", names=header, usecols=keep, dtype=str, on_bad_lines="skip")
    df["chr"] = df["chr"].map(normalize_chr)
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
    df = df.dropna(subset=["pos"])
    df["pos"] = df["pos"].astype(int)
    df = df[(df["chr"] == chrom) & (df["pos"] >= args.start) & (df["pos"] <= args.end)]
    df.to_parquet(out_pq, index=False)
    slice_path.unlink(missing_ok=True)
    print(f"Saved {out_pq} ({len(df):,} rows)")


if __name__ == "__main__":
    main()
