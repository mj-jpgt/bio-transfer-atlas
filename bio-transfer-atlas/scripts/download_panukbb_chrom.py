"""
Fetch Pan-UKBB multi-ancestry summary stats for a target chromosome using tabix byte ranges.
"""
from __future__ import annotations

import argparse
import gzip
import struct
import subprocess
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
S3 = "https://pan-ukb-us-east-1.s3.amazonaws.com"
TBI_DIR = ROOT / "data/raw/panukbb/tbi"
TBI_DIR.mkdir(parents=True, exist_ok=True)

TRAITS = {
    "T2D": "sumstats_flat_files/phecode-250.2-both_sexes.tsv.bgz",
    "CAD": "sumstats_flat_files/phecode-411.4-both_sexes.tsv.bgz",
    "BMI": "sumstats_flat_files/continuous-21001-both_sexes-irnt.tsv.bgz",
    "LDL": "sumstats_flat_files/continuous-LDLC-both_sexes-medadj_irnt.tsv.bgz",
    # UKB field 30000 — white blood cell (leukocyte) count; Duffy-null positive control
    "WBC": "sumstats_flat_files/continuous-30000-both_sexes-irnt.tsv.bgz",
}
ANCS = ["AFR", "AMR", "CSA", "EAS", "EUR", "MID"]


def normalize_chr(value: str) -> str:
    s = str(value).strip()
    return s[3:] if s.lower().startswith("chr") else s


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--chrom", default="22", help="Chromosome number/name (e.g. 22 or chr22)")
    p.add_argument(
        "--traits",
        default="",
        help="Comma-separated trait keys (default: all). Example: WBC",
    )
    p.add_argument(
        "--force-full-chrom",
        action="store_true",
        help="Override lean-machine guard that blocks chr1-7 full downloads.",
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


def parse_tbi_target_coffset(tbi_path: Path, target_chrom: str) -> int:
    raw = gzip.open(tbi_path, "rb").read()
    p = 0
    magic = raw[p : p + 4]
    p += 4
    if magic != b"TBI\x01":
        raise RuntimeError(f"bad tbi magic {magic!r}")
    (n_ref, _fmt, _col_seq, _col_beg, _col_end, _meta, _skip, l_nm) = struct.unpack("<8i", raw[p : p + 32])
    p += 32
    names_blob = raw[p : p + l_nm]
    p += l_nm
    names = [n.decode() for n in names_blob.split(b"\x00") if n]

    target_idx = None
    for i, nm in enumerate(names):
        if normalize_chr(nm) == target_chrom:
            target_idx = i
            break
    if target_idx is None:
        raise RuntimeError(f"chr{target_chrom} not found in tbi contigs: {names}")

    min_coffset = None
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
        if ref == target_idx:
            nz = [(o >> 16) for o in ioffs if o > 0]
            min_coffset = min(nz) if nz else 0
    if min_coffset is None:
        raise RuntimeError(f"could not derive chr{target_chrom} coffset")
    return int(min_coffset)


def bgzf_decompress(data: bytes) -> bytes:
    import zlib

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
    tmp = out_dir / "_hdr.bgz"
    curl_range(url, 0, 200_000, tmp)
    txt = bgzf_decompress(tmp.read_bytes()).decode("utf-8", errors="replace")
    tmp.unlink(missing_ok=True)
    first_line = txt.split("\n", 1)[0]
    return first_line.rstrip("\r").split("\t")


def fetch_trait(trait: str, key: str, chrom: str, out_dir: Path) -> Path:
    out_pq = out_dir / f"{trait}.chr{chrom}.parquet"
    if out_pq.exists():
        print(f"  {trait}: cached -> {out_pq.name}")
        return out_pq

    url = f"{S3}/{key}"
    tbi_url = f"{S3}/sumstats_flat_files_tabix/{Path(key).name}.tbi"
    tbi_path = TBI_DIR / f"{trait}.chr{chrom}.tbi"

    print(f"  {trait}: downloading tabix index ...", flush=True)
    curl_full(tbi_url, tbi_path)
    coffset = parse_tbi_target_coffset(tbi_path, chrom)
    print(f"  {trait}: chr{chrom} starts at compressed byte {coffset:,}")

    print(f"  {trait}: fetching header ...", flush=True)
    header = read_header(url, out_dir)

    print(f"  {trait}: range-downloading chr{chrom} slice (coffset..EOF) ...", flush=True)
    tail = out_dir / f"{trait}.chr{chrom}.bgz"
    curl_range(url, coffset, None, tail)
    print(f"  {trait}: downloaded {tail.stat().st_size / 1e6:.1f} MB slice")

    keep = ["chr", "pos", "ref", "alt"]
    for a in ANCS:
        for pre in ("beta", "se", "pval", "low_confidence"):
            c = f"{pre}_{a}"
            if c in header:
                keep.append(c)
    for c in ("beta_meta", "se_meta", "pval_meta", "af_meta"):
        if c in header:
            keep.append(c)
    keep = list(dict.fromkeys(keep))

    # Stream through BGZF blocks one chunk at a time to avoid loading the
    # entire decompressed slice (can be 10+ GB for chr1) into RAM.
    print(f"  {trait}: streaming decompression ...", flush=True)
    chunks = []
    with gzip.open(str(tail), "rt", encoding="utf-8", errors="replace") as gz:
        gz.readline()  # skip first partial line (tail of previous chr's block)
        reader = pd.read_csv(
            gz,
            sep="\t",
            names=header,
            usecols=keep,
            dtype=str,
            engine="python",
            on_bad_lines="skip",
            chunksize=100_000,
        )
        for chunk in reader:
            sub = chunk[chunk["chr"].map(normalize_chr) == chrom]
            if not sub.empty:
                chunks.append(sub)
    df = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame(columns=header)
    df = df.copy()
    df = df[[c for c in keep if c in df.columns]].copy()
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
    df = df.dropna(subset=["pos"])
    df["pos"] = df["pos"].astype(int)

    df.to_parquet(out_pq, index=False)
    tail.unlink(missing_ok=True)
    print(f"  {trait}: saved {len(df):,} variants -> {out_pq.name}")
    return out_pq


def main() -> None:
    args = parse_args()
    chrom = normalize_chr(args.chrom)
    # Safety: full-chr downloads of large autosomes crash 16 GB machines.
    # Prefer scripts/download_panukbb_region.py for windows; allow override.
    large = chrom.isdigit() and int(chrom) <= 7
    if large and not getattr(args, "force_full_chrom", False):
        raise SystemExit(
            f"Refusing full-chr{chrom} download on lean local budget "
            f"(~2GB compressed for large chroms). Use:\n"
            f"  python scripts/download_panukbb_region.py --trait WBC --chrom {chrom} "
            f"--start START --end END\n"
            f"Or pass --force-full-chrom if you truly have RAM/disk headroom."
        )
    out_dir = ROOT / "data/raw/panukbb" / f"chr{chrom}"
    out_dir.mkdir(parents=True, exist_ok=True)

    want = {t.strip().upper() for t in args.traits.split(",") if t.strip()}
    traits = {k: v for k, v in TRAITS.items() if not want or k.upper() in want}
    if not traits:
        raise SystemExit(f"No traits matched --traits={args.traits!r}; known={list(TRAITS)}")

    print(f"Fetching Pan-UKBB chr{chrom} slices (multi-ancestry): {list(traits)}\n")
    for trait, key in traits.items():
        fetch_trait(trait, key, chrom, out_dir)
    print("\nDone.")


if __name__ == "__main__":
    main()
