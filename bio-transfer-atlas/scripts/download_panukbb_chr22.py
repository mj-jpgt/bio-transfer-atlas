"""
FAIRGEN-Open Stage 8.1: Fetch Pan-UKBB chr22 multi-ancestry summary stats
==========================================================================
Pan-UKBB flat files are ~2.4 GB bgzipped, sorted by chrom/pos, with a tabix
index. We only need chr22, so we:
  1. download the small .tbi index (~2.2 MB),
  2. parse it to find the compressed byte offset where chr22 begins,
  3. HTTP Range-GET from that offset to EOF (~30-40 MB, complete BGZF blocks),
  4. decompress (gzip reads concatenated BGZF members), keep chr==22 rows,
  5. fetch the header separately from the first block,
  6. save per-ancestry betas/SEs as parquet.

No aws CLI / tabix binary required — only curl + Python stdlib.

Output:
  data/raw/panukbb/chr22/<trait>.chr22.parquet
"""
import gzip
import io
import struct
import subprocess
from pathlib import Path

import pandas as pd

root    = Path(__file__).resolve().parents[1]
out_dir = root / "data/raw/panukbb/chr22"
tbi_dir = root / "data/raw/panukbb/tbi"
out_dir.mkdir(parents=True, exist_ok=True)
tbi_dir.mkdir(parents=True, exist_ok=True)

S3 = "https://pan-ukb-us-east-1.s3.amazonaws.com"

TRAITS = {
    "T2D": "sumstats_flat_files/phecode-250.2-both_sexes.tsv.bgz",
    "CAD": "sumstats_flat_files/phecode-411.4-both_sexes.tsv.bgz",
    "BMI": "sumstats_flat_files/continuous-21001-both_sexes-irnt.tsv.bgz",
    "LDL": "sumstats_flat_files/continuous-LDLC-both_sexes-medadj_irnt.tsv.bgz",
}
ANCS = ["AFR", "AMR", "CSA", "EAS", "EUR", "MID"]
CHR22 = "22"


def curl_range(url: str, start: int, end, out_path: Path):
    """HTTP Range GET [start, end] (end inclusive, or '' for EOF) -> file."""
    from bta_curl import curl_bin

    rng = f"{start}-{end if end is not None else ''}"
    cmd = [curl_bin(), "-sL", "-r", rng, "-o", str(out_path), url]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"curl range {rng} failed: {r.stderr[-300:]}")
    return out_path


def curl_full(url: str, out_path: Path):
    from bta_curl import curl_bin

    cmd = [curl_bin(), "-sL", "-o", str(out_path), url]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"curl failed: {r.stderr[-300:]}")
    return out_path


def parse_tbi_chr22_coffset(tbi_path: Path) -> int:
    """Parse tabix index, return min compressed offset (coffset) for chr22."""
    raw = gzip.open(tbi_path, "rb").read()
    p = 0
    magic = raw[p:p+4]; p += 4
    assert magic == b"TBI\x01", f"bad tbi magic {magic!r}"
    (n_ref, fmt, col_seq, col_beg, col_end, meta, skip, l_nm) = \
        struct.unpack("<8i", raw[p:p+32]); p += 32
    names_blob = raw[p:p+l_nm]; p += l_nm
    names = names_blob.split(b"\x00")
    names = [n.decode() for n in names if n]

    # find index of chr22 contig (accept '22' or 'chr22')
    target_idx = None
    for i, nm in enumerate(names):
        if nm in ("22", "chr22"):
            target_idx = i; break
    if target_idx is None:
        raise RuntimeError(f"chr22 not found in tbi contigs: {names}")

    # Use ONLY the linear index (intervals). Bin chunks include a metadata
    # pseudo-bin (37450) whose "chunks" store record counts, not offsets, which
    # would corrupt a min over chunk offsets. The linear index ioff>>16 gives
    # the compressed block offset of the first record in each 16kb window.
    min_coffset = None
    for ref in range(n_ref):
        (n_bin,) = struct.unpack("<i", raw[p:p+4]); p += 4
        for _ in range(n_bin):
            (bin_id, n_chunk) = struct.unpack("<Ii", raw[p:p+8]); p += 8
            p += 16 * n_chunk  # skip chunks
        (n_intv,) = struct.unpack("<i", raw[p:p+4]); p += 4
        ioffs = struct.unpack(f"<{n_intv}Q", raw[p:p+8*n_intv]); p += 8 * n_intv
        if ref == target_idx:
            nz = [(o >> 16) for o in ioffs if o > 0]
            min_coffset = min(nz) if nz else 0
    if min_coffset is None:
        raise RuntimeError("could not derive chr22 coffset")
    return int(min_coffset)


def bgzf_decompress(data: bytes) -> bytes:
    """Decompress concatenated BGZF blocks, stopping at the last COMPLETE block.
    Tolerates a truncated final block (e.g. when only a byte-prefix was fetched).
    """
    import zlib
    out = bytearray()
    n = len(data)
    p = 0
    while p + 18 <= n:
        if data[p] != 0x1f or data[p+1] != 0x8b:
            break
        xlen = struct.unpack("<H", data[p+10:p+12])[0]
        extra = data[p+12:p+12+xlen]
        bsize = None
        q = 0
        while q + 4 <= len(extra):
            si1, si2, slen = extra[q], extra[q+1], struct.unpack("<H", extra[q+2:q+4])[0]
            if si1 == 66 and si2 == 67:  # 'BC'
                bsize = struct.unpack("<H", extra[q+4:q+6])[0]
            q += 4 + slen
        if bsize is None:
            break
        block_len = bsize + 1
        if p + block_len > n:
            break  # truncated final block
        cdata = data[p+12+xlen : p+block_len-8]
        try:
            out.extend(zlib.decompress(cdata, -15))
        except Exception:
            break
        p += block_len
    return bytes(out)


def read_header(url: str) -> list:
    """Fetch first ~200 KB, decompress complete blocks, return header columns."""
    tmp = out_dir / "_hdr.bgz"
    curl_range(url, 0, 200_000, tmp)
    txt = bgzf_decompress(tmp.read_bytes()).decode("utf-8", errors="replace")
    tmp.unlink(missing_ok=True)
    first_line = txt.split("\n", 1)[0]
    return first_line.rstrip("\r").split("\t")


def fetch_trait(trait: str, key: str) -> Path:
    out_pq = out_dir / f"{trait}.chr22.parquet"
    if out_pq.exists():
        print(f"  {trait}: cached -> {out_pq.name}")
        return out_pq
    url = f"{S3}/{key}"
    tbi_url = f"{S3}/sumstats_flat_files_tabix/{Path(key).name}.tbi"

    print(f"  {trait}: downloading tabix index ...", flush=True)
    tbi_path = tbi_dir / f"{trait}.tbi"
    curl_full(tbi_url, tbi_path)
    coffset = parse_tbi_chr22_coffset(tbi_path)
    print(f"  {trait}: chr22 starts at compressed byte {coffset:,}")

    print(f"  {trait}: fetching header ...", flush=True)
    header = read_header(url)
    print(f"  {trait}: {len(header)} columns")

    print(f"  {trait}: range-downloading chr22 slice (coffset..EOF) ...", flush=True)
    tail = out_dir / f"{trait}.chr22.bgz"
    curl_range(url, coffset, None, tail)
    size_mb = tail.stat().st_size / 1e6
    print(f"  {trait}: downloaded {size_mb:.1f} MB slice")

    # decompress slice (range ends at real EOF -> complete BGZF blocks)
    text = bgzf_decompress(tail.read_bytes()).decode("utf-8", errors="replace")
    # drop a possibly-partial first line (block boundary may split a record)
    lines = text.split("\n")
    if lines:
        lines = lines[1:]
    df = pd.read_csv(
        io.StringIO("\n".join(lines)), sep="\t", names=header,
        dtype=str, low_memory=False, on_bad_lines="skip",
    )
    df = df[df["chr"] == CHR22].copy()

    keep = ["chr", "pos", "ref", "alt"]
    for a in ANCS:
        for pre in ("beta", "se", "pval", "low_confidence"):
            c = f"{pre}_{a}"
            if c in df.columns:
                keep.append(c)
    for c in ("beta_meta", "se_meta", "pval_meta", "af_meta"):
        if c in df.columns:
            keep.append(c)
    df = df[[c for c in keep if c in df.columns]].copy()
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
    df = df.dropna(subset=["pos"]); df["pos"] = df["pos"].astype(int)

    df.to_parquet(out_pq, index=False)
    tail.unlink(missing_ok=True)
    print(f"  {trait}: saved {len(df):,} chr22 variants -> {out_pq.name}")
    return out_pq


if __name__ == "__main__":
    print("Fetching Pan-UKBB chr22 slices (multi-ancestry)\n")
    for trait, key in TRAITS.items():
        fetch_trait(trait, key)
    print("\nDone.")
