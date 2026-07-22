"""
Shared GWAS I/O utilities used by download/label scripts.
"""
from __future__ import annotations

import io
import zipfile
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional

import pandas as pd
import requests


def normalize_chr(value) -> str:
    s = str(value).strip()
    if s.lower().startswith("chr"):
        s = s[3:]
    return s


def read_tsv_gz_chunks(
    source: str,
    chunksize: int = 400_000,
    usecols: Optional[List[str]] = None,
) -> Iterator[pd.DataFrame]:
    """
    Stream-read a gzip-compressed TSV from URL or local path.
    """
    yield from pd.read_csv(
        source,
        sep="\t",
        compression="gzip",
        dtype=str,
        chunksize=chunksize,
        low_memory=False,
        usecols=usecols,
    )


def read_zip_member_chunks(
    zip_url: str,
    member_name: Optional[str] = None,
    chunksize: int = 250_000,
) -> Iterator[pd.DataFrame]:
    """
    Download a ZIP and stream-read one TSV-like member in chunks.
    """
    resp = requests.get(zip_url, timeout=180)
    resp.raise_for_status()
    with zipfile.ZipFile(io.BytesIO(resp.content)) as zf:
        names = [n for n in zf.namelist() if not n.endswith("/")]
        if not names:
            raise RuntimeError(f"ZIP had no files: {zip_url}")
        chosen = member_name
        if chosen is None:
            # Prefer explicit text-like members.
            txt = [n for n in names if n.lower().endswith((".txt", ".tsv", ".txt.gz", ".tsv.gz"))]
            chosen = txt[0] if txt else names[0]
        with zf.open(chosen) as f:
            if chosen.lower().endswith(".gz"):
                import gzip

                with gzip.open(f, mode="rt", encoding="utf-8", errors="replace") as g:
                    yield from pd.read_csv(
                        g,
                        sep="\t",
                        dtype=str,
                        chunksize=chunksize,
                        low_memory=False,
                    )
            else:
                text = io.TextIOWrapper(f, encoding="utf-8", errors="replace")
                yield from pd.read_csv(
                    text,
                    sep="\t",
                    dtype=str,
                    chunksize=chunksize,
                    low_memory=False,
                )


def filter_chr_and_project(
    chunks: Iterable[pd.DataFrame],
    chrom: str,
    chrom_col: str,
    col_map: Dict[str, str],
) -> pd.DataFrame:
    """
    Keep rows on a target chromosome and project to canonical columns.
    """
    out: List[pd.DataFrame] = []
    target = normalize_chr(chrom)
    required = list(col_map)
    for chunk in chunks:
        missing = [c for c in required if c not in chunk.columns]
        if missing:
            raise RuntimeError(f"Missing required columns: {missing}")
        c = chunk[required].copy()
        c = c[c[chrom_col].map(normalize_chr) == target]
        if len(c) == 0:
            continue
        c = c.rename(columns=col_map)
        out.append(c)
    if not out:
        return pd.DataFrame(columns=list(col_map.values()))
    return pd.concat(out, ignore_index=True)


def finalize_gwas_frame(df: pd.DataFrame) -> pd.DataFrame:
    """
    Apply canonical typing/cleanup for shared GWAS columns.
    """
    for c in ("chr", "ref", "alt"):
        if c in df.columns:
            df[c] = df[c].astype(str)
    if "chr" in df.columns:
        df["chr"] = df["chr"].map(normalize_chr)
    if "pos" in df.columns:
        df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
        df = df.dropna(subset=["pos"]).copy()
        df["pos"] = df["pos"].astype(int)
    for c in ("beta", "se", "pval", "af", "af_alt", "sebeta"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    for c in ("ref", "alt"):
        if c in df.columns:
            df[c] = df[c].str.upper()
    return df
