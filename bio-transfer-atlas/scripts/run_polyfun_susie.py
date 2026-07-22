"""
Phase B1: PolyFun+SuSiE runner (real CS/PIP when R+susieR available).

Runs susieR::susie_rss per LD-block on associated variants using z_meta from the
master/labels table. Writes data/labels/susie/susie_{trait}_{anc}.parquet with
columns variant_id, pip, in_cs.

If susieR is missing, writes heuristic_fallback status (legacy behavior).
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--trait", default="T2D")
    p.add_argument("--chrom", default="22", help="Restrict to chrom prefix in variant_id if set")
    p.add_argument("--anc", default="EUR", help="Label only; z_meta is meta-analytic")
    p.add_argument("--out-dir", default=str(ROOT / "data/labels/susie"))
    p.add_argument(
        "--master",
        default=str(ROOT / "data/modeling/master_variant_table_genomewide_genomewide.parquet"),
    )
    p.add_argument(
        "--ld-blocks",
        default=str(ROOT / "data/modeling/ld_block_assignments_genomewide.parquet"),
    )
    p.add_argument("--max-blocks", type=int, default=200, help="Cap blocks for runtime")
    p.add_argument("--jobs", type=int, default=8)
    p.add_argument("--L", type=int, default=10, help="SuSiE L")
    p.add_argument(
        "--prefer-real-ld",
        action="store_true",
        help="Extract PLINK LD correlation matrix from 1000G score pfiles when possible",
    )
    p.add_argument(
        "--pfile-root",
        default=str(ROOT / "data/interim/1000g_grch38"),
        help="Directory with chr{N}.score pfiles",
    )
    return p.parse_args()


def have_susier() -> bool:
    if not shutil.which("Rscript"):
        return False
    r = subprocess.run(
        ["Rscript", "-e", "cat(requireNamespace('susieR', quietly=TRUE))"],
        capture_output=True,
        text=True,
    )
    return "TRUE" in (r.stdout or "").upper()


def load_associated(master: Path, trait: str, chrom: str | None) -> pd.DataFrame:
    slim = ROOT / "data/labels/_tmp_associated_labels.parquet"
    df = None
    if slim.exists():
        df = pd.read_parquet(slim)
        if "z_meta" not in df.columns:
            df = None
    if df is None:
        import pyarrow.dataset as ds

        dataset = ds.dataset(str(master), format="parquet")
        cols = [c for c in ["variant_id", "trait", "z_meta", "I2", "associated"] if c in set(dataset.schema.names)]
        filt = ds.field("associated") == True  # noqa: E712
        chunks = []
        for batch in dataset.scanner(columns=cols, filter=filt, batch_size=200_000).to_batches():
            chunks.append(batch.to_pandas())
        df = pd.concat(chunks, ignore_index=True)
    df = df[df["trait"].astype(str) == trait].copy()
    if chrom:
        # variant_id like chr22:pos:ref:alt or 22:pos:...
        pref = (df["variant_id"].astype(str).str.startswith(f"chr{chrom}:")) | (
            df["variant_id"].astype(str).str.startswith(f"{chrom}:")
        )
        df = df[pref]
    return df.drop_duplicates("variant_id")


def _plink2() -> str | None:
    candidates = [
        str(ROOT / "tools/plink2/plink2"),
        str(ROOT / "bin/plink2"),
        "plink2",
    ]
    for cand in candidates:
        p = Path(cand)
        if p.is_file() or (cand == "plink2" and shutil.which(cand)):
            return str(p) if p.is_file() else cand
    try:
        from bta_plink import plink2_bin

        return plink2_bin()
    except Exception:
        return None


def extract_ld_matrix(
    block_df: pd.DataFrame,
    work: Path,
    pfile_root: Path,
    chrom: str,
    max_snps: int = 80,
) -> tuple[np.ndarray, list[str], str, dict]:
    """
    Build correlation matrix R via PLINK --r-unphased (signed) with r2 fallback.
    Returns (R, kept_vids, mode_tag, ld_meta).
    """
    empty_meta = {"ld_type": "identity", "fallback_mode": "identity_ld", "ridge_added": 0.0, "plink_flag": "none"}
    vids = block_df["variant_id"].astype(str).tolist()
    if len(vids) > max_snps:
        # keep highest |z|
        order = np.argsort(-np.abs(block_df["z_meta"].to_numpy(dtype=float)))
        keep_idx = order[:max_snps]
        block_df = block_df.iloc[keep_idx].reset_index(drop=True)
        vids = block_df["variant_id"].astype(str).tolist()

    plink = _plink2()
    pfile = pfile_root / f"chr{chrom}.score"
    if plink is None or not (Path(str(pfile) + ".pgen").exists() or Path(str(pfile) + ".bed").exists()):
        n = len(vids)
        return (np.eye(n), vids, "identity", empty_meta)

    snplist = work / "snps.txt"
    # PLINK IDs may be chr:pos:ref:alt — write as-is
    snplist.write_text("\n".join(vids) + "\n", encoding="utf-8")
    out_prefix = work / "ld"
    # Primary: signed unphased correlations (--r-unphased) required by susie_rss
    cmd = [
        plink,
        "--pfile",
        str(pfile),
        "--extract",
        str(snplist),
        "--r-unphased",
        "square",
        "--out",
        str(out_prefix),
    ]
    subprocess.run(cmd, capture_output=True, text=True)
    signed_mats = [
        Path(str(out_prefix) + ".unphased.vcor1"),  # PLINK2 a.7+ signed --r-unphased
        Path(str(out_prefix) + ".unphased.vcor"),
        Path(str(out_prefix) + ".vcor1"),
        Path(str(out_prefix) + ".vcor"),
    ]
    unsigned_mats = [
        Path(str(out_prefix) + ".unphased.vcor2"),
        Path(str(out_prefix) + ".vcor2"),
        Path(str(out_prefix) + ".ld"),
    ]
    signed_vars = [
        Path(str(out_prefix) + ".unphased.vcor1.vars"),
        Path(str(out_prefix) + ".unphased.vcor.vars"),
        Path(str(out_prefix) + ".vcor1.vars"),
        Path(str(out_prefix) + ".vcor.vars"),
    ]
    unsigned_vars = [
        Path(str(out_prefix) + ".unphased.vcor2.vars"),
        Path(str(out_prefix) + ".vcor2.vars"),
    ]
    mat_path = next((p for p in signed_mats if p.exists()), None)
    vars_path = next((p for p in signed_vars if p.exists()), None)
    mode = "signed_ld"
    if mat_path is None:
        # r2 fallback is NOT primary (unsigned |r|)
        cmd2 = [
            plink,
            "--pfile",
            str(pfile),
            "--extract",
            str(snplist),
            "--r2-unphased",
            "square",
            "--out",
            str(out_prefix),
        ]
        subprocess.run(cmd2, capture_output=True, text=True)
        mat_path = next((p for p in unsigned_mats if p.exists()), None)
        vars_path = next((p for p in unsigned_vars if p.exists()), None)
        mode = "unsigned_abs_r"
        if mat_path is None:
            n = len(vids)
            return (np.eye(n), vids, "identity", empty_meta)

    try:
        R = np.loadtxt(mat_path)
        if R.ndim != 2 or R.shape[0] != R.shape[1]:
            raise ValueError("bad shape")
        if vars_path is not None and vars_path.exists():
            kept = [ln.strip() for ln in vars_path.read_text(encoding="utf-8").splitlines() if ln.strip()]
            if kept and len(kept) == R.shape[0]:
                vids = kept
        R = np.nan_to_num(R, nan=0.0, posinf=0.0, neginf=0.0)
        R = 0.5 * (R + R.T)
        np.fill_diagonal(R, 1.0)
        sym_ok = bool(np.allclose(R, R.T, atol=1e-8))
        diag_ok = bool(np.allclose(np.diag(R), 1.0, atol=1e-6))
        min_before = float(np.min(np.linalg.eigvalsh(R))) if R.size else float("nan")
        has_neg = bool(np.nanmin(R) < -1e-8) if mode == "signed_ld" else False
        if mode == "unsigned_abs_r":
            R = np.sqrt(np.clip(R, 0, 1))
            np.fill_diagonal(R, 1.0)
        ridge = 0.0
        evals = np.linalg.eigvalsh(R)
        min_ev = float(np.min(evals))
        if min_ev < 1e-6:
            ridge = float(1e-6 - min_ev)
            R = R + ridge * np.eye(R.shape[0])
            np.fill_diagonal(R, 1.0)
        else:
            R = R + 1e-3 * np.eye(R.shape[0])
            ridge = 1e-3
            np.fill_diagonal(R, 1.0)
        if R.shape[0] != len(vids):
            n = len(vids)
            return (
                np.eye(n),
                vids,
                "identity",
                {"ld_type": "identity", "fallback_mode": "identity_ld", "ridge_added": 0.0},
            )
        meta = {
            "ld_type": "signed_r" if mode == "signed_ld" else ("unsigned_abs_r" if mode == "unsigned_abs_r" else "identity"),
            "plink_flag": "--r-unphased square" if mode == "signed_ld" else ("--r2-unphased square" if mode == "unsigned_abs_r" else "none"),
            "min_eigenvalue_before": min_before,
            "ridge_added": ridge,
            "variant_order_verified": True,
            "symmetric_ok": sym_ok,
            "diag_ok": diag_ok,
            "has_negative_entry_heuristic": has_neg,
            "fallback_mode": "signed_ld" if mode == "signed_ld" else ("unsigned_abs_r" if mode == "unsigned_abs_r" else "identity_ld"),
        }
        return R, vids, mode, meta
    except Exception:
        n = len(vids)
        return (
            np.eye(n),
            vids,
            "identity",
            {"ld_type": "identity", "fallback_mode": "identity_ld", "ridge_added": 0.0},
        )

def run_susie_block(
    block_df: pd.DataFrame,
    L: int,
    work: Path,
    *,
    prefer_real_ld: bool = False,
    pfile_root: Path | None = None,
    chrom: str = "22",
) -> pd.DataFrame:
    """Run susie_rss with real LD when available; label fallback_mode explicitly."""
    z = block_df["z_meta"].to_numpy(dtype=float)
    vids = block_df["variant_id"].astype(str).tolist()
    if len(vids) < 2 or not np.isfinite(z).any():
        return pd.DataFrame(columns=["variant_id", "pip", "in_cs", "fallback_mode"])

    mode = "identity"
    R = np.eye(len(vids))
    ld_meta: dict = {
        "ld_type": "identity",
        "plink_flag": "none",
        "min_eigenvalue_before": 1.0,
        "ridge_added": 0.0,
        "variant_order_verified": True,
        "fallback_mode": "identity_ld",
    }
    if prefer_real_ld and pfile_root is not None:
        R, vids, mode, ld_meta = extract_ld_matrix(block_df, work, pfile_root, chrom)
        # realign z to kept vids
        zmap = dict(zip(block_df["variant_id"].astype(str), z))
        z = np.array([zmap[v] for v in vids], dtype=float)

    zpath = work / "z.csv"
    rpath = work / "R.csv"
    opath = work / "out.csv"
    # Symmetrize and lightly ridge for numerical stability
    R = 0.5 * (R + R.T)
    np.fill_diagonal(R, 1.0)
    R = R + 1e-3 * np.eye(len(vids))
    np.fill_diagonal(R, 1.0)
    pd.DataFrame({"z": z}).to_csv(zpath, index=False)
    pd.DataFrame(R).to_csv(rpath, index=False, header=False)
    # Reference LD: do NOT estimate residual variance
    est_res = "FALSE" if mode in ("signed_ld", "unsigned_abs_r") else "TRUE"
    rscript = f"""
suppressPackageStartupMessages(library(susieR))
z <- as.numeric(read.csv('{zpath.as_posix()}')$z)
R <- as.matrix(read.csv('{rpath.as_posix()}', header=FALSE))
R <- (R + t(R)) / 2
diag(R) <- 1
n <- length(z)
res <- tryCatch(
  susie_rss(z = z, R = R, L = {L}, estimate_residual_variance = {est_res}, check_prior = FALSE),
  error = function(e) NULL
)
if (is.null(res)) {{
  pip <- rep(0, n)
  in_cs <- rep(FALSE, n)
}} else {{
  pip <- as.numeric(res$pip)
  in_cs <- rep(FALSE, n)
  if (!is.null(res$sets) && !is.null(res$sets$cs) && length(res$sets$cs) > 0) {{
    for (s in res$sets$cs) {{
      in_cs[s] <- TRUE
    }}
  }} else {{
    in_cs <- pip >= 0.5
  }}
}}
write.csv(data.frame(idx=seq_len(n), pip=pip, in_cs=in_cs), '{opath.as_posix()}', row.names=FALSE)
"""
    rfile = work / "run_susie.R"
    rfile.write_text(rscript, encoding="utf-8")
    r = subprocess.run(["Rscript", str(rfile)], capture_output=True, text=True)
    if mode == "signed_ld":
        fallback = "signed_ld"
    elif mode == "unsigned_abs_r":
        fallback = "unsigned_abs_r"
    else:
        fallback = "identity_ld"
    if r.returncode != 0 or not opath.exists():
        # Fallback: soft PIP from |z| ranks when susieR fails
        az = np.abs(z)
        az = np.nan_to_num(az, nan=0.0)
        if az.sum() <= 0:
            pip = np.zeros(len(vids))
        else:
            pip = az / az.sum()
        in_cs = pip >= max(float(np.quantile(pip, 0.9)), 1e-12)
        return pd.DataFrame(
            {
                "variant_id": vids,
                "pip": pip,
                "in_cs": in_cs,
                "fallback_mode": "absz_weight",
            }
        )
    out = pd.read_csv(opath)
    out["variant_id"] = [vids[i - 1] for i in out["idx"]]
    if float(out["pip"].fillna(0).max()) <= 0:
        az = np.abs(z)
        az = np.nan_to_num(az, nan=0.0)
        pip = az / az.sum() if az.sum() > 0 else np.zeros(len(vids))
        out["pip"] = pip
        out["in_cs"] = pip >= max(float(np.quantile(pip, 0.9)), 1e-12)
        fallback = "absz_weight"
    out["fallback_mode"] = fallback
    for k, v in ld_meta.items():
        if k == "fallback_mode":
            continue
        out[k] = v
    # Prefer LD-extract fallback_mode when susie succeeded with that matrix
    if fallback in ("signed_ld", "unsigned_abs_r", "identity_ld"):
        out["fallback_mode"] = ld_meta.get("fallback_mode", fallback)
    keep = [
        "variant_id",
        "pip",
        "in_cs",
        "fallback_mode",
        "ld_type",
        "plink_flag",
        "min_eigenvalue_before",
        "ridge_added",
        "variant_order_verified",
        "symmetric_ok",
        "diag_ok",
        "has_negative_entry_heuristic",
    ]
    keep = [c for c in keep if c in out.columns]
    return out[keep]


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    status = {
        "trait": args.trait,
        "chrom": args.chrom,
        "anc": args.anc,
        "polyfun": (ROOT / "tools/polyfun").exists(),
        "Rscript": shutil.which("Rscript") is not None,
        "susieR": False,
        "mode": "heuristic_fallback",
    }
    status["susieR"] = have_susier()
    status_path = out_dir / f"susie_status_{args.trait}_chr{args.chrom}_{args.anc}.json"

    if not status["susieR"]:
        status_path.write_text(json.dumps(status, indent=2), encoding="utf-8")
        print(json.dumps(status, indent=2))
        print("susieR unavailable — heuristic_fallback")
        return

    print("Loading associated + LD blocks ...", flush=True)
    labels = load_associated(Path(args.master), args.trait, args.chrom if args.chrom else None)
    blocks = pd.read_parquet(args.ld_blocks, columns=["variant_id", "ld_block"])
    df = labels.merge(blocks, on="variant_id", how="inner")
    df = df.dropna(subset=["z_meta"])
    print(f"{args.trait} chr{args.chrom}: {len(df):,} associated with blocks", flush=True)

    # Largest blocks first
    sizes = df.groupby("ld_block").size().sort_values(ascending=False)
    keep_blocks = sizes.head(args.max_blocks).index.tolist()
    df = df[df["ld_block"].isin(keep_blocks)]

    parts = []
    block_groups = list(df.groupby("ld_block"))
    pfile_root = Path(args.pfile_root)

    def _one(item):
        bname, bdf = item
        with tempfile.TemporaryDirectory(prefix="susie_") as td:
            return run_susie_block(
                bdf.reset_index(drop=True),
                args.L,
                Path(td),
                prefer_real_ld=bool(args.prefer_real_ld),
                pfile_root=pfile_root,
                chrom=str(args.chrom or "22"),
            )

    jobs = max(1, args.jobs)
    if jobs == 1:
        for item in block_groups:
            parts.append(_one(item))
    else:
        with ThreadPoolExecutor(max_workers=jobs) as ex:
            futs = [ex.submit(_one, it) for it in block_groups]
            for fut in as_completed(futs):
                parts.append(fut.result())

    if not parts:
        status["mode"] = "susie_no_blocks"
        status_path.write_text(json.dumps(status, indent=2), encoding="utf-8")
        raise SystemExit("No SuSiE outputs")

    out = pd.concat(parts, ignore_index=True)
    mode_counts = out["fallback_mode"].value_counts().to_dict() if "fallback_mode" in out.columns else {}
    if "fallback_mode" in out.columns:
        def _mode_pref(s):
            s = set(s)
            for pref in ("signed_ld", "unsigned_abs_r", "identity_ld", "absz_weight"):
                if pref in s:
                    return pref
            return list(s)[0] if s else "unknown"

        out = out.groupby("variant_id", as_index=False).agg(
            pip=("pip", "max"),
            in_cs=("in_cs", "max"),
            fallback_mode=("fallback_mode", _mode_pref),
        )
    else:
        out = out.groupby("variant_id", as_index=False).agg(pip=("pip", "max"), in_cs=("in_cs", "max"))
    out_path = out_dir / f"susie_{args.trait}_{args.anc}.parquet"
    if out_path.exists() and args.chrom:
        prev = pd.read_parquet(out_path)
        out = pd.concat([prev, out], ignore_index=True)
        if "fallback_mode" in out.columns:
            out = out.groupby("variant_id", as_index=False).agg(
                pip=("pip", "max"),
                in_cs=("in_cs", "max"),
                fallback_mode=("fallback_mode", lambda s: "signed_ld" if "signed_ld" in set(s) else list(s)[0]),
            )
        else:
            out = out.groupby("variant_id", as_index=False).agg(pip=("pip", "max"), in_cs=("in_cs", "max"))
    out.to_parquet(out_path, index=False)
    n_signed = int((out.get("fallback_mode", pd.Series(dtype=str)) == "signed_ld").sum()) if "fallback_mode" in out.columns else 0
    status["mode"] = (
        "susie_cs_signed_ld"
        if n_signed > 0.5 * len(out)
        else ("susie_cs_mixed" if n_signed > 0 else "susie_cs_nonprimary")
    )
    status["n_variants"] = int(len(out))
    status["n_in_cs"] = int(out["in_cs"].sum())
    status["n_signed_ld"] = n_signed
    status["n_real_ld"] = n_signed  # legacy alias
    status["fallback_mode_counts"] = {str(k): int(v) for k, v in mode_counts.items()}
    status["primary_requires"] = "signed_ld"
    status_path.write_text(json.dumps(status, indent=2), encoding="utf-8")
    summary_path = out_dir / "susie_real_ld_status_summary.csv"
    row = {
        "trait": args.trait,
        "anc": args.anc,
        "chrom": args.chrom,
        "mode": status["mode"],
        "n_variants": status["n_variants"],
        "n_signed_ld": n_signed,
        "n_real_ld": n_signed,
    }
    if summary_path.exists():
        prev_s = pd.read_csv(summary_path)
        prev_s = pd.concat([prev_s, pd.DataFrame([row])], ignore_index=True)
        prev_s.to_csv(summary_path, index=False)
    else:
        pd.DataFrame([row]).to_csv(summary_path, index=False)
    print(json.dumps(status, indent=2))
    print(f"Saved {out_path}")


if __name__ == "__main__":
    main()
