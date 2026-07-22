"""
Science deepen: pathway enrichment with BH-FDR, LD-block permutation, leave-one-locus-out.

Outputs:
  results/tables/pathway_enrichment_fdr_genomewide.csv
  results/tables/pathway_lolo_sensitivity_genomewide.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

ROOT = Path(__file__).resolve().parents[1]
SEED = 719
RNG = np.random.default_rng(SEED)
TRAITS = ["CAD", "T2D", "BMI", "LDL"]
ANN = ROOT / "data/annotations"
TAB = ROOT / "results/tables"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Pathway FDR + LD-block null + LOLO.")
    p.add_argument(
        "--master",
        default=str(ROOT / "data/modeling/master_variant_table_genomewide_genomewide.parquet"),
    )
    p.add_argument("--tag", default="genomewide")
    p.add_argument("--n-perm", type=int, default=200, help="LD-block permutations per trait")
    p.add_argument("--block-mb", type=float, default=1.0, help="Genomic bin size (Mb) if no LDetect")
    p.add_argument("--high-i2-thresh", type=float, default=0.25)
    p.add_argument("--lolo-window-bp", type=int, default=250_000)
    p.add_argument("--min-assoc", type=int, default=5)
    return p.parse_args()


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(1, n + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out = np.empty(n, dtype=float)
    out[order] = q
    return out


def load_variant_gene() -> pd.DataFrame:
    frames = []
    for chrom in range(1, 23):
        path = ANN / f"variant_to_gene.chr{chrom}.parquet"
        if path.exists():
            frames.append(pd.read_parquet(path))
    if not frames:
        legacy = ANN / "variant_to_gene.parquet"
        if legacy.exists():
            frames.append(pd.read_parquet(legacy))
    if not frames:
        raise SystemExit("No variant_to_gene annotations found")
    return pd.concat(frames, ignore_index=True).drop_duplicates(["variant_id", "ensg"])


def load_gene_pathway() -> pd.DataFrame:
    # gene_to_pathway files appear identical across chroms; load once
    path = ANN / "gene_to_pathway.parquet"
    if not path.exists():
        path = ANN / "gene_to_pathway.chr22.parquet"
    g = pd.read_parquet(path)
    return g.drop_duplicates(["ensg", "rid"])


def load_associated_labels(master: Path) -> pd.DataFrame:
    import pyarrow.dataset as ds

    dataset = ds.dataset(str(master), format="parquet")
    cols = ["variant_id", "trait", "I2", "y_high_I2", "associated"]
    use = [c for c in cols if c in set(dataset.schema.names)]
    filt = ds.field("associated") == True  # noqa: E712
    scanner = dataset.scanner(columns=use, filter=filt, batch_size=1_000_000)
    chunks = []
    n = 0
    for batch in scanner.to_batches():
        chunks.append(batch.to_pandas())
        n += len(chunks[-1])
        if n % 2_000_000 < 1_000_000:
            print(f"  labels loaded {n:,} ...", flush=True)
    out = pd.concat(chunks, ignore_index=True).drop_duplicates(["variant_id", "trait"])
    print(f"  associated labels: {len(out):,}", flush=True)
    return out


def parse_variant_id(vid: str) -> tuple[str, int] | None:
    parts = str(vid).split(":")
    if len(parts) < 2:
        return None
    try:
        return parts[0], int(parts[1])
    except ValueError:
        return None


def assign_blocks(variant_ids: pd.Series, block_mb: float) -> pd.Series:
    block_bp = int(block_mb * 1_000_000)
    blocks = []
    for vid in variant_ids:
        parsed = parse_variant_id(vid)
        if parsed is None:
            blocks.append("NA")
            continue
        chrom, pos = parsed
        blocks.append(f"{chrom}:{pos // block_bp}")
    return pd.Series(blocks, index=variant_ids.index)


def fisher_enrichment(n_sig_in: int, n_sig: int, n_pw: int, n_bg: int) -> float:
    # contingency: [[in_pw & sig, not_in_pw & sig], [in_pw & not_sig, not_in_pw & not_sig]]
    a = n_sig_in
    b = n_sig - n_sig_in
    c = n_pw - n_sig_in
    d = n_bg - n_sig - c
    if min(a + b, a + c, b + d, c + d) <= 0:
        return 1.0
    _, p = fisher_exact([[a, b], [c, d]], alternative="greater")
    return float(p)


def main() -> None:
    args = parse_args()
    TAB.mkdir(parents=True, exist_ok=True)

    print("Loading annotations ...", flush=True)
    v2g = load_variant_gene()
    g2p = load_gene_pathway()
    print(f"  variant_to_gene={len(v2g):,}  gene_to_pathway={len(g2p):,}", flush=True)

    print("Loading associated labels ...", flush=True)
    labels = load_associated_labels(Path(args.master))
    labels["high_i2"] = (labels["I2"] >= args.high_i2_thresh) | (labels["y_high_I2"] == 1)

    # Map variants -> genes -> pathways at gene level per trait
    vg = labels.merge(v2g, on="variant_id", how="inner")
    # gene-level: gene is high-I2 if any associated variant is high-I2
    gene_trait = (
        vg.groupby(["trait", "ensg"], as_index=False)
        .agg(high_i2=("high_i2", "max"), n_var=("variant_id", "nunique"))
    )
    gene_pw = gene_trait.merge(g2p, on="ensg", how="inner")

    fdr_rows = []
    lolo_rows = []

    for trait in TRAITS:
        print(f"\n=== {trait} ===", flush=True)
        gt = gene_trait[gene_trait["trait"] == trait].copy()
        gp = gene_pw[gene_pw["trait"] == trait].copy()
        if gt.empty or gp.empty:
            print(f"  skip {trait}: empty gene/pathway map")
            continue

        bg_genes = set(gt["ensg"])
        sig_genes = set(gt.loc[gt["high_i2"] == 1, "ensg"])
        n_bg = len(bg_genes)
        n_sig = len(sig_genes)
        print(f"  genes bg={n_bg} high_I2={n_sig}", flush=True)

        # Observed enrichment per pathway
        pw_stats = []
        for rid, sub in gp.groupby("rid"):
            pw_genes = set(sub["ensg"])
            n_pw = len(pw_genes)
            n_sig_in = len(pw_genes & sig_genes)
            if n_pw < 3:
                continue
            p = fisher_enrichment(n_sig_in, n_sig, n_pw, n_bg)
            pathway = str(sub["pathway"].iloc[0])
            pw_stats.append(
                {
                    "trait": trait,
                    "rid": rid,
                    "pathway": pathway,
                    "n_genes": n_pw,
                    "n_high_i2_genes": n_sig_in,
                    "frac_high_i2": n_sig_in / n_pw if n_pw else np.nan,
                    "fisher_p": p,
                }
            )
        obs = pd.DataFrame(pw_stats)
        if obs.empty:
            continue
        obs["bh_q"] = bh_fdr(obs["fisher_p"].to_numpy())

        # LD-block permutation null on gene high_i2 labels via variant blocks
        # Assign each gene a representative block from its variants
        vt = vg[vg["trait"] == trait][["variant_id", "ensg", "high_i2"]].copy()
        vt["block"] = assign_blocks(vt["variant_id"], args.block_mb)
        gene_block = (
            vt.groupby("ensg")
            .agg(block=("block", lambda s: s.mode().iloc[0] if len(s) else "NA"), high_i2=("high_i2", "max"))
            .reset_index()
        )
        blocks = gene_block["block"].to_numpy()
        labels_arr = gene_block["high_i2"].astype(int).to_numpy()
        ensg_list = gene_block["ensg"].tolist()
        ensg_to_idx = {g: i for i, g in enumerate(ensg_list)}

        # Precompute pathway gene index lists
        pw_idx = {}
        for rid, sub in gp.groupby("rid"):
            idxs = [ensg_to_idx[g] for g in sub["ensg"].unique() if g in ensg_to_idx]
            if len(idxs) >= 3:
                pw_idx[rid] = np.array(idxs, dtype=int)

        # Unique blocks for permutation
        uniq_blocks = pd.unique(blocks)
        block_to_genes = {b: np.where(blocks == b)[0] for b in uniq_blocks}

        # Observed max -log10 p for empirical; store per-pathway perm p
        perm_ge = {rid: 0 for rid in pw_idx}
        for pi in range(args.n_perm):
            # shuffle high_i2 labels across blocks (keep within-block structure by swapping block labels)
            block_labels = {}
            for b, idxs in block_to_genes.items():
                block_labels[b] = int(labels_arr[idxs].max()) if len(idxs) else 0
            perm_block_vals = RNG.permutation(list(block_labels.values()))
            new_block_label = dict(zip(block_labels.keys(), perm_block_vals))
            perm_lab = np.array([new_block_label[b] for b in blocks], dtype=int)
            n_sig_p = int(perm_lab.sum())
            if n_sig_p == 0:
                continue
            for rid, idxs in pw_idx.items():
                n_pw = len(idxs)
                n_sig_in = int(perm_lab[idxs].sum())
                p_perm = fisher_enrichment(n_sig_in, n_sig_p, n_pw, n_bg)
                # compare to observed fisher_p for this rid
                obs_p = float(obs.loc[obs["rid"] == rid, "fisher_p"].iloc[0])
                if p_perm <= obs_p:
                    perm_ge[rid] += 1
            if (pi + 1) % 50 == 0:
                print(f"  perm {pi + 1}/{args.n_perm}", flush=True)

        obs["ld_block_perm_p"] = obs["rid"].map(
            lambda r: (perm_ge.get(r, 0) + 1) / (args.n_perm + 1)
        )
        obs["ld_block_perm_q"] = bh_fdr(obs["ld_block_perm_p"].fillna(1.0).to_numpy())
        obs["fdr_significant"] = (obs["bh_q"] < 0.05) & (obs["ld_block_perm_q"] < 0.05)
        fdr_rows.append(obs)

        n_sig_pw = int(obs["fdr_significant"].sum())
        print(f"  pathways tested={len(obs)} FDR+LD significant={n_sig_pw}", flush=True)

        # Leave-one-locus-out on top pathways by mean frac / fisher
        top = obs.sort_values("fisher_p").head(20)
        # Need variant positions for top locus per pathway
        for _, row in top.iterrows():
            rid = row["rid"]
            pw_genes = set(gp.loc[gp["rid"] == rid, "ensg"])
            vars_pw = vt[vt["ensg"].isin(pw_genes)].copy()
            if vars_pw.empty:
                continue
            # top locus = variant with high_i2, pick max count block or first high_i2
            high = vars_pw[vars_pw["high_i2"] == 1]
            if high.empty:
                high = vars_pw
            # choose locus center as mode position of high_i2 variants
            coords = [parse_variant_id(v) for v in high["variant_id"]]
            coords = [c for c in coords if c is not None]
            if not coords:
                continue
            # pick the chrom:pos of the densest 1 variant — use first high_i2 by I2 proxy (any)
            chrom, pos = coords[0]
            # drop variants within window of this locus
            keep_mask = []
            for vid in vars_pw["variant_id"]:
                parsed = parse_variant_id(vid)
                if parsed is None:
                    keep_mask.append(True)
                    continue
                c, p = parsed
                if c == chrom and abs(p - pos) <= args.lolo_window_bp:
                    keep_mask.append(False)
                else:
                    keep_mask.append(True)
            vars_kept = vars_pw.loc[keep_mask]
            genes_kept = set(vars_kept["ensg"])
            # recompute gene high_i2 without dropped locus variants
            gene_high = (
                vars_kept.groupby("ensg")["high_i2"].max().reindex(list(pw_genes)).fillna(0)
            )
            n_pw = len(pw_genes)
            n_sig_in_lolo = int((gene_high == 1).sum())
            # background unchanged
            p_lolo = fisher_enrichment(n_sig_in_lolo, n_sig, n_pw, n_bg)
            fragile = (row["fisher_p"] < 0.05) and (p_lolo > 0.05 or n_sig_in_lolo < row["n_high_i2_genes"] * 0.5)
            lolo_rows.append(
                {
                    "trait": trait,
                    "rid": rid,
                    "pathway": row["pathway"],
                    "fisher_p": row["fisher_p"],
                    "bh_q": row["bh_q"],
                    "ld_block_perm_q": row["ld_block_perm_q"],
                    "n_high_i2_genes": row["n_high_i2_genes"],
                    "n_high_i2_genes_lolo": n_sig_in_lolo,
                    "fisher_p_lolo": p_lolo,
                    "dropped_locus": f"{chrom}:{pos}",
                    "fragile": fragile,
                    "fdr_significant": bool(row["fdr_significant"]),
                }
            )

    if fdr_rows:
        fdr = pd.concat(fdr_rows, ignore_index=True)
        fdr_path = TAB / f"pathway_enrichment_fdr_{args.tag}.csv"
        fdr.sort_values(["trait", "bh_q", "fisher_p"]).to_csv(fdr_path, index=False, float_format="%.6g")
        print(f"\nSaved {fdr_path} ({len(fdr)} rows)")
        for trait in TRAITS:
            sub = fdr[fdr["trait"] == trait]
            n = int(sub["fdr_significant"].sum()) if len(sub) else 0
            print(f"  {trait}: {n} pathways survive BH-FDR q<0.05 AND LD-block perm q<0.05")
    else:
        print("No FDR rows produced")

    if lolo_rows:
        lolo = pd.DataFrame(lolo_rows)
        lolo_path = TAB / f"pathway_lolo_sensitivity_{args.tag}.csv"
        lolo.to_csv(lolo_path, index=False, float_format="%.6g")
        print(f"Saved {lolo_path} ({len(lolo)} rows)")


if __name__ == "__main__":
    main()
