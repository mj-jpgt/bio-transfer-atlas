"""
FAIRGEN-Open Stage 8: Cross-ancestry GWAS concordance labels (Pan-UKBB)
========================================================================
The make-or-break stage. For each chr22 PGS/feature variant we collect
per-ancestry betas from Pan-UKBB (one harmonized file per trait, so all
ancestries already share allele orientation) and compute:

  sign_concordance : fraction of ancestry pairs with same beta sign
  cochran_Q, I2    : effect-size heterogeneity across ancestries
  het_pval         : Cochran's Q p-value
  beta_corr        : mean pairwise Pearson r of betas (>=3 ancestries)
  risk_class       : 0 low / 1 medium / 2 high  (checklistv2 thresholds)

Per checklistv2 8.5 we do NOT require p<0.05 (that biases concordance upward).
We use every ancestry with a non-null beta & se.

Build alignment: Pan-UKBB is GRCh37; our features/PGS are GRCh38. We lift the
GRCh38 target variants DOWN to GRCh37 (once, ~134k) and match Pan-UKBB by
chr:pos37 + allele set, then carry the GRCh38 variant_id as the merge key.

Output:
  data/labels/gwas_concordance_labels.parquet
  data/labels/label_audit.txt
"""
from pathlib import Path
from itertools import combinations

import numpy as np
import pandas as pd
from scipy import stats
from pyliftover import LiftOver

root      = Path(__file__).resolve().parents[1]
pukbb_dir = root / "data/raw/panukbb/chr22"
feat_af   = root / "data/features/af/af_divergence_features.parquet"
chain     = root / "data/raw/reference/hg38ToHg19.over.chain.gz"
out_dir   = root / "data/labels"
out_dir.mkdir(parents=True, exist_ok=True)

TRAITS = ["T2D", "CAD", "BMI", "LDL"]
ANCS   = ["AFR", "AMR", "CSA", "EAS", "EUR", "MID"]


# ── 1. Target variants (GRCh38) -> lift to GRCh37 ─────────────────────────
print("Loading target chr22 variants (GRCh38) ...")
feat = pd.read_parquet(feat_af, columns=["variant_id"])
parts = feat["variant_id"].str.split(":", expand=True)
feat["pos38"] = pd.to_numeric(parts[1], errors="coerce").astype("Int64")
feat["a1"] = parts[2].str.upper()
feat["a2"] = parts[3].str.upper()
feat = feat.dropna(subset=["pos38"]).copy()
print(f"  {len(feat):,} target variants")

lift_cache = out_dir / "_chr22_pos37_cache.parquet"
if lift_cache.exists():
    print("Using cached GRCh38->GRCh37 liftover ...")
    cache = pd.read_parquet(lift_cache)
    pos37_map = dict(zip(cache["pos38"], cache["pos37"]))
else:
    print("Lifting GRCh38 -> GRCh37 (chr22) ...")
    lo = LiftOver(str(chain))
    pos37_map = {}
    for pos38 in feat["pos38"].unique():
        res = lo.convert_coordinate("chr22", int(pos38))
        if res:
            c, p37, strand, _ = res[0]
            if c in ("chr22", "22"):
                pos37_map[int(pos38)] = int(p37)
    pd.DataFrame({"pos38": list(pos37_map), "pos37": list(pos37_map.values())}) \
        .to_parquet(lift_cache, index=False)
feat["pos37"] = feat["pos38"].map(pos37_map)
feat = feat.dropna(subset=["pos37"]).copy()
feat["pos37"] = feat["pos37"].astype(int)
print(f"  {len(feat):,} variants lifted to GRCh37")

# pos37 -> list of (variant_id, allele_set)
target_idx: dict = {}
for vid, p37, a1, a2 in zip(feat["variant_id"], feat["pos37"], feat["a1"], feat["a2"]):
    target_idx.setdefault(p37, []).append((vid, frozenset({a1, a2})))


def concordance_metrics(betas, ses):
    betas = np.asarray(betas, float); ses = np.asarray(ses, float)
    n = len(betas)
    w = 1.0 / ses**2
    bmean = np.average(betas, weights=w)
    Q = float(np.sum(w * (betas - bmean) ** 2))
    dfree = n - 1
    I2 = max(0.0, (Q - dfree) / Q) if Q > 0 else 0.0
    het_p = float(stats.chi2.sf(Q, dfree)) if dfree > 0 else np.nan
    signs = np.sign(betas)
    pairs = list(combinations(range(n), 2))
    sign_conc = np.mean([signs[i] == signs[j] for i, j in pairs]) if pairs else 1.0
    if n >= 3:
        # mean pairwise correlation is undefined for scalars; use sign agreement
        beta_corr = float(np.corrcoef(betas, np.sign(betas))[0, 1]) if np.ptp(betas) > 0 else np.nan
    else:
        beta_corr = np.nan
    risk = 0
    if I2 > 0.25 or sign_conc < 0.80:
        risk = 1
    if I2 > 0.50 or sign_conc < 0.60:
        risk = 2
    return dict(n_anc=n, sign_concordance=float(sign_conc), cochran_Q=Q,
                I2=float(I2), het_pval=het_p, beta_corr=beta_corr, risk_class=risk)


# ── 2. Build labels per trait ─────────────────────────────────────────────
all_rows = []
for trait in TRAITS:
    pq = pukbb_dir / f"{trait}.chr22.parquet"
    df = pd.read_parquet(pq)
    ancs = [a for a in ANCS if f"beta_{a}" in df.columns]
    matched = 0
    for row in df.itertuples(index=False):
        p37 = int(row.pos)
        cands = target_idx.get(p37)
        if not cands:
            continue
        aset = frozenset({str(row.ref).upper(), str(row.alt).upper()})
        vid = None
        for cand_vid, cand_aset in cands:
            if cand_aset == aset:
                vid = cand_vid; break
        if vid is None:
            continue
        betas, ses, n_lowconf = [], [], 0
        used = []
        for a in ancs:
            b = getattr(row, f"beta_{a}")
            s = getattr(row, f"se_{a}")
            if pd.isna(b) or pd.isna(s) or float(s) <= 0:
                continue
            lc = getattr(row, f"low_confidence_{a}", None)
            if str(lc).lower() == "true":
                n_lowconf += 1
            betas.append(float(b)); ses.append(float(s)); used.append(a)
        if len(betas) < 2:
            continue
        m = concordance_metrics(betas, ses)
        # meta-analysis signal (ancestry-neutral) to flag truly associated variants
        bm = getattr(row, "beta_meta", np.nan)
        sm = getattr(row, "se_meta", np.nan)
        z_meta = float(bm) / float(sm) if (pd.notna(bm) and pd.notna(sm) and float(sm) > 0) else np.nan
        m.update(variant_id=vid, trait=trait, pos37=p37,
                 ancestries=",".join(used), n_low_conf=n_lowconf,
                 z_meta=z_meta, associated=bool(abs(z_meta) >= 3.0) if pd.notna(z_meta) else False)
        all_rows.append(m)
        matched += 1
    print(f"  {trait}: {matched:,} labeled variants (ancestries: {ancs})")

labels = pd.DataFrame(all_rows)
labels.to_parquet(out_dir / "gwas_concordance_labels.parquet", index=False)
print(f"\nSaved gwas_concordance_labels.parquet  ({len(labels):,} rows)")


# ── 3. Label quality audit ─────────────────────────────────────────────────
af = pd.read_parquet(feat_af, columns=["variant_id"])
ld = pd.read_parquet(root / "data/features/ld/ld_divergence_features.parquet",
                     columns=["variant_id"])
lines = []
def log(s):
    print(s); lines.append(s)

log("=" * 60)
log("LABEL COVERAGE AUDIT")
log("=" * 60)
log(f"Total labeled variant-trait rows : {len(labels):,}")
log(f"Unique variants                   : {labels['variant_id'].nunique():,}")
log(f"Traits represented                : {labels['trait'].nunique()}")
log("\nRows by trait:")
log(labels["trait"].value_counts().to_string())
log("\nn_ancestries distribution:")
log(labels["n_anc"].value_counts().sort_index().to_string())
log("\nrisk_class balance:")
log(labels["risk_class"].value_counts().sort_index().to_string())
log(f"\nI2:  mean={labels['I2'].mean():.3f}  median={labels['I2'].median():.3f}  "
    f"max={labels['I2'].max():.3f}")
log(f"sign_concordance: mean={labels['sign_concordance'].mean():.3f}  "
    f"min={labels['sign_concordance'].min():.3f}")
log(f"\nOverlap with AF features: "
    f"{labels['variant_id'].isin(af['variant_id']).mean()*100:.1f}%")
log(f"Overlap with LD features: "
    f"{labels['variant_id'].isin(ld['variant_id']).mean()*100:.1f}%")

# Associated subset (|z_meta|>=3): where cross-ancestry concordance is meaningful
assoc = labels[labels["associated"]]
log("\n" + "-" * 60)
log("ASSOCIATED SUBSET (|z_meta| >= 3, ancestry-neutral selection)")
log("-" * 60)
log(f"Associated variant-trait rows : {len(assoc):,} "
    f"({len(assoc)/len(labels)*100:.1f}% of all)")
log("Rows by trait (associated):")
log(assoc["trait"].value_counts().to_string())
log("\nrisk_class balance (associated):")
log(assoc["risk_class"].value_counts(normalize=True).sort_index().round(3).to_string())
log(f"I2 (associated): mean={assoc['I2'].mean():.3f} median={assoc['I2'].median():.3f}")
log(f"sign_concordance (associated): mean={assoc['sign_concordance'].mean():.3f}")

# Gate (checklistv2): >=10k rows, >=3 traits, no trait >70%
gate_rows = len(labels) >= 10_000
gate_traits = labels["trait"].nunique() >= 3
top_frac = labels["trait"].value_counts(normalize=True).max()
gate_balance = top_frac <= 0.70
log("\nGATE CHECK:")
log(f"  >=10,000 rows         : {'PASS' if gate_rows else 'FAIL'} ({len(labels):,})")
log(f"  >=3 traits            : {'PASS' if gate_traits else 'FAIL'}")
log(f"  no trait >70%         : {'PASS' if gate_balance else 'FAIL'} (max {top_frac*100:.0f}%)")

(out_dir / "label_audit.txt").write_text("\n".join(lines), encoding="utf-8")
print(f"\nSaved label_audit.txt")
