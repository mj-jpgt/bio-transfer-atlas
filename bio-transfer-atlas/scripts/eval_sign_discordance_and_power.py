#!/usr/bin/env python3
"""Sign-discordance endpoint + power-confound sensitivity (SE/N covariates)."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
SEED = 719


def auroc_safe(y, p):
    y = np.asarray(y).astype(int)
    if len(np.unique(y)) < 2:
        return float("nan")
    return float(roc_auc_score(y, p))


def main() -> None:
    sample = ROOT / "data/modeling/_tmp_ldblock_associated_sample.parquet"
    labels = ROOT / "data/labels/_tmp_associated_labels.parquet"
    groups = json.loads(
        (ROOT / "data/modeling/feature_groups_genomewide_genomewide.json").read_text(encoding="utf-8")
    )
    out = ROOT / "results/tables/sign_discordance_endpoint.csv"
    sens = ROOT / "results/tables/power_confound_sensitivity.csv"
    cont = ROOT / "results/tables/continuous_i2_endpoint.csv"
    strata_out = ROOT / "results/tables/power_quartile_strata_auroc.csv"
    if not sample.exists():
        pd.DataFrame([{"status": "missing_sample"}]).to_csv(out, index=False)
        return
    df = pd.read_parquet(sample)
    if labels.exists():
        lab = pd.read_parquet(labels)
        extra = [c for c in lab.columns if c != "variant_id" and c not in df.columns]
        if extra:
            df = df.merge(lab[["variant_id"] + extra].drop_duplicates("variant_id"), on="variant_id", how="left")
    blocks = ROOT / "data/modeling/ld_block_assignments_genomewide.parquet"
    if "split_ld_block" not in df.columns and blocks.exists():
        ld = pd.read_parquet(blocks, columns=["variant_id", "split_ld_block"])
        df = df.merge(ld.drop_duplicates("variant_id"), on="variant_id", how="inner")
    if "y_sign_discordant" not in df.columns and "sign_concordance" in df.columns:
        df["y_sign_discordant"] = (pd.to_numeric(df["sign_concordance"], errors="coerce") < 1.0).astype(int)
    if "y_sign_discordant" not in df.columns:
        pd.DataFrame(
            [{"status": "no_sign_label", "note": "master lacks y_sign_discordant"}]
        ).to_csv(out, index=False)
        pd.DataFrame([{"status": "skipped"}]).to_csv(sens, index=False)
        print("no sign label")
        return

    af = [f for f in groups.get("AF_LD_SEL", []) if f in df.columns]
    split = "split_ld_block" if "split_ld_block" in df.columns else None
    if split:
        tr = df[df[split] == "train"]
        te = df[df[split] == "test"]
    else:
        rng = np.random.default_rng(SEED)
        mask = rng.random(len(df)) < 0.2
        te, tr = df[mask], df[~mask]

    rows = []
    for endpoint, ycol in [("high_I2", "y_high_I2"), ("sign_discordant", "y_sign_discordant")]:
        if ycol not in df.columns:
            continue
        ytr = tr[ycol].astype(int)
        yte = te[ycol].astype(int)
        if yte.nunique() < 2:
            continue
        imp = SimpleImputer(strategy="median")
        Xtr = imp.fit_transform(tr[af])
        Xte = imp.transform(te[af])
        clf = HistGradientBoostingClassifier(
            max_depth=6, learning_rate=0.05, max_iter=200, random_state=SEED
        )
        clf.fit(Xtr, ytr)
        p = clf.predict_proba(Xte)[:, 1]
        rows.append(
            {
                "endpoint": endpoint,
                "feature_group": "AF_LD_SEL",
                "AUROC": auroc_safe(yte, p),
                "pos_rate": float(yte.mean()),
                "n_test": len(yte),
                "status": "ok",
            }
        )
    pd.DataFrame(rows).to_csv(out, index=False, float_format="%.4f")

    # Continuous I2 / Δβ tables when joinable (do not block freeze)
    cont_rows = []
    y_cont_cols = [c for c in ("I2", "i2", "abs_delta_beta", "std_delta_beta") if c in df.columns]
    if not y_cont_cols:
        cont_rows.append({"status": "skipped", "reason": "no_continuous_I2_or_delta_beta"})
    else:
        from scipy.stats import spearmanr

        for yc in y_cont_cols:
            y = pd.to_numeric(te[yc], errors="coerce")
            for feat in af[:8]:
                x = pd.to_numeric(te[feat], errors="coerce")
                m = x.notna() & y.notna()
                if m.sum() < 200:
                    cont_rows.append(
                        {"status": "skipped", "outcome": yc, "feature": feat, "reason": "n<200"}
                    )
                    continue
                rho, _ = spearmanr(x[m], y[m])
                mae = float(np.mean(np.abs(x[m] - y[m])))
                cont_rows.append(
                    {
                        "status": "ok",
                        "outcome": yc,
                        "feature": feat,
                        "spearman_r": float(rho),
                        "mae": mae,
                        "n": int(m.sum()),
                    }
                )
    pd.DataFrame(cont_rows).to_csv(cont, index=False, float_format="%.4f")

    # Power confound: add SE/N covariates or filter min MAC
    sens_rows = []
    ytr = tr["y_high_I2"].astype(int)
    yte = te["y_high_I2"].astype(int)
    base_feats = af
    power_cols = [c for c in ("se_meta", "n_meta", "N_meta", "mac", "MAC") if c in df.columns]
    for label, feats in [
        ("AF_LD_SEL", base_feats),
        ("AF_LD_SEL_plus_power", base_feats + power_cols),
    ]:
        feats = [f for f in feats if f in df.columns]
        if not feats:
            continue
        imp = SimpleImputer(strategy="median")
        clf = HistGradientBoostingClassifier(
            max_depth=6, learning_rate=0.05, max_iter=200, random_state=SEED
        )
        clf.fit(imp.fit_transform(tr[feats]), ytr)
        p = clf.predict_proba(imp.transform(te[feats]))[:, 1]
        sens_rows.append(
            {"setting": label, "AUROC": auroc_safe(yte, p), "n_feats": len(feats), "status": "ok"}
        )
    if not power_cols:
        sens_rows.append(
            {"setting": "AF_LD_SEL_plus_power", "status": "skipped", "reason": "no_se_n_mac"}
        )
    if "mac" in df.columns or "MAC" in df.columns:
        mac_c = "mac" if "mac" in df.columns else "MAC"
        thr = float(np.nanpercentile(df[mac_c], 25))
        tr2 = tr[tr[mac_c] >= thr]
        te2 = te[te[mac_c] >= thr]
        if len(te2) > 200 and te2["y_high_I2"].nunique() > 1:
            imp = SimpleImputer(strategy="median")
            clf = HistGradientBoostingClassifier(
                max_depth=6, learning_rate=0.05, max_iter=200, random_state=SEED
            )
            clf.fit(imp.fit_transform(tr2[base_feats]), tr2["y_high_I2"].astype(int))
            p = clf.predict_proba(imp.transform(te2[base_feats]))[:, 1]
            sens_rows.append(
                {
                    "setting": f"AF_LD_SEL_min_mac_p25_{thr:.1f}",
                    "AUROC": auroc_safe(te2["y_high_I2"].astype(int), p),
                    "n_test": len(te2),
                    "status": "ok",
                }
            )
    pd.DataFrame(sens_rows).to_csv(sens, index=False, float_format="%.4f")

    # Quartile strata by se_meta or MAC
    strata_rows = []
    strat_col = next((c for c in ("se_meta", "mac", "MAC") if c in df.columns), None)
    if strat_col is None:
        strata_rows.append({"status": "skipped", "reason": "no_se_meta_or_mac"})
    else:
        te_s = te.copy()
        te_s["_q"] = pd.qcut(
            pd.to_numeric(te_s[strat_col], errors="coerce"), 4, labels=False, duplicates="drop"
        )
        for q, g in te_s.groupby("_q"):
            if len(g) < 100 or g["y_high_I2"].nunique() < 2:
                strata_rows.append(
                    {
                        "status": "skipped",
                        "stratum_col": strat_col,
                        "quartile": int(q),
                        "n": len(g),
                        "reason": "underpowered",
                    }
                )
                continue
            imp = SimpleImputer(strategy="median")
            clf = HistGradientBoostingClassifier(
                max_depth=6, learning_rate=0.05, max_iter=200, random_state=SEED
            )
            clf.fit(imp.fit_transform(tr[base_feats]), ytr)
            p = clf.predict_proba(imp.transform(g[base_feats]))[:, 1]
            strata_rows.append(
                {
                    "status": "ok",
                    "stratum_col": strat_col,
                    "quartile": int(q),
                    "AUROC": auroc_safe(g["y_high_I2"].astype(int), p),
                    "n_test": len(g),
                    "feature_group": "AF_LD_SEL",
                }
            )
    pd.DataFrame(strata_rows).to_csv(strata_out, index=False, float_format="%.4f")
    print(f"Saved {out}, {sens}, {cont}, {strata_out}")


if __name__ == "__main__":
    main()
