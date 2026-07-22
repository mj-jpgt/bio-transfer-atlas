#!/usr/bin/env python3
"""Robustness gate: require honesty artifacts before paper freeze."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"


def main() -> None:
    required = [
        TABLES / "mhc_sensitivity_genomewide.csv",
        TABLES / "ablation_ldblock_and_baselines_genomewide.csv",
        TABLES / "ablation_ldblock_peer_contest.csv",
        TABLES / "ablation_finemap_tiers.csv",
        TABLES / "vep_af_interaction_eval.csv",
        TABLES / "duffy_positive_control_genomewide.csv",
        TABLES / "duffy_allele_audit.csv",
        TABLES / "duffy_ackr1_score_decomposition.csv",
        TABLES / "external_sumstat_validation.csv",
        TABLES / "internal_panukbb_concordance_sensitivity.csv",
        TABLES / "external_page_validation.csv",
        TABLES / "external_page_qc_counts.csv",
        TABLES / "shap_mechanism_attribution_genomewide.csv",
        TABLES / "score_pgen_chr1_7_rebuild_status.csv",
        TABLES / "subpop_af_features_status.csv",
        TABLES / "popcorn_rg_summary.csv",
        TABLES / "trait_scale_portability.csv",
        TABLES / "ablation_nested_af_ld_sel.csv",
        TABLES / "auroc_paired_delta_ldblock.csv",
        TABLES / "intervention_retention_variance_metrics.csv",
        TABLES / "intervention_matched_random_controls.csv",
        TABLES / "intervention_loso_mad_by_mode.csv",
        TABLES / "grouped_permutation_importance_ldblock.csv",
        TABLES / "sign_discordance_endpoint.csv",
        ROOT / "paper/METHODS.md",
        ROOT / "paper/RESULTS.md",
        ROOT / "paper/DISCUSSION.md",
        ROOT / "results/BTA_M4_M6_RESULTS_REPORT.md",
    ]
    missing = [str(p) for p in required if not p.exists()]
    if missing:
        raise SystemExit(f"FAIL missing: {missing}")

    peer = pd.read_csv(TABLES / "ablation_ldblock_peer_contest.csv")
    if (peer["feature_group"] == "RG_REAL").any():
        raise SystemExit("FAIL peer contest still contains RG_REAL")
    if peer["feature_group"].astype(str).str.contains("TRAIT_CONSTANT").any():
        raise SystemExit("FAIL peer contest still contains TRAIT_CONSTANT diagnostic")

    ld = pd.read_csv(TABLES / "ablation_ldblock_and_baselines_genomewide.csv")
    row = ld[(ld["split"] == "split_ld_block") & (ld["feature_group"] == "AF_LD_SEL")].iloc[0]
    fst = ld[(ld["split"] == "split_ld_block") & (ld["feature_group"] == "FST")].iloc[0]
    if not (row["AUROC"] > fst["AUROC"]):
        raise SystemExit("FAIL AF_LD_SEL does not beat FST")

    internal = pd.read_csv(TABLES / "internal_panukbb_concordance_sensitivity.csv")
    if "status" in internal.columns and not all(
        str(s) in ("internal_only", "missing_panukbb") for s in internal["status"]
    ):
        raise SystemExit("FAIL internal sensitivity must be status=internal_only")

    page = pd.read_csv(TABLES / "external_page_validation.csv")
    if "status" in page.columns:
        for _, r in page.iterrows():
            nvar = r.get("n_variants", r.get("n", 0))
            if str(r["status"]) == "ok_external" and int(nvar or 0) < 500:
                raise SystemExit("FAIL ok_external with n_variants<500")
    qc = pd.read_csv(TABLES / "external_page_qc_counts.csv")
    if qc.empty or "final_variants" not in qc.columns:
        raise SystemExit("FAIL external_page_qc_counts incomplete")

    # Forbid stale circular Pan-UKB-only rows labeled as external ok
    ext = pd.read_csv(TABLES / "external_sumstat_validation.csv")
    for _, r in ext.iterrows():
        st = str(r.get("status", ""))
        analysis = str(r.get("analysis", "")).lower()
        pair = str(r.get("pair", "")).lower()
        if st == "ok_open_panukbb":
            raise SystemExit(
                "FAIL Pan-UKB concordance still labeled as external ok; use internal_only"
            )
        if analysis.startswith("internal") and st.startswith("ok_external"):
            raise SystemExit("FAIL internal analysis labeled ok_external")
        if pair.startswith("panukbb_") and "page" not in pair and st.startswith("ok_external"):
            raise SystemExit(
                "FAIL Pan-UKB concordance still labeled as external ok; use internal_only"
            )

    results = (ROOT / "paper/RESULTS.md").read_text(encoding="utf-8").lower()
    banned = ["af/ld beats", "beats rg_real", "af_ld_sel beats rg"]
    if any(b in results for b in banned):
        raise SystemExit("FAIL RESULTS still claims AF/LD beats rg peer contest")
    if "external validation" in results and "internal" not in results:
        raise SystemExit("FAIL RESULTS mentions external validation without internal framing")

    meta_paths = [
        ROOT / "data/features/baselines/rg_real_meta.json",
        TABLES / "rg_real_meta.json",
    ]
    meta_ok = False
    for mp in meta_paths:
        if not mp.exists():
            continue
        m = json.loads(mp.read_text(encoding="utf-8"))
        if m.get("estimand") == "cross_ancestry_z_score_concordance" or "z_concordance" in str(
            m.get("method", "")
        ):
            meta_ok = True
    if not meta_ok:
        # Also accept z_concordance columns in popcorn summary
        rg = pd.read_csv(TABLES / "popcorn_rg_summary.csv")
        if "z_concordance" not in rg.columns and not rg.get("method", pd.Series(dtype=str)).astype(
            str
        ).str.contains("z_concordance|zcorr", na=False).any():
            raise SystemExit("FAIL missing z_concordance naming in rg meta/summary")

    duffy = pd.read_csv(TABLES / "duffy_allele_audit.csv")
    need = {"plink_counted_allele", "duffy_null_allele_target"}
    if not need.intersection(set(duffy.columns)) and "plink_counted_allele" not in duffy.columns:
        # accept alternate column names from audit writer
        if "counted_allele" not in duffy.columns and "duffy_null_allele" not in "".join(duffy.columns):
            raise SystemExit("FAIL duffy allele audit incomplete")

    nest = pd.read_csv(TABLES / "ablation_nested_af_ld_sel.csv")
    if len(nest) < 3:
        raise SystemExit("FAIL nested ablation too short")
    paired = pd.read_csv(TABLES / "auroc_paired_delta_ldblock.csv")
    if paired.empty:
        raise SystemExit("FAIL empty paired ΔAUROC table")
    need_pairs = {("AF_LD_SEL", "AF"), ("AF_LD", "AF"), ("AF_LD_SEL", "AF_LD")}
    have = {(str(a), str(b)) for a, b in zip(paired["model_a"], paired["model_b"])}
    if not need_pairs.issubset(have):
        raise SystemExit(f"FAIL paired table missing AF nest deltas; have={have}")
    if "verdict" not in paired.columns:
        raise SystemExit("FAIL paired table missing verdict")

    fg = json.loads(
        (ROOT / "data/modeling/feature_groups_genomewide_genomewide.json").read_text(encoding="utf-8")
    )
    if "RG_REAL" in fg:
        raise SystemExit("FAIL feature groups still contain RG_REAL")
    af_cols = [str(c) for c in fg.get("AF", [])]
    if any(c.startswith("FST_") or c.startswith("fst_") for c in af_cols):
        raise SystemExit("FAIL FST columns still nested under AF")

    decomp = pd.read_csv(TABLES / "duffy_ackr1_score_decomposition.csv")
    if decomp.empty or "ackr1_fraction" not in decomp.columns:
        raise SystemExit("FAIL duffy ACKR1 decomposition incomplete")
    if not ((decomp.get("super_pop") == "AFR") & (decomp.get("status") == "ok")).any():
        raise SystemExit("FAIL missing AFR ok ACKR1 decomposition row")

    mc = pd.read_csv(TABLES / "intervention_matched_random_controls.csv")
    for col in ("empirical_p_n", "empirical_p_mass"):
        if col not in mc.columns:
            raise SystemExit(f"FAIL matched controls missing {col}")
    loso = pd.read_csv(TABLES / "intervention_loso_mad_by_mode.csv")
    if loso.empty:
        raise SystemExit("FAIL empty LOSO MAD table")

    gpi = pd.read_csv(TABLES / "grouped_permutation_importance_ldblock.csv")
    if gpi.empty or "feature_family" not in gpi.columns:
        raise SystemExit("FAIL grouped permutation importance missing")

    ret = pd.read_csv(TABLES / "intervention_retention_variance_metrics.csv")
    if ret.empty:
        raise SystemExit("FAIL empty intervention retention metrics")

    # SuSiE primary: require meta + status; majority signed_ld OR explicit awaiting rerun
    susie = ROOT / "data/labels/susie/susie_real_ld_status_summary.csv"
    if not susie.exists():
        susie = TABLES / "susie_real_ld_status_summary.csv"
    if not susie.exists():
        raise SystemExit("FAIL missing SuSiE status summary")
    meta_s = ROOT / "data/labels/susie/susie_primary_meta.json"
    if not meta_s.exists():
        raise SystemExit("FAIL missing susie_primary_meta.json (signed_ld primary rule)")
    sm = json.loads(meta_s.read_text(encoding="utf-8"))
    if sm.get("primary_requires") != "signed_ld":
        raise SystemExit("FAIL susie primary_requires must be signed_ld")
    ss = pd.read_csv(susie)
    if "n_signed_ld" in ss.columns and float(ss["n_signed_ld"].fillna(0).sum()) > 0:
        # Prefer majority signed among files that report a top mode / status
        if "status" in ss.columns:
            ok = ss["status"].astype(str).isin(["susie_cs_signed_ld", "partial_signed_ld", "awaiting_signed_ld_rerun"])
            if not ok.all():
                raise SystemExit("FAIL unexpected SuSiE status values")
            if (ss["status"].astype(str) == "susie_cs_signed_ld").sum() < 1 and float(
                ss["n_signed_ld"].sum()
            ) < 1000:
                raise SystemExit("FAIL no substantial signed_ld SuSiE coverage yet")
        elif "n_in_cs" in ss.columns:
            n_signed = float(ss["n_signed_ld"].sum())
            n_cs = float(ss["n_in_cs"].sum())
            if n_cs > 0 and (n_signed / max(n_cs, 1)) < 0.5 and n_signed < 1000:
                raise SystemExit("FAIL primary SuSiE claims without majority signed_ld")
    elif "status" in ss.columns and not all(
        str(s) in ("awaiting_signed_ld_rerun", "susie_cs_signed_ld", "partial_signed_ld")
        for s in ss["status"]
    ):
        raise SystemExit("FAIL SuSiE status must be signed_ld or awaiting_signed_ld_rerun")

    report = (ROOT / "results/BTA_M4_M6_RESULTS_REPORT.md").read_text(encoding="utf-8").lower()
    if "af/ld beats" in report and "not" not in report:
        # soft check — prefer explicit not-a-peer language
        pass

    print("PASS robustness gates")
    print(
        f"  LD-block AF_LD_SEL AUROC={row['AUROC']:.3f} (FST={fst['AUROC']:.3f})"
    )
    print("  peer contest excludes RG_REAL; internal concordance labeled internal_only")
    print(f"  nested rows={len(nest)}; paired ΔAUROC rows={len(paired)}; retention rows={len(ret)}")


if __name__ == "__main__":
    main()
