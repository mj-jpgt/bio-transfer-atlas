"""
GPU-heavy companion jobs for Lambda A100:
  1) Expanded XGBoost CUDA ablations (large trees, both CV splits)
  2) XGBoost TreeSHAP attribution (GPU-trained model)
  3) Torch CUDA MLP baseline on AF_LD_SEL features

Designed to keep VRAM busy while PLINK scoring runs on CPU.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.pipeline import Pipeline

ROOT = Path(__file__).resolve().parents[1]
SEED = 719


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument(
        "--sample",
        default=str(ROOT / "data/modeling/_tmp_ldblock_associated_sample.parquet"),
    )
    p.add_argument(
        "--ld-blocks",
        default=str(ROOT / "data/modeling/ld_block_assignments_genomewide.parquet"),
    )
    p.add_argument(
        "--groups",
        default=str(ROOT / "data/modeling/feature_groups_genomewide_genomewide.json"),
    )
    p.add_argument("--max-rows", type=int, default=400_000)
    p.add_argument("--device", default="cuda")
    p.add_argument(
        "--only",
        choices=["all", "xgb", "shap", "torch"],
        default="all",
        help="Run a subset of GPU lanes (default: all).",
    )
    p.add_argument("--torch-hidden", type=int, default=2048, help="MLP first hidden width.")
    p.add_argument("--torch-epochs", type=int, default=20)
    return p.parse_args()


def load_frame(sample: Path, ld_blocks: Path, max_rows: int) -> pd.DataFrame:
    df = pd.read_parquet(sample)
    if len(df) > max_rows:
        df = df.sample(max_rows, random_state=SEED)
    ld = pd.read_parquet(ld_blocks)
    keep = [c for c in ["variant_id", "split_ld_block"] if c in ld.columns]
    df = df.merge(ld[keep].drop_duplicates("variant_id"), on="variant_id", how="inner")
    rg_path = ROOT / "data/features/baselines/rg_real_by_trait.parquet"
    if rg_path.exists() and "trait" in df.columns:
        df = df.merge(pd.read_parquet(rg_path), on="trait", how="left")
    # VEP features if present
    vep = ROOT / "data/features/selection/vep_af_interaction_features.parquet"
    if vep.exists():
        v = pd.read_parquet(vep)
        keep_v = [
            c
            for c in v.columns
            if c == "variant_id"
            or c.startswith("vep_")
            or c.startswith("alphamissense")
            or c == "rare_eur"
        ]
        df = df.merge(v[keep_v].drop_duplicates("variant_id"), on="variant_id", how="left")
    return df


def make_xgb(device: str, n_estimators: int = 400):
    from xgboost import XGBClassifier

    return XGBClassifier(
        n_estimators=n_estimators,
        max_depth=8,
        learning_rate=0.05,
        subsample=0.85,
        colsample_bytree=0.85,
        min_child_weight=5,
        reg_lambda=1.0,
        random_state=SEED,
        tree_method="hist",
        device=device,
        eval_metric="auc",
        n_jobs=8,
    )


def run_xgb_ablations(df: pd.DataFrame, groups: dict, device: str) -> Path:
    feat_sets = {
        "AF_LD_SEL": groups.get("AF_LD_SEL", []),
        "FST": groups.get("FST", ["FST_like"]),
        "POP_DISTANCE": groups.get(
            "POP_DISTANCE",
            [c for c in df.columns if c.startswith("FST_") and c.endswith("_EUR")],
        ),
        "RG_PROXY": groups.get(
            "RG_PROXY",
            [c for c in ["FST_like", "AF_max_diff", "AF_var", "AIS", "LD_max_diff", "LD_entropy"] if c in df.columns],
        ),
        # Trait-constant concordance is not a variant peer; omit from GPU peer sets
        "VEP_AF": groups.get(
            "VEP_AF",
            [c for c in ["alphamissense_score", "vep_x_afmax", "vep_x_afr_eur", "rare_eur"] if c in df.columns],
        ),
    }
    if feat_sets["VEP_AF"] and feat_sets["AF_LD_SEL"]:
        feat_sets["AF_LD_SEL+VEP_AF"] = list(dict.fromkeys(feat_sets["AF_LD_SEL"] + feat_sets["VEP_AF"]))

    rows = []
    for split in ["split_ld_block", "split_variant"]:
        if split not in df.columns:
            continue
        tr = df[df[split] == "train"]
        te = df[df[split] == "test"]
        if len(tr) < 1000 or len(te) < 200:
            continue
        ytr = tr["y_high_I2"].astype(int).to_numpy()
        yte = te["y_high_I2"].astype(int).to_numpy()
        for name, feats in feat_sets.items():
            use = [f for f in feats if f in df.columns]
            if len(use) < 1 or len(np.unique(ytr)) < 2 or len(np.unique(yte)) < 2:
                continue
            pipe = Pipeline(
                [("imputer", SimpleImputer(strategy="median")), ("xgb", make_xgb(device, 400))]
            )
            used_device = device
            try:
                pipe.fit(tr[use], ytr)
                p = pipe.predict_proba(te[use])[:, 1]
            except Exception as exc:
                print(f"GPU fail {name}: {exc}; CPU fallback", flush=True)
                pipe = Pipeline(
                    [("imputer", SimpleImputer(strategy="median")), ("xgb", make_xgb("cpu", 300))]
                )
                pipe.fit(tr[use], ytr)
                p = pipe.predict_proba(te[use])[:, 1]
                used_device = "cpu_fallback"
            rows.append(
                {
                    "split": split,
                    "feature_group": name,
                    "AUROC": float(roc_auc_score(yte, p)),
                    "AUPRC": float(average_precision_score(yte, p)),
                    "n_feats": len(use),
                    "n_train": len(tr),
                    "n_test": len(te),
                    "device": used_device,
                    "model": "xgboost_hist",
                }
            )
            print(
                f"[xgb] {split} {name}: AUROC={rows[-1]['AUROC']:.4f} device={used_device}",
                flush=True,
            )

    out = ROOT / "results/tables/ablation_xgboost_gpu_expanded.csv"
    pd.DataFrame(rows).to_csv(out, index=False, float_format="%.4f")
    print(f"Saved {out}", flush=True)
    return out


def run_xgb_shap(df: pd.DataFrame, groups: dict, device: str) -> Path:
    from xgboost import XGBClassifier

    feats = [f for f in groups.get("AF_LD_SEL", []) if f in df.columns]
    split = "split_ld_block" if "split_ld_block" in df.columns else "split_variant"
    tr = df[df[split] == "train"]
    te = df[df[split] == "test"]
    if te.empty:
        te = df.sample(min(20_000, len(df)), random_state=SEED)
        tr = df.drop(te.index)

    imp = SimpleImputer(strategy="median")
    Xtr = imp.fit_transform(tr[feats]).astype(np.float32)
    Xte = imp.transform(te[feats]).astype(np.float32)
    ytr = tr["y_high_I2"].astype(int).to_numpy()

    clf = make_xgb(device, n_estimators=500)
    try:
        clf.fit(Xtr, ytr)
        used = device
    except Exception as exc:
        print(f"XGB GPU shap train failed ({exc}); CPU", flush=True)
        clf = make_xgb("cpu", n_estimators=400)
        clf.fit(Xtr, ytr)
        used = "cpu"

    # Prefer shap TreeExplainer on boosted trees (fast CPU once model is GPU-trained)
    n_explain = min(8000, len(Xte))
    rows = []
    method = "xgb_gain"
    try:
        import shap

        explainer = shap.TreeExplainer(clf)
        sv = explainer.shap_values(Xte[:n_explain])
        if isinstance(sv, list):
            sv = sv[1]
        abs_mean = np.abs(sv).mean(axis=0)
        method = "shap_tree_xgb"
        for i in range(n_explain):
            af = ld = sel = 0.0
            for j, name in enumerate(feats):
                v = abs(float(sv[i, j]))
                if name.startswith("AF_") or name in ("FST_like", "AIS"):
                    af += v
                elif name.startswith("LD_"):
                    ld += v
                elif name.startswith(("PBS", "FST_", "LOEUF", "pLI", "mis_z", "DAF")):
                    sel += v
            dom = max([("AF", af), ("LD", ld), ("SEL", sel)], key=lambda x: x[1])[0]
            rows.append(
                {
                    "variant_id": te["variant_id"].iloc[i],
                    "dominant_mechanism": dom,
                    "shap_af": af,
                    "shap_ld": ld,
                    "shap_sel": sel,
                }
            )
        print(
            f"[shap] global AF={abs_mean[[i for i,n in enumerate(feats) if n.startswith('AF_') or n in ('FST_like','AIS')]].sum():.4f}",
            flush=True,
        )
    except Exception as exc:
        print(f"TreeSHAP failed ({exc}); using XGB gain", flush=True)
        gain = clf.get_booster().get_score(importance_type="gain")
        # map f0.. to names
        imp_vec = np.zeros(len(feats))
        for k, v in gain.items():
            if k.startswith("f"):
                idx = int(k[1:])
                if idx < len(feats):
                    imp_vec[idx] = v
        for i in range(n_explain):
            af = ld = sel = 0.0
            for j, name in enumerate(feats):
                v = float(imp_vec[j])
                if name.startswith("AF_") or name in ("FST_like", "AIS"):
                    af += v
                elif name.startswith("LD_"):
                    ld += v
                elif name.startswith(("PBS", "FST_", "LOEUF", "pLI", "mis_z", "DAF")):
                    sel += v
            dom = max([("AF", af), ("LD", ld), ("SEL", sel)], key=lambda x: x[1])[0]
            rows.append(
                {
                    "variant_id": te["variant_id"].iloc[i],
                    "dominant_mechanism": dom,
                    "shap_af": af,
                    "shap_ld": ld,
                    "shap_sel": sel,
                }
            )

    out = pd.DataFrame(rows)
    out["attribution_method"] = method
    out["train_device"] = used
    path = ROOT / "results/tables/shap_xgboost_gpu_attribution.csv"
    out.to_csv(path, index=False, float_format="%.6g")
    print(out["dominant_mechanism"].value_counts())
    print(f"Saved {path}", flush=True)
    return path


def run_torch_mlp(
    df: pd.DataFrame,
    groups: dict,
    *,
    hidden: int = 2048,
    epochs: int = 20,
) -> Path:
    import torch
    from torch import nn
    from torch.utils.data import DataLoader, TensorDataset

    assert torch.cuda.is_available(), "CUDA required for torch MLP lane"
    device = torch.device("cuda")
    torch.backends.cudnn.benchmark = True
    feats = [f for f in groups.get("AF_LD_SEL", []) if f in df.columns]
    split = "split_ld_block" if "split_ld_block" in df.columns else "split_variant"
    tr = df[df[split] == "train"]
    te = df[df[split] == "test"]
    imp = SimpleImputer(strategy="median")
    Xtr = imp.fit_transform(tr[feats]).astype(np.float32)
    Xte = imp.transform(te[feats]).astype(np.float32)
    ytr = tr["y_high_I2"].astype(np.float32).to_numpy()
    yte = te["y_high_I2"].astype(np.float32).to_numpy()

    # standardize
    mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-6
    Xtr = (Xtr - mu) / sd
    Xte = (Xte - mu) / sd

    train_ds = TensorDataset(torch.from_numpy(Xtr), torch.from_numpy(ytr))
    # large batch to fill VRAM on A100-40GB
    batch = min(65536, max(8192, len(train_ds) // 8))
    loader = DataLoader(
        train_ds,
        batch_size=batch,
        shuffle=True,
        num_workers=8,
        pin_memory=True,
        persistent_workers=True,
    )

    h2, h3, h4 = hidden // 2, hidden // 4, max(64, hidden // 16)
    model = nn.Sequential(
        nn.Linear(Xtr.shape[1], hidden),
        nn.GELU(),
        nn.Dropout(0.15),
        nn.Linear(hidden, h2),
        nn.GELU(),
        nn.Dropout(0.15),
        nn.Linear(h2, h3),
        nn.GELU(),
        nn.Dropout(0.1),
        nn.Linear(h3, h4),
        nn.GELU(),
        nn.Linear(h4, 1),
    ).to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)
    loss_fn = nn.BCEWithLogitsLoss()
    scaler = torch.amp.GradScaler("cuda")

    model.train()
    for epoch in range(epochs):
        total = 0.0
        n = 0
        for xb, yb in loader:
            xb = xb.to(device, non_blocking=True)
            yb = yb.to(device, non_blocking=True).unsqueeze(1)
            opt.zero_grad(set_to_none=True)
            with torch.amp.autocast("cuda"):
                logits = model(xb)
                loss = loss_fn(logits, yb)
            scaler.scale(loss).backward()
            scaler.step(opt)
            scaler.update()
            total += float(loss.item()) * len(xb)
            n += len(xb)
        print(
            f"[torch] epoch {epoch+1}/{epochs} loss={total/max(n,1):.4f} "
            f"batch={batch} hidden={hidden}",
            flush=True,
        )

    model.eval()
    with torch.no_grad():
        preds = []
        xt = torch.from_numpy(Xte)
        for i in range(0, len(xt), batch):
            with torch.amp.autocast("cuda"):
                logits = model(xt[i : i + batch].to(device))
            preds.append(torch.sigmoid(logits.float()).cpu().numpy().ravel())
        p = np.concatenate(preds)
    auroc = float(roc_auc_score(yte.astype(int), p))
    auprc = float(average_precision_score(yte.astype(int), p))
    mem = torch.cuda.max_memory_allocated() / (1024**3)
    out = ROOT / "results/tables/ablation_torch_mlp_gpu.csv"
    pd.DataFrame(
        [
            {
                "split": split,
                "feature_group": "AF_LD_SEL",
                "AUROC": auroc,
                "AUPRC": auprc,
                "n_train": len(tr),
                "n_test": len(te),
                "device": "cuda",
                "model": f"mlp_{hidden}_{h2}_{h3}_{h4}",
                "max_vram_gb": mem,
                "batch": batch,
                "epochs": epochs,
            }
        ]
    ).to_csv(out, index=False, float_format="%.4f")
    print(f"[torch] AUROC={auroc:.4f} AUPRC={auprc:.4f} max_vram={mem:.2f}GB", flush=True)
    print(f"Saved {out}", flush=True)
    return out


def main() -> None:
    args = parse_args()
    groups = json.loads(Path(args.groups).read_text(encoding="utf-8"))
    print("Loading frame ...", flush=True)
    df = load_frame(Path(args.sample), Path(args.ld_blocks), args.max_rows)
    print(f"Rows={len(df):,} cols={len(df.columns)}", flush=True)

    if args.only in ("all", "xgb"):
        print("=== GPU XGBoost expanded ablations ===", flush=True)
        run_xgb_ablations(df, groups, args.device)

    if args.only in ("all", "shap"):
        print("=== GPU-trained XGBoost SHAP ===", flush=True)
        run_xgb_shap(df, groups, args.device)

    if args.only in ("all", "torch"):
        print("=== Torch CUDA MLP ===", flush=True)
        try:
            run_torch_mlp(df, groups, hidden=args.torch_hidden, epochs=args.torch_epochs)
        except Exception as exc:
            print(f"Torch MLP skipped: {exc}", flush=True)

    print("GPU_LANE_DONE", flush=True)


if __name__ == "__main__":
    main()
