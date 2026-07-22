#!/usr/bin/env python3
"""Burn A100 VRAM productively while CPU atlas/intervention runs.

1) Large-depth XGBoost CUDA on AF_LD_SEL (QuantileDMatrix on device)
2) Torch MLP with random Fourier features to inflate width and VRAM use
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.metrics import average_precision_score, roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
SEED = 719


def load_df(max_rows: int = 500_000) -> tuple[pd.DataFrame, list[str]]:
    import sys

    sys.path.insert(0, str(ROOT / "scripts"))
    from run_gpu_lane import load_frame

    groups = json.loads((ROOT / "data/modeling/feature_groups_genomewide_genomewide.json").read_text())
    df = load_frame(
        ROOT / "data/modeling/_tmp_ldblock_associated_sample.parquet",
        ROOT / "data/modeling/ld_block_assignments_genomewide.parquet",
        max_rows,
    )
    feats = [f for f in groups.get("AF_LD_SEL", []) if f in df.columns]
    return df, feats


def run_xgb_heavy(df: pd.DataFrame, feats: list[str]) -> Path:
    import xgboost as xgb

    split = "split_ld_block"
    tr = df[df[split] == "train"]
    te = df[df[split] == "test"]
    imp = SimpleImputer(strategy="median")
    Xtr = imp.fit_transform(tr[feats]).astype(np.float32)
    Xte = imp.transform(te[feats]).astype(np.float32)
    ytr = tr["y_high_I2"].astype(np.float32).to_numpy()
    yte = te["y_high_I2"].astype(int).to_numpy()

    dtrain = xgb.QuantileDMatrix(Xtr, label=ytr, max_bin=256)
    dtest = xgb.QuantileDMatrix(Xte, ref=dtrain)
    params = {
        "device": "cuda",
        "tree_method": "hist",
        "max_depth": 12,
        "eta": 0.03,
        "subsample": 0.9,
        "colsample_bytree": 0.9,
        "min_child_weight": 3,
        "lambda": 1.0,
        "objective": "binary:logistic",
        "eval_metric": "auc",
        "seed": SEED,
    }
    print(f"[xgb-heavy] train n={len(tr):,} trees=1500 depth=12", flush=True)
    booster = xgb.train(params, dtrain, num_boost_round=1500, evals=[(dtest, "test")], verbose_eval=100)
    p = booster.predict(dtest)
    auroc = float(roc_auc_score(yte, p))
    auprc = float(average_precision_score(yte, p))
    out = ROOT / "results/tables/ablation_xgboost_gpu_heavy.csv"
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
                "model": "xgb_hist_depth12_1500",
            }
        ]
    ).to_csv(out, index=False, float_format="%.4f")
    print(f"[xgb-heavy] AUROC={auroc:.4f} AUPRC={auprc:.4f} -> {out}", flush=True)
    return out


def run_torch_rff(df: pd.DataFrame, feats: list[str], rff_dim: int = 8192) -> Path:
    import torch
    from torch import nn
    from torch.utils.data import DataLoader, TensorDataset

    assert torch.cuda.is_available()
    device = torch.device("cuda")
    torch.backends.cudnn.benchmark = True
    split = "split_ld_block"
    tr = df[df[split] == "train"]
    te = df[df[split] == "test"]
    imp = SimpleImputer(strategy="median")
    Xtr = imp.fit_transform(tr[feats]).astype(np.float32)
    Xte = imp.transform(te[feats]).astype(np.float32)
    ytr = tr["y_high_I2"].astype(np.float32).to_numpy()
    yte = te["y_high_I2"].astype(int).to_numpy()
    mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-6
    Xtr = (Xtr - mu) / sd
    Xte = (Xte - mu) / sd

    rng = np.random.default_rng(SEED)
    W = rng.normal(size=(Xtr.shape[1], rff_dim)).astype(np.float32) * 0.5
    b = rng.uniform(0, 2 * np.pi, size=(rff_dim,)).astype(np.float32)

    def rff(x: np.ndarray) -> np.ndarray:
        z = x @ W + b
        return np.concatenate([np.cos(z), np.sin(z)], axis=1).astype(np.float32)

    print(f"[torch-rff] expanding {Xtr.shape[1]} -> {2 * rff_dim} features", flush=True)
    Ztr = rff(Xtr)
    Zte = rff(Xte)
    train_ds = TensorDataset(torch.from_numpy(Ztr), torch.from_numpy(ytr))
    batch = min(16384, max(4096, len(train_ds) // 16))
    loader = DataLoader(train_ds, batch_size=batch, shuffle=True, num_workers=8, pin_memory=True)

    model = nn.Sequential(
        nn.Linear(Ztr.shape[1], 4096),
        nn.GELU(),
        nn.Dropout(0.1),
        nn.Linear(4096, 2048),
        nn.GELU(),
        nn.Dropout(0.1),
        nn.Linear(2048, 512),
        nn.GELU(),
        nn.Linear(512, 1),
    ).to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=8e-4, weight_decay=1e-4)
    loss_fn = nn.BCEWithLogitsLoss()
    scaler = torch.amp.GradScaler("cuda")
    epochs = 15
    model.train()
    for epoch in range(epochs):
        total = 0.0
        n = 0
        for xb, yb in loader:
            xb = xb.to(device, non_blocking=True)
            yb = yb.to(device, non_blocking=True).unsqueeze(1)
            opt.zero_grad(set_to_none=True)
            with torch.amp.autocast("cuda"):
                loss = loss_fn(model(xb), yb)
            scaler.scale(loss).backward()
            scaler.step(opt)
            scaler.update()
            total += float(loss.item()) * len(xb)
            n += len(xb)
        print(f"[torch-rff] epoch {epoch+1}/{epochs} loss={total/max(n,1):.4f} batch={batch}", flush=True)

    model.eval()
    with torch.no_grad():
        preds = []
        xt = torch.from_numpy(Zte)
        for i in range(0, len(xt), batch):
            with torch.amp.autocast("cuda"):
                logits = model(xt[i : i + batch].to(device))
            preds.append(torch.sigmoid(logits.float()).cpu().numpy().ravel())
        p = np.concatenate(preds)
    auroc = float(roc_auc_score(yte, p))
    auprc = float(average_precision_score(yte, p))
    mem = torch.cuda.max_memory_allocated() / (1024**3)
    out = ROOT / "results/tables/ablation_torch_rff_gpu.csv"
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
                "model": f"rff_{rff_dim}_mlp_4096",
                "max_vram_gb": mem,
                "batch": batch,
            }
        ]
    ).to_csv(out, index=False, float_format="%.4f")
    print(f"[torch-rff] AUROC={auroc:.4f} AUPRC={auprc:.4f} max_vram={mem:.2f}GB -> {out}", flush=True)
    return out


def main() -> None:
    print("Loading ...", flush=True)
    df, feats = load_df(500_000)
    print(f"rows={len(df):,} feats={len(feats)}", flush=True)
    run_xgb_heavy(df, feats)
    run_torch_rff(df, feats, rff_dim=8192)
    print("GPU_HEAVY_DONE", flush=True)


if __name__ == "__main__":
    main()
