#!/usr/bin/env python3
"""
M5 LD-graph GAT: per LD-block graphs, AF/LD node features, LD-block CV vs trees.

Edges: r2 proxy via same-block adjacency within 250kb (feature-space neighbor graph
when genotype LD matrix is unavailable). Uses torch_geometric GATConv on CUDA.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.metrics import average_precision_score, roc_auc_score

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
    p.add_argument("--max-blocks", type=int, default=400)
    p.add_argument("--epochs", type=int, default=25)
    p.add_argument("--hidden", type=int, default=64)
    p.add_argument("--device", default="cuda")
    return p.parse_args()


def build_block_graphs(df: pd.DataFrame, feats: list[str], max_nodes: int = 80):
    import torch
    from torch_geometric.data import Data

    graphs = []
    for bid, g in df.groupby("ld_block", sort=False):
        if len(g) < 8:
            continue
        if len(g) > max_nodes:
            g = g.sample(max_nodes, random_state=SEED)
        g = g.reset_index(drop=True)
        X = g[feats].to_numpy(dtype=np.float32)
        y = g["y_high_I2"].astype(np.float32).to_numpy()
        # Position-based edges within 250kb if pos available, else kNN in feature space
        edges = []
        if "pos" in g.columns or "position" in g.columns:
            pos = pd.to_numeric(g.get("pos", g.get("position")), errors="coerce").to_numpy()
            for i in range(len(g)):
                for j in range(i + 1, len(g)):
                    if np.isfinite(pos[i]) and np.isfinite(pos[j]) and abs(pos[i] - pos[j]) <= 250_000:
                        edges.append((i, j))
                        edges.append((j, i))
        if len(edges) < len(g):
            # fallback: connect each node to 4 nearest in feature space
            from sklearn.neighbors import NearestNeighbors

            nn = NearestNeighbors(n_neighbors=min(5, len(g))).fit(X)
            idx = nn.kneighbors(return_distance=False)
            for i, nbrs in enumerate(idx):
                for j in nbrs:
                    if i != j:
                        edges.append((i, j))
        if not edges:
            continue
        edge_index = torch.tensor(np.array(edges).T, dtype=torch.long)
        data = Data(
            x=torch.from_numpy(X),
            edge_index=edge_index,
            y=torch.from_numpy(y),
            split=str(g["split_ld_block"].iloc[0]) if "split_ld_block" in g.columns else "train",
            block_id=str(bid),
        )
        graphs.append(data)
    return graphs


def main() -> None:
    args = parse_args()
    import torch
    import torch.nn.functional as F
    from torch_geometric.nn import GATConv, global_mean_pool

    assert torch.cuda.is_available() or args.device == "cpu"
    device = torch.device(args.device if torch.cuda.is_available() else "cpu")
    torch.backends.cudnn.benchmark = True

    groups = json.loads(Path(args.groups).read_text(encoding="utf-8"))
    feats = [f for f in groups.get("AF_LD_SEL", []) if True]
    df = pd.read_parquet(args.sample)
    ld = pd.read_parquet(args.ld_blocks)
    keep = [c for c in ["variant_id", "split_ld_block", "ld_block"] if c in ld.columns]
    if "ld_block" not in keep and "block_id" in ld.columns:
        ld = ld.rename(columns={"block_id": "ld_block"})
        keep = [c for c in ["variant_id", "split_ld_block", "ld_block"] if c in ld.columns]
    df = df.merge(ld[keep].drop_duplicates("variant_id"), on="variant_id", how="inner")
    if "ld_block" not in df.columns:
        # synthesize blocks from chrom+mb
        parts = df["variant_id"].astype(str).str.split(":", n=2, expand=True)
        df["ld_block"] = parts[0] + ":" + (pd.to_numeric(parts[1], errors="coerce") // 1_000_000).astype(str)

    feats = [f for f in feats if f in df.columns]
    imp = SimpleImputer(strategy="median")
    df[feats] = imp.fit_transform(df[feats])

    # subsample blocks for tractability
    blocks = df["ld_block"].drop_duplicates().tolist()
    rng = np.random.default_rng(SEED)
    if len(blocks) > args.max_blocks:
        blocks = list(rng.choice(blocks, size=args.max_blocks, replace=False))
    df = df[df["ld_block"].isin(blocks)].copy()
    print(f"rows={len(df):,} feats={len(feats)} blocks={df['ld_block'].nunique()}", flush=True)

    graphs = build_block_graphs(df, feats)
    print(f"graphs={len(graphs)}", flush=True)
    train_g = [g for g in graphs if g.split == "train"]
    test_g = [g for g in graphs if g.split == "test"]
    if len(test_g) < 5:
        # random block holdout
        rng.shuffle(graphs)
        n_te = max(1, len(graphs) // 5)
        test_g = graphs[:n_te]
        train_g = graphs[n_te:]

    class GATNet(torch.nn.Module):
        def __init__(self, in_dim: int, hidden: int):
            super().__init__()
            self.conv1 = GATConv(in_dim, hidden, heads=4, concat=True, dropout=0.1)
            self.conv2 = GATConv(hidden * 4, hidden, heads=1, concat=True, dropout=0.1)
            self.lin = torch.nn.Linear(hidden, 1)

        def forward(self, data):
            x, edge_index = data.x, data.edge_index
            x = F.elu(self.conv1(x, edge_index))
            x = F.elu(self.conv2(x, edge_index))
            return self.lin(x).squeeze(-1)

    model = GATNet(len(feats), args.hidden).to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)
    scaler = torch.amp.GradScaler("cuda", enabled=device.type == "cuda")

    def run_epoch(graphs_list, train: bool):
        model.train(train)
        total = 0.0
        n = 0
        ys, ps = [], []
        for data in graphs_list:
            data = data.to(device)
            with torch.set_grad_enabled(train):
                with torch.amp.autocast("cuda", enabled=device.type == "cuda"):
                    logits = model(data)
                    loss = F.binary_cross_entropy_with_logits(logits, data.y)
                if train:
                    opt.zero_grad(set_to_none=True)
                    scaler.scale(loss).backward()
                    scaler.step(opt)
                    scaler.update()
            total += float(loss.item()) * len(data.y)
            n += len(data.y)
            ys.append(data.y.detach().cpu().numpy())
            ps.append(torch.sigmoid(logits).detach().cpu().numpy())
        y = np.concatenate(ys)
        p = np.concatenate(ps)
        auroc = float(roc_auc_score(y.astype(int), p)) if len(np.unique(y)) > 1 else float("nan")
        auprc = float(average_precision_score(y.astype(int), p)) if len(np.unique(y)) > 1 else float("nan")
        return total / max(n, 1), auroc, auprc

    best = 0.0
    for epoch in range(args.epochs):
        tr_loss, tr_auc, _ = run_epoch(train_g, True)
        te_loss, te_auc, te_auprc = run_epoch(test_g, False)
        best = max(best, te_auc if np.isfinite(te_auc) else 0.0)
        print(
            f"[gat] epoch {epoch+1}/{args.epochs} train_loss={tr_loss:.4f} "
            f"train_auc={tr_auc:.4f} test_auc={te_auc:.4f} test_auprc={te_auprc:.4f}",
            flush=True,
        )

    # baseline XGB on same rows
    from xgboost import XGBClassifier

    tr = df[df.get("split_ld_block", "train") == "train"] if "split_ld_block" in df.columns else df.sample(frac=0.8, random_state=SEED)
    te = df[df.get("split_ld_block", "test") == "test"] if "split_ld_block" in df.columns else df.drop(tr.index)
    if te.empty:
        te = df.sample(frac=0.2, random_state=SEED + 1)
        tr = df.drop(te.index)
    clf = XGBClassifier(
        n_estimators=300,
        max_depth=6,
        learning_rate=0.05,
        tree_method="hist",
        device="cuda" if device.type == "cuda" else "cpu",
        eval_metric="auc",
        random_state=SEED,
    )
    clf.fit(tr[feats], tr["y_high_I2"].astype(int))
    p_xgb = clf.predict_proba(te[feats])[:, 1]
    xgb_auc = float(roc_auc_score(te["y_high_I2"].astype(int), p_xgb))
    xgb_auprc = float(average_precision_score(te["y_high_I2"].astype(int), p_xgb))

    _, gat_auc, gat_auprc = run_epoch(test_g, False)
    out = ROOT / "results/tables/ablation_gat_ldblock.csv"
    pd.DataFrame(
        [
            {
                "model": "gat",
                "split": "split_ld_block",
                "AUROC": gat_auc,
                "AUPRC": gat_auprc,
                "device": str(device),
                "n_graphs_train": len(train_g),
                "n_graphs_test": len(test_g),
                "n_feats": len(feats),
            },
            {
                "model": "xgboost_hist",
                "split": "split_ld_block",
                "AUROC": xgb_auc,
                "AUPRC": xgb_auprc,
                "device": "cuda" if device.type == "cuda" else "cpu",
                "n_graphs_train": len(train_g),
                "n_graphs_test": len(test_g),
                "n_feats": len(feats),
            },
        ]
    ).to_csv(out, index=False, float_format="%.4f")
    ckpt = ROOT / "data/modeling/gat_ldblock.pt"
    ckpt.parent.mkdir(parents=True, exist_ok=True)
    torch.save({"model": model.state_dict(), "feats": feats, "hidden": args.hidden}, ckpt)
    print(f"Saved {out} and {ckpt}", flush=True)
    print("GAT_DONE", flush=True)


if __name__ == "__main__":
    main()
