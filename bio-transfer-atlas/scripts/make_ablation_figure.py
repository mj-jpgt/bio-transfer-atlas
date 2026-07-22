"""
FAIRGEN-Open Stage 11 (cont.): Mechanism-ablation figure
=========================================================
Renders the headline result: AF/LD/selection vs FST-only for predicting
cross-ancestry concordance, plus the variant- vs trait-holdout gap.

Outputs:
  results/figures/fig_ablation_auroc.png
  results/figures/fig_ablation_generalization.png
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root   = Path(__file__).resolve().parents[1]
tables = root / "results/tables"
figs   = root / "results/figures"
figs.mkdir(parents=True, exist_ok=True)

ORDER = ["FST", "AF", "LD", "SEL", "AF_LD", "AF_SEL", "LD_SEL", "AF_LD_SEL"]
COL = {"FST": "#999999", "AF": "#1f77b4", "LD": "#ff7f0e", "SEL": "#2ca02c",
       "AF_LD": "#9467bd", "AF_SEL": "#8c564b", "LD_SEL": "#e377c2",
       "AF_LD_SEL": "#d62728"}

cls = pd.read_csv(tables / "ablation_classification.csv")
reg = pd.read_csv(tables / "ablation_regression.csv")

# ── Figure 1: AUROC + I2 Spearman, associated subset, variant split ─────────
c = cls[(cls.subset == "associated") & (cls.split == "split_variant") & (cls.model == "hgb")] \
    .set_index("feature_group").reindex(ORDER)
r = reg[(reg.subset == "associated") & (reg.split == "split_variant") & (reg.model == "hgb")] \
    .set_index("feature_group").reindex(ORDER)

fig, ax = plt.subplots(1, 2, figsize=(13, 5))
x = np.arange(len(ORDER)); colors = [COL[g] for g in ORDER]
ax[0].bar(x, c["AUROC"], color=colors)
ax[0].axhline(0.5, ls="--", c="k", lw=0.8, label="chance")
ax[0].axhline(c.loc["FST", "AUROC"], ls=":", c="#999", lw=1.2, label="FST-only")
ax[0].set_xticks(x); ax[0].set_xticklabels(ORDER, rotation=45, ha="right")
ax[0].set_ylabel("AUROC"); ax[0].set_ylim(0.45, 0.95)
ax[0].set_title("Predicting cross-ancestry heterogeneity (I²>0.25)")
for i, v in enumerate(c["AUROC"]):
    ax[0].text(i, v + 0.005, f"{v:.2f}", ha="center", fontsize=8)
ax[0].legend()

ax[1].bar(x, r["spearman_rho"], color=colors)
ax[1].axhline(r.loc["FST", "spearman_rho"], ls=":", c="#999", lw=1.2, label="FST-only")
ax[1].set_xticks(x); ax[1].set_xticklabels(ORDER, rotation=45, ha="right")
ax[1].set_ylabel("Spearman ρ"); ax[1].set_title("I² regression (rank accuracy)")
for i, v in enumerate(r["spearman_rho"]):
    ax[1].text(i, v + 0.005, f"{v:.2f}", ha="center", fontsize=8)
ax[1].legend()
fig.suptitle("Mechanism ablation (chr22, associated variants, unseen-variant split)",
             fontsize=13)
plt.tight_layout()
fig.savefig(figs / "fig_ablation_auroc.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("Saved fig_ablation_auroc.png")

# ── Figure 2: generalization gap (variant vs trait holdout), AF_LD_SEL ──────
fig, ax = plt.subplots(figsize=(7, 5))
width = 0.35
splits = ["split_variant", "split_trait"]
labels = ["unseen variants", "unseen trait (LDL)"]
for i, grp in enumerate(["FST", "AF_LD", "AF_LD_SEL"]):
    vals = []
    for sp in splits:
        row = cls[(cls.subset == "associated") & (cls.split == sp) &
                  (cls.model == "hgb") & (cls.feature_group == grp)]
        vals.append(row["AUROC"].iloc[0] if len(row) else np.nan)
    ax.bar(np.arange(2) + i * width, vals, width, label=grp,
           color=COL.get(grp, "#555"))
ax.axhline(0.5, ls="--", c="k", lw=0.8)
ax.set_xticks(np.arange(2) + width)
ax.set_xticklabels(labels)
ax.set_ylabel("AUROC (I²>0.25)")
ax.set_title("Generalization: unseen variants vs unseen trait (chr22)")
ax.legend(title="features")
plt.tight_layout()
fig.savefig(figs / "fig_ablation_generalization.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("Saved fig_ablation_generalization.png")
