#!/usr/bin/env python3
import pandas as pd

r = pd.read_csv("results/tables/intervention_results.genomewide.csv")
m = r[r.metric == "mean_abs_delta_EUR"][["mode", "pgs_id", "reduction"]].copy()
print(m.groupby("mode")["reduction"].agg(["mean", "median"]).round(4).sort_values("mean", ascending=False).to_string())
print("--- best mode/pgs ---")
print(m.sort_values("reduction", ascending=False).head(10).to_string(index=False))
