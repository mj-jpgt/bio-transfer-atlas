import pandas as pd


def main():
    ci = pd.read_csv("results/tables/headline_metrics_ci.csv")
    neg = ci[ci["feature_group"] == "PERMUTED"]
    assert len(neg) > 0, "missing permuted controls"
    assert (neg["AUROC"] >= 0.45).all() and (neg["AUROC"] <= 0.57).all(), "permuted AUROC not near chance"
    real = ci[(ci["feature_group"] == "AF_LD") & (ci["split"] == "split_variant")]
    assert len(real) > 0
    assert real["AUROC_lo"].max() > 0.5, "AF_LD should beat chance"
    print("GATE PASS: CIs + negative control valid")


if __name__ == "__main__":
    main()
