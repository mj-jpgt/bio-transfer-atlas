import pandas as pd


def main():
    d = pd.read_parquet("data/labels/gwas_concordance_labels_multisource.parquet")
    assert len(d) > 10_000, "too few rows in multisource labels"
    assert d["n_sources"].max() >= 3, "need >=3 sources for at least some variants"
    assoc = d[d["associated_multisource"]] if "associated_multisource" in d.columns else d[d["associated"]]
    assert len(assoc) > 0, "associated multisource subset empty"
    assert assoc["I2"].median() > 0.0, "associated multisource median I2 should be > 0"
    print("GATE PASS: multi-source labels with positive heterogeneity signal")


if __name__ == "__main__":
    main()
