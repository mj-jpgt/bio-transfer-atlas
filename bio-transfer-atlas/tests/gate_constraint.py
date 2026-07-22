import pandas as pd


def main():
    c = pd.read_parquet("data/features/selection/constraint_features.parquet")
    required = {"variant_id", "LOEUF", "pLI", "mis_z"}
    assert required.issubset(c.columns), f"missing columns: {required - set(c.columns)}"
    assert c["LOEUF"].notna().mean() > 0.5, "LOEUF coverage too low"
    print("GATE PASS: constraint features attached")


if __name__ == "__main__":
    main()
