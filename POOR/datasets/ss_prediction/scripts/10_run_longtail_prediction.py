#!/usr/bin/env python3
"""
Compute OOD-LongTail based on 8-state DSSP label distribution.

Definition:
- Long-tail states: I (π-helix), P (Polyproline II), B (β-bridge)
  These are the three rarest states in the global 8-state distribution.
- For each protein, compute: LT_ratio = (I_count + P_count + B_count) / valid_residues
- OOD-LongTail: proteins in the top N% by LT_ratio.

Default: top 5% (configurable via --percentile).
"""

import argparse
import csv
from pathlib import Path

import pandas as pd

TASK_DIR = Path(__file__).resolve().parent.parent  # datasets/ss_prediction

LONGTAIL_STATES = {2, 4, 7}  # I, B, P
DEFAULT_PERCENTILE = 5


def compute_lt_ratio(labels_str: str) -> tuple:
    """Parse space-separated labels and compute LT_ratio."""
    labels = [int(x) for x in labels_str.strip().strip("[]").split() if x]
    valid = [l for l in labels if l != -1]
    valid_len = len(valid)
    if valid_len == 0:
        return 0, 0, valid_len
    lt_count = sum(1 for l in valid if l in LONGTAIL_STATES)
    return lt_count / valid_len, lt_count, valid_len


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--percentile", type=float, default=DEFAULT_PERCENTILE,
                        help="Top N%% by LT_ratio marked as OOD-LongTail (default: 5)")
    parser.add_argument("--full-data", default=str(TASK_DIR / "output/ssp_full_data.csv"))
    parser.add_argument("--test-csv", default=str(TASK_DIR / "splits/ssp_test.csv"))
    parser.add_argument("--output", default=str(TASK_DIR / "output/longtail_predictions.csv"))
    args = parser.parse_args()

    print(f"Loading full data: {args.full_data}")
    df = pd.read_csv(args.full_data)
    print(f"Total proteins: {len(df)}")

    # Compute LT_ratio for all proteins
    results = []
    for _, row in df.iterrows():
        lt_ratio, lt_count, valid_len = compute_lt_ratio(str(row["labels"]))
        results.append({
            "unique_id": row["unique_id"],
            "valid_len": valid_len,
            "lt_count": lt_count,
            "lt_ratio": round(lt_ratio, 6),
        })

    rdf = pd.DataFrame(results)

    # Determine threshold: top N% by LT_ratio
    threshold = rdf["lt_ratio"].quantile(1 - args.percentile / 100)
    rdf["ood_longtail"] = rdf["lt_ratio"] >= threshold

    print(f"\nLongTail threshold (top {args.percentile}%): {threshold:.6f}")
    print(f"OOD-LongTail count: {rdf['ood_longtail'].sum()} ({rdf['ood_longtail'].mean()*100:.2f}%)")
    print(f"LT_ratio stats:")
    print(rdf["lt_ratio"].describe(percentiles=[0.5, 0.75, 0.9, 0.95, 0.99]))

    # Save predictions
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    rdf.to_csv(args.output, index=False)
    print(f"\nSaved: {args.output}")

    # Merge into test CSV
    test_df = pd.read_csv(args.test_csv)
    test_df = test_df.merge(
        rdf[["unique_id", "ood_longtail"]].rename(columns={"ood_longtail": "OOD_LongTail"}),
        on="unique_id",
        how="left",
    )
    test_df["OOD_LongTail"] = test_df["OOD_LongTail"].fillna(False).astype(bool)

    # Default=False（极端长度）行直接置 NaN（"未计算"语义，全项目统一惯例，
    # 流水线位置写入；finalize_ood_columns.py --mode nan 仅作兜底校验）
    if "Default" in test_df.columns:
        nd = ~test_df["Default"].astype(bool)
        if nd.any():
            test_df["OOD_LongTail"] = test_df["OOD_LongTail"].astype(object)
            test_df.loc[nd, "OOD_LongTail"] = pd.NA
            print(f"Default=False rows -> NaN: {int(nd.sum())}")

    # Move to end
    cols = [c for c in test_df.columns if c != "OOD_LongTail"] + ["OOD_LongTail"]
    test_df = test_df[cols]
    test_df.to_csv(args.test_csv, index=False)

    n_test = int((test_df["OOD_LongTail"] == True).sum())  # noqa: E712
    print(f"Test set OOD-LongTail: {n_test} / {len(test_df)} ({n_test/len(test_df)*100:.2f}%)")


if __name__ == "__main__":
    main()
