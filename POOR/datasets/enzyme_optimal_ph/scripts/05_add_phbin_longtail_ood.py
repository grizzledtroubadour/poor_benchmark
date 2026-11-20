#!/usr/bin/env python3
"""
Compute pH-bin long-tail OOD flags:
- OOD_LongTail_pHbin_le50 / OOD_LongTail_pHbin_le100: Default test samples whose
  pH bin has low training frequency (width 0.5, range [1.5,13.5), counts <=50 / <=100).

These flags apply only to Default test samples. Default=False extreme-length rows
receive NaN (empty cells: not computed), following project convention.
finalize_ood_columns.py --mode nan is only a fallback check.

Historical note: extracted from the pHbin section of archived compute_additional_ood.py.
Its legacy EC columns (OOD_NewFunction / OOD_LongTailFunction_*) were replaced by
ec-function-ood-annotation skill columns (OOD_NewEC_* / OOD_LongTail_EC_*).
Canonical recomputation tool: .skills/ood-annotation-toolkit/scripts/value_bin_longtail_ood.py.
"""
import argparse
import json
from pathlib import Path

import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SPLITS_DIR = ROOT / "splits"
OUTPUT_DIR = ROOT / "output"

TRAIN_PATH = SPLITS_DIR / "optimal_ph_prediction_train.csv"
TEST_PATH = SPLITS_DIR / "optimal_ph_prediction_test.csv"
SUMMARY_PATH = OUTPUT_DIR / "ood_markers_summary.json"

BIN_EDGES = np.arange(1.5, 13.5, 0.5)
THRESHOLDS = (50, 100)


def main():
    parser = argparse.ArgumentParser(description="Compute OOD_LongTail_pHbin_* flags")
    parser.add_argument("--train", type=Path, default=TRAIN_PATH)
    parser.add_argument("--test", type=Path, default=TEST_PATH)
    parser.add_argument("--summary", type=Path, default=SUMMARY_PATH)
    args = parser.parse_args()

    print("[INFO] Loading splits...")
    train = pd.read_csv(args.train, dtype=str, keep_default_na=False)
    test = pd.read_csv(args.test, dtype=str, keep_default_na=False)

    train_ph = pd.to_numeric(train["label"], errors="coerce")
    test_ph = pd.to_numeric(test["label"], errors="coerce")

    # Use pH intervals of width 0.5.
    train_bin = pd.cut(train_ph, bins=BIN_EDGES)
    test["_ph_bin"] = pd.cut(test_ph, bins=BIN_EDGES)
    bin_counts = train_bin.value_counts().to_dict()

    mask_default = test["Default"] == "True"

    for thr in THRESHOLDS:
        col = f"OOD_LongTail_pHbin_le{thr}"
        # Default=False rows receive empty/NaN; Default rows receive True/False.
        test[col] = ""
        test.loc[mask_default, col] = test.loc[mask_default, "_ph_bin"].apply(
            lambda x: "True" if (pd.notna(x) and bin_counts.get(x, 0) <= thr) else "False"
        )

    # Drop temporary ph_bin intervals, which must not remain in final splits.
    test = test.drop(columns=["_ph_bin"])

    # Save.
    test.to_csv(args.test, index=False, lineterminator="\n")
    print(f"[INFO] Updated {args.test}")

    # Update summary.
    if args.summary.exists():
        with open(args.summary) as f:
            summary = json.load(f)
    else:
        summary = {"n_test": int(len(test))}

    for thr in THRESHOLDS:
        col = f"OOD_LongTail_pHbin_le{thr}"
        summary[col] = int((test[col] == "True").sum())

    with open(args.summary, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[INFO] Summary updated: {args.summary}")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
