#!/usr/bin/env python3
"""Add value-bin long-tail OOD columns to the kcat test set.

Definition (aligned with optimal_ph OOD_LongTail_pHbin_*):
- label = log10(kcat), clipped/censored upstream by CatPred-DB to [-6, 6].
  Original BRENDA range: 2.24e-11 to 9.3e8; label==6.0 represents true values >=1e6 s^-1.
- Bin width 0.5, edges [-6.5, 7.0). Keep clipped boundary values in separate bins (fix 2026-08-27):
  label==-6.0 and label==6.0 each form independent bins, with training counts 3 / 92.
  This prevents 92 censored values at 6.0 inflating the (5.5,6.0] count to 168 and hiding the right tail.
- Count training samples per bin.
- Default=True test rows: flag True when training bin count <=50 / <=100.
- Default=False extreme-length rows receive NaN; missing/unbinnable Default labels receive False.

I/O: read all columns as strings; update only the two KcatBin columns and preserve other cells bytewise.
After writing, compare cell by cell against output/backup_before_kcatbin_clipfix/.
"""
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent.parent
TRAIN = BASE / "splits" / "kcat_train.csv"
TEST = BASE / "splits" / "kcat_test.csv"
BACKUP = BASE / "output" / "backup_before_kcatbin_clipfix" / "kcat_test.csv"
OUT_SUMMARY = BASE / "output" / "kcat_longtail_kcatbin_summary.json"

BINS = np.arange(-6.5, 7.0, 0.5)
CLIP_LO, CLIP_HI = -6.0, 6.0
CLIP_LO_KEY, CLIP_HI_KEY = "clip_eq_-6.0", "clip_eq_6.0"
THRESHOLDS = (50, 100)
TARGET_COLS = [f"OOD_LongTail_KcatBin_le{t}" for t in THRESHOLDS]


def bin_keys(values: pd.Series) -> pd.Series:
    """Map each value to a bin key; clipped boundary values use separate bins."""
    keys = pd.cut(values, bins=BINS).astype(object)
    keys[values == CLIP_LO] = CLIP_LO_KEY
    keys[values == CLIP_HI] = CLIP_HI_KEY
    return keys


train = pd.read_csv(TRAIN)
test = pd.read_csv(TEST, dtype=str, keep_default_na=False)
test_label = pd.to_numeric(test["label"], errors="coerce")

train_keys = bin_keys(train["label"].astype(float))
bin_counts = train_keys.value_counts().to_dict()
print(f"train bins: {len(bin_counts)}, "
      f"le50={sum(c <= 50 for c in bin_counts.values())}, "
      f"le100={sum(c <= 100 for c in bin_counts.values())}")
print(f"  clip bins: {CLIP_LO_KEY}={bin_counts.get(CLIP_LO_KEY, 0)}, "
      f"{CLIP_HI_KEY}={bin_counts.get(CLIP_HI_KEY, 0)}")

test_keys = bin_keys(test_label)
mask_default = test["Default"] == "True"

for thr, col in zip(THRESHOLDS, TARGET_COLS):
    new = pd.Series("", index=test.index, dtype=object)
    flagged = test_keys[mask_default].apply(
        lambda k: bin_counts.get(k, 0) <= thr if pd.notna(k) else False
    )
    new[mask_default] = flagged.map({True: "True", False: "False"})
    test[col] = new

# ---- Byte-level validation: non-target cells must match the backup. -----------------------
old = pd.read_csv(BACKUP, dtype=str, keep_default_na=False)
assert len(old) == len(test) and list(old.columns) == list(test.columns)
cmp_new = test.drop(columns=TARGET_COLS)
cmp_old = old.drop(columns=TARGET_COLS)
assert cmp_new.equals(cmp_old), "Non-target cells changed!"

n_na = int((~mask_default).sum())
summary = {
    "train_rows": int(len(train)),
    "test_rows": int(len(test)),
    "test_default": int(mask_default.sum()),
    "bin_width": 0.5,
    "bin_range": [-6.5, 7.0],
    "clip_bins_separate": {CLIP_LO_KEY: bin_counts.get(CLIP_LO_KEY, 0),
                           CLIP_HI_KEY: bin_counts.get(CLIP_HI_KEY, 0)},
    "non_default_rows_set_to_NA": n_na,
}
for col in TARGET_COLS:
    summary[col] = int((test[col] == "True").sum())

test.to_csv(TEST, index=False, lineterminator="\n")
OUT_SUMMARY.write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
print(f"updated {TEST}")
print(f"summary -> {OUT_SUMMARY}")
