#!/usr/bin/env python3
"""为 ligand_binding_affinity 测试集追加值域分 bin 的 LongTail OOD 列。

口径（对齐 optimal_ph 的 OOD_LongTail_pHbin_*）：
- label = -log10(K)（pK 型亲和力），范围约 0.4–15.2；bin 宽 0.5，边 [-0.5, 16.5)
- 统计 **train** 各 bin 频次
- test Default=True 行：所属 bin 的 train 计数 <= 阈值（50 / 100）→ True
- 本任务 test 全部 Default=True；label 缺失/无法分箱行 → False
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent.parent
TRAIN = BASE / "splits" / "ligand_binding_affinity_train.csv"
TEST = BASE / "splits" / "ligand_binding_affinity_test.csv"
OUT_SUMMARY = BASE / "output" / "ligand_longtail_affbin_summary.json"

BINS = np.arange(-0.5, 16.5, 0.5)
THRESHOLDS = (50, 100)

train = pd.read_csv(TRAIN)
test = pd.read_csv(TEST)

bin_counts = pd.cut(train["label"], bins=BINS).value_counts().to_dict()
print(f"train bins: {len(bin_counts)}, min={min(bin_counts.values())}, "
      f"le50={sum(c <= 50 for c in bin_counts.values())}, le100={sum(c <= 100 for c in bin_counts.values())}")

test["_bin"] = pd.cut(test["label"], bins=BINS)
mask_default = test["Default"] == True  # noqa: E712

for thr in THRESHOLDS:
    col = f"OOD_LongTail_AffBin_le{thr}"
    test[col] = False
    test.loc[mask_default, col] = (
        test.loc[mask_default, "_bin"]
        .apply(lambda x: bin_counts.get(x, 0) <= thr if pd.notna(x) else False)
    )

test = test.drop(columns=["_bin"])

summary = {
    "train_rows": int(len(train)),
    "test_rows": int(len(test)),
    "test_default": int(mask_default.sum()),
    "bin_width": 0.5,
    "bin_range": [-0.5, 16.5],
    "non_default_rows": int((~mask_default).sum()),
}
for thr in THRESHOLDS:
    col = f"OOD_LongTail_AffBin_le{thr}"
    summary[col] = int((test[col] == True).sum())  # noqa: E712
print(json.dumps(summary, indent=2))

test.to_csv(TEST, index=False)
OUT_SUMMARY.write_text(json.dumps(summary, indent=2))
print(f"updated {TEST}")
print(f"summary -> {OUT_SUMMARY}")
