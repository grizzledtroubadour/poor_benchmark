#!/usr/bin/env python3
"""为 kcat 测试集追加值域分 bin 的 LongTail OOD 列。

口径（对齐 optimal_ph 的 OOD_LongTail_pHbin_*）：
- label = log10(kcat)；上游 CatPred-DB 已将其截断（删失）到 [-6, 6]
  （BRENDA 原始范围 2.24e-11 ~ 9.3e8，==6.0 的真实值 ≥ 1e6 s^-1）
- bin 宽 0.5，边 [-6.5, 7.0)；**clip 边界值单独成箱**（2026-08-27 修复）：
  label==-6.0、label==6.0 各作为独立 bin 统计 train 频次（3 / 92），
  避免 ==6.0 的删失堆积（92 个）把 (5.5,6.0] bin 计数抬到 168 而遮蔽右尾
- 统计 **train** 各 bin 频次
- test Default=True 行：所属 bin 的 train 计数 <= 阈值（50 / 100）→ True
- Default=False（极端长度）行置 NaN（同本任务其他 OOD 列惯例）；label 缺失/无法分箱的 Default 行 → False

I/O：整表 dtype=str 读入，仅覆写两个 KcatBin 列，其余单元格字节级保留；
写回后与 output/backup_before_kcatbin_clipfix/ 的备份做逐单元格 diff 验证。
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
    """每个值 -> bin 键；clip 边界值单独成箱。"""
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

# ---- 字节级验证：除两个目标列外所有单元格与备份一致 -----------------------
old = pd.read_csv(BACKUP, dtype=str, keep_default_na=False)
assert len(old) == len(test) and list(old.columns) == list(test.columns)
cmp_new = test.drop(columns=TARGET_COLS)
cmp_old = old.drop(columns=TARGET_COLS)
assert cmp_new.equals(cmp_old), "非目标单元格发生变化!"

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
