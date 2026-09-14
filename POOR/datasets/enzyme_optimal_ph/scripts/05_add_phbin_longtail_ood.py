#!/usr/bin/env python3
"""
计算 pH 值域长尾 OOD 标记：
- OOD_LongTail_pHbin_le50 / OOD_LongTail_pHbin_le100：Default 测试样本的
  pH 落在训练集低频区间（bin 宽 0.5，范围 [1.5, 13.5)，阈值 <= 50 / <= 100）

注意：这些 OOD 标记仅针对 Default 测试样本计算；Default=False（极端长度）行
在本脚本内直接写 NaN（空单元格，"未计算"语义，全项目统一惯例）。
finalize_ood_columns.py --mode nan 仅作兜底校验，不再承担 NaN 化职责。

历史注记：本脚本抽取自已归档的 compute_additional_ood.py 的 pHbin 部分；
原脚本的 EC 口径列（OOD_NewFunction / OOD_LongTailFunction_*）已被
ec-function-ood-annotation skill 的统一口径（OOD_NewEC_* / OOD_LongTail_EC_*）取代。
统一复算工具见 .skills/ood-annotation-toolkit/scripts/value_bin_longtail_ood.py。
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
    parser = argparse.ArgumentParser(description="计算 OOD_LongTail_pHbin_* 标记")
    parser.add_argument("--train", type=Path, default=TRAIN_PATH)
    parser.add_argument("--test", type=Path, default=TEST_PATH)
    parser.add_argument("--summary", type=Path, default=SUMMARY_PATH)
    args = parser.parse_args()

    print("[INFO] Loading splits...")
    train = pd.read_csv(args.train, dtype=str, keep_default_na=False)
    test = pd.read_csv(args.test, dtype=str, keep_default_na=False)

    train_ph = pd.to_numeric(train["label"], errors="coerce")
    test_ph = pd.to_numeric(test["label"], errors="coerce")

    # 使用 0.5 为间隔的 pH 区间
    train_bin = pd.cut(train_ph, bins=BIN_EDGES)
    test["_ph_bin"] = pd.cut(test_ph, bins=BIN_EDGES)
    bin_counts = train_bin.value_counts().to_dict()

    mask_default = test["Default"] == "True"

    for thr in THRESHOLDS:
        col = f"OOD_LongTail_pHbin_le{thr}"
        # Default=False 行直接写 NaN（空），Default 行写 True/False
        test[col] = ""
        test.loc[mask_default, col] = test.loc[mask_default, "_ph_bin"].apply(
            lambda x: "True" if (pd.notna(x) and bin_counts.get(x, 0) <= thr) else "False"
        )

    # 删除临时列（ph_bin 为区间类型，不能保留在最终 splits 中）
    test = test.drop(columns=["_ph_bin"])

    # 保存
    test.to_csv(args.test, index=False, lineterminator="\n")
    print(f"[INFO] Updated {args.test}")

    # 更新摘要
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
