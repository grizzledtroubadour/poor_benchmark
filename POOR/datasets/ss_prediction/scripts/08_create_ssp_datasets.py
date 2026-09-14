#!/usr/bin/env python3
"""Stage 3.1：SSP 数据集划分（ssp_train/val/test.csv）。

从 Stage 1/2 的中间产物一次性重建划分：

1. 输入 ``output/ssp_full_data.csv``（代表链：unique_id/aa_seq/labels），
   链长取自 ``output/stage4_cluster_representatives.csv``，沉积日期取自
   ``output/pdb_metadata_api.csv``。
2. 单链长度 > 2000 的样本直接剔除（全库仅 4O9X_A，2,082 残基）。
3. QC 剔除 ``QC_EXCLUDE`` 中的低质量链（7UYL_K：59% UNK 残基，序列提取
   跳过 UNK 而结构保留导致标注错帧，见 DATA_PROCESS.md 注意事项）。
4. 极长链（长度 > 1000，82 条）不参与时间切分，直接进入 test，
   标记 ``Default=False`` / ``OOD_ExtremeLong=True``，按 unique_id 排序附在
   test 末尾；其余 OOD 列由后续脚本/skill 标注（Default=False 行置 NaN）。
5. 其余链按 deposition_date 升序（稳定排序）后按 **72:8:20 分位切分**：
   ``n_train = int(n*0.72)``、``n_val = round(n*0.08)``、剩余进 test。
   时间切分样本 ``Default=True`` / ``OOD_ExtremeLong=False``。
   注：cutoff 为运行时对当前池动态计算的分位数，不钉死；当前发布数据的
   实际 cutoff 为 2019-06-06 / 2021-02-09（系历史路线下对 ≥60 池分位
   计算、短链后按该边界补入，见 DATA_PROCESS.md 历史沿革），后续复现
   的 cutoff 随池组成变化，不要求与之一致。

输出列：train/val 为 ``unique_id,aa_seq,struct_file,label``；test 追加
``Default,OOD_ExtremeLong``。``struct_file = unique_id + ".pdb"``。

用法：
    python scripts/08_create_ssp_datasets.py                    # 写入 splits/
    python scripts/08_create_ssp_datasets.py --out-dir output/verify
"""
import argparse
from pathlib import Path

import pandas as pd

TASK_DIR = Path(__file__).resolve().parent.parent
OUTPUT_DIR = TASK_DIR / "output"

FULL_DATA = OUTPUT_DIR / "ssp_full_data.csv"
REPRESENTATIVES = OUTPUT_DIR / "stage4_cluster_representatives.csv"
METADATA = OUTPUT_DIR / "pdb_metadata_api.csv"

HARD_MAX_LEN = 2000   # 单链 > 2000 直接剔除
EXTREME_LONG = 1000   # > 1000 直入 test

# 72:8:20 时间切分（旧版口径，cutoff 运行时动态计算）
TRAIN_FRAC, VAL_FRAC = 0.72, 0.08

# QC 剔除：7UYL_K 的 120 个结构残基中 71 个为 UNK（59%），序列提取按
# AA_MAP 跳过 UNK（aa_seq 仅 49 残基）而结构提取保留 UNK，导致 label 与
# struct_label 契约破坏、标注错帧。
QC_EXCLUDE = {"7UYL_K"}

BASE_COLS = ["unique_id", "aa_seq", "struct_file", "label"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", type=Path, default=TASK_DIR / "splits",
                    help="输出目录（默认 splits/，试跑用 output/verify）")
    args = ap.parse_args()

    full = pd.read_csv(FULL_DATA, dtype=str, keep_default_na=False)
    rep = pd.read_csv(REPRESENTATIVES, dtype=str, keep_default_na=False)
    meta = pd.read_csv(METADATA, dtype=str, keep_default_na=False)

    rep["num_residues"] = rep["num_residues"].astype(int)
    df = full.merge(rep[["chain_key", "pdb_id", "num_residues"]],
                    left_on="unique_id", right_on="chain_key", how="left")
    df = df.merge(meta[["pdb_id", "deposition_date"]], on="pdb_id", how="left")
    assert df["num_residues"].notna().all() and (df["deposition_date"] != "").all()

    # 单链 > 2000 硬性剔除（4O9X_A，2,082 残基）
    dropped = df[df["num_residues"] > HARD_MAX_LEN]
    df = df[df["num_residues"] <= HARD_MAX_LEN].copy()
    print(f"dropped len>{HARD_MAX_LEN}: {len(dropped)} "
          f"({', '.join(dropped['unique_id'])})")

    # QC 剔除低质量链
    qc = df[df["unique_id"].isin(QC_EXCLUDE)]
    df = df[~df["unique_id"].isin(QC_EXCLUDE)].copy()
    print(f"dropped QC_EXCLUDE: {len(qc)} ({', '.join(qc['unique_id'])})")

    extreme = df[df["num_residues"] > EXTREME_LONG].sort_values("unique_id")
    rest = df[df["num_residues"] <= EXTREME_LONG].sort_values(
        "deposition_date", kind="stable")
    print(f"extreme-long (>{EXTREME_LONG}): {len(extreme)}, time-split pool: {len(rest)}")

    # 72:8:20 时间切分（deposition_date 升序稳定排序后按数量分位）
    n = len(rest)
    n_train = int(n * TRAIN_FRAC)
    n_val = round(n * VAL_FRAC)
    train = rest.iloc[:n_train]
    val = rest.iloc[n_train:n_train + n_val]
    test_time = rest.iloc[n_train + n_val:]

    def base(frame):
        out = frame[["unique_id", "aa_seq"]].copy()
        out["struct_file"] = out["unique_id"] + ".pdb"
        out["label"] = frame["labels"]
        return out[BASE_COLS]

    test = pd.concat([base(test_time), base(extreme)])
    test["Default"] = ["True"] * len(test_time) + ["False"] * len(extreme)
    test["OOD_ExtremeLong"] = ["False"] * len(test_time) + ["True"] * len(extreme)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    for name, frame in [("ssp_train", base(train)), ("ssp_val", base(val)),
                        ("ssp_test", test)]:
        path = args.out_dir / f"{name}.csv"
        frame.to_csv(path, index=False, lineterminator="\n")
        print(f"  {path}  rows={len(frame)}")
    print(f"cutoffs (computed): train <= {train['deposition_date'].max()}, "
          f"val <= {val['deposition_date'].max()}")
    print(f"date range: train {train['deposition_date'].min()} ~ {train['deposition_date'].max()}, "
          f"val {val['deposition_date'].min()} ~ {val['deposition_date'].max()}, "
          f"test {test_time['deposition_date'].min()} ~ {test_time['deposition_date'].max()}")


if __name__ == "__main__":
    main()
