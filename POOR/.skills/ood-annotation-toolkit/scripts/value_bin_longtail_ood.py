#!/usr/bin/env python3
"""回归任务值域分 bin 的 LongTail OOD 计算（ood-annotation-toolkit）。

统一口径（kcat `08_add_value_bin_longtail_ood.py` / affinity `12_add_value_bin_longtail_ood.py` /
optimal_ph `05_add_phbin_longtail_ood.py`，均由历史任务脚本演进而来）：
- 把 label（回归值）按固定宽度分 bin；**clip 边界值单独成箱**（可选，
  kcat 的 log10(kcat) 被删失到 [-6, 6]，==±6.0 各作为独立 bin 统计 train
  频次，避免删失堆积遮蔽尾部）
- 统计 **train** 各 bin 频次
- test Default=True 行：所属 bin 的 train 计数 <= 阈值 → True；
  label 缺失/超出 bin 范围 → False
- Default=False 行：默认置 NaN（--nondefault nan，全项目统一惯例），
  可选 false（optimal_ph 历史惯例）

各任务历史参数（重算时对齐，见 SKILL.md）：
- kcat：label=log10(kcat)，bin [-6.5,7.0) 宽 0.5，clip ±6.0，le50/le100，前缀 OOD_LongTail_KcatBin
- affinity：label=-log10(K)，bin [-0.5,16.5) 宽 0.5，无 clip，le50/le100，前缀 OOD_LongTail_AffBin
- optimal_ph：label=pH，bin [1.5,13.5) 宽 0.5，无 clip，le50/le100，前缀 OOD_LongTail_pHbin

用法：
    python value_bin_longtail_ood.py --train kcat_train.csv --test kcat_test.csv \
        --out kcat_test_flagged.csv --bin-start -6.5 --bin-end 7.0 --bin-width 0.5 \
        --clip -6,6 --thresholds 50,100 --col-prefix OOD_LongTail_KcatBin
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--train", required=True, type=Path)
    ap.add_argument("--test", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path,
                    help="输出 csv（可与 --test 同路径就地写回）")
    ap.add_argument("--label-col", default="label")
    ap.add_argument("--bin-start", type=float, required=True)
    ap.add_argument("--bin-end", type=float, required=True)
    ap.add_argument("--bin-width", type=float, default=0.5)
    ap.add_argument("--clip", default=None,
                    help="可选 lo,hi：==lo 与 ==hi 的值各自单独成箱（kcat 口径）")
    ap.add_argument("--thresholds", default="50,100",
                    help="train bin 频次阈值，逗号分隔（<= 阈值 → True）")
    ap.add_argument("--col-prefix", required=True,
                    help="输出列前缀，如 OOD_LongTail_KcatBin → OOD_LongTail_KcatBin_le50")
    ap.add_argument("--default-col", default="Default")
    ap.add_argument("--nondefault", choices=["nan", "false"], default="nan",
                    help="Default=False 行的取值（nan=全项目统一惯例，流水线位置直接写入；"
                         "false=旧惯例，仅为兼容保留）")
    return ap.parse_args()


def main():
    args = parse_args()
    bins = np.arange(args.bin_start, args.bin_end, args.bin_width)
    thresholds = [int(t) for t in args.thresholds.split(",")]
    clip = None
    if args.clip:
        lo, hi = (float(x) for x in args.clip.split(","))
        clip = (lo, hi, f"clip_eq_{lo}", f"clip_eq_{hi}")

    def bin_keys(values):
        """每个值 -> bin 键；clip 边界值单独成箱；范围外 -> NaN。"""
        keys = pd.cut(values, bins=bins).astype(object)
        if clip:
            keys[values == clip[0]] = clip[2]
            keys[values == clip[1]] = clip[3]
        return keys

    train = pd.read_csv(args.train)
    test = pd.read_csv(args.test, dtype=str, keep_default_na=False)
    if args.label_col not in train.columns or args.label_col not in test.columns:
        raise SystemExit(f"错误：train/test 缺少 label 列 {args.label_col!r}")

    train_keys = bin_keys(train[args.label_col].astype(float))
    bin_counts = train_keys.value_counts().to_dict()
    print(f"train bins: {len(bin_counts)}, "
          + ", ".join(f"le{t}={sum(c <= t for c in bin_counts.values())}" for t in thresholds))
    if clip:
        print(f"  clip bins: {clip[2]}={bin_counts.get(clip[2], 0)}, "
              f"{clip[3]}={bin_counts.get(clip[3], 0)}")

    test_keys = bin_keys(pd.to_numeric(test[args.label_col], errors="coerce"))
    mask_default = test[args.default_col] == "True"
    n_nd = int((~mask_default).sum())

    for thr in thresholds:
        col = f"{args.col_prefix}_le{thr}"
        if args.nondefault == "nan":
            new = pd.Series("", index=test.index, dtype=object)
        else:
            new = pd.Series("False", index=test.index, dtype=object)
        flagged = test_keys[mask_default].apply(
            lambda k: bin_counts.get(k, 0) <= thr if pd.notna(k) else False)
        new[mask_default] = flagged.map({True: "True", False: "False"})
        test[col] = new
        print(f"  {col}: True={int(flagged.sum())}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    test.to_csv(args.out, index=False, lineterminator="\n")
    print(f"-> {args.out}（{len(test)} 行，Default=False {n_nd} 行置 "
          f"{'NaN' if args.nondefault == 'nan' else 'False'}）")


if __name__ == "__main__":
    main()
