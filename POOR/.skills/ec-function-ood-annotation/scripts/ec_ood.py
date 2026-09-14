#!/usr/bin/env python3
"""NewEC / LongTailEC OOD 计算（ec-function-ood-annotation skill 第 2 步）。

统一口径（2026-08-27 定稿）：
- 层级解析 levels()：EC 段数 ≥3 → L3=前三段；段数 ==4 → L4=完整编号；允许部分 "-"
- **NewEC**：参考集 = train+val 出现的 EC 集合；**all-not-in** 语义——
  样本有 EC 注释且其所有 EC（对应层级）都不在参考集中 → True；无 EC 注释 → False
- **LongTailEC**：参考集 = train（逐样本计数：携带该类别的 train 样本数）；
  **all-not-in** 语义——样本有 EC 注释且其所有 EC（对应层级）的 train 频次都 < 5 / < 10（**含 0**）→ True

输入两种形式（与 cath_ted_holdout.py 一致）：
- 单文件：含 `split` 列（train/val/test，`;` 分隔多归属）+ `ec_numbers` 列
- 双文件：--labels 为 test，--ref-labels 为 train+val（需含 `split` 列区分）

输出：仅 test 行，追加 6 列：
  OOD_NewEC_L3, OOD_NewEC_L4,
  OOD_LongTail_EC_L3_le5, OOD_LongTail_EC_L3_le10,
  OOD_LongTail_EC_L4_le5, OOD_LongTail_EC_L4_le10

用法：
    python ec_ood.py --labels all_ec_labels.csv --split-col split --out test_flagged.csv
    python ec_ood.py --labels test_ec.csv --ref-labels trainval_ec.csv --out test_flagged.csv
"""
import argparse
import csv
from collections import Counter
from pathlib import Path

NEW_COLS = ["OOD_NewEC_L3", "OOD_NewEC_L4",
            "OOD_LongTail_EC_L3_le5", "OOD_LongTail_EC_L3_le10",
            "OOD_LongTail_EC_L4_le5", "OOD_LongTail_EC_L4_le10"]
THRS = (5, 10)


def levels(ec):
    parts = ec.split(".")
    return (".".join(parts[:3]) if len(parts) >= 3 else None,
            ec if len(parts) == 4 else None)


def parse_ecs(s):
    return {x.strip() for x in (s or "").split(";") if x.strip() and x.strip() != "-"}


def level_sets(ecs):
    return ({levels(e)[0] for e in ecs} - {None},
            {levels(e)[1] for e in ecs} - {None})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", required=True, type=Path, help="test 行（或全量，配 --split-col）")
    ap.add_argument("--ref-labels", type=Path, help="train+val 参考文件（需 split 列）")
    ap.add_argument("--split-col", default="split")
    ap.add_argument("--ec-col", default="ec_numbers")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    with open(args.labels) as f:
        rows = list(csv.DictReader(f))
    fieldnames = list(rows[0].keys())

    if args.ref_labels:
        with open(args.ref_labels) as f:
            ref_rows = list(csv.DictReader(f))
        test_rows = rows
    else:
        sc = args.split_col
        ref_rows = [r for r in rows if "train" in r[sc].split(";") or "val" in r[sc].split(";")]
        test_rows = [r for r in rows if "test" in r[sc].split(";")]

    # train+val 集合（NewEC 参考）与 train 频次（LongTailEC 参考）
    tv3, tv4 = set(), set()
    f3, f4 = Counter(), Counter()
    for r in ref_rows:
        ecs = parse_ecs(r[args.ec_col])
        l3s, l4s = level_sets(ecs)
        tv3 |= l3s
        tv4 |= l4s
        is_train = ("train" in r[args.split_col].split(";")) if args.split_col in r else True
        if is_train:
            for x in l3s:
                f3[x] += 1
            for x in l4s:
                f4[x] += 1
    print(f"ref rows: {len(ref_rows)}, train+val L3={len(tv3)} L4={len(tv4)}; "
          f"train L3={len(f3)} L4={len(f4)}")

    for col in NEW_COLS:
        if col not in fieldnames:
            fieldnames.append(col)
    counts = Counter()
    for r in test_rows:
        ecs = parse_ecs(r[args.ec_col])
        l3s, l4s = level_sets(ecs)
        vals = {
            "OOD_NewEC_L3": bool(l3s) and all(x not in tv3 for x in l3s),
            "OOD_NewEC_L4": bool(l4s) and all(x not in tv4 for x in l4s),
            "OOD_LongTail_EC_L3_le5": bool(l3s) and all(f3.get(x, 0) < 5 for x in l3s),
            "OOD_LongTail_EC_L3_le10": bool(l3s) and all(f3.get(x, 0) < 10 for x in l3s),
            "OOD_LongTail_EC_L4_le5": bool(l4s) and all(f4.get(x, 0) < 5 for x in l4s),
            "OOD_LongTail_EC_L4_le10": bool(l4s) and all(f4.get(x, 0) < 10 for x in l4s),
        }
        for col, v in vals.items():
            r[col] = str(v)
            counts[col] += v

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(test_rows)
    print(f"test rows: {len(test_rows)}")
    for col in NEW_COLS:
        print(f"  {col}: {counts[col]}")
    print(f"-> {args.out}")


if __name__ == "__main__":
    main()
