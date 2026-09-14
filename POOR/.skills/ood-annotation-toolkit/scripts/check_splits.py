#!/usr/bin/env python3
"""划分文件一致性检查（ood-annotation-toolkit，只读不写）。

吸收 analysis/scripts/format_consistency.py 的核心检查（A–G）与
survey_columns.py 的列清点，作用于显式给定的三个划分文件：

- 行数 / 列清单（含 train/val/test 列集对比：train==val，test = 核心列 + OOD 块）
- unique_id split 内唯一性
- **F. unique_id 泄露检查**：train∩val、train∩test、val∩test
- 核心列完整性（unique_id / aa_seq(1,2) / struct_file(1,2) / label）
- 废弃列检查（OOD-* 连字符、labels、pdb_file*、protein_*_seq 等）
- label 空值；train/val 不应携带 OOD/Default 列，test 应有 Default + OOD 列
- Default / InD 统计；Default=False 行非极端 OOD_* 列的 NaN 惯例违例计数
- OOD_Orphan ⊆ seq_Redundancy_30 子集关系（若两列都存在）

用法：
    python check_splits.py --train kcat_train.csv --val kcat_val.csv --test kcat_test.csv
    python check_splits.py --train ppi_train.csv --val ppi_val.csv --test ppi_test.csv \
        --seq-cols aa_seq1,aa_seq2 --struct-cols struct_file1,struct_file2
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

DEPRECATED_EXACT = {
    "labels", "file_type", "pdb_file", "pdb_file_A", "pdb_file_B",
    "file_type_A", "file_type_B", "protein_A_seq", "protein_B_seq",
    "is_time_cutoff", "OOD_NewFunction", "OOD_NewFunction_EC",
    "OOD_LongTailFunction", "OOD_LongTail_EC",
    "OOD_CATH_FoldHoldout", "OOD_CATH_SuperfamilyHoldout",
    "OOD_LongTailFunction_L3_le5", "OOD_LongTailFunction_L3_le10",
}
KEEP_NONDEFAULT = {"OOD_ExtremeShort", "OOD_ExtremeLong"}
# test 特有但非 OOD 标记的合法标注列（OOD_IDR 的原始指标）
ANNOTATION_COLS = {"idr_ratio"}


def is_deprecated(col):
    return col in DEPRECATED_EXACT or col.startswith("OOD-")


def ood_like(col):
    return (col.startswith("OOD_") or col.startswith("seq_Redundancy_")
            or col.startswith("TM-score_") or col in ("Default", "InD"))


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--train", required=True, type=Path)
    ap.add_argument("--val", required=True, type=Path)
    ap.add_argument("--test", required=True, type=Path)
    ap.add_argument("--seq-cols", default="aa_seq",
                    help="序列核心列，逗号分隔（ppi 用 aa_seq1,aa_seq2）")
    ap.add_argument("--struct-cols", default="struct_file",
                    help="结构核心列，逗号分隔（ppi 用 struct_file1,struct_file2）")
    ap.add_argument("--label-col", default="label")
    return ap.parse_args()


def main():
    args = parse_args()
    issues = []

    dfs = {}
    for name, p in (("train", args.train), ("val", args.val), ("test", args.test)):
        # label 强制按字符串读，避免 1.10 -> 1.1 之类的浮点折叠影响检查
        df = pd.read_csv(p, dtype=str, keep_default_na=False)
        dfs[name] = df
        print(f"{name}: {p}  rows={len(df)}  cols={len(df.columns)}")

    seq_cols = [c.strip() for c in args.seq_cols.split(",")]
    struct_cols = [c.strip() for c in args.struct_cols.split(",")]
    core = ["unique_id"] + seq_cols + struct_cols + [args.label_col]

    # A/B/C/D：逐文件检查
    for name, df in dfs.items():
        miss = [c for c in core if c not in df.columns]
        if miss:
            issues.append(f"{name}: 缺核心列 {miss}")
        if "unique_id" in df.columns and not df["unique_id"].is_unique:
            n_dup = int(df["unique_id"].duplicated().sum())
            issues.append(f"{name}: unique_id split 内重复 {n_dup} 行")
        dep = [c for c in df.columns if is_deprecated(c)]
        if dep:
            issues.append(f"{name}: 含废弃列 {dep}")
        if args.label_col in df.columns:
            n_null = int((df[args.label_col] == "").sum())
            if n_null:
                issues.append(f"{name}: label 空值 {n_null} 行")
        has_ood = [c for c in df.columns if ood_like(c)]
        if name in ("train", "val") and has_ood:
            issues.append(f"{name}: 不应携带 OOD/Default 列 {has_ood}")
        if name == "test" and ("Default" not in df.columns or len(has_ood) < 2):
            issues.append(f"test: 应有 Default + OOD 列，实际 {has_ood}")

    # F. unique_id 泄露
    idsets = {n: set(df["unique_id"]) for n, df in dfs.items() if "unique_id" in df.columns}
    if len(idsets) == 3:
        for a, b in (("train", "val"), ("train", "test"), ("val", "test")):
            inter = idsets[a] & idsets[b]
            status = f"{len(inter)} 个" if inter else "无"
            print(f"泄露检查 {a}∩{b}: {status}")
            if inter:
                issues.append(f"unique_id 泄露 {a}∩{b}: {len(inter)} 个，"
                              f"如 {sorted(inter)[:5]}")

    # G. 列集对比
    col_train, col_val, col_test = (list(dfs[n].columns) for n in ("train", "val", "test"))
    if set(col_train) != set(col_val):
        issues.append(f"train/val 列集不一致: "
                      f"train-only={sorted(set(col_train)-set(col_val))}, "
                      f"val-only={sorted(set(col_val)-set(col_train))}")
    extra_test = [c for c in col_test if c not in col_train]
    non_ood_extra = [c for c in extra_test
                     if not ood_like(c) and c not in ANNOTATION_COLS
                     and not c.startswith("idr_ratio")]
    missing_test = [c for c in col_train if c not in col_test]
    if non_ood_extra:
        issues.append(f"test 含 train 没有的非 OOD 列 {non_ood_extra}")
    if missing_test:
        issues.append(f"test 缺少 train 有的列 {missing_test}")

    # 列清单
    print("\n=== 列清单 ===")
    for n in ("train", "val", "test"):
        print(f"{n} ({len(dfs[n].columns)} 列): {list(dfs[n].columns)}")

    # Default / InD 统计 + NaN 惯例
    test = dfs["test"]
    if "Default" in test.columns:
        n_false = int((test["Default"] == "False").sum())
        print(f"\nDefault: True={int((test['Default'] == 'True').sum())}, "
              f"False={n_false}")
        if n_false:
            viol = {}
            for c in test.columns:
                if c.startswith("OOD_") and c not in KEEP_NONDEFAULT:
                    n = int((test.loc[test["Default"] == "False", c] != "").sum())
                    if n:
                        viol[c] = n
            if viol:
                issues.append(f"Default=False 行非极端 OOD_* 列应为 NaN，违例: {viol}")
            else:
                print("NaN 惯例：Default=False 行非极端 OOD_* 列全部为空（合规）")
    if "InD" in test.columns:
        print(f"InD: True={int((test['InD'] == 'True').sum())} / {len(test)}")
    else:
        print("InD: 列不存在（可用 finalize_ood_columns.py --mode ind 补算）")

    # Orphan ⊆ seq_Redundancy_30
    if "OOD_Orphan" in test.columns and "seq_Redundancy_30" in test.columns:
        bad = int(((test["OOD_Orphan"] == "True")
                   & (test["seq_Redundancy_30"] != "True")).sum())
        print(f"Orphan ⊆ seq_Redundancy_30: {'成立' if bad == 0 else f'违例 {bad} 行'}")
        if bad:
            issues.append(f"OOD_Orphan 不是 seq_Redundancy_30 的子集（{bad} 行）")

    print("\n=== 检查结果 ===")
    if issues:
        for i in issues:
            print(f"  [FAIL] {i}")
        print(f"共 {len(issues)} 个问题")
        return 1
    print("全部通过")
    return 0


if __name__ == "__main__":
    sys.exit(main())
