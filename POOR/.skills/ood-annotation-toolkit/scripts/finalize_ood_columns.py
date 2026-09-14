#!/usr/bin/env python3
"""OOD 列收尾：Default=False 行 NaN 化 + InD 列重算（ood-annotation-toolkit）。

角色说明（2026-08-28 起）：各任务流水线脚本在生成各 OOD 列的位置直接为
Default=False 行写入 NaN（"未计算"语义）；本脚本的 `--mode nan` 转为
**兜底校验**（幂等，期望为 no-op 或仅修正遗漏），不再承担主要的 NaN 化职责。

吸收并通用化两个 analysis 脚本的口径：

1. `--mode nan`（原 analysis/scripts/nondefault_ood_nan.py）：
   Default==False 的行，除 OOD_ExtremeShort / OOD_ExtremeLong 外的所有
   OOD_* 列置为 NaN（空单元格），表示"未计算"，与 False 区分。
   特殊文件如需跳过 nan（整表 Default=False 属刻意设计时），传
   `--exempt-nondefault`（现有任务均不需要：PDA 两个 design 集合
   已整表 Default=True）。
   注意覆盖 OOD_ 前缀列与 seq_Redundancy_* / TM-score_*；idr_ratio 等
   连续标注列不在本模式范围（其 NaN 由各任务生成脚本在流水线位置负责）。

2. `--mode ind`（原 analysis/scripts/add_ind_column.py）：
   InD = 所有 OOD_* + seq_Redundancy_* + TM-score_* 标记均为 False（NaN
   视为 False；Default=False 行的 ExtremeShort/Long 为 True，自然 InD=False）。
   已有 InD 列 → 原位覆写（幂等）；有 Pure_ID 列（fold / PDA 旧名）→
   原位更名 InD 并重算；都没有 → 追加到末尾。

`--mode both` = 先 nan 后 ind。`--dry-run` 只报告差异不写回。

实现：整表 dtype=str, keep_default_na=False 读入，仅覆写目标列，
to_csv(index=False, lineterminator="\\n")，其余单元格字节级保留。

用法：
    python finalize_ood_columns.py --test-csv kcat_test.csv --mode both
    python finalize_ood_columns.py --test-csv cath_test.csv --mode both --dry-run
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

FLAG_PREFIXES = ("OOD_", "seq_Redundancy_", "TM-score_")
KEEP_NONDEFAULT = {"OOD_ExtremeShort", "OOD_ExtremeLong"}


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--test-csv", required=True, type=Path)
    ap.add_argument("--mode", required=True, choices=["nan", "ind", "both"])
    ap.add_argument("--default-col", default="Default")
    ap.add_argument("--exempt-nondefault", action="store_true",
                    help="豁免 nan 模式（仅用于整表 Default=False 属刻意设计的特殊文件；"
                         "现有任务均不需要）")
    ap.add_argument("--dry-run", action="store_true", help="只报告差异，不写回")
    return ap.parse_args()


def flag_columns(df):
    return [c for c in df.columns if c.startswith(FLAG_PREFIXES)]


def mode_nan(df, default_col, dry_run):
    """Default=False 行的非极端长度 OOD 标记列置 NaN（兜底校验，幂等）。

    覆盖 OOD_*（除 OOD_ExtremeShort/Long）与 seq_Redundancy_* / TM-score_*；
    idr_ratio 等连续标注列不在范围。返回 (df, 报告)。"""
    if default_col not in df.columns:
        print(f"  nan: 无 {default_col} 列，跳过")
        return df, {}
    nd = df[default_col] == "False"
    n_nd = int(nd.sum())
    ood_cols = [c for c in df.columns if c.startswith(FLAG_PREFIXES)
                and c not in KEEP_NONDEFAULT]
    if n_nd == 0 or not ood_cols:
        print(f"  nan: no-op（Default=False 行数={n_nd}，标记列数={len(ood_cols)}）")
        return df, {"nondefault_rows": n_nd, "cleared": {}}
    before = {c: int((df.loc[nd, c] != "").sum()) for c in ood_cols}
    cleared = {c: n for c, n in before.items() if n}
    print(f"  nan: {len(ood_cols)} 列 x {n_nd} 行置 NaN；"
          f"清除非空单元格：{cleared or '无（本就全空）'}")
    if not dry_run:
        for c in ood_cols:
            # 不用 astype(object)：pandas 3 的 dtype=str 读入产生 str dtype 列，
            # astype(object) 会改变列 dtype 导致下方防御性 equals 检查误报。
            ser = df[c].copy()
            ser.loc[nd] = ""
            df[c] = ser
    return df, {"nondefault_rows": n_nd, "cleared": cleared}


def mode_ind(df, dry_run):
    """重算 InD 列。返回 (df, 报告)。"""
    renamed = False
    recompute = "InD" in df.columns
    if not recompute and "Pure_ID" in df.columns:
        df = df.rename(columns={"Pure_ID": "InD"})
        renamed = True
    cols = flag_columns(df)
    if not cols:
        print("  ind: 无任何 OOD 标记列，跳过")
        return df, {}
    flags = df[cols].apply(lambda s: s.eq("True"))
    ind = ~flags.any(axis=1)
    new_vals = ind.map({True: "True", False: "False"})
    action = "覆写" if recompute else ("Pure_ID 原位更名" if renamed else "末尾追加")
    if recompute or renamed:
        diff = int((df["InD"] != new_vals).sum())
        print(f"  ind: {action}重算，{len(cols)} 个标记列；"
              f"InD=True {int(ind.sum())}/{len(df)} ({ind.sum()/len(df)*100:.1f}%)；"
              f"与现有值差异 {diff} 行")
    else:
        diff = None
        print(f"  ind: {action}，{len(cols)} 个标记列；"
              f"InD=True {int(ind.sum())}/{len(df)} ({ind.sum()/len(df)*100:.1f}%)")
    if not dry_run:
        df["InD"] = new_vals
    return df, {"n_flag_cols": len(cols), "ind_true": int(ind.sum()),
                "renamed": renamed, "recomputed": recompute, "diff_rows": diff}


def main():
    args = parse_args()
    df = pd.read_csv(args.test_csv, dtype=str, keep_default_na=False)
    print(f"{args.test_csv}: {len(df)} 行, {len(df.columns)} 列, "
          f"mode={args.mode}{' (dry-run)' if args.dry_run else ''}")
    orig = df.copy()

    if args.mode in ("nan", "both"):
        if args.exempt_nondefault:
            print(f"  nan: --exempt-nondefault，豁免跳过（{args.test_csv.name}）")
        else:
            df, _ = mode_nan(df, args.default_col, args.dry_run)
    if args.mode in ("ind", "both"):
        df, _ = mode_ind(df, args.dry_run)

    if args.dry_run:
        changed_cols = [c for c in df.columns
                        if c in orig.columns and not df[c].equals(orig[c])]
        new_cols = [c for c in df.columns if c not in orig.columns]
        print(f"dry-run：将改动的列={changed_cols}，新增列={new_cols}；未写回。")
        return

    # 防御：除目标列外不得有变化（nan 只动 OOD_* + seq_Redundancy_* / TM-score_*，
    # ind 只动 InD/Pure_ID）
    touched = {"InD", "Pure_ID"} | {
        c for c in df.columns if c.startswith(FLAG_PREFIXES) and c not in KEEP_NONDEFAULT}
    for c in orig.columns:
        if c not in touched and c in df.columns:
            assert df[c].equals(orig[c]), f"非目标列 {c} 发生变化!"
    df.to_csv(args.test_csv, index=False, lineterminator="\n")
    print(f"-> 已写回 {args.test_csv}")


if __name__ == "__main__":
    sys.exit(main())
