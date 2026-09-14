#!/usr/bin/env python3
"""为 func_prediction 的 go_mf / go_bp / go_cc 测试集添加 EC 功能 OOD 六列。

背景：GO 版 OOD（NewGO/LongTailGO）经分析不可行（标签词表封闭 → NewGO 恒 0；
全量 GO 注释 all-not-in 也恒 0）。按用户决定，改为将统一的 EC 功能 OOD 口径
（NewEC / LongTailEC）扩展到 go_* 三个测试集。

口径（与 add_unified_newec.py / add_unified_longtail_ec.py 完全一致）：
- EC 来源：SIFTS 链级 EC（data/sifts/sifts_chain_ec.tsv.gz），
  uid 映射 (u[:4].lower(), u[5:])，与 func_ec 相同
- levels()：段数>=3 → L3=前三段；段数==4 → L4（允许部分 "-"）
- OOD_NewEC_L3/L4：参考 train+val，有 EC 且所有 EC（对应层级）∉ 参考集 → True
  （all-not-in）；无 EC → False
- OOD_LongTail_EC_{L3,L4}_le{5,10}：参考 train only（逐样本频次，含 0），
  有 EC 且所有 EC（对应层级）train 频次 < 阈值 → True（all-not-in）；无 EC → False
- Default=False 行（极端长度预挑出样本）：六列均为 NaN（未计算，区别于 False）
- 六列追加到 csv 末尾

用法：python3 .skills/ec-function-ood-annotation/scripts/add_go_ec_ood.py
"""
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

# 仓库根：.skills/ec-function-ood-annotation/scripts/ 上溯 3 级
POOR = Path(__file__).resolve().parents[3]
SIFTS_EC_TSV = POOR / "data" / "sifts" / "sifts_chain_ec.tsv.gz"
FUNC = POOR / "datasets" / "func_prediction"

L3_COL, L4_COL = "OOD_NewEC_L3", "OOD_NewEC_L4"
LT_COLS = ["OOD_LongTail_EC_L3_le5", "OOD_LongTail_EC_L3_le10",
           "OOD_LongTail_EC_L4_le5", "OOD_LongTail_EC_L4_le10"]


def levels(ec):
    parts = ec.split(".")
    return (".".join(parts[:3]) if len(parts) >= 3 else None,
            ec if len(parts) == 4 else None)


def load_sifts_chain_ec():
    m = defaultdict(set)
    with gzip.open(SIFTS_EC_TSV, "rt") as f:
        next(f)
        for line in f:
            p = line.strip().split("\t")
            if len(p) >= 4 and p[3] and p[3] != "-":
                m[(p[0].lower(), p[1])].add(p[3])
    return m


def main():
    sifts = load_sifts_chain_ec()
    ec_of = lambda u: sifts.get((u[:4].lower(), u[5:]), set())  # noqa: E731

    summary = {}
    for task in ["go_mf", "go_bp", "go_cc"]:
        test_path = FUNC / "splits" / f"{task}_test.csv"
        test = pd.read_csv(test_path)
        train_ids = pd.read_csv(FUNC / "splits" / f"{task}_train.csv",
                                usecols=["unique_id"])["unique_id"]
        val_ids = pd.read_csv(FUNC / "splits" / f"{task}_val.csv",
                              usecols=["unique_id"])["unique_id"]

        # 参考集：NewEC 用 train+val；LongTailEC 用 train-only 逐样本频次
        tv3, tv4 = set(), set()
        f3, f4 = Counter(), Counter()
        for u in train_ids:
            s = ec_of(u)
            l3s = {levels(e)[0] for e in s} - {None}
            l4s = {levels(e)[1] for e in s} - {None}
            tv3 |= l3s
            tv4 |= l4s
            for x in l3s:
                f3[x] += 1
            for x in l4s:
                f4[x] += 1
        for u in val_ids:
            s = ec_of(u)
            tv3 |= {levels(e)[0] for e in s} - {None}
            tv4 |= {levels(e)[1] for e in s} - {None}

        def flags(u):
            s = ec_of(u)
            l3s = {levels(e)[0] for e in s} - {None}
            l4s = {levels(e)[1] for e in s} - {None}
            return {
                L3_COL: bool(l3s) and all(x not in tv3 for x in l3s),
                L4_COL: bool(l4s) and all(x not in tv4 for x in l4s),
                LT_COLS[0]: bool(l3s) and all(f3.get(x, 0) < 5 for x in l3s),
                LT_COLS[1]: bool(l3s) and all(f3.get(x, 0) < 10 for x in l3s),
                LT_COLS[2]: bool(l4s) and all(f4.get(x, 0) < 5 for x in l4s),
                LT_COLS[3]: bool(l4s) and all(f4.get(x, 0) < 10 for x in l4s),
            }

        mask = test["Default"] == True  # noqa: E712
        cols = [L3_COL, L4_COL] + LT_COLS
        lut = {u: flags(u) for u in test.loc[mask, "unique_id"]}
        for col in cols:
            ser = pd.Series(pd.NA, index=test.index, dtype=object)
            ser[mask] = test.loc[mask, "unique_id"].map(
                {u: f[col] for u, f in lut.items()}).to_numpy()
            test[col] = ser
        test.to_csv(test_path, index=False)

        counts = {col: int((test[col] == True).sum()) for col in cols}  # noqa: E712
        summary[task] = {"test_rows": int(len(test)),
                         "computed_rows": int(mask.sum()),
                         "no_ec_rows": int(sum(1 for u in test.loc[mask, "unique_id"]
                                               if not ec_of(u))),
                         **counts}
        print(f"{task}: computed={int(mask.sum())}, " +
              " ".join(f"{c.replace('OOD_', '')}={n}" for c, n in counts.items()))

    (FUNC / "output" / "go_ec_ood_summary.json").write_text(
        json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
