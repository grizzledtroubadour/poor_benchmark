#!/usr/bin/env python3
"""为 enzyme_kinetics (kcat) test 集追加 OOD_FoldHoldout / OOD_SuperfamilyHoldout（CATH+TED 口径）。

样本键为 UniProt accession（unique_id = AF-<acc>-F1-model_v6，09 修复后
可能带配体哈希后缀 AF-<acc>-F1-model_v6-<hash8>，extract_acc 两者兼容）：
skill 的 uniprot 模式（CATH 经 SIFTS 反查并集 + TED 全长标注），
train-only 对称参考、any 语义、无注释 -> False、极端长度行（Default=False）-> NaN。

流程：
  0. build_acc_list：从 splits/kcat_{train,val,test}.csv 生成 accession 清单
     output/kcat_ted_aug_acc_list.csv（补缺口：原为历史一次性生成、脚本已删除）
  1. cath_ted_annotate.py  生成 accession 级 CATH+TED 标签（--offline 仅用缓存）
  2. cath_ted_holdout.py   按统一口径重算 holdout（参考集 = train 行并集，any 语义）
两段均调用 .skills/cath-ted-domain-annotation 的统一脚本；
本脚本只负责把 holdout 结果按 accession 映射回 splits/kcat_test.csv 逐行写回
（同一 accession 多底物多行），并执行 Default=False 行置 NaN 的项目统一惯例。
"""
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd

TASK_DIR = Path(__file__).resolve().parent.parent
ROOT = TASK_DIR.parent.parent
SKILL_DIR = ROOT / ".skills/cath-ted-domain-annotation/scripts"
ANNOTATE = SKILL_DIR / "cath_ted_annotate.py"
HOLDOUT = SKILL_DIR / "cath_ted_holdout.py"
ACC_LIST = TASK_DIR / "output/kcat_ted_aug_acc_list.csv"
LABELS = TASK_DIR / "output/kcat_acc_cath_ted_labels.csv"
HOLDOUT_OUT = TASK_DIR / "output/kcat_cath_ted_holdout.csv"
TEST_CSV = TASK_DIR / "splits/kcat_test.csv"

OOD_COLS = ["OOD_FoldHoldout", "OOD_SuperfamilyHoldout"]

# unique_id 两种形态：AF-<acc>-F1-model_v6 或带配体哈希后缀
# AF-<acc>-F1-model_v6-<hash8>（09_fix_unique_id_uniqueness.py 引入）
ACC_RE = re.compile(r"^AF-(.+?)-F1-model_v6(?:-[0-9a-f]{8})?$")


def extract_acc(unique_id: pd.Series) -> pd.Series:
    """从 unique_id 提取 UniProt accession（兼容配体哈希后缀）。"""
    return unique_id.str.extract(ACC_RE, expand=False)


def build_acc_list(out_path):
    """从 splits 生成 accession 清单（列：uniprot, splits；每 accession 一行）。"""
    rows = []
    for split in ["train", "val", "test"]:
        df = pd.read_csv(TASK_DIR / f"splits/kcat_{split}.csv",
                         dtype=str, keep_default_na=False)
        accs = extract_acc(df["unique_id"])
        for acc in sorted(accs.unique()):
            rows.append({"uniprot": acc, "splits": split})
    out = pd.DataFrame(rows)
    out.to_csv(out_path, index=False, lineterminator="\n")
    print(f"acc list: {len(out)} rows -> {out_path.name}")


# 0. 生成 accession 清单
build_acc_list(ACC_LIST)

# 1. 标注（train/val/test 对称，--offline 仅用 TED 缓存）
subprocess.run([sys.executable, str(ANNOTATE), "--chains", str(ACC_LIST),
                "--out", str(LABELS), "--offline"], check=True)

# 2. holdout 重算（统一实现：train-only 参考、any 语义、无注释 -> False）
subprocess.run([sys.executable, str(HOLDOUT), "--labels", str(LABELS),
                "--split-col", "splits", "--out", str(HOLDOUT_OUT)], check=True)

# 3. 按 accession 映射回 test csv 逐行写回
flags = pd.read_csv(HOLDOUT_OUT, dtype=str, keep_default_na=False)
flag_by_acc = flags.groupby("uniprot")[OOD_COLS].first().to_dict("index")

test = pd.read_csv(TEST_CSV, dtype=str, keep_default_na=False)
mask_default = test["Default"] == "True"
accs = extract_acc(test["unique_id"])

for col in OOD_COLS:
    if col not in test.columns:
        test[col] = ""
    mapped = accs.map(lambda a: flag_by_acc.get(a, {}).get(col, "False"))
    test.loc[mask_default, col] = mapped[mask_default]
    test.loc[~mask_default, col] = ""  # 极端长度行不参与计算，置 NaN（项目统一惯例）

test.to_csv(TEST_CSV, index=False, lineterminator="\n")

d = test[mask_default]
print(f"updated {TEST_CSV} (rows={len(test)}, Default={mask_default.sum()})")
for col in OOD_COLS:
    print(f"  {col}=True: {(d[col] == 'True').sum()}")
unann_accs = set(flags.loc[flags["topos"].eq("") & flags["sfs"].eq(""), "uniprot"])
n_unann = accs[mask_default].isin(unann_accs).sum()
print(f"  unannotated (-> False): {n_unann} rows")
