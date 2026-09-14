#!/usr/bin/env python3
"""为 ligand_binding_site test 集追加 OOD_FoldHoldout / OOD_SuperfamilyHoldout（CATH+TED 口径）。

链清单为 output/lbs_ted_aug_chain_list.csv（历史一次性生成，生成脚本已删除；
unique_id 前两段为 pdb/chain）。
样本粒度：链级（同一链的不同配体样本共享同一标注）。
口径：CATH 优先 + TED 并集补缺、50% 域重叠过滤、train-only 对称参考、any 语义、
无注释 -> False、Default=False（极端长度）行 -> NaN。

Holdout 判定统一调用 `.skills/cath-ted-domain-annotation/scripts/cath_ted_holdout.py`
（参考集范围、any 语义、无注释 -> False 等口径以该脚本为唯一实现）。
"""
import subprocess
import sys
from pathlib import Path

import pandas as pd

TASK_DIR = Path(__file__).resolve().parent.parent
ROOT = TASK_DIR.parent.parent
ANNOTATE = ROOT / ".skills/cath-ted-domain-annotation/scripts/cath_ted_annotate.py"
HOLDOUT = ROOT / ".skills/cath-ted-domain-annotation/scripts/cath_ted_holdout.py"
CHAIN_LIST = TASK_DIR / "output/lbs_ted_aug_chain_list.csv"
LABELS = TASK_DIR / "output/lbs_chain_cath_ted_labels.csv"
DOMAINS = TASK_DIR / "output/lbs_chain_domains.csv"
FLAGS = TASK_DIR / "output/lbs_chain_holdout_flags.csv"
TEST_CSV = TASK_DIR / "splits/ligand_binding_site_test.csv"

FLAG_COLS = ["OOD_FoldHoldout", "OOD_SuperfamilyHoldout"]

# 1) 链级 CATH+TED 标注（train/val/test 对称增强）
subprocess.run([sys.executable, str(ANNOTATE), "--chains", str(CHAIN_LIST),
                "--out", str(LABELS), "--domains-out", str(DOMAINS), "--offline"],
               check=True)

# 2) Holdout 判定（skill 唯一实现：train-only 参考、any 语义、无注释 -> False）
subprocess.run([sys.executable, str(HOLDOUT), "--labels", str(LABELS),
                "--split-col", "splits", "--out", str(FLAGS)],
               check=True)

# 3) 写回 splits：仅覆写两列，Default=False（极端长度）行置 NaN
flags = pd.read_csv(FLAGS, dtype=str, keep_default_na=False)
flag_map = flags.set_index("unique_id")[FLAG_COLS].to_dict("index")

test = pd.read_csv(TEST_CSV, dtype=str, keep_default_na=False)
for col in FLAG_COLS:
    if col not in test.columns:
        test[col] = ""

n_f = n_s = n_un = n_miss = 0
for idx, r in test.iterrows():
    if r.get("Default") != "True":
        test.at[idx, "OOD_FoldHoldout"] = ""
        test.at[idx, "OOD_SuperfamilyHoldout"] = ""
        continue
    p = r["unique_id"].split("_")
    lab = flag_map.get(f"{p[0].lower()}_{p[1]}")
    if lab is None:
        n_miss += 1
        fh = sh = "False"
    else:
        fh, sh = lab["OOD_FoldHoldout"], lab["OOD_SuperfamilyHoldout"]
    test.at[idx, "OOD_FoldHoldout"] = fh
    test.at[idx, "OOD_SuperfamilyHoldout"] = sh
    n_f += fh == "True"
    n_s += sh == "True"

# 无注释链数（-> False）从 holdout 输出旁证统计
labels = pd.read_csv(LABELS, dtype=str, keep_default_na=False)
test_labels = labels[labels["splits"].str.split(";").apply(lambda xs: "test" in xs)]
n_un = int(((test_labels["topos"] == "") & (test_labels["sfs"] == "")).sum())

test.to_csv(TEST_CSV, index=False, lineterminator="\n")
n_default = int((test["Default"] == "True").sum())
print(f"updated {TEST_CSV} (rows={len(test)}, Default={n_default})")
print(f"  OOD_FoldHoldout=True:        {n_f}")
print(f"  OOD_SuperfamilyHoldout=True: {n_s}")
print(f"  unannotated (-> False):      {n_un}")
print(f"  rows w/o chain in list:      {n_miss}")
