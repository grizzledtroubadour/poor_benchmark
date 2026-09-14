#!/usr/bin/env python3
"""为 func_prediction 四个子任务（ec/go_bp/go_cc/go_mf）test 集追加
OOD_FoldHoldout / OOD_SuperfamilyHoldout（CATH+TED 口径）。

流程（每个子任务，train/val/test 对称增强、独立 train-only 参考集）：
  1. build_chain_list：从 splits/{sub}_{train,val,test}.csv 生成链清单
     output/func_{sub}_ted_aug_chain_list.csv（unique_id 大写 PDB 转小写键）
  2. 调 .skills/cath-ted-domain-annotation/scripts/cath_ted_annotate.py 生成链级标注
  3. 调同 skill 的 cath_ted_holdout.py 做 holdout 判定（train-only 参考集、any 语义、无注释 -> False）
  4. 本脚本只负责把判定结果写回 splits

Default=False（极端长度）行此处直接写 NaN（空单元格，"未计算"语义，
全项目统一惯例）；finalize_ood_columns.py --mode nan 仅作兜底校验
（见 DATA_PROCESS.md 阶段 4.1）。
"""
import subprocess
import sys
from pathlib import Path

import pandas as pd

TASK_DIR = Path(__file__).resolve().parent.parent
ROOT = TASK_DIR.parent.parent
SKILL_DIR = ROOT / ".skills/cath-ted-domain-annotation/scripts"
ANNOTATE = SKILL_DIR / "cath_ted_annotate.py"
HOLDOUT = SKILL_DIR / "cath_ted_holdout.py"

SUBS = ["ec", "go_bp", "go_cc", "go_mf"]


def uid_key(uid):
    """test csv 的 unique_id 为大写 PDB（如 8OZ1_A），labels 键为小写 pdb_chain"""
    pdb, chain = uid.rsplit("_", 1)
    return f"{pdb.lower()}_{chain}"


def build_chain_list(sub, out_path):
    """从 splits 生成 CATH+TED 标注链清单（补缺口：原为历史一次性生成、脚本已删除）。

    输出列：unique_id（小写 pdb_chain 键）、pdb、chain、splits（train/val/test）。
    """
    rows = []
    for split in ["train", "val", "test"]:
        df = pd.read_csv(TASK_DIR / f"splits/{sub}_{split}.csv",
                         dtype=str, keep_default_na=False)
        for uid in df["unique_id"]:
            pdb, chain = uid.rsplit("_", 1)
            rows.append({"unique_id": f"{pdb.lower()}_{chain}",
                         "pdb": pdb.lower(), "chain": chain, "splits": split})
    out = pd.DataFrame(rows)
    out.to_csv(out_path, index=False, lineterminator="\n")
    print(f"  chain list: {len(out)} rows -> {out_path.name}")


for sub in SUBS:
    chain_list = TASK_DIR / f"output/func_{sub}_ted_aug_chain_list.csv"
    labels_f = TASK_DIR / f"output/func_{sub}_chain_cath_ted_labels.csv"
    domains_f = TASK_DIR / f"output/func_{sub}_chain_domains.csv"
    flagged_f = TASK_DIR / f"output/func_{sub}_test_holdout_flags.csv"
    test_csv = TASK_DIR / f"splits/{sub}_test.csv"
    print(f"===== {sub} =====")

    build_chain_list(sub, chain_list)

    subprocess.run([sys.executable, str(ANNOTATE), "--chains", str(chain_list),
                    "--out", str(labels_f), "--domains-out", str(domains_f),
                    "--offline"], check=True)

    # holdout 判定由 skill 脚本统一实现（train-only 参考集、any 语义、无注释 -> False）
    subprocess.run([sys.executable, str(HOLDOUT), "--labels", str(labels_f),
                    "--split-col", "splits", "--out", str(flagged_f)], check=True)

    flags = pd.read_csv(flagged_f, dtype=str, keep_default_na=False)
    flag_map = {r["unique_id"]: r for r in flags.to_dict("records")}

    test = pd.read_csv(test_csv, dtype=str, keep_default_na=False)
    for col in ["OOD_FoldHoldout", "OOD_SuperfamilyHoldout"]:
        if col not in test.columns:
            test[col] = ""

    n_f = n_s = n_un = n_miss = 0
    for i, r in test.iterrows():
        if r.get("Default") != "True":
            # 极端长度行不参与本维度计算，直接写 NaN（空）
            test.at[i, "OOD_FoldHoldout"] = ""
            test.at[i, "OOD_SuperfamilyHoldout"] = ""
            continue
        lab = flag_map.get(uid_key(r["unique_id"]))
        if lab is None:
            n_miss += 1
            fh = sh = "False"
        else:
            fh, sh = lab["OOD_FoldHoldout"], lab["OOD_SuperfamilyHoldout"]
            n_un += not (lab.get("topos") or lab.get("sfs"))
        test.at[i, "OOD_FoldHoldout"] = fh
        test.at[i, "OOD_SuperfamilyHoldout"] = sh
        n_f += fh == "True"
        n_s += sh == "True"

    test.to_csv(test_csv, index=False, lineterminator="\n")
    n_default = int((test["Default"] == "True").sum())
    print(f"updated {test_csv.name} (rows={len(test)}, Default={n_default})")
    print(f"  OOD_FoldHoldout=True: {n_f}, OOD_SuperfamilyHoldout=True: {n_s}, unannotated: {n_un}")
    if n_miss:
        print(f"  WARNING: {n_miss} Default=True rows missing from holdout flags (-> False)")
