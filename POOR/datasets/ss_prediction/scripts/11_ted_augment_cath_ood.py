#!/usr/bin/env python3
"""为 ss_prediction test 集追加 OOD_FoldHoldout / OOD_SuperfamilyHoldout（CATH+TED 口径）。

链清单 output/ss_ted_aug_chain_list.csv 由本脚本生成（步骤 0）：取
splits/ssp_{train,val,test}.csv 的 unique_id（PDB_chain，PDB 大写），
拆为小写 pdb + chain 并附 splits 列（40,066 条，即全库减去划分时剔除的
4O9X_A）。标注与 holdout 计算均走 .skills/cath-ted-domain-annotation：
- cath_ted_annotate.py：CATH 优先 + TED 并集补缺、50% 域重叠过滤
- cath_ted_holdout.py：train-only 对称参考、any 语义、无注释 -> False
写回时 Default=False（极端长度）行两列置 NaN（ood-annotation-toolkit 惯例，
不参与计算）。
"""
import csv
import subprocess
import sys
from pathlib import Path

TASK_DIR = Path(__file__).resolve().parent.parent
ROOT = TASK_DIR.parent.parent
SKILL_DIR = ROOT / ".skills/cath-ted-domain-annotation/scripts"
ANNOTATE = SKILL_DIR / "cath_ted_annotate.py"
HOLDOUT = SKILL_DIR / "cath_ted_holdout.py"
CHAIN_LIST = TASK_DIR / "output/ss_ted_aug_chain_list.csv"
LABELS = TASK_DIR / "output/ss_chain_cath_ted_labels.csv"
DOMAINS = TASK_DIR / "output/ss_chain_domains.csv"
HOLDOUT_CSV = TASK_DIR / "output/ss_chain_cath_ted_holdout.csv"
TEST_CSV = TASK_DIR / "splits/ssp_test.csv"


def build_chain_list():
    """从 splits 生成链清单（unique_id 为小写 pdb_chain，附 splits 列）。"""
    rows = []
    for split in ["train", "val", "test"]:
        with open(TASK_DIR / f"splits/ssp_{split}.csv") as f:
            for r in csv.DictReader(f):
                pdb, chain = r["unique_id"].rsplit("_", 1)
                rows.append({
                    "unique_id": f"{pdb.lower()}_{chain}",
                    "pdb": pdb.lower(),
                    "chain": chain,
                    "splits": split,
                })
    CHAIN_LIST.parent.mkdir(parents=True, exist_ok=True)
    with open(CHAIN_LIST, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["unique_id", "pdb", "chain", "splits"])
        w.writeheader()
        w.writerows(rows)
    print(f"chain list: {CHAIN_LIST} ({len(rows)} chains)")


# 0) 链清单生成
build_chain_list()

# 1) 链级 CATH+TED 标注（离线，全部命中项目级 TED 缓存）
subprocess.run([sys.executable, str(ANNOTATE), "--chains", str(CHAIN_LIST),
                "--out", str(LABELS), "--domains-out", str(DOMAINS), "--offline"],
               check=True)

# 2) holdout 重算（skill 统一口径；labels 的划分列为 splits）
subprocess.run([sys.executable, str(HOLDOUT), "--labels", str(LABELS),
                "--split-col", "splits", "--out", str(HOLDOUT_CSV)],
               check=True)

# 3) 写回 ssp_test.csv


def uid_key(uid):
    """test csv 的 unique_id 为大写 PDB（如 7E3T_A），labels 键为小写 pdb_chain"""
    pdb, chain = uid.rsplit("_", 1)
    return f"{pdb.lower()}_{chain}"


flags = {r["unique_id"]: r for r in csv.DictReader(open(HOLDOUT_CSV))}

with open(TEST_CSV) as f:
    test = list(csv.DictReader(f))
fieldnames = list(test[0].keys())
for col in ["OOD_FoldHoldout", "OOD_SuperfamilyHoldout"]:
    if col not in fieldnames:
        fieldnames.append(col)

n_f = n_s = n_un = 0
for r in test:
    if r.get("Default") != "True":
        # 极端长度行不参与计算，置 NaN（空单元格）
        r["OOD_FoldHoldout"] = ""
        r["OOD_SuperfamilyHoldout"] = ""
        continue
    flag = flags.get(uid_key(r["unique_id"]))
    fh = flag["OOD_FoldHoldout"] if flag else "False"
    sh = flag["OOD_SuperfamilyHoldout"] if flag else "False"
    r["OOD_FoldHoldout"] = fh
    r["OOD_SuperfamilyHoldout"] = sh
    n_f += fh == "True"
    n_s += sh == "True"
    n_un += flag is None

with open(TEST_CSV, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(test)
n_default = sum(1 for r in test if r.get("Default") == "True")
print(f"updated {TEST_CSV} (rows={len(test)}, Default={n_default})")
print(f"  OOD_FoldHoldout=True:        {n_f}")
print(f"  OOD_SuperfamilyHoldout=True: {n_s}")
print(f"  test rows missing in holdout output (-> False): {n_un}")
