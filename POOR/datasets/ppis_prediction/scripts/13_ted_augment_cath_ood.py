#!/usr/bin/env python3
"""用 CATH+TED 统一标注重算 PPIS test 的 OOD_FoldHoldout / OOD_SuperfamilyHoldout。

流程：
1. seq_hash -> (pdb, chain)（复用已归档 add_cath_holdout_ood.py 的映射逻辑），生成链清单；
2. 调用 .skills/cath-ted-domain-annotation/scripts/cath_ted_annotate.py（--offline，共享缓存）；
3. 调用 .skills/cath-ted-domain-annotation/scripts/cath_ted_holdout.py 计算 Holdout
   （train-only 参考、any 语义、无注释 -> False，口径以该脚本为唯一实现）；
4. 就地更新 splits/ppis_test.csv 的两列（其余列不动）；
   Default=False（极端长度）行直接写 NaN（空单元格，"未计算"语义，
   全项目统一惯例）；finalize_ood_columns.py --mode nan 仅作兜底校验。
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
OUT = TASK_DIR / "output"
CHAIN_LIST = OUT / "ppis_ted_aug_chain_list.csv"
LABELS = OUT / "ppis_chain_cath_ted_labels.csv"
DOMAINS = OUT / "ppis_chain_domains.csv"
HOLDOUT_OUT = OUT / "ppis_ted_aug_holdout.csv"
TEST_CSV = TASK_DIR / "splits/ppis_test.csv"

# ---------- 1. chain list ----------
pair_meta = {}
with open(OUT / "dips_plus_metadata.csv") as f:
    for r in csv.DictReader(f):
        pair_meta[r["unique_id"]] = r
meta_by_hash = {}
with open(OUT / "chain_centric_dedup_metadata.csv") as f:
    for r in csv.DictReader(f):
        meta_by_hash[r["seq_hash"]] = r


def src_chain(h):
    m = meta_by_hash.get(h, {})
    uid = m.get("structure_source_unique_id", "")
    role = m.get("structure_source_role", "").split("->")[0].strip()
    if not uid or not role:
        return None
    pm = pair_meta.get(uid)
    if not pm:
        return None
    pdb = uid.split(".")[0].lower()
    ch = pm.get("receptor_chain" if role == "receptor" else "ligand_chain", "")
    return (pdb, ch) if ch else None


rows = {}
for split in ["train", "val", "test"]:
    with open(TASK_DIR / f"splits/ppis_{split}.csv") as f:
        for r in csv.DictReader(f):
            h = r["unique_id"]
            c = src_chain(h)
            if c is None:
                continue
            if h in rows:
                rows[h]["split"] += f";{split}"
                continue
            rows[h] = {"seq_hash": h, "pdb": c[0], "chain": c[1], "split": split}
with open(CHAIN_LIST, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["seq_hash", "pdb", "chain", "split"])
    w.writeheader()
    w.writerows(rows.values())
print(f"chain list: {len(rows)} hashes")

# ---------- 2. annotate via skill ----------
subprocess.run([sys.executable, str(ANNOTATE), "--chains", str(CHAIN_LIST),
                "--out", str(LABELS), "--domains-out", str(DOMAINS), "--offline"],
               check=True)

# ---------- 3. holdout via skill ----------
subprocess.run([sys.executable, str(HOLDOUT), "--labels", str(LABELS),
                "--split-col", "split", "--out", str(HOLDOUT_OUT)],
               check=True)

# ---------- 4. write back to splits ----------
holdout = {r["seq_hash"]: r for r in csv.DictReader(open(HOLDOUT_OUT))}

with open(TEST_CSV) as f:
    test = list(csv.DictReader(f))
fieldnames = list(test[0].keys())

n_f = n_s = n_un = 0
for r in test:
    if r.get("Default") != "True":
        r["OOD_FoldHoldout"] = ""
        r["OOD_SuperfamilyHoldout"] = ""
        continue
    h = holdout.get(r["unique_id"])
    fh = h["OOD_FoldHoldout"] if h else "False"
    sh = h["OOD_SuperfamilyHoldout"] if h else "False"
    r["OOD_FoldHoldout"] = fh
    r["OOD_SuperfamilyHoldout"] = sh
    n_f += fh == "True"
    n_s += sh == "True"
    n_un += h is None

with open(TEST_CSV, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(test)

n_default = sum(1 for r in test if r.get("Default") == "True")
print(f"updated {TEST_CSV}")
print(f"  Default rows: {n_default}")
print(f"  OOD_FoldHoldout=True:        {n_f}")
print(f"  OOD_SuperfamilyHoldout=True: {n_s}")
print(f"  no holdout record (-> False): {n_un}")
