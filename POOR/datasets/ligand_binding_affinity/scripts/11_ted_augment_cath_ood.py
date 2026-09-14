#!/usr/bin/env python3
"""用 CATH+TED 统一标注重算 ligand_binding_affinity test 的
OOD_FoldHoldout / OOD_SuperfamilyHoldout。

流程：
1. 链清单：11,663 个 PDB ×（SIFTS 两文件 ∪ cath-domain-list.txt 中出现的链）——
   仅 SIFTS 播种会漏 944 条 CATH 已标注链（含 87 个 test PDB 的 278 链）；
2. 调用 .skills/cath-ted-domain-annotation/scripts/cath_ted_annotate.py（--offline）；
3. 链级标签按 PDB 聚合取并集（沿用已归档 16_cath_ec_ood.py 的 PDB 级口径），
   holdout 判定调用 .skills/cath-ted-domain-annotation/scripts/cath_ted_holdout.py
   （train-only 参考、any 语义、无注释 -> False，口径以该脚本为唯一实现）；
4. 就地更新 splits/ligand_binding_affinity_test.csv 的两列（其余列不动）。
"""
import csv
import gzip
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

TASK_DIR = Path(__file__).resolve().parent.parent
ROOT = TASK_DIR.parent.parent
SKILL = ROOT / ".skills/cath-ted-domain-annotation/scripts/cath_ted_annotate.py"
HOLDOUT = ROOT / ".skills/cath-ted-domain-annotation/scripts/cath_ted_holdout.py"
OUT = TASK_DIR / "output"
CHAIN_LIST = OUT / "ligand_ted_aug_chain_list.csv"
LABELS = OUT / "ligand_chain_cath_ted_labels.csv"
DOMAINS = OUT / "ligand_chain_domains.csv"
PDB_LABELS = OUT / "ligand_pdb_cath_ted_labels.csv"
FLAGGED = OUT / "ligand_pdb_cath_ted_holdout.csv"
TEST_CSV = TASK_DIR / "splits/ligand_binding_affinity_test.csv"

# ---------- 1. chain list: SIFTS (both files) U cath-domain-list ----------
pdb_splits = defaultdict(set)
for split in ["train", "val", "test"]:
    with open(TASK_DIR / f"splits/ligand_binding_affinity_{split}.csv") as f:
        for r in csv.DictReader(f):
            pdb_splits[r["unique_id"].lower()].add(split)
pdbs = set(pdb_splits)

chains = defaultdict(set)  # pdb -> {chain}
for fn in ["pdb_chain_cath.tsv.gz", "sifts_chain_uniprot.tsv.gz"]:
    with gzip.open(ROOT / f"data/sifts/{fn}", "rt") as f:
        for line in f:
            if line.startswith("#") or line.startswith("PDB"):
                continue
            p = line.split("\t")
            if len(p) >= 2 and p[0] in pdbs:
                chains[p[0]].add(p[1])
with open(ROOT / "data/CATHv44/cath-classification-data/cath-domain-list.txt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        p = line.split()
        if len(p) >= 5 and p[0][:4].lower() in pdbs:
            chains[p[0][:4].lower()].add(p[0][4])

with open(CHAIN_LIST, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["pdb", "chain", "splits"])
    for pdb in sorted(chains):
        for ch in sorted(chains[pdb]):
            w.writerow([pdb, ch, ";".join(sorted(pdb_splits[pdb]))])
n_rows = sum(len(v) for v in chains.values())
print(f"chain list: {len(chains)} PDBs, {n_rows} chains")

# ---------- 2. annotate via skill ----------
subprocess.run([sys.executable, str(SKILL), "--chains", str(CHAIN_LIST),
                "--out", str(LABELS), "--domains-out", str(DOMAINS), "--offline"],
               check=True)

# ---------- 3. PDB-level union + holdout (via skill cath_ted_holdout.py) ----------
def toset(s):
    return {x for x in str(s).split(";") if x}


pdb_labels = defaultdict(lambda: [set(), set()])
pdb_split_map = {}
for r in csv.DictReader(open(LABELS)):
    pdb_labels[r["pdb"]][0] |= toset(r["topos"])
    pdb_labels[r["pdb"]][1] |= toset(r["sfs"])
    pdb_split_map[r["pdb"]] = r["splits"]

with open(PDB_LABELS, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["pdb", "splits", "topos", "sfs"])
    for pdb in sorted(pdb_labels):
        tt, ss = pdb_labels[pdb]
        w.writerow([pdb, pdb_split_map[pdb], ";".join(sorted(tt)), ";".join(sorted(ss))])

# holdout 判定统一由 skill 脚本实现（train-only 参考、any 语义、无注释 -> False）
subprocess.run([sys.executable, str(HOLDOUT), "--labels", str(PDB_LABELS),
                "--split-col", "splits", "--out", str(FLAGGED)],
               check=True)

# ---------- 4. 写回 splits（仅两列，其余列不动） ----------
flagged = pd.read_csv(FLAGGED, dtype=str, keep_default_na=False).set_index("pdb")
test = pd.read_csv(TEST_CSV, dtype=str, keep_default_na=False)
pdb_key = test["unique_id"].str.lower()
# 链清单中无任何链记录的 test PDB（无注释）按口径 -> False
test["OOD_FoldHoldout"] = pdb_key.map(flagged["OOD_FoldHoldout"]).fillna("False")
test["OOD_SuperfamilyHoldout"] = pdb_key.map(flagged["OOD_SuperfamilyHoldout"]).fillna("False")
assert test[["OOD_FoldHoldout", "OOD_SuperfamilyHoldout"]].isin(["True", "False"]).all().all()
test.to_csv(TEST_CSV, index=False, lineterminator="\n")

print(f"updated {TEST_CSV} ({len(test)} rows)")
print(f"  OOD_FoldHoldout=True:        {(test['OOD_FoldHoldout'] == 'True').sum()}")
print(f"  OOD_SuperfamilyHoldout=True: {(test['OOD_SuperfamilyHoldout'] == 'True').sum()}")
