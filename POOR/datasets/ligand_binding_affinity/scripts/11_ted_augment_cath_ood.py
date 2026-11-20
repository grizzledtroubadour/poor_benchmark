#!/usr/bin/env python3
"""Recompute ligand_binding_affinity test holdouts using unified CATH+TED annotations:
OOD_FoldHoldout / OOD_SuperfamilyHoldout。

Workflow:
1. Chain list: 11,663 PDBs × chains in the union of both SIFTS files and cath-domain-list.txt.
   SIFTS-only seeding misses 944 CATH-annotated chains (278 chains in 87 test PDBs).
2. Run .skills/cath-ted-domain-annotation/scripts/cath_ted_annotate.py (--offline).
3. Union chain-level labels by PDB, retaining the archived 16_cath_ec_ood.py aggregation unit.
   Compute holdouts with .skills/cath-ted-domain-annotation/scripts/cath_ted_holdout.py
   (training-only reference, any semantics, unannotated -> False; sole canonical implementation).
4. Update only the two holdout columns in splits/ligand_binding_affinity_test.csv.
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

# Canonical skill holdout implementation: training-only reference, any semantics, unannotated -> False.
subprocess.run([sys.executable, str(HOLDOUT), "--labels", str(PDB_LABELS),
                "--split-col", "splits", "--out", str(FLAGGED)],
               check=True)

# ---------- 4. Write back only the two holdout columns; retain all others. ----------
flagged = pd.read_csv(FLAGGED, dtype=str, keep_default_na=False).set_index("pdb")
test = pd.read_csv(TEST_CSV, dtype=str, keep_default_na=False)
pdb_key = test["unique_id"].str.lower()
# Test PDBs without any listed chains are unannotated and receive False.
test["OOD_FoldHoldout"] = pdb_key.map(flagged["OOD_FoldHoldout"]).fillna("False")
test["OOD_SuperfamilyHoldout"] = pdb_key.map(flagged["OOD_SuperfamilyHoldout"]).fillna("False")
assert test[["OOD_FoldHoldout", "OOD_SuperfamilyHoldout"]].isin(["True", "False"]).all().all()
test.to_csv(TEST_CSV, index=False, lineterminator="\n")

print(f"updated {TEST_CSV} ({len(test)} rows)")
print(f"  OOD_FoldHoldout=True:        {(test['OOD_FoldHoldout'] == 'True').sum()}")
print(f"  OOD_SuperfamilyHoldout=True: {(test['OOD_SuperfamilyHoldout'] == 'True').sum()}")
