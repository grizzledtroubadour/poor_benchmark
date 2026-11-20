#!/usr/bin/env python3
"""Add CATH+TED OOD_FoldHoldout / OOD_SuperfamilyHoldout to the SSP test set.

Step 0 generates output/ss_ted_aug_chain_list.csv from
unique_id in splits/ssp_{train,val,test}.csv (uppercase PDB_chain),
normalizing PDB to lowercase and adding splits. It has 40,066 entries, excluding
4O9X_A removed during splitting. Annotation and holdouts use .skills/cath-ted-domain-annotation:
- cath_ted_annotate.py: prioritize CATH, supplement with TED, require 50% domain overlap.
- cath_ted_holdout.py: symmetric training-only reference, any semantics, unannotated -> False.
Set both columns to NaN on Default=False extreme-length rows, per toolkit convention;
these rows are not evaluated.
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
    """Generate a chain list from splits: lowercase pdb_chain unique_id and splits column."""
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


# 0) Generate chain list.
build_chain_list()

# 1) Chain-level CATH+TED annotation offline; all entries hit the project TED cache.
subprocess.run([sys.executable, str(ANNOTATE), "--chains", str(CHAIN_LIST),
                "--out", str(LABELS), "--domains-out", str(DOMAINS), "--offline"],
               check=True)

# 2) Canonical holdout recomputation; labels use the splits membership column.
subprocess.run([sys.executable, str(HOLDOUT), "--labels", str(LABELS),
                "--split-col", "splits", "--out", str(HOLDOUT_CSV)],
               check=True)

# 3) Write back ssp_test.csv.


def uid_key(uid):
    """Test IDs use uppercase PDB (e.g. 7E3T_A); label keys use lowercase pdb_chain."""
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
        # Extreme-length rows are not evaluated; write NaN (empty cells).
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
