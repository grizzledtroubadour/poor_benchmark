#!/usr/bin/env python3
"""Add CATH+TED OOD_FoldHoldout / OOD_SuperfamilyHoldout to ligand_binding_site test.

Chain list: output/lbs_ted_aug_chain_list.csv, historically generated once by a deleted script;
the first two unique_id components identify pdb/chain.
Annotation unit: chain; different ligand samples of one chain share annotations.
CATH priority + TED union gap filling, 50% domain overlap, symmetric training-only reference, any semantics;
unannotated -> False; Default=False extreme-length rows -> NaN.

Holdout decisions use .skills/cath-ted-domain-annotation/scripts/cath_ted_holdout.py,
the sole canonical implementation of reference scope, any semantics, and unannotated -> False.
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

# 1. Chain-level CATH+TED annotation with symmetric train/val/test augmentation.
subprocess.run([sys.executable, str(ANNOTATE), "--chains", str(CHAIN_LIST),
                "--out", str(LABELS), "--domains-out", str(DOMAINS), "--offline"],
               check=True)

# 2. Canonical holdout implementation: training-only reference, any semantics, unannotated -> False.
subprocess.run([sys.executable, str(HOLDOUT), "--labels", str(LABELS),
                "--split-col", "splits", "--out", str(FLAGS)],
               check=True)

# 3. Update only two split columns; Default=False extreme-length rows receive NaN.
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

# Count unannotated chains (False) from holdout outputs as supporting evidence.
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
