#!/usr/bin/env python3
"""Add CATH+TED OOD_FoldHoldout / OOD_SuperfamilyHoldout to PPI test.

Samples are protein pairs; protein_A_id/protein_B_id are UniProt accessions:
- Annotate all accessions in uniprot mode: reverse-map CATH through SIFTS, union with full-length TED.
- Union both proteins' labels per sample; any semantics; symmetric training-only reference.
- Reuse .skills/cath-ted-domain-annotation/scripts/cath_ted_holdout.py:
      first write pair labels (unique_id, topos, sfs, split), then compute canonical holdouts.
- Accession list: output/ppi_ted_aug_acc_list.csv, historically generated once; original script deleted.
"""
import csv
import subprocess
import sys
from pathlib import Path

import pandas as pd

TASK_DIR = Path(__file__).resolve().parent.parent
ROOT = TASK_DIR.parent.parent
SKILL_DIR = ROOT / ".skills/cath-ted-domain-annotation/scripts"
ANNOTATE = SKILL_DIR / "cath_ted_annotate.py"
HOLDOUT = SKILL_DIR / "cath_ted_holdout.py"
ACC_LIST = TASK_DIR / "output/ppi_ted_aug_acc_list.csv"
LABELS = TASK_DIR / "output/ppi_acc_cath_ted_labels.csv"
PAIR_LABELS = TASK_DIR / "output/ppi_pair_cath_ted_labels.csv"
PAIR_FLAGS = TASK_DIR / "output/ppi_pair_holdout_flags.csv"
TEST_CSV = TASK_DIR / "splits/ppi_test.csv"

HOLDOUT_COLS = ["OOD_FoldHoldout", "OOD_SuperfamilyHoldout"]


def toset(s):
    return {x for x in str(s).split(";") if x}


# 1. Annotate all accessions with reverse-mapped CATH + full-length TED, using offline cache.
subprocess.run([sys.executable, str(ANNOTATE), "--chains", str(ACC_LIST),
                "--out", str(LABELS), "--offline"], check=True)

labels = {r["uniprot"]: r for r in csv.DictReader(open(LABELS))}

# 2. Aggregate pair labels as the union of both accessions' labels.
frames = []
for split in ["train", "val", "test"]:
    df = pd.read_csv(TASK_DIR / f"splits/ppi_{split}.csv",
                     dtype=str, keep_default_na=False)
    df["split"] = split
    frames.append(df)
all_rows = pd.concat(frames, ignore_index=True)

pair_rows = []
for r in all_rows.itertuples(index=False):
    topos, sfs = set(), set()
    for acc in (r.protein_A_id, r.protein_B_id):
        lab = labels.get(acc)
        if lab:
            topos |= toset(lab["topos"])
            sfs |= toset(lab["sfs"])
    pair_rows.append({
        "unique_id": r.unique_id,
        "topos": ";".join(sorted(topos)),
        "sfs": ";".join(sorted(sfs)),
        "split": r.split,
    })
pd.DataFrame(pair_rows).to_csv(PAIR_LABELS, index=False, lineterminator="\n")
print(f"pair-level labels -> {PAIR_LABELS} (rows={len(pair_rows)})")

# 3. Canonical holdouts: training-only reference, any semantics, unannotated -> False.
subprocess.run([sys.executable, str(HOLDOUT), "--labels", str(PAIR_LABELS),
                "--split-col", "split", "--out", str(PAIR_FLAGS)], check=True)

# 4. Write test splits, setting Default=False rows to NaN by project convention.
flags = pd.read_csv(PAIR_FLAGS, dtype=str, keep_default_na=False).set_index("unique_id")
test = pd.read_csv(TEST_CSV, dtype=str, keep_default_na=False)
for col in HOLDOUT_COLS:
    if col not in test.columns:
        test[col] = ""

n_f = n_s = n_un = 0
for i, r in test.iterrows():
    if r["Default"] != "True":
        test.at[i, "OOD_FoldHoldout"] = ""
        test.at[i, "OOD_SuperfamilyHoldout"] = ""
        continue
    f = flags.loc[r["unique_id"]]
    test.at[i, "OOD_FoldHoldout"] = f["OOD_FoldHoldout"]
    test.at[i, "OOD_SuperfamilyHoldout"] = f["OOD_SuperfamilyHoldout"]
    n_f += f["OOD_FoldHoldout"] == "True"
    n_s += f["OOD_SuperfamilyHoldout"] == "True"
    n_un += not (f["topos"] or f["sfs"])

test.to_csv(TEST_CSV, index=False, lineterminator="\n")
n_default = int((test["Default"] == "True").sum())
print(f"updated {TEST_CSV} (rows={len(test)}, Default={n_default})")
print(f"  OOD_FoldHoldout=True:        {n_f}")
print(f"  OOD_SuperfamilyHoldout=True: {n_s}")
print(f"  unannotated (-> False):      {n_un}")
