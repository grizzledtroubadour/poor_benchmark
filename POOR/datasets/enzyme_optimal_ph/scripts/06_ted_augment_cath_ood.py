#!/usr/bin/env python3
"""Add CATH+TED OOD_FoldHoldout / OOD_SuperfamilyHoldout to enzyme_optimal_ph test.

Sample key: UniProt accession (unique_id). Use the skill in uniprot mode:
CATH union via reverse SIFTS mapping + full-length TED; symmetric training-only reference, any semantics,
unannotated -> False. Default=False extreme-length rows receive empty cells/NaN (not computed),
following project convention. finalize_ood_columns.py --mode nan is only a fallback check.

Use .skills/cath-ted-domain-annotation/scripts/cath_ted_holdout.py for all holdout decisions;
this is the sole project-wide implementation. This wrapper handles annotation, invocation, and split updates.
"""
import csv
import subprocess
import sys
import tempfile
from pathlib import Path

TASK_DIR = Path(__file__).resolve().parent.parent
ROOT = TASK_DIR.parent.parent
SKILL_DIR = ROOT / ".skills/cath-ted-domain-annotation/scripts"
ANNOTATE = SKILL_DIR / "cath_ted_annotate.py"
HOLDOUT = SKILL_DIR / "cath_ted_holdout.py"
ACC_LIST = TASK_DIR / "output/optimal_ph_ted_aug_acc_list.csv"
LABELS = TASK_DIR / "output/optimal_ph_acc_cath_ted_labels.csv"
TEST_CSV = TASK_DIR / "splits/optimal_ph_prediction_test.csv"


subprocess.run([sys.executable, str(ANNOTATE), "--chains", str(ACC_LIST),
                "--out", str(LABELS), "--offline"], check=True)

# Recompute holdouts with the canonical skill: union labels from rows whose splits contain train,
# any semantics, unannotated -> False; the membership column is named splits.
with tempfile.TemporaryDirectory(prefix="holdout_", dir=TASK_DIR / "output") as tmp:
    flagged = Path(tmp) / "test_flagged.csv"
    subprocess.run([sys.executable, str(HOLDOUT), "--labels", str(LABELS),
                    "--split-col", "splits", "--out", str(flagged)], check=True)
    with open(flagged) as f:
        flags = {r["uniprot"]: r for r in csv.DictReader(f)}

with open(TEST_CSV) as f:
    test = list(csv.DictReader(f))
fieldnames = list(test[0].keys())
for col in ["OOD_FoldHoldout", "OOD_SuperfamilyHoldout"]:
    if col not in fieldnames:
        fieldnames.append(col)

n_f = n_s = n_un = 0
for r in test:
    if r.get("Default") != "True":
        r["OOD_FoldHoldout"] = ""
        r["OOD_SuperfamilyHoldout"] = ""
        continue
    fl = flags.get(r["unique_id"])
    fh = fl is not None and fl["OOD_FoldHoldout"] == "True"
    sh = fl is not None and fl["OOD_SuperfamilyHoldout"] == "True"
    r["OOD_FoldHoldout"] = str(fh)
    r["OOD_SuperfamilyHoldout"] = str(sh)
    n_f += fh
    n_s += sh
    n_un += fl is None or (not fl["topos"] and not fl["sfs"])

with open(TEST_CSV, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(test)
n_default = sum(1 for r in test if r.get("Default") == "True")
print(f"updated {TEST_CSV} (rows={len(test)}, Default={n_default})")
print(f"  OOD_FoldHoldout=True:        {n_f}")
print(f"  OOD_SuperfamilyHoldout=True: {n_s}")
print(f"  unannotated (-> False):      {n_un}")
