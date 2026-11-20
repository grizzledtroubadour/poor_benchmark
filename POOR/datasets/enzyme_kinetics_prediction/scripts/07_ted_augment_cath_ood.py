#!/usr/bin/env python3
"""Add CATH+TED OOD_FoldHoldout / OOD_SuperfamilyHoldout flags to the kcat test set.

Samples are keyed by UniProt accession: unique_id = AF-<acc>-F1-model_v6; after the step-09
fix, IDs may include a ligand hash suffix AF-<acc>-F1-model_v6-<hash8>. extract_acc supports both.
Use the skill in uniprot mode: union reverse-mapped SIFTS CATH labels and full-length TED labels.
Symmetric training-only reference; any semantics; unannotated -> False; Default=False -> NaN.

Workflow:
  0. build_acc_list: generate an accession list from splits/kcat_{train,val,test}.csv.
     output/kcat_ted_aug_acc_list.csv replaces a historical one-off output whose script was deleted.
  1. cath_ted_annotate.py: accession-level CATH+TED labels (--offline uses cached results only).
  2. cath_ted_holdout.py: canonical holdout computation (union of training rows; any semantics).
Both annotation steps use canonical .skills/cath-ted-domain-annotation scripts.
This wrapper maps holdouts by accession back to every splits/kcat_test.csv row
(multiple substrates per accession), setting Default=False rows to NaN by project convention.
"""
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd

TASK_DIR = Path(__file__).resolve().parent.parent
ROOT = TASK_DIR.parent.parent
SKILL_DIR = ROOT / ".skills/cath-ted-domain-annotation/scripts"
ANNOTATE = SKILL_DIR / "cath_ted_annotate.py"
HOLDOUT = SKILL_DIR / "cath_ted_holdout.py"
ACC_LIST = TASK_DIR / "output/kcat_ted_aug_acc_list.csv"
LABELS = TASK_DIR / "output/kcat_acc_cath_ted_labels.csv"
HOLDOUT_OUT = TASK_DIR / "output/kcat_cath_ted_holdout.csv"
TEST_CSV = TASK_DIR / "splits/kcat_test.csv"

OOD_COLS = ["OOD_FoldHoldout", "OOD_SuperfamilyHoldout"]

# unique_id formats: AF-<acc>-F1-model_v6 or the ligand-hash-suffixed form
# AF-<acc>-F1-model_v6-<hash8> introduced by 09_fix_unique_id_uniqueness.py.
ACC_RE = re.compile(r"^AF-(.+?)-F1-model_v6(?:-[0-9a-f]{8})?$")


def extract_acc(unique_id: pd.Series) -> pd.Series:
    """Extract UniProt accession from unique_id, supporting ligand-hash suffixes."""
    return unique_id.str.extract(ACC_RE, expand=False)


def build_acc_list(out_path):
    """Build one accession-list row per protein from splits (columns: uniprot, splits)."""
    rows = []
    for split in ["train", "val", "test"]:
        df = pd.read_csv(TASK_DIR / f"splits/kcat_{split}.csv",
                         dtype=str, keep_default_na=False)
        accs = extract_acc(df["unique_id"])
        for acc in sorted(accs.unique()):
            rows.append({"uniprot": acc, "splits": split})
    out = pd.DataFrame(rows)
    out.to_csv(out_path, index=False, lineterminator="\n")
    print(f"acc list: {len(out)} rows -> {out_path.name}")


# 0. Generate the accession list.
build_acc_list(ACC_LIST)

# 1. Symmetric train/val/test annotation (--offline uses TED cache only).
subprocess.run([sys.executable, str(ANNOTATE), "--chains", str(ACC_LIST),
                "--out", str(LABELS), "--offline"], check=True)

# 2. Canonical holdouts: training-only reference, any semantics, unannotated -> False.
subprocess.run([sys.executable, str(HOLDOUT), "--labels", str(LABELS),
                "--split-col", "splits", "--out", str(HOLDOUT_OUT)], check=True)

# 3. Map results by accession back to each test CSV row.
flags = pd.read_csv(HOLDOUT_OUT, dtype=str, keep_default_na=False)
flag_by_acc = flags.groupby("uniprot")[OOD_COLS].first().to_dict("index")

test = pd.read_csv(TEST_CSV, dtype=str, keep_default_na=False)
mask_default = test["Default"] == "True"
accs = extract_acc(test["unique_id"])

for col in OOD_COLS:
    if col not in test.columns:
        test[col] = ""
    mapped = accs.map(lambda a: flag_by_acc.get(a, {}).get(col, "False"))
    test.loc[mask_default, col] = mapped[mask_default]
    test.loc[~mask_default, col] = ""  # Extreme-length rows are not evaluated: NaN by project convention.

test.to_csv(TEST_CSV, index=False, lineterminator="\n")

d = test[mask_default]
print(f"updated {TEST_CSV} (rows={len(test)}, Default={mask_default.sum()})")
for col in OOD_COLS:
    print(f"  {col}=True: {(d[col] == 'True').sum()}")
unann_accs = set(flags.loc[flags["topos"].eq("") & flags["sfs"].eq(""), "uniprot"])
n_unann = accs[mask_default].isin(unann_accs).sum()
print(f"  unannotated (-> False): {n_unann} rows")
