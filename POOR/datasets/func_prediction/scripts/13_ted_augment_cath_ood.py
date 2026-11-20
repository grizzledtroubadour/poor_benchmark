#!/usr/bin/env python3
"""Add test-set holdout flags to four func_prediction tasks (ec/go_bp/go_cc/go_mf):
OOD_FoldHoldout / OOD_SuperfamilyHoldout using CATH+TED.

For each task, augment train/val/test symmetrically and use its own training-only reference:
  1. build_chain_list: generate a chain list from splits/{sub}_{train,val,test}.csv.
     Write output/func_{sub}_ted_aug_chain_list.csv, lowercasing PDB IDs in unique_id.
  2. Call .skills/cath-ted-domain-annotation/scripts/cath_ted_annotate.py for chain labels.
  3. Call cath_ted_holdout.py (training-only reference, any semantics, unannotated -> False).
  4. Write the resulting flags back to splits.

Default=False extreme-length rows receive empty cells (NaN: not computed), following
project convention. finalize_ood_columns.py --mode nan is only a fallback check
(see DATA_PROCESS.md Stage 4.1).
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
    """Test unique_id uses uppercase PDB (e.g. 8OZ1_A); label keys use lowercase pdb_chain."""
    pdb, chain = uid.rsplit("_", 1)
    return f"{pdb.lower()}_{chain}"


def build_chain_list(sub, out_path):
    """Generate CATH+TED chain lists from splits, replacing a deleted historical one-off script.

    Columns: unique_id (lowercase pdb_chain), pdb, chain, splits (train/val/test).
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

    # Use the canonical skill holdout script: training-only reference, any semantics, unannotated -> False.
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
            # Extreme-length rows are not evaluated on this dimension: write empty/NaN.
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
