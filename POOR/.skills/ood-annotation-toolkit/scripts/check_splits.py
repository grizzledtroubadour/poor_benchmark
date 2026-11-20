#!/usr/bin/env python3
"""Read-only split consistency checks (ood-annotation-toolkit).

Combines core checks A–G from analysis/scripts/format_consistency.py with
survey_columns.py column inventories for three explicitly supplied split files:

- Row counts and columns: train==val; test adds OOD annotations to core columns.
- Within-split unique_id uniqueness.
- **F. unique_id leakage checks**：train∩val、train∩test、val∩test
- Required-column completeness（unique_id / aa_seq(1,2) / struct_file(1,2) / label）
- Deprecated columns: hyphenated OOD-*, labels, pdb_file*, protein_*_seq, etc.
- Missing labels; train/val must lack OOD/Default columns, test must contain them.
- Default/InD counts; non-extreme OOD_* NaN violations on Default=False rows.
- OOD_Orphan ⊆ seq_Redundancy_30 subset relation when both columns exist.

Usage:
    python check_splits.py --train kcat_train.csv --val kcat_val.csv --test kcat_test.csv
    python check_splits.py --train ppi_train.csv --val ppi_val.csv --test ppi_test.csv \
        --seq-cols aa_seq1,aa_seq2 --struct-cols struct_file1,struct_file2
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

DEPRECATED_EXACT = {
    "labels", "file_type", "pdb_file", "pdb_file_A", "pdb_file_B",
    "file_type_A", "file_type_B", "protein_A_seq", "protein_B_seq",
    "is_time_cutoff", "OOD_NewFunction", "OOD_NewFunction_EC",
    "OOD_LongTailFunction", "OOD_LongTail_EC",
    "OOD_CATH_FoldHoldout", "OOD_CATH_SuperfamilyHoldout",
    "OOD_LongTailFunction_L3_le5", "OOD_LongTailFunction_L3_le10",
}
KEEP_NONDEFAULT = {"OOD_ExtremeShort", "OOD_ExtremeLong"}
# Valid test-only non-flag annotation columns (raw OOD_IDR metrics).
ANNOTATION_COLS = {"idr_ratio"}


def is_deprecated(col):
    return col in DEPRECATED_EXACT or col.startswith("OOD-")


def ood_like(col):
    return (col.startswith("OOD_") or col.startswith("seq_Redundancy_")
            or col.startswith("TM-score_") or col in ("Default", "InD"))


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--train", required=True, type=Path)
    ap.add_argument("--val", required=True, type=Path)
    ap.add_argument("--test", required=True, type=Path)
    ap.add_argument("--seq-cols", default="aa_seq",
                    help="Comma-separated sequence columns; ppi uses aa_seq1,aa_seq2.")
    ap.add_argument("--struct-cols", default="struct_file",
                    help="Comma-separated structure columns; ppi uses struct_file1,struct_file2.")
    ap.add_argument("--label-col", default="label")
    return ap.parse_args()


def main():
    args = parse_args()
    issues = []

    dfs = {}
    for name, p in (("train", args.train), ("val", args.val), ("test", args.test)):
        # Read labels as strings to avoid floating-point collapsing such as 1.10 -> 1.1.
        df = pd.read_csv(p, dtype=str, keep_default_na=False)
        dfs[name] = df
        print(f"{name}: {p}  rows={len(df)}  cols={len(df.columns)}")

    seq_cols = [c.strip() for c in args.seq_cols.split(",")]
    struct_cols = [c.strip() for c in args.struct_cols.split(",")]
    core = ["unique_id"] + seq_cols + struct_cols + [args.label_col]

    # A/B/C/D: per-file checks.
    for name, df in dfs.items():
        miss = [c for c in core if c not in df.columns]
        if miss:
            issues.append(f"{name}: missing core columns {miss}")
        if "unique_id" in df.columns and not df["unique_id"].is_unique:
            n_dup = int(df["unique_id"].duplicated().sum())
            issues.append(f"{name}: duplicate within-split unique_id: {n_dup} rows")
        dep = [c for c in df.columns if is_deprecated(c)]
        if dep:
            issues.append(f"{name}: deprecated columns {dep}")
        if args.label_col in df.columns:
            n_null = int((df[args.label_col] == "").sum())
            if n_null:
                issues.append(f"{name}: missing labels: {n_null} rows")
        has_ood = [c for c in df.columns if ood_like(c)]
        if name in ("train", "val") and has_ood:
            issues.append(f"{name}: unexpected OOD/Default columns {has_ood}")
        if name == "test" and ("Default" not in df.columns or len(has_ood) < 2):
            issues.append(f"test: expected Default + OOD columns; found {has_ood}")

    # F. unique_id leakage
    idsets = {n: set(df["unique_id"]) for n, df in dfs.items() if "unique_id" in df.columns}
    if len(idsets) == 3:
        for a, b in (("train", "val"), ("train", "test"), ("val", "test")):
            inter = idsets[a] & idsets[b]
            status = f"{len(inter)} items" if inter else "none"
            print(f"Leakage check {a}∩{b}: {status}")
            if inter:
                issues.append(f"unique_id leakage {a}∩{b}: {len(inter)} items, "
                              f"e.g. {sorted(inter)[:5]}")

    # G. Column-set comparison
    col_train, col_val, col_test = (list(dfs[n].columns) for n in ("train", "val", "test"))
    if set(col_train) != set(col_val):
        issues.append(f"train/val column sets differ: "
                      f"train-only={sorted(set(col_train)-set(col_val))}, "
                      f"val-only={sorted(set(col_val)-set(col_train))}")
    extra_test = [c for c in col_test if c not in col_train]
    non_ood_extra = [c for c in extra_test
                     if not ood_like(c) and c not in ANNOTATION_COLS
                     and not c.startswith("idr_ratio")]
    missing_test = [c for c in col_train if c not in col_test]
    if non_ood_extra:
        issues.append(f"test has non-OOD columns absent from train: {non_ood_extra}")
    if missing_test:
        issues.append(f"test lacks training columns: {missing_test}")

    # Column inventory
    print("\n=== Column inventory ===")
    for n in ("train", "val", "test"):
        print(f"{n} ({len(dfs[n].columns)} columns): {list(dfs[n].columns)}")

    # Default/InD statistics and NaN convention.
    test = dfs["test"]
    if "Default" in test.columns:
        n_false = int((test["Default"] == "False").sum())
        print(f"\nDefault: True={int((test['Default'] == 'True').sum())}, "
              f"False={n_false}")
        if n_false:
            viol = {}
            for c in test.columns:
                if c.startswith("OOD_") and c not in KEEP_NONDEFAULT:
                    n = int((test.loc[test["Default"] == "False", c] != "").sum())
                    if n:
                        viol[c] = n
            if viol:
                issues.append(f"Default=False non-extreme OOD_* must be NaN; violations: {viol}")
            else:
                print("NaN convention satisfied: Default=False non-extreme OOD_* cells are all empty.")
    if "InD" in test.columns:
        print(f"InD: True={int((test['InD'] == 'True').sum())} / {len(test)}")
    else:
        print("InD: column absent; compute using finalize_ood_columns.py --mode ind.")

    # Orphan ⊆ seq_Redundancy_30
    if "OOD_Orphan" in test.columns and "seq_Redundancy_30" in test.columns:
        bad = int(((test["OOD_Orphan"] == "True")
                   & (test["seq_Redundancy_30"] != "True")).sum())
        print(f"Orphan ⊆ seq_Redundancy_30: {'satisfied' if bad == 0 else f'violations: {bad} rows'}")
        if bad:
            issues.append(f"OOD_Orphan is not a subset of seq_Redundancy_30（{bad} rows）")

    print("\n=== Check results ===")
    if issues:
        for i in issues:
            print(f"  [FAIL] {i}")
        print(f"Total:{len(issues)} issues")
        return 1
    print("All checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
