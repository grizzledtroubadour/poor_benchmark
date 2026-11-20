#!/usr/bin/env python3
"""Finalize OOD columns: Default=False NaN masking and InD recomputation.

Since 2026-08-28, task pipelines write NaN directly when generating OOD columns
for Default=False rows, meaning not computed. Here --mode nan provides
idempotent fallback validation: normally a no-op, repairing only missed values.

Generalizes two historical analysis scripts:

1. --mode nan (formerly analysis/scripts/nondefault_ood_nan.py):
   On Default=False rows, set all OOD_* flags except OOD_ExtremeShort/Long
   to NaN (empty cells), meaning not computed rather than False.
   For intentionally all-Default=False special files, skip masking with
   --exempt-nondefault. No current task needs this: both PDA design sets
   now have Default=True throughout.
   Includes OOD_ flags plus seq_Redundancy_* / TM-score_*; continuous annotations
   such as idr_ratio are handled by their generation scripts, not this mode.

2. --mode ind (formerly analysis/scripts/add_ind_column.py):
   InD means all OOD_*, seq_Redundancy_*, and TM-score_* flags are False, treating
   NaN as False. Default=False extreme rows have ExtremeShort/Long=True, hence InD=False.
   Overwrite existing InD idempotently; rename legacy Pure_ID (fold/PDA) in place
   and recompute; otherwise append InD.

--mode both applies nan then ind. --dry-run reports changes without writing.

Read all columns as strings with keep_default_na=False; update target columns only.
Write to_csv(index=False, lineterminator="\n"), preserving all other cells bytewise.

Usage:
    python finalize_ood_columns.py --test-csv kcat_test.csv --mode both
    python finalize_ood_columns.py --test-csv cath_test.csv --mode both --dry-run
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

FLAG_PREFIXES = ("OOD_", "seq_Redundancy_", "TM-score_")
KEEP_NONDEFAULT = {"OOD_ExtremeShort", "OOD_ExtremeLong"}


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--test-csv", required=True, type=Path)
    ap.add_argument("--mode", required=True, choices=["nan", "ind", "both"])
    ap.add_argument("--default-col", default="Default")
    ap.add_argument("--exempt-nondefault", action="store_true",
                    help="Skip nan for intentionally all-Default=False special files; "
                         "not needed for current tasks")
    ap.add_argument("--dry-run", action="store_true", help="Report differences without writing")
    return ap.parse_args()


def flag_columns(df):
    return [c for c in df.columns if c.startswith(FLAG_PREFIXES)]


def mode_nan(df, default_col, dry_run):
    """Idempotent fallback: NaN-mask non-length OOD flags on Default=False rows.

    Includes OOD_* except OOD_ExtremeShort/Long, plus seq_Redundancy_* / TM-score_*.
    Excludes continuous annotations such as idr_ratio. Return (df, report)."""
    if default_col not in df.columns:
        print(f"  nan: no {default_col} column; skipping")
        return df, {}
    nd = df[default_col] == "False"
    n_nd = int(nd.sum())
    ood_cols = [c for c in df.columns if c.startswith(FLAG_PREFIXES)
                and c not in KEEP_NONDEFAULT]
    if n_nd == 0 or not ood_cols:
        print(f"  nan: no-op (Default=False rows={n_nd}, flag columns={len(ood_cols)})")
        return df, {"nondefault_rows": n_nd, "cleared": {}}
    before = {c: int((df.loc[nd, c] != "").sum()) for c in ood_cols}
    cleared = {c: n for c, n in before.items() if n}
    print(f"  nan: set {len(ood_cols)} columns x {n_nd} rows to NaN; "
          f"nonempty cells cleared: {cleared or 'none (already empty)'}")
    if not dry_run:
        for c in ood_cols:
            # Do not use astype(object): pandas 3 dtype=str reads string-dtype columns.
            # Changing their dtype would falsely fail the defensive equals check below.
            ser = df[c].copy()
            ser.loc[nd] = ""
            df[c] = ser
    return df, {"nondefault_rows": n_nd, "cleared": cleared}


def mode_ind(df, dry_run):
    """Recompute InD. Return (df, report)."""
    renamed = False
    recompute = "InD" in df.columns
    if not recompute and "Pure_ID" in df.columns:
        df = df.rename(columns={"Pure_ID": "InD"})
        renamed = True
    cols = flag_columns(df)
    if not cols:
        print("  ind: no OOD flag columns; skipping")
        return df, {}
    flags = df[cols].apply(lambda s: s.eq("True"))
    ind = ~flags.any(axis=1)
    new_vals = ind.map({True: "True", False: "False"})
    action = "overwrite" if recompute else ("rename Pure_ID in place" if renamed else "append")
    if recompute or renamed:
        diff = int((df["InD"] != new_vals).sum())
        print(f"  ind: {action}, recompute over {len(cols)} flag columns; "
              f"InD=True {int(ind.sum())}/{len(df)} ({ind.sum()/len(df)*100:.1f}%)；"
              f"{diff} rows differ from existing values")
    else:
        diff = None
        print(f"  ind: {action}, {len(cols)} flag columns; "
              f"InD=True {int(ind.sum())}/{len(df)} ({ind.sum()/len(df)*100:.1f}%)")
    if not dry_run:
        df["InD"] = new_vals
    return df, {"n_flag_cols": len(cols), "ind_true": int(ind.sum()),
                "renamed": renamed, "recomputed": recompute, "diff_rows": diff}


def main():
    args = parse_args()
    df = pd.read_csv(args.test_csv, dtype=str, keep_default_na=False)
    print(f"{args.test_csv}: {len(df)} rows, {len(df.columns)} columns, "
          f"mode={args.mode}{' (dry-run)' if args.dry_run else ''}")
    orig = df.copy()

    if args.mode in ("nan", "both"):
        if args.exempt_nondefault:
            print(f"  nan: --exempt-nondefault; skipping ({args.test_csv.name})")
        else:
            df, _ = mode_nan(df, args.default_col, args.dry_run)
    if args.mode in ("ind", "both"):
        df, _ = mode_ind(df, args.dry_run)

    if args.dry_run:
        changed_cols = [c for c in df.columns
                        if c in orig.columns and not df[c].equals(orig[c])]
        new_cols = [c for c in df.columns if c not in orig.columns]
        print(f"dry-run: changed columns={changed_cols}, added columns={new_cols}; not written.")
        return

    # Defensive check: change target columns only (nan: OOD_* + seq_Redundancy_* / TM-score_*;
    # ind: InD/Pure_ID only).
    touched = {"InD", "Pure_ID"} | {
        c for c in df.columns if c.startswith(FLAG_PREFIXES) and c not in KEEP_NONDEFAULT}
    for c in orig.columns:
        if c not in touched and c in df.columns:
            assert df[c].equals(orig[c]), f"Non-target column {c} changed!"
    df.to_csv(args.test_csv, index=False, lineterminator="\n")
    print(f"-> Written {args.test_csv}")


if __name__ == "__main__":
    sys.exit(main())
