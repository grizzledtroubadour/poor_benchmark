#!/usr/bin/env python3
"""Compute FoldHoldout / SuperfamilyHoldout OOD flags from CATH+TED label tables.

Semantics (POOD-Benchmark convention):
  - reference set = union of labels over `train` rows (train-only)
  - a test row is Holdout if ANY of its labels at that level is absent from the reference
  - unannotated rows -> False (never True)
  - train/val/test are taken from a split column; rows may carry several splits
    separated by ';' (a row in both train and test still contributes to the reference)

IMPORTANT: annotate train/val/test symmetrically (run cath_ted_annotate.py on all
splits together) before computing holdouts, otherwise TED-only test labels appear
spuriously "unseen".

Usage:
  python cath_ted_holdout.py --labels my_labels.csv --split-col split --out test_flagged.csv
  # or evaluate a test-only file against a separate reference file:
  python cath_ted_holdout.py --labels test_labels.csv --ref-labels trainval_labels.csv --out out.csv
"""
import argparse
import csv
from pathlib import Path


def toset(s):
    return {x for x in str(s).split(";") if x}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", type=Path, required=True,
                    help="CSV produced by cath_ted_annotate.py")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--split-col", default="split",
                    help="column with train/val/test membership (';' separated). "
                         "If absent, all rows are treated as test")
    ap.add_argument("--ref-labels", type=Path, default=None,
                    help="optional separate CSV whose rows define the reference sets")
    ap.add_argument("--topo-col", default="topos")
    ap.add_argument("--sf-col", default="sfs")
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.labels)))

    def has_split(r, name):
        return name in str(r.get(args.split_col, "")).split(";")

    if args.ref_labels:
        ref_rows = list(csv.DictReader(open(args.ref_labels)))
        test_rows = rows
    elif any(args.split_col in r for r in rows):
        ref_rows = [r for r in rows if has_split(r, "train")]
        test_rows = [r for r in rows if has_split(r, "test")]
    else:
        ref_rows, test_rows = [], rows

    ref_topo = set().union(*[toset(r[args.topo_col]) for r in ref_rows]) if ref_rows else set()
    ref_sf = set().union(*[toset(r[args.sf_col]) for r in ref_rows]) if ref_rows else set()
    print(f"reference: {len(ref_rows)} train rows, {len(ref_topo)} topos, {len(ref_sf)} sfs")
    if not ref_rows:
        print("WARNING: empty reference set — every annotated row will be Holdout=True")

    n_f = n_s = n_un = 0
    for r in test_rows:
        topos, sfs = toset(r[args.topo_col]), toset(r[args.sf_col])
        fh = bool(topos) and any(t not in ref_topo for t in topos)
        sh = bool(sfs) and any(s not in ref_sf for s in sfs)
        r["OOD_FoldHoldout"] = str(fh)
        r["OOD_SuperfamilyHoldout"] = str(sh)
        n_f += fh
        n_s += sh
        n_un += not (topos or sfs)

    fieldnames = list(rows[0].keys()) + ["OOD_FoldHoldout", "OOD_SuperfamilyHoldout"]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(test_rows)
    print(f"test rows: {len(test_rows)}")
    print(f"  OOD_FoldHoldout=True:        {n_f} ({n_f/len(test_rows)*100:.1f}%)")
    print(f"  OOD_SuperfamilyHoldout=True: {n_s} ({n_s/len(test_rows)*100:.1f}%)")
    print(f"  unannotated (-> False):      {n_un}")
    print(f"-> {args.out}")


if __name__ == "__main__":
    main()
