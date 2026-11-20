#!/usr/bin/env python3
"""Regression value-bin LongTail OOD computation (ood-annotation-toolkit).

Canonical definition shared by kcat 08_add_value_bin_longtail_ood.py, affinity
12_add_value_bin_longtail_ood.py, and optimal_ph 05_add_phbin_longtail_ood.py:
- Bin regression labels at fixed width; optionally give clipped boundary values separate bins.
  kcat log10(kcat) is censored to [-6, 6]; values exactly ±6.0 each form separate bins
  for training counts, preventing censored boundary mass from hiding tails.
- Count training samples per bin.
- Default=True test rows: training count for their bin <=threshold -> True;
  missing labels or values outside bin range -> False.
- Default=False rows: NaN by default (--nondefault nan, project-wide convention);
  false is available for historical optimal_ph compatibility.

Historical task parameters, to preserve on recomputation; see SKILL.md:
- kcat: label=log10(kcat), bins [-6.5,7.0), width 0.5, clip ±6.0, le50/le100, prefix OOD_LongTail_KcatBin.
- affinity: label=-log10(K), bins [-0.5,16.5), width 0.5, no clip, le50/le100, prefix OOD_LongTail_AffBin.
- optimal_ph: label=pH, bins [1.5,13.5), width 0.5, no clip, le50/le100, prefix OOD_LongTail_pHbin.

Usage:
    python value_bin_longtail_ood.py --train kcat_train.csv --test kcat_test.csv \
        --out kcat_test_flagged.csv --bin-start -6.5 --bin-end 7.0 --bin-width 0.5 \
        --clip -6,6 --thresholds 50,100 --col-prefix OOD_LongTail_KcatBin
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--train", required=True, type=Path)
    ap.add_argument("--test", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path,
                    help="Output CSV; may match --test for in-place updates")
    ap.add_argument("--label-col", default="label")
    ap.add_argument("--bin-start", type=float, required=True)
    ap.add_argument("--bin-end", type=float, required=True)
    ap.add_argument("--bin-width", type=float, default=0.5)
    ap.add_argument("--clip", default=None,
                    help="Optional lo,hi: boundary values equal to lo/hi form separate bins, as in kcat")
    ap.add_argument("--thresholds", default="50,100",
                    help="Comma-separated training-bin count thresholds; <= threshold means True")
    ap.add_argument("--col-prefix", required=True,
                    help="Output prefix, e.g. OOD_LongTail_KcatBin -> OOD_LongTail_KcatBin_le50")
    ap.add_argument("--default-col", default="Default")
    ap.add_argument("--nondefault", choices=["nan", "false"], default="nan",
                    help="Default=False values: nan follows the project-wide generation-time convention; "
                         "false is retained only for legacy compatibility")
    return ap.parse_args()


def main():
    args = parse_args()
    bins = np.arange(args.bin_start, args.bin_end, args.bin_width)
    thresholds = [int(t) for t in args.thresholds.split(",")]
    clip = None
    if args.clip:
        lo, hi = (float(x) for x in args.clip.split(","))
        clip = (lo, hi, f"clip_eq_{lo}", f"clip_eq_{hi}")

    def bin_keys(values):
        """Map values to bin keys; clipped boundaries get separate bins; out-of-range -> NaN."""
        keys = pd.cut(values, bins=bins).astype(object)
        if clip:
            keys[values == clip[0]] = clip[2]
            keys[values == clip[1]] = clip[3]
        return keys

    train = pd.read_csv(args.train)
    test = pd.read_csv(args.test, dtype=str, keep_default_na=False)
    if args.label_col not in train.columns or args.label_col not in test.columns:
        raise SystemExit(f"Error: train/test lacks label column {args.label_col!r}")

    train_keys = bin_keys(train[args.label_col].astype(float))
    bin_counts = train_keys.value_counts().to_dict()
    print(f"train bins: {len(bin_counts)}, "
          + ", ".join(f"le{t}={sum(c <= t for c in bin_counts.values())}" for t in thresholds))
    if clip:
        print(f"  clip bins: {clip[2]}={bin_counts.get(clip[2], 0)}, "
              f"{clip[3]}={bin_counts.get(clip[3], 0)}")

    test_keys = bin_keys(pd.to_numeric(test[args.label_col], errors="coerce"))
    mask_default = test[args.default_col] == "True"
    n_nd = int((~mask_default).sum())

    for thr in thresholds:
        col = f"{args.col_prefix}_le{thr}"
        if args.nondefault == "nan":
            new = pd.Series("", index=test.index, dtype=object)
        else:
            new = pd.Series("False", index=test.index, dtype=object)
        flagged = test_keys[mask_default].apply(
            lambda k: bin_counts.get(k, 0) <= thr if pd.notna(k) else False)
        new[mask_default] = flagged.map({True: "True", False: "False"})
        test[col] = new
        print(f"  {col}: True={int(flagged.sum())}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    test.to_csv(args.out, index=False, lineterminator="\n")
    print(f"-> {args.out} ({len(test)} rows; {n_nd} Default=False rows set to "
          f"{'NaN' if args.nondefault == 'nan' else 'False'}）")


if __name__ == "__main__":
    main()
