#!/usr/bin/env python3
"""Intrinsic-disorder (IDR) OOD computation (ood-annotation-toolkit).

Unified tool with task-specific historical definitions: select mode explicitly.
- residue-ratio: fraction of residues with disorder >0.5 (ss/func/fold, threshold 0.3).
- region-ratio: longest continuous disordered region / sequence length (ppi, threshold 0.3).
Set --threshold explicitly for affinity/optimal_ph, historically 0.1.

Sequence cleaning also differs by task; select --clean:
- strip (default): delete characters outside the standard 20 amino acids (fold/ss).
- map: B→N, Z→Q, J→L, U→C, O→K; other nonstandard characters→A (ppi).

Cache cleaned sequences incrementally with --cache. CSV columns:
sequence, residue_ratio, max_region_ratio. Cache both definitions so changes in
mode/threshold require no predictions; save per batch to support resuming.

For multiple sequence columns (--seq-col aa_seq1,aa_seq2), compute each chain separately.
OOD_IDR is True if either chain exceeds threshold (ppi any semantics).

Usage:
    python idr_ood.py --test ssp_test.csv --out test_idr.csv \
        --mode residue-ratio --threshold 0.3 --cache output/idr_cache.csv
    python idr_ood.py --test ppi_test.csv --seq-col aa_seq1,aa_seq2 \
        --mode region-ratio --clean map --out ppi_idr.csv --cache output/idr_cache.csv
"""
import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

STD_AA = "ACDEFGHIKLMNPQRSTVWY"


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--test", required=True, type=Path, help="Test CSV containing sequence columns")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--seq-col", default="aa_seq",
                    help="Comma-separated sequence columns, e.g. aa_seq1,aa_seq2")
    ap.add_argument("--mode", required=True, choices=["residue-ratio", "region-ratio"],
                    help="residue-ratio=disordered-residue fraction (ss/func/fold); "
                         "region-ratio=longest continuous disordered-region fraction (ppi)")
    ap.add_argument("--threshold", type=float, default=0.3,
                    help="OOD_IDR threshold: default 0.3; affinity/pH historically use 0.1")
    ap.add_argument("--clean", choices=["strip", "map"], default="strip",
                    help="Nonstandard residue cleaning: strip for fold/ss, map for ppi")
    ap.add_argument("--cache", type=Path, default=None,
                    help="Incremental sequence cache CSV; reuse existing entries and append predictions")
    ap.add_argument("--batch-size", type=int, default=512)
    ap.add_argument("--default-col", default="Default",
                    help="Default column; apply --nondefault to OOD_IDR on Default=False rows")
    ap.add_argument("--nondefault", choices=["nan", "false"], default="nan",
                    help="Default=False values: nan follows the project-wide generation-time convention; "
                         "false is retained only for legacy compatibility")
    return ap.parse_args()


def clean_sequence(seq, mode):
    seq = seq.upper()
    if mode == "strip":
        return re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", seq)
    seq = seq.replace("B", "N").replace("Z", "Q").replace("J", "L")
    seq = seq.replace("U", "C").replace("O", "K")
    return re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "A", seq)


def metrics_from_scores(scores):
    """Compute both definitions from per-residue disorder scores."""
    arr = np.asarray(scores, dtype=float)
    if len(arr) == 0:
        return np.nan, np.nan
    binary = arr > 0.5
    residue_ratio = float(binary.mean())
    # Longest continuous True region divided by sequence length.
    max_run = run = 0
    for b in binary:
        run = run + 1 if b else 0
        if run > max_run:
            max_run = run
    return residue_ratio, max_run / len(arr)


def main():
    args = parse_args()
    seq_cols = [c.strip() for c in args.seq_col.split(",") if c.strip()]
    multi = len(seq_cols) > 1

    test = pd.read_csv(args.test, dtype=str, keep_default_na=False)
    for col in seq_cols:
        if col not in test.columns:
            sys.exit(f"Error: test lacks sequence column {col!r}")
    print(f"test={len(test)}，seq_cols={seq_cols}，mode={args.mode}，"
          f"threshold={args.threshold}，clean={args.clean}")

    # Collect deduplicated cleaned sequences.
    clean_seqs = set()
    for col in seq_cols:
        clean_seqs.update(clean_sequence(s, args.clean) for s in test[col] if s)
    clean_seqs = sorted(clean_seqs)
    print(f"Unique cleaned sequences: {len(clean_seqs)}")

    cache = {}
    if args.cache is not None and args.cache.exists():
        cdf = pd.read_csv(args.cache, dtype=str, keep_default_na=False)
        for r in cdf.itertuples(index=False):
            cache[r.sequence] = (float(r.residue_ratio), float(r.max_region_ratio))
        print(f"Cache hits: {len(cache)} ({args.cache})")

    missing = [s for s in clean_seqs if s not in cache]
    print(f"Pending predictions: {len(missing)}")
    if missing:
        try:
            import metapredict as meta
        except ImportError:
            sys.exit("Error: metapredict v3 is not installed in this environment. "
                     "Switch to a Python with metapredict, e.g. /root/miniconda3/bin/python, "
                     "or supply --cache covering these sequences for classification only.")
        print(f"metapredict version: {meta.__version__}")
        empty = [s for s in missing if len(s) == 0]
        for s in empty:
            cache[s] = (np.nan, np.nan)
        todo = [s for s in missing if len(s) > 0]
        for i in range(0, len(todo), args.batch_size):
            chunk = todo[i:i + args.batch_size]
            scores = meta.predict_disorder_batch(chunk, show_progress_bar=False)
            for seq, sc in zip(chunk, scores):
                cache[seq] = metrics_from_scores(sc[1])
            print(f"  {min(i + args.batch_size, len(todo))}/{len(todo)}")
            if args.cache is not None:  # Save each batch to support resuming.
                save_cache(cache, args.cache)
        if args.cache is not None:
            save_cache(cache, args.cache)

    metric_idx = 0 if args.mode == "residue-ratio" else 1

    def metric_of(seq):
        if not seq:
            return np.nan
        return cache.get(clean_sequence(seq, args.clean), (np.nan, np.nan))[metric_idx]

    flags = []
    ratio_cols = []
    for col in seq_cols:
        ratio_col = "idr_ratio" if not multi else f"idr_ratio_{col}"
        test[ratio_col] = test[col].map(metric_of)
        ratio_cols.append(ratio_col)
        flags.append(test[ratio_col] > args.threshold)
    ood = flags[0]
    for f in flags[1:]:
        ood = ood | f
    test["OOD_IDR"] = ood.map({True: "True", False: "False"})

    # Default=False extreme-length rows receive NaN (not computed) at generation time,
    # by project convention; finalize_ood_columns.py --mode nan is fallback validation only.
    if args.default_col in test.columns:
        nd = test[args.default_col] != "True"
        if nd.any():
            fill = "" if args.nondefault == "nan" else "False"
            test.loc[nd, "OOD_IDR"] = fill
            print(f"  Default=False: {int(nd.sum())} rows OOD_IDR set to {'NaN' if fill == '' else 'False'}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    test.to_csv(args.out, index=False, lineterminator="\n")
    n_nan = int(sum(test[c].isna().sum() for c in ratio_cols))
    print(f"-> {args.out} ({len(test)} rows)")
    print(f"  OOD_IDR=True: {int(ood.sum())}; chains with NaN metrics (flagged False): {n_nan}")


def save_cache(cache, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(
        [(s, v[0], v[1]) for s, v in cache.items()],
        columns=["sequence", "residue_ratio", "max_region_ratio"])
    df.to_csv(path, index=False, lineterminator="\n")


if __name__ == "__main__":
    main()
