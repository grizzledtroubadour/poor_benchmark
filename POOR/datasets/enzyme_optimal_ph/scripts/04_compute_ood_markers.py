#!/usr/bin/env python3
"""
Compute OOD flags for optimal_ph_prediction_test.csv.
- seq_Redundancy_90 ~ seq_Redundancy_30（MMseqs2 easy-search）
- TM-score_0.9 ~ TM-score_0.3（Foldseek easy-search）
- OOD_Orphan: no significant MMseqs2 easy-search hit.
- OOD_IDR: metapredict v3 batch predictions, IDR_ratio >0.1, Default samples only.
  The historical threshold for this task is 0.1, unlike 0.3 for other tasks; see
  the idr_ood.py definitions table in .skills/ood-annotation-toolkit/SKILL.md.
  This logic incorporates the archived standalone compute_ood_idr.py implementation.

Set all these OOD columns to NaN for Default=False extreme-length rows before saving.
NaN means not computed, unlike False; retain only the OOD_ExtremeShort /
OOD_ExtremeLong reason flags. finalize_ood_columns.py --mode nan is only a
fallback check, no longer responsible for NaN conversion.
"""
import os
import re
import json
import subprocess
import tempfile
from pathlib import Path
from collections import defaultdict

import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SPLITS_DIR = ROOT / "splits"
PDBS_DIR = ROOT / "pdbs"
OUTPUT_DIR = ROOT / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

TRAIN_PATH = SPLITS_DIR / "optimal_ph_prediction_train.csv"
VAL_PATH = SPLITS_DIR / "optimal_ph_prediction_val.csv"
TEST_PATH = SPLITS_DIR / "optimal_ph_prediction_test.csv"

SEQ_REDUNDANCY_THRESHOLDS = [90, 80, 70, 60, 50, 40, 30]
TM_SCORE_THRESHOLDS = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3]


def clean_seq(seq):
    return re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", seq.upper())


def write_fasta(df, fasta_path):
    with open(fasta_path, "w") as f:
        for _, row in df.iterrows():
            seq = clean_seq(row["aa_seq"])
            f.write(f">{row['unique_id']}\n{seq}\n")


def run_mmseqs_search(query_fasta, target_fasta, output_tsv, tmp_dir, sensitivity=7):
    """Run MMseqs2 easy-search and return the result TSV."""
    tmp_dir = Path(tmp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "mmseqs", "easy-search",
        str(query_fasta), str(target_fasta), str(output_tsv), str(tmp_dir),
        "-s", str(sensitivity),
        "--format-output", "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits",
    ]
    print(f"[INFO] Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def run_foldseek_search(query_dir, target_dir, output_tsv, tmp_dir):
    """Run Foldseek easy-search and return the result TSV."""
    tmp_dir = Path(tmp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "foldseek", "easy-search",
        str(query_dir), str(target_dir), str(output_tsv), str(tmp_dir),
        "--format-output", "query,target,alntmscore",
    ]
    print(f"[INFO] Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def parse_mmseqs_max_identity(tsv_path, test_ids):
    """Parse MMseqs2 results and return maximum pident per query."""
    max_identity = {uid: 0.0 for uid in test_ids}
    if not Path(tsv_path).exists():
        return max_identity
    with open(tsv_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            query, pident = parts[0], float(parts[2])
            if query in max_identity:
                max_identity[query] = max(max_identity[query], pident)
    return max_identity


def parse_foldseek_max_tmscore(tsv_path, test_ids):
    """Parse Foldseek results and return maximum alntmscore per query."""
    max_tmscore = {uid: 0.0 for uid in test_ids}
    if not Path(tsv_path).exists():
        return max_tmscore
    with open(tsv_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            query, tmscore = parts[0], float(parts[2])
            if query in max_tmscore:
                max_tmscore[query] = max(max_tmscore[query], tmscore)
    return max_tmscore


def compute_seq_redundancy(test_df, train_val_fasta, tmp_dir):
    """Compute seq_Redundancy flags."""
    print("[INFO] Computing seq_Redundancy...")
    test_fasta = tmp_dir / "test.fa"
    result_tsv = tmp_dir / "mmseqs_search.tsv"
    write_fasta(test_df, test_fasta)
    run_mmseqs_search(test_fasta, train_val_fasta, result_tsv, tmp_dir / "mmseqs_tmp")

    max_identity = parse_mmseqs_max_identity(result_tsv, test_df["unique_id"].tolist())

    for thr in SEQ_REDUNDANCY_THRESHOLDS:
        col = f"seq_Redundancy_{thr}"
        # Since 2026-08-29, use strict less-than consistently with the project;
        # the former <= flagged 35 threshold-equal rows, which were recomputed.
        test_df[col] = test_df["unique_id"].map(lambda uid: max_identity.get(uid, 0.0) < thr)

    return test_df


def compute_tm_score(test_df, train_val_dir, tmp_dir):
    """Compute TM-score flags."""
    print("[INFO] Computing TM-score...")
    # Create a temporary query directory containing all test structures.
    query_dir = tmp_dir / "test_structures"
    query_dir.mkdir(parents=True, exist_ok=True)
    for uid in test_df["unique_id"]:
        src = PDBS_DIR / f"{uid}.cif"
        if not src.exists():
            src = PDBS_DIR / f"{uid}.pdb"
        if src.exists():
            dst = query_dir / f"{uid}{src.suffix}"
            if not dst.exists():
                os.link(src, dst)

    result_tsv = tmp_dir / "foldseek_search.tsv"
    run_foldseek_search(query_dir, train_val_dir, result_tsv, tmp_dir / "foldseek_tmp")

    max_tmscore = parse_foldseek_max_tmscore(result_tsv, test_df["unique_id"].tolist())

    for thr in TM_SCORE_THRESHOLDS:
        col = f"TM-score_{thr:.1f}".replace(".0", "")
        # Use standardized column names such as TM-score_0.9.
        col = f"TM-score_{thr}"
        # Strict less-than since 2026-08-29; no threshold-equal samples, so values are unchanged.
        test_df[col] = test_df["unique_id"].map(lambda uid: max_tmscore.get(uid, 0.0) < thr)

    return test_df


def compute_orphan(test_df, train_val_fasta, tmp_dir):
    """Compute OOD_Orphan: no significant test-to-train+val hit."""
    print("[INFO] Computing OOD_Orphan...")
    test_fasta = tmp_dir / "test_orphan.fa"
    result_tsv = tmp_dir / "mmseqs_orphan.tsv"
    write_fasta(test_df, test_fasta)
    run_mmseqs_search(test_fasta, train_val_fasta, result_tsv, tmp_dir / "mmseqs_orphan_tmp", sensitivity=7)

    has_hit = set()
    if Path(result_tsv).exists():
        with open(result_tsv) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 1:
                    has_hit.add(parts[0])

    test_df["OOD_Orphan"] = ~test_df["unique_id"].isin(has_hit)
    return test_df


def compute_idr(test_df, idr_threshold=0.1):
    """Compute OOD_IDR using metapredict v3 batches for Default samples only; other rows -> False."""
    print("[INFO] Computing OOD_IDR...")
    try:
        from metapredict import predict_disorder_fasta
    except ImportError as e:
        print(f"[ERROR] metapredict not available: {e}, skipping OOD_IDR")
        test_df["OOD_IDR"] = False
        return test_df

    # Write FASTA, retained as an intermediate under output/.
    fasta_path = OUTPUT_DIR / "test_sequences_for_idr.fasta"
    print(f"[INFO] Writing FASTA: {fasta_path}")
    with open(fasta_path, "w") as f:
        for _, row in test_df.iterrows():
            seq = clean_seq(row["aa_seq"])
            f.write(f">{row['unique_id']}\n{seq}\n")

    # Batch prediction.
    print("[INFO] Computing IDR ratios with metapredict v3 (batch)...")
    disorder_dict = predict_disorder_fasta(str(fasta_path), version="V3")

    idr_ratios = {}
    for uid, item in disorder_dict.items():
        # predict_disorder_fasta returns [sequence, scores].
        if isinstance(item, (list, tuple)) and len(item) == 2:
            scores = np.array(item[1])
        else:
            scores = np.array(item)
        idr_len = int((scores > 0.5).sum())
        idr_ratio = idr_len / len(scores) if len(scores) > 0 else 0.0
        idr_ratios[uid] = idr_ratio

    test_df["OOD_IDR"] = test_df.apply(
        lambda row: idr_ratios.get(row["unique_id"], 0.0) > idr_threshold if row["Default"] else False,
        axis=1,
    )
    return test_df


def main():
    print("[INFO] Loading splits...")
    train_df = pd.read_csv(TRAIN_PATH)
    val_df = pd.read_csv(VAL_PATH)
    test_df = pd.read_csv(TEST_PATH)

    print(f"[INFO] Train: {len(train_df)}, Val: {len(val_df)}, Test: {len(test_df)}")

    train_val_df = pd.concat([train_df, val_df], ignore_index=True)
    print(f"[INFO] Train+Val: {len(train_val_df)}")

    with tempfile.TemporaryDirectory(prefix="ood_", dir=str(OUTPUT_DIR)) as tmp:
        tmp_dir = Path(tmp)

        # 1. Prepare train+val FASTA.
        train_val_fasta = tmp_dir / "train_val.fa"
        write_fasta(train_val_df, train_val_fasta)

        # 2. Prepare train+val structures for Foldseek.
        train_val_dir = tmp_dir / "train_val_structures"
        train_val_dir.mkdir(parents=True, exist_ok=True)
        for uid in train_val_df["unique_id"]:
            src = PDBS_DIR / f"{uid}.cif"
            if not src.exists():
                src = PDBS_DIR / f"{uid}.pdb"
            if src.exists():
                dst = train_val_dir / f"{uid}{src.suffix}"
                if not dst.exists():
                    os.link(src, dst)

        # 3. Compute OOD flags.
        test_df = compute_seq_redundancy(test_df, train_val_fasta, tmp_dir)
        test_df = compute_tm_score(test_df, train_val_dir, tmp_dir)
        test_df = compute_orphan(test_df, train_val_fasta, tmp_dir)
        test_df = compute_idr(test_df)

    # 4. Ensure OOD columns are boolean.
    ood_cols = [c for c in test_df.columns if c.startswith(("seq_Redundancy_", "TM-score_", "OOD_"))]
    for col in ood_cols:
        test_df[col] = test_df[col].astype(bool)

    # 4b. Default=False extreme-length rows: set all OOD columns except length flags
    #     to NaN (not computed), following project convention at this pipeline stage.
    nd_mask = ~test_df["Default"]
    keep_cols = {"OOD_ExtremeShort", "OOD_ExtremeLong"}
    nan_cols = [c for c in ood_cols if c not in keep_cols]
    for col in nan_cols:
        test_df[col] = test_df[col].astype(object)
        test_df.loc[nd_mask, col] = pd.NA
    print(f"[INFO] Default=False rows -> NaN: {int(nd_mask.sum())} rows x {len(nan_cols)} cols")

    # 5. Save.
    test_df.to_csv(TEST_PATH, index=False)
    print(f"[INFO] Updated {TEST_PATH}")

    # 6. Summary statistics.
    summary = {
        "n_test": int(len(test_df)),
        "n_default": int(test_df["Default"].sum()),
        "n_extreme_long": int(test_df["OOD_ExtremeLong"].sum()),
        "n_extreme_short": int(test_df["OOD_ExtremeShort"].sum()),
        "OOD_IDR_threshold": 0.1,
        "OOD_IDR_default_only": True,
        "nondefault_nan_cols": len(nan_cols),
    }

    for col in ood_cols:
        summary[col] = int((test_df[col] == True).sum())  # noqa: E712

    summary_path = OUTPUT_DIR / "ood_markers_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[INFO] Summary saved to {summary_path}")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
