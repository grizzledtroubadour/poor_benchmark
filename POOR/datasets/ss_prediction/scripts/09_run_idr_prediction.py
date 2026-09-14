#!/usr/bin/env python3
"""
Run metapredict v3 IDR prediction on SSP representative sequences.
Computes IDR_ratio and marks OOD-IDR (IDR_ratio > 0.3).

口径说明：与 .skills/ood-annotation-toolkit/scripts/idr_ood.py 的
`--mode residue-ratio --threshold 0.3 --clean strip` 口径一致
（residue-ratio、阈值 0.3、预测前剔除非标准氨基酸）。

Usage:
    python 09_run_idr_prediction.py
"""

import os
import re
import csv
import time
from pathlib import Path

from Bio import SeqIO
import metapredict

# ------------------------------------------------------------------
# Config
# ------------------------------------------------------------------
TASK_DIR = Path(__file__).resolve().parent.parent  # datasets/ss_prediction
FASTA_PATH = TASK_DIR / "output" / "stage4_nonredundant_sequences.fa"
OUTPUT_CSV = TASK_DIR / "output" / "idr_predictions.csv"
TEST_CSV = TASK_DIR / "splits" / "ssp_test.csv"
NONSTANDARD_AA = set("XBZJU")
IDR_THRESHOLD = 0.3   # IDR_ratio > 0.3 => OOD-IDR
BATCH_SIZE = 2000


def clean_sequence(seq: str) -> str:
    """Remove non-standard amino acids before metapredict prediction."""
    return "".join(aa for aa in seq if aa not in NONSTANDARD_AA)


def main():
    print(f"metapredict version: {metapredict.__version__}")
    print(f"Input FASTA: {FASTA_PATH}")
    print(f"Output CSV:  {OUTPUT_CSV}")
    print()

    records = list(SeqIO.parse(FASTA_PATH, "fasta"))
    total = len(records)
    print(f"Total sequences: {total}")

    # ------------------------------------------------------------------
    # Run predictions in batches
    # ------------------------------------------------------------------
    results = []
    start_time = time.time()

    for batch_start in range(0, total, BATCH_SIZE):
        batch_end = min(batch_start + BATCH_SIZE, total)
        batch_recs = records[batch_start:batch_end]

        # Prepare cleaned sequences and IDs
        ids = []
        clean_seqs = []
        for rec in batch_recs:
            raw_seq = str(rec.seq).upper()
            clean_seq = clean_sequence(raw_seq)
            ids.append(rec.id)
            clean_seqs.append(clean_seq)

        # Batch predict with domains
        disorder_objects = metapredict.predict_disorder_batch(
            clean_seqs,
            return_domains=True,
            show_progress_bar=False,
        )

        for uid, clean_seq, obj in zip(ids, clean_seqs, disorder_objects):
            seq_len = len(clean_seq)
            if seq_len == 0:
                idr_ratio = 0.0
            else:
                # IDR_ratio = residues with disorder score > 0.5 / total length
                idr_count = int((obj.disorder > 0.5).sum())
                idr_ratio = idr_count / seq_len

            results.append({
                "unique_id": uid,
                "seq_len": seq_len,
                "idr_residues": idr_count if seq_len > 0 else 0,
                "idr_ratio": round(idr_ratio, 6),
                "ood_idr": idr_ratio > IDR_THRESHOLD,
            })

        elapsed = time.time() - start_time
        rate = (batch_end) / elapsed if elapsed > 0 else 0
        eta = (total - batch_end) / rate if rate > 0 else 0
        print(f"  Processed {batch_end}/{total}  ({rate:.1f} seq/s, ETA {eta/60:.1f} min)")

    total_time = time.time() - start_time
    print(f"\nPrediction done in {total_time/60:.2f} minutes")

    # ------------------------------------------------------------------
    # Save IDR predictions
    # ------------------------------------------------------------------
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["unique_id", "seq_len", "idr_residues", "idr_ratio", "ood_idr"],
        )
        writer.writeheader()
        writer.writerows(results)

    # Summary stats
    ood_count = sum(1 for r in results if r["ood_idr"])
    avg_ratio = sum(r["idr_ratio"] for r in results) / len(results)
    print(f"\nIDR prediction summary:")
    print(f"  Total sequences: {len(results)}")
    print(f"  Mean IDR_ratio:  {avg_ratio:.4f}")
    print(f"  OOD-IDR (>0.3):  {ood_count} ({ood_count/len(results)*100:.2f}%)")

    # ------------------------------------------------------------------
    # Merge OOD_IDR into ssp_test.csv
    # ------------------------------------------------------------------
    if TEST_CSV.exists():
        import pandas as pd

        df = pd.read_csv(TEST_CSV)
        idr_df = pd.DataFrame(results)

        # Merge on unique_id
        df = df.merge(
            idr_df[["unique_id", "ood_idr"]].rename(columns={"ood_idr": "OOD_IDR"}),
            on="unique_id",
            how="left",
        )

        # Ensure boolean
        df["OOD_IDR"] = df["OOD_IDR"].fillna(False).astype(bool)

        # Default=False（极端长度）行直接置 NaN（"未计算"语义，全项目统一惯例，
        # 流水线位置写入；finalize_ood_columns.py --mode nan 仅作兜底校验）
        if "Default" in df.columns:
            nd = ~df["Default"].astype(bool)
            if nd.any():
                df["OOD_IDR"] = df["OOD_IDR"].astype(object)
                df.loc[nd, "OOD_IDR"] = pd.NA
                print(f"  Default=False rows -> NaN: {int(nd.sum())}")

        # Move OOD_IDR to end (after OOD_ExtremeLong)
        cols = list(df.columns)
        cols.remove("OOD_IDR")
        cols.append("OOD_IDR")
        df = df[cols]

        df.to_csv(TEST_CSV, index=False)

        test_ood_count = int((df["OOD_IDR"] == True).sum())  # noqa: E712
        print(f"\nMerged into {TEST_CSV}:")
        print(f"  Test set OOD-IDR: {test_ood_count} / {len(df)} ({test_ood_count/len(df)*100:.2f}%)")
    else:
        print(f"\nWarning: {TEST_CSV} not found, skipping merge.")


if __name__ == "__main__":
    main()
