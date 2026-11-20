#!/usr/bin/env python3
"""
EnzyBase12k dataset splitting
- Reserve samples with seq_length <60 or >1000 as OOD-ExtremeLength test samples.
- Assign all remaining structure_source == "PDB" samples to test.
- Randomly split the remaining AF samples 9:1 into train and validation.
"""
import json
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
SPLITS_DIR = ROOT / "splits"
SPLITS_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR = ROOT / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

META_PATH = DATA_DIR / "metadata.csv"


def main(random_seed=42):
    print(f"[INFO] Loading metadata: {META_PATH}")
    df = pd.read_csv(META_PATH)
    n_total = len(df)
    print(f"[INFO] Total samples: {n_total}")

    # 1. OOD-ExtremeLength: length <60 or >1000.
    mask_extreme = (df["seq_length"] < 60) | (df["seq_length"] > 1000)
    df_extreme = df[mask_extreme].copy()
    df_rest = df[~mask_extreme].copy()

    print(f"[INFO] OOD-ExtremeLength samples: {len(df_extreme)} ({len(df_extreme)/n_total*100:.2f}%)")
    print(f"[INFO] Remaining samples: {len(df_rest)} ({len(df_rest)/n_total*100:.2f}%)")

    # 2. Assign remaining PDB samples to test.
    mask_pdb = df_rest["structure_source"] == "PDB"
    df_test = df_rest[mask_pdb].copy()
    df_af = df_rest[~mask_pdb].copy()

    print(f"[INFO] Test (PDB) samples: {len(df_test)} ({len(df_test)/n_total*100:.2f}%)")
    print(f"[INFO] AF samples for train/val split: {len(df_af)} ({len(df_af)/n_total*100:.2f}%)")

    # 3. Split AF samples 9:1 into train/validation.
    rng = np.random.default_rng(random_seed)
    n_af = len(df_af)
    indices = np.arange(n_af)
    rng.shuffle(indices)

    n_val = max(1, int(n_af * 0.1))
    val_indices = indices[:n_val]
    train_indices = indices[n_val:]

    df_val = df_af.iloc[val_indices].copy()
    df_train = df_af.iloc[train_indices].copy()

    print(f"[INFO] Train samples: {len(df_train)} ({len(df_train)/n_total*100:.2f}%)")
    print(f"[INFO] Val samples: {len(df_val)} ({len(df_val)/n_total*100:.2f}%)")

    # 4. Save split results.
    # Retain core fields.
    columns_to_keep = [
        "uniprot_id",
        "ph_optimum",
        "structure_source",
        "structure_mode",
        "pdb_id_final",
        "resolution",
        "pdb_chain",
        "ec_id",
        "cut_to_uniprot_domain",
        "systematic_name",
        "organism",
        "seq_length",
        "uniprot_seq_cut",
    ]

    df_train[columns_to_keep].to_csv(SPLITS_DIR / "train.csv", index=False)
    df_val[columns_to_keep].to_csv(SPLITS_DIR / "val.csv", index=False)
    df_test[columns_to_keep].to_csv(SPLITS_DIR / "test.csv", index=False)
    df_extreme[columns_to_keep].to_csv(SPLITS_DIR / "ood_extremelength.csv", index=False)

    print(f"[INFO] Splits saved to {SPLITS_DIR}")

    # 5. Generate summary statistics.
    def split_stats(sub_df, name):
        return {
            "split": name,
            "n_samples": int(len(sub_df)),
            "pct_of_total": round(len(sub_df) / n_total * 100, 2),
            "ph_mean": round(float(sub_df["ph_optimum"].mean()), 3),
            "ph_median": round(float(sub_df["ph_optimum"].median()), 3),
            "ph_std": round(float(sub_df["ph_optimum"].std()), 3),
            "seq_length_mean": round(float(sub_df["seq_length"].mean()), 1),
            "seq_length_median": int(sub_df["seq_length"].median()),
            "seq_length_min": int(sub_df["seq_length"].min()),
            "seq_length_max": int(sub_df["seq_length"].max()),
            "n_pdb": int((sub_df["structure_source"] == "PDB").sum()),
            "n_af": int((sub_df["structure_source"] == "AF").sum()),
            "n_with_ec": int(sub_df["ec_id"].notna().sum()),
        }

    summary = {
        "random_seed": random_seed,
        "n_total": n_total,
        "splits": [
            split_stats(df_train, "train"),
            split_stats(df_val, "val"),
            split_stats(df_test, "test"),
            split_stats(df_extreme, "ood_extremelength"),
        ],
    }

    summary_path = OUTPUT_DIR / "split_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[INFO] Summary saved to {summary_path}")

    print("\n===== Split Summary =====")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--random_seed", type=int, default=42)
    args = parser.parse_args()
    main(args.random_seed)
