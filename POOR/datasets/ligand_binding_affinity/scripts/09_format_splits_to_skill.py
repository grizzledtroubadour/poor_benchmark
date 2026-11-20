#!/usr/bin/env python3
"""
Script: 09_format_splits_to_skill.py
Purpose: Export splits using the unified POOD-Benchmark schema:
      - Filenames: {task_name}_train.csv / {task_name}_val.csv / {task_name}_test.csv
      - Core columns: unique_id, aa_seq, struct_file, ligand_smiles, ligand_ecfp4, label
      - Additional test OOD columns: Default, OOD_ExtremeShort, OOD_ExtremeLong,
        seq_Redundancy_90..30, OOD_Orphan
      - Store label (pK) as a string with six decimal places.

Note (superseded): the mmseqs2 homology code here is the historical first-pass implementation.
The canonical implementation and replay tool are in .skills/homology-ood-annotation:
scripts/seq_homology_ood.py (supports easy-search and --m8-cache replay).

Usage:
    cd ligand_binding_affinity/
    python scripts/09_format_splits_to_skill.py
"""

import argparse
import json
import shutil
import subprocess
from pathlib import Path

import pandas as pd

TASK_NAME = "ligand_binding_affinity"
SEQ_THRESHOLDS = [90, 80, 70, 60, 50, 40, 30]


def write_fasta(df, path: Path, id_col="unique_id", seq_col="aa_seq"):
    """Write FASTA; remove chain separators '|' because mmseqs2 accepts only standard amino acids."""
    with open(path, "w", encoding="utf-8") as f:
        for _, row in df.iterrows():
            uid = row[id_col]
            seq = row[seq_col]
            if pd.isna(seq):
                seq = ""
            f.write(f">{uid}\n{seq.replace('|', '')}\n")


def parse_max_identity(out_tsv):
    """Parse convertalis output (query,target,pident,evalue) and return {query: max_identity}."""
    max_identity = {}
    if out_tsv.exists():
        with open(out_tsv, "r", encoding="utf-8") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 3:
                    continue
                q, _, pident = parts[0], parts[1], parts[2]
                pid = float(pident)
                max_identity[q] = max(max_identity.get(q, 0.0), pid)
    return max_identity


def run_mmseqs2(query_fasta, target_fasta, out_tsv, tmp_dir):
    """Run mmseqs2 search and return the best-hit identity mapping {query: max_identity}.

    Superseded: use seq_homology_ood.py in the homology-ood-annotation skill as the canonical implementation.
    """
    tmp_dir.mkdir(parents=True, exist_ok=True)
    query_db = tmp_dir / "query_db"
    target_db = tmp_dir / "target_db"
    result_db = tmp_dir / "result_db"

    subprocess.run(["mmseqs", "createdb", str(query_fasta), str(query_db)], check=True)
    subprocess.run(["mmseqs", "createdb", str(target_fasta), str(target_db)], check=True)
    subprocess.run(
        ["mmseqs", "search", str(query_db), str(target_db), str(result_db), str(tmp_dir / "tmp")],
        check=True,
    )
    subprocess.run(
        [
            "mmseqs", "convertalis",
            str(query_db), str(target_db), str(result_db), str(out_tsv),
            "--format-output", "query,target,pident,evalue",
        ],
        check=True,
    )

    return parse_max_identity(out_tsv)


def build_base_df(meta_df):
    """Build a unified-schema base DataFrame from the complete metadata."""
    out = pd.DataFrame()
    out["unique_id"] = meta_df["pdb_id"].astype(str)
    out["aa_seq"] = meta_df["protein_sequence"]
    out["struct_file"] = meta_df["pdb_id"].astype(str) + ".pdb"
    out["ligand_smiles"] = meta_df["ligand_smiles"]
    out["ligand_ecfp4"] = meta_df["ligand_ecfp4"]
    # Store labels as strings with six decimal places.
    out["label"] = meta_df["pK"].apply(lambda x: f"{x:.6f}" if pd.notna(x) else "")
    # Retain auxiliary columns needed for subsequent OOD annotation.
    out["_protein_length"] = meta_df["protein_length"]
    return out


def main():
    parser = argparse.ArgumentParser(description="Format split CSVs to POOD-Benchmark skill spec")
    parser.add_argument(
        "--dataset",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "output" / "pdbbind_v2020_dataset.csv",
        help="Full dataset CSV with metadata",
    )
    parser.add_argument(
        "--splits-dir",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "splits",
        help="Directory containing split CSVs",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "output",
        help="Directory for intermediate files",
    )
    parser.add_argument(
        "--mmseqs-cache",
        type=Path,
        default=None,
        help="Optional existing convertalis TSV (query,target,pident,evalue); "
             "if supplied, parse the cache instead of rerunning mmseqs2.",
    )
    args = parser.parse_args()

    # Read complete metadata.
    meta_df = pd.read_csv(args.dataset)
    base_df = build_base_df(meta_df)

    # Read unique_id lists from current split files to determine membership.
    def load_split_ids(name):
        path = args.splits_dir / f"{TASK_NAME}_{name}.csv"
        if not path.exists():
            # Support legacy filenames during the initial migration.
            path = args.splits_dir / f"{name}.csv"
        if not path.exists():
            raise FileNotFoundError(f"Split file not found for {name}")
        df = pd.read_csv(path)
        return set(df["unique_id"].astype(str))

    train_ids = load_split_ids("train")
    val_ids = load_split_ids("val")
    test_ids = load_split_ids("test")

    train_base = base_df[base_df["unique_id"].isin(train_ids)].copy()
    val_base = base_df[base_df["unique_id"].isin(val_ids)].copy()
    test_base = base_df[base_df["unique_id"].isin(test_ids)].copy()

    # Compute test sequence-redundancy OOD (superseded by the homology-ood-annotation skill).
    print("[OOD] Computing test-to-train+val sequence homology with mmseqs2...")
    train_val_df = pd.concat([train_base, val_base], ignore_index=True)
    query_fasta = args.output_dir / "test_for_mmseqs.fasta"
    target_fasta = args.output_dir / "train_val_for_mmseqs.fasta"
    result_tsv = args.output_dir / "mmseqs_test_vs_trainval.tsv"
    mmseqs_tmp = args.output_dir / "mmseqs_tmp"

    write_fasta(train_val_df, target_fasta)
    write_fasta(test_base, query_fasta)
    if args.mmseqs_cache is not None:
        print(f"[OOD] Using cached mmseqs2 results: {args.mmseqs_cache}")
        max_identity = parse_max_identity(args.mmseqs_cache)
    else:
        max_identity = run_mmseqs2(query_fasta, target_fasta, result_tsv, mmseqs_tmp)

    # OOD columns
    test_lengths = test_base["_protein_length"].values
    test_base["OOD_ExtremeShort"] = test_lengths < 60
    test_base["OOD_ExtremeLong"] = test_lengths > 1000

    for thr in SEQ_THRESHOLDS:
        col = f"seq_Redundancy_{thr}"
        test_base[col] = test_base["unique_id"].apply(
            lambda uid: (max_identity.get(uid, -1.0) < thr) if max_identity.get(uid, -1.0) >= 0 else True
        )

    test_base["OOD_Orphan"] = test_base["unique_id"].apply(
        lambda uid: max_identity.get(uid, -1.0) < 0
    )

    # Default: this task does not exclude test samples beforehand; set all to True.
    test_base["Default"] = True

    # Remove auxiliary columns.
    train_base = train_base.drop(columns=["_protein_length"])
    val_base = val_base.drop(columns=["_protein_length"])
    test_base = test_base.drop(columns=["_protein_length"])

    # Reorder columns.
    base_cols = ["unique_id", "aa_seq", "struct_file", "ligand_smiles", "ligand_ecfp4", "label"]
    ood_cols_ordered = ["Default", "OOD_ExtremeShort", "OOD_ExtremeLong"] + \
                       [f"seq_Redundancy_{thr}" for thr in SEQ_THRESHOLDS] + \
                       ["OOD_Orphan"]
    test_base = test_base[base_cols + ood_cols_ordered]

    # Save
    args.splits_dir.mkdir(parents=True, exist_ok=True)
    train_out_path = args.splits_dir / f"{TASK_NAME}_train.csv"
    val_out_path = args.splits_dir / f"{TASK_NAME}_val.csv"
    test_out_path = args.splits_dir / f"{TASK_NAME}_test.csv"

    train_base.to_csv(train_out_path, index=False)
    val_base.to_csv(val_out_path, index=False)
    test_base.to_csv(test_out_path, index=False)

    # Remove legacy filenames if present.
    for old in ["train.csv", "val.csv", "test.csv"]:
        old_path = args.splits_dir / old
        if old_path.exists():
            old_path.unlink()

    # Remove temporary mmseqs files.
    if mmseqs_tmp.exists():
        shutil.rmtree(mmseqs_tmp)

    print("[Output] Created:")
    print(f"  {train_out_path}: {len(train_base)} rows")
    print(f"  {val_out_path}: {len(val_base)} rows")
    print(f"  {test_out_path}: {len(test_base)} rows")
    print(f"\n[Test OOD statistics]")
    print(f"  Default=True: {test_base['Default'].sum()}")
    print(f"  OOD_ExtremeShort=True: {test_base['OOD_ExtremeShort'].sum()}")
    print(f"  OOD_ExtremeLong=True: {test_base['OOD_ExtremeLong'].sum()}")
    print(f"  OOD_Orphan=True: {test_base['OOD_Orphan'].sum()}")
    for thr in SEQ_THRESHOLDS:
        col = f"seq_Redundancy_{thr}"
        print(f"  {col}=True: {test_base[col].sum()}")


if __name__ == "__main__":
    main()
