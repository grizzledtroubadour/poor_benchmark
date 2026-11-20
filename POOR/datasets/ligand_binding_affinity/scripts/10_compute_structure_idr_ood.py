#!/usr/bin/env python3
"""
Script: 10_compute_structure_idr_ood.py
Purpose: Compute test structural OOD using Foldseek-Multimer complex-level symmetric TM-score
      and intrinsic-disorder OOD using metapredict with threshold 0.1.

Output columns (appended/reordered after step 09):
  - TM-score_0.9, ..., TM-score_0.3
  - idr_ratio
  - OOD_IDR
  (Default / OOD_ExtremeShort / OOD_ExtremeLong / seq_Redundancy_* / OOD_Orphan come from 09.)

Structural OOD definition (since 2026-08-29): foldseek `easy-multimersearch` complex alignment
(Foldseek-Multimer, Kim et al., Nat Methods 2025), taking the maximum across complex pairs
of their symmetric mean complexqtmscore/complexttmscore. Complex-level alignment matches
the complex-ligand affinity task. The previous single-chain `search -a` + qtmscore/ttmscore
symmetric mean was replaced on 2026-08-29 (backup: output/backup_before_multimer_tm_unify/;
comparison: analysis/output/foldseek_easysearch_audit/lba/multimer/).
Canonical IDR tool: .skills/ood-annotation-toolkit/idr_ood.py (threshold 0.1 for this task).

Usage:
    cd ligand_binding_affinity/
    python scripts/10_compute_structure_idr_ood.py
"""

import argparse
import json
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import metapredict

TASK_NAME = "ligand_binding_affinity"
SEQ_THRESHOLDS = [90, 80, 70, 60, 50, 40, 30]
TM_THRESHOLDS = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3]
IDR_THRESHOLD = 0.1


def write_file_list(file_paths, out_path: Path):
    with open(out_path, "w", encoding="utf-8") as f:
        for p in file_paths:
            f.write(f"{Path(p).resolve()}\n")


def run_foldseek(trainval_list: Path, test_list: Path, out_dir: Path):
    """Run Foldseek-Multimer easy-multimersearch at complex level and return the result TSV."""
    out_dir.mkdir(parents=True, exist_ok=True)
    result_tsv = out_dir / "foldseek_multimer_test_vs_trainval.tsv"
    tmp_dir = out_dir / "tmp_multimer"
    q_dir = out_dir / "query_structures"
    t_dir = out_dir / "target_structures"
    # multimersearch requires structure directories; use symlinks to avoid copies.
    for link_dir, list_file in ((q_dir, test_list), (t_dir, trainval_list)):
        link_dir.mkdir(exist_ok=True)
        with open(list_file) as f:
            for line in f:
                src = Path(line.strip())
                if not src.name:
                    continue
                dst = link_dir / src.name
                if not dst.exists():
                    dst.symlink_to(src)
    subprocess.run(
        ["foldseek", "easy-multimersearch", str(q_dir), str(t_dir),
         str(result_tsv), str(tmp_dir),
         "--format-output", "query,target,complexqtmscore,complexttmscore"],
        check=True,
    )
    return result_tsv


def complex_stem(chain_name: str, known_stems: set):
    """foldseek createdb splits multichain PDBs into <file>_<chain> records; single-chain
    files retain bare filenames. Normalize m8 query/target names back to complex stems."""
    if chain_name in known_stems:
        return chain_name
    return chain_name.rsplit("_", 1)[0]


def compute_max_tmscore(result_tsv: Path, test_stems: set, ref_stems: set):
    """Aggregate maximum symmetric complex TM-score (mean of complexqtm/complexttm) per test complex."""
    max_tm = {}
    with open(result_tsv, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            q_pdb = complex_stem(parts[0], test_stems)
            t_pdb = complex_stem(parts[1], ref_stems)
            if q_pdb == t_pdb:
                continue
            sym_score = (float(parts[2]) + float(parts[3])) / 2.0
            max_tm[q_pdb] = max(max_tm.get(q_pdb, 0.0), sym_score)
    return max_tm


def main():
    parser = argparse.ArgumentParser(description="Compute structure and IDR OOD for test set")
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
        "--pdbs-dir",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "pdbs",
        help="Directory containing whole-protein PDB files referenced by struct_file",
    )
    args = parser.parse_args()

    train_df = pd.read_csv(args.splits_dir / f"{TASK_NAME}_train.csv")
    val_df = pd.read_csv(args.splits_dir / f"{TASK_NAME}_val.csv")
    test_df = pd.read_csv(args.splits_dir / f"{TASK_NAME}_test.csv")

    print(f"[Input] train: {len(train_df)}, val: {len(val_df)}, test: {len(test_df)}")

    # Build Foldseek file lists by resolving struct_file filenames against pdbs/.
    fs_dir = args.output_dir / "foldseek_ood"
    fs_dir.mkdir(parents=True, exist_ok=True)
    trainval_list = fs_dir / "trainval_files.tsv"
    test_list = fs_dir / "test_files.tsv"

    trainval_files = [args.pdbs_dir / f for f in pd.concat([train_df, val_df])["struct_file"]]
    test_files = [args.pdbs_dir / f for f in test_df["struct_file"]]
    write_file_list(trainval_files, trainval_list)
    write_file_list(test_files, test_list)

    # Run Foldseek-Multimer.
    print("[Foldseek-Multimer] Computing test-to-train+val complex structural similarity...")
    result_tsv = run_foldseek(trainval_list, test_list, fs_dir)
    test_stems = {Path(f).stem for f in test_df["struct_file"]}
    ref_stems = {Path(f).stem for f in pd.concat([train_df, val_df])["struct_file"]}
    max_tm = compute_max_tmscore(result_tsv, test_stems, ref_stems)

    # Add TM-score OOD columns.
    for thr in TM_THRESHOLDS:
        col = f"TM-score_{thr:.1f}"
        test_df[col] = test_df["unique_id"].apply(
            lambda uid: max_tm.get(uid, 0.0) < thr
        )

    # Compute IDR in batches.
    print("[metapredict] Computing test sequence disorder fractions...")
    # metapredict accepts standard amino-acid sequences only; remove chain separators '|'.
    seqs = test_df["aa_seq"].fillna("").str.replace("|", "").tolist()
    batch_size = 500
    idr_ratios = []
    for i in range(0, len(seqs), batch_size):
        batch = seqs[i:i + batch_size]
        outs = metapredict.predict_disorder_batch(batch, show_progress_bar=False)
        for out in outs:
            scores = np.asarray(out[1])
            mask = scores > 0.5
            if not mask.any():
                idr_ratios.append(0.0)
                continue
            diff = np.diff(np.concatenate([[0], mask.astype(int), [0]]))
            starts = np.where(diff == 1)[0]
            ends = np.where(diff == -1)[0]
            idr_ratios.append(float(sum(ends - starts) / len(scores)))
        if (i + batch_size) % 500 == 0 or (i + batch_size) >= len(seqs):
            print(f"  Completed {min(i + batch_size, len(seqs))}/{len(seqs)}")

    test_df["idr_ratio"] = idr_ratios
    test_df["OOD_IDR"] = test_df["idr_ratio"] > IDR_THRESHOLD

    # Default: set True for every regular test sample.
    test_df["Default"] = True

    # Extreme length: amino-acid count excluding chain separators '|'.
    test_lengths = test_df["aa_seq"].str.replace("|", "").str.len().values
    test_df["OOD_ExtremeShort"] = test_lengths < 60
    test_df["OOD_ExtremeLong"] = test_lengths > 1000

    # Sequence-redundancy OOD: reuse step 09 outputs if available.
    if "OOD_Orphan" not in test_df.columns:
        # If not yet computed, fail here; normal split files already contain these columns.
        raise RuntimeError("Run 09_format_splits_to_skill.py first to generate sequence-homology OOD columns")

    # Order known columns by the unified schema and retain later-added columns unchanged,
    # preventing loss of downstream OOD annotations when rerunning on complete splits.
    base_cols = ["unique_id", "aa_seq", "struct_file", "ligand_smiles", "ligand_ecfp4", "label"]
    ood_cols_ordered = (
        ["Default", "OOD_ExtremeShort", "OOD_ExtremeLong"]
        + [f"seq_Redundancy_{thr}" for thr in SEQ_THRESHOLDS]
        + [f"TM-score_{thr:.1f}" for thr in TM_THRESHOLDS]
        + ["idr_ratio", "OOD_IDR", "OOD_Orphan"]
    )
    known_cols = base_cols + ood_cols_ordered
    extra_cols = [c for c in test_df.columns if c not in known_cols]
    test_df = test_df[known_cols + extra_cols]

    # Save.
    test_out_path = args.splits_dir / f"{TASK_NAME}_test.csv"
    test_df.to_csv(test_out_path, index=False)

    # Update statistics JSON.
    stats_path = args.output_dir / "cluster_split_stats.json"
    if stats_path.exists():
        with open(stats_path, "r", encoding="utf-8") as f:
            stats = json.load(f)
    else:
        stats = {}

    tm_counts = {f"TM-score_{thr:.1f}": int(test_df[f"TM-score_{thr:.1f}"].sum()) for thr in TM_THRESHOLDS}
    stats["OOD_stats"] = {
        "idr_threshold": IDR_THRESHOLD,
        "idr_ratio": test_df["idr_ratio"].describe().to_dict(),
        "OOD_IDR": int(test_df["OOD_IDR"].sum()),
        "TM_score_counts": tm_counts,
        "Default": int(test_df["Default"].sum()),
    }
    with open(stats_path, "w", encoding="utf-8") as f:
        json.dump(stats, f, indent=2)

    # Remove temporary Foldseek files.
    if fs_dir.exists():
        shutil.rmtree(fs_dir)

    print(f"\n[Output] Updated {test_out_path}")
    print(f"[OOD] OOD_IDR=True (idr>{IDR_THRESHOLD}): {test_df['OOD_IDR'].sum()}")
    for thr in TM_THRESHOLDS:
        col = f"TM-score_{thr:.1f}"
        print(f"[OOD] {col}=True: {test_df[col].sum()}")
    print(f"[OOD] Default=True: {test_df['Default'].sum()}")


if __name__ == "__main__":
    main()
