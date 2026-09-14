#!/usr/bin/env python3
"""Rebuild final splits so that all extreme-length proteins are in the test set.

- ExtremeShort: protein_length < 60
- ExtremeLong:  protein_length > 1000
- Default:      60 <= protein_length <= 1000

All extreme-length rows from train/val/test are moved to test.
Other OOD annotations are recomputed only on Default test proteins.

NOTE (2026-08 统一整理): 列名已适配统一 schema（labels->label, pdb_file->struct_file,
file_type 列取消）；protein_length 缺失时由 aa_seq 派生；修复了原脚本中未定义的
ROOT 路径（缓存清理/备份均指向本任务 output/）。输出的 ExtremeShort/ExtremeLong、
OOD_NewFold_*、is_time_cutoff 等为历史列名，现已被 OOD_ExtremeShort/OOD_ExtremeLong、
TM-score_* 等统一口径取代（见 DATA_PROCESS.md 阶段四）。
"""
import argparse
import logging
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s')
logger = logging.getLogger(__name__)

TASK_DIR = Path(__file__).resolve().parent.parent  # datasets/ligand_binding_site
DATASET_DIR = TASK_DIR / "splits"
OUTPUT_DIR = TASK_DIR / "output"
OOD_SCRIPT = TASK_DIR / "scripts" / "11_build_ood_splits.py"

CORE_COLS = ["unique_id", "aa_seq", "label", "struct_file"]
SEQ_THRESHOLDS = [90, 80, 70, 60, 50, 40, 30]
TM_THRESHOLDS = [0.9, 0.8, 0.7, 0.6, 0.5]


def is_extreme(length):
    return (length < 60) | (length > 1000)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out_dir", type=str, default=str(DATASET_DIR))
    args = parser.parse_args()

    out_dir = Path(args.out_dir)

    logger.info("Loading current final splits ...")
    train = pd.read_csv(DATASET_DIR / "ligand_binding_site_train.csv", keep_default_na=False, na_values=[''])
    val = pd.read_csv(DATASET_DIR / "ligand_binding_site_val.csv", keep_default_na=False, na_values=[''])
    test = pd.read_csv(DATASET_DIR / "ligand_binding_site_test.csv", keep_default_na=False, na_values=[''])

    for df, name in [(train, 'train'), (val, 'val'), (test, 'test')]:
        if 'protein_length' not in df.columns:
            df['protein_length'] = df['aa_seq'].str.len()
        df['protein_length'] = pd.to_numeric(df['protein_length'], errors='coerce')
        df['is_extreme'] = is_extreme(df['protein_length'])

    # Split into normal and extreme
    train_normal = train[~train['is_extreme']].copy()
    val_normal = val[~val['is_extreme']].copy()
    test_normal = test[~test['is_extreme']].copy()

    train_extreme = train[train['is_extreme']].copy()
    val_extreme = val[val['is_extreme']].copy()
    test_extreme = test[test['is_extreme']].copy()

    # Mark origin for time_cutoff
    train_extreme['_origin'] = 'trainval'
    val_extreme['_origin'] = 'trainval'
    test_extreme['_origin'] = 'test'
    extreme_all = pd.concat([train_extreme, val_extreme, test_extreme], ignore_index=True)

    logger.info(f"Normal lengths -> train={len(train_normal)}, val={len(val_normal)}, test_default={len(test_normal)}")
    logger.info(f"Extreme lengths -> moved to test: {len(extreme_all)} (train {len(train_extreme)}, val {len(val_extreme)}, test {len(test_extreme)})")

    # Recompute OOD only on Default test proteins
    # Force full OOD recomputation for the new Default-only reference.
    for cache in ["mmseqs_max_identity.csv", "foldseek_max_tmscore.csv", "idr_predictions.csv"]:
        cp = OUTPUT_DIR / "ood" / cache
        if cp.exists():
            cp.unlink()
            logger.info(f"Removed stale OOD cache: {cp}")

    # Avoid duplicate is_time_cutoff columns during merge; let OOD script provide it.
    test_normal = test_normal.drop(columns=["is_time_cutoff"], errors="ignore")

    with tempfile.TemporaryDirectory(prefix="lbs_ood_default_", dir=OUTPUT_DIR) as tmp:
        tmp = Path(tmp)
        train_normal.to_csv(tmp / "train_normal.csv", index=False)
        val_normal.to_csv(tmp / "val_normal.csv", index=False)
        test_normal.to_csv(tmp / "test_default.csv", index=False)
        (tmp / "ood_out").mkdir(parents=True, exist_ok=True)

        logger.info("Recomputing OOD for Default test proteins ...")
        subprocess.run(
            [
                sys.executable, str(OOD_SCRIPT),
                "--train", str(tmp / "train_normal.csv"),
                "--val", str(tmp / "val_normal.csv"),
                "--test", str(tmp / "test_default.csv"),
                "--out_dir", str(tmp / "ood_out"),
            ],
            check=True,
        )

        ood_test = pd.read_csv(tmp / "ood_out" / "ligand_binding_site_test.csv", keep_default_na=False, na_values=[''])

    # Collect OOD columns from the recomputed Default test
    ood_cols = ["is_time_cutoff"]
    ood_cols += [f"seq_Redundancy_{thr}" for thr in SEQ_THRESHOLDS]
    ood_cols += [f"OOD_NewFold_{thr:.1f}" for thr in TM_THRESHOLDS]
    ood_cols += ["OOD_Orphan", "idr_ratio", "OOD_IDR"]
    ood_cols = [c for c in ood_cols if c in ood_test.columns]
    ood_df = ood_test[["unique_id"] + ood_cols].copy()

    # Merge OOD onto default test rows
    test_default = test_normal.merge(ood_df, on="unique_id", how="left")
    test_default["Default"] = True
    test_default["ExtremeShort"] = False
    test_default["ExtremeLong"] = False

    # Prepare extreme rows: Default=False；其余 OOD 列与 idr_ratio 直接置 NaN
    # （"未计算"语义，全项目统一惯例，流水线位置写入；
    #   finalize_ood_columns.py --mode nan 仅作兜底校验）
    extreme_all["Default"] = False
    extreme_all["ExtremeShort"] = extreme_all["protein_length"] < 60
    extreme_all["ExtremeLong"] = extreme_all["protein_length"] > 1000
    extreme_all["is_time_cutoff"] = extreme_all["_origin"] == "test"
    for c in ood_cols:
        if c == "is_time_cutoff":
            continue
        extreme_all[c] = pd.NA

    # Combine test
    test_final = pd.concat([test_default, extreme_all], ignore_index=True)

    # Ensure boolean dtype for OOD columns (except idr_ratio)；extreme 行保留 NaN
    dmask = test_final["Default"].astype(bool)
    bool_ood_cols = [c for c in ood_cols if c != "idr_ratio"] + ["Default", "ExtremeShort", "ExtremeLong"]
    for c in bool_ood_cols:
        if c in ("Default", "ExtremeShort", "ExtremeLong"):
            test_final[c] = test_final[c].fillna(False).astype(bool)
        else:
            col = test_final[c].astype(object)
            col.loc[dmask] = col.loc[dmask].fillna(False).astype(bool)
            col.loc[~dmask] = pd.NA
            test_final[c] = col

    # Build final column orders
    test_final_cols = CORE_COLS + ["is_time_cutoff", "Default", "ExtremeShort", "ExtremeLong"]
    test_final_cols += [f"seq_Redundancy_{thr}" for thr in SEQ_THRESHOLDS]
    test_final_cols += [f"OOD_NewFold_{thr:.1f}" for thr in TM_THRESHOLDS]
    test_final_cols += ["OOD_Orphan", "idr_ratio", "OOD_IDR"]
    test_final_cols = [c for c in test_final_cols if c in test_final.columns]

    train_final = train_normal[CORE_COLS]
    val_final = val_normal[CORE_COLS]
    test_final = test_final[test_final_cols]

    # Backup old final files
    backup_dir = OUTPUT_DIR / "intermediate_splits"
    backup_dir.mkdir(parents=True, exist_ok=True)
    for fname in ["ligand_binding_site_train.csv", "ligand_binding_site_val.csv", "ligand_binding_site_test.csv"]:
        src = DATASET_DIR / fname
        if src.exists():
            shutil.move(str(src), str(backup_dir / fname))

    # Write new final files
    train_final.to_csv(out_dir / "ligand_binding_site_train.csv", index=False)
    val_final.to_csv(out_dir / "ligand_binding_site_val.csv", index=False)
    test_final.to_csv(out_dir / "ligand_binding_site_test.csv", index=False)

    logger.info(f"Saved final splits: train={len(train_final)}, val={len(val_final)}, test={len(test_final)}")
    logger.info(f"  test Default={test_final['Default'].sum()}, ExtremeShort={test_final['ExtremeShort'].sum()}, ExtremeLong={test_final['ExtremeLong'].sum()}")

    # Sanity checks
    assert len(train_final) == len(train_normal)
    assert len(val_final) == len(val_normal)
    assert len(test_final) == len(test_default) + len(extreme_all)
    assert test_final["Default"].sum() + test_final["ExtremeShort"].sum() + test_final["ExtremeLong"].sum() == len(test_final)


if __name__ == "__main__":
    main()
