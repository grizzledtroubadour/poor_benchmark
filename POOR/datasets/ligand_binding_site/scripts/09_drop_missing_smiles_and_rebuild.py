#!/usr/bin/env python3
"""Drop samples without ligand_smiles and rebuild final splits.

By default this script only incrementally adds (or refreshes) the
OOD-NewFunction and OOD-LongTail labels.  The expensive mmseqs2 / Foldseek /
metapredict OOD computations are NOT re-run; their cached values are carried
over from the existing test split.

Use --recompute_full_ood only when you explicitly want to regenerate all OOD
columns from scratch (will take a long time).

NOTE (2026-08 统一整理): 列名已适配统一 schema（labels->label, pdb_file->struct_file,
file_type 列取消；ExtremeShort/ExtremeLong -> OOD_ExtremeShort/OOD_ExtremeLong）。
has_ec_annotation / OOD_NewFunction / OOD_NewFold_* / is_time_cutoff 等为历史列名，
现已被 OOD_NewEC_* / TM-score_* 等统一口径取代（由对应 skill 脚本落盘，
见 DATA_PROCESS.md 阶段四）。
"""
import argparse
import ast
import gzip
import logging
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s')
logger = logging.getLogger(__name__)

TASK_DIR = Path(__file__).resolve().parent.parent  # datasets/ligand_binding_site
ROOT = TASK_DIR.parent.parent                      # 项目根（POOR/）
DATASET_DIR = TASK_DIR / "splits"
OOD_SCRIPT = TASK_DIR / "scripts" / "11_build_ood_splits.py"
SIFTS_EC_FILE = ROOT / "data" / "sifts" / "sifts_chain_ec.tsv.gz"

CORE_COLS = ["unique_id", "ligand_smiles", "ligand_ecfp4", "aa_seq", "label", "struct_file"]
SEQ_THRESHOLDS = [90, 80, 70, 60, 50, 40, 30]
TM_THRESHOLDS = [0.9, 0.8, 0.7, 0.6, 0.5]


def chain_key_from_unique_id(unique_id):
    parts = str(unique_id).split("_")
    return f"{parts[0].upper()}_{parts[1]}"


def load_ec_annotations(ec_file):
    if not ec_file.exists():
        logger.warning(f"EC annotation file not found: {ec_file}")
        return {}
    logger.info(f"Loading EC annotations from {ec_file}")
    df = pd.read_csv(ec_file, sep="\t", comment="#", keep_default_na=False, na_values=[""])
    df["PDB"] = df["PDB"].str.upper()
    df["chain_key"] = df["PDB"] + "_" + df["CHAIN"]
    ec_map = (
        df.groupby("chain_key")["EC_NUMBER"]
        .apply(lambda x: set(v for v in x if v and v.strip()))
        .to_dict()
    )
    logger.info(f"Loaded EC annotations for {len(ec_map)} unique chains")
    return ec_map


def binding_site_ratio(labels_str, seq_len):
    try:
        labels = ast.literal_eval(str(labels_str))
    except Exception:
        return 0.0
    if not isinstance(labels, (list, tuple)) or len(labels) == 0:
        return 0.0
    if seq_len <= 0:
        return 0.0
    return float(sum(1 for v in labels if v)) / seq_len


def add_meta_for_ood(df):
    df = df.copy()
    df["pdb_id"] = df["unique_id"].str.split("_").str[0]
    df["chain_id"] = df["unique_id"].str.split("_").str[1]
    return df


def recompute_full_ood(train_clean, val_clean, test_default):
    """Run the expensive 11_build_ood_splits.py pipeline and return OOD df."""
    for cache in ["mmseqs_max_identity.csv", "foldseek_max_tmscore.csv", "idr_predictions.csv"]:
        cp = TASK_DIR / "output" / "ood" / cache
        if cp.exists():
            cp.unlink()
            logger.info(f"Removed stale OOD cache: {cp}")

    test_default_ood_input = test_default.drop(columns=["is_time_cutoff"], errors="ignore")
    test_default_ood_input["protein_length"] = 100

    with tempfile.TemporaryDirectory(prefix="lbs_clean_ood_", dir=TASK_DIR / "output") as tmp:
        tmp = Path(tmp)
        add_meta_for_ood(train_clean).to_csv(tmp / "train_clean.csv", index=False)
        add_meta_for_ood(val_clean).to_csv(tmp / "val_clean.csv", index=False)
        add_meta_for_ood(test_default_ood_input).to_csv(tmp / "test_default.csv", index=False)
        (tmp / "ood_out").mkdir(parents=True, exist_ok=True)

        logger.info("Recomputing OOD for cleaned Default test proteins ...")
        subprocess.run(
            [
                sys.executable, str(OOD_SCRIPT),
                "--train", str(tmp / "train_clean.csv"),
                "--val", str(tmp / "val_clean.csv"),
                "--test", str(tmp / "test_default.csv"),
                "--out_dir", str(tmp / "ood_out"),
            ],
            check=True,
        )
        ood_test = pd.read_csv(
            tmp / "ood_out" / "ligand_binding_site_test.csv",
            keep_default_na=False, na_values=['']
        )

    ood_cols = ["is_time_cutoff"]
    ood_cols += [f"seq_Redundancy_{thr}" for thr in SEQ_THRESHOLDS]
    ood_cols += [f"OOD_NewFold_{thr:.1f}" for thr in TM_THRESHOLDS]
    ood_cols += ["OOD_Orphan", "idr_ratio", "OOD_IDR"]
    ood_cols += ["has_ec_annotation", "OOD_NewFunction", "binding_site_ratio", "OOD_LongTail"]
    ood_cols = [c for c in ood_cols if c in ood_test.columns]
    return ood_test[["unique_id"] + ood_cols].copy()


def add_newfunction_longtail_ood(train_clean, val_clean, test_default):
    """Incrementally compute only OOD-NewFunction and OOD-LongTail."""
    ec_map = load_ec_annotations(SIFTS_EC_FILE)

    trainval_ec = set()
    for df in [train_clean, val_clean]:
        df["chain_key"] = df["unique_id"].apply(chain_key_from_unique_id)
        for ck in df["chain_key"]:
            trainval_ec.update(ec_map.get(ck, set()))
    logger.info(f"Train+Val unique EC numbers: {len(trainval_ec)}")

    test_default = test_default.copy()
    test_default["chain_key"] = test_default["unique_id"].apply(chain_key_from_unique_id)
    test_default["ec_set"] = test_default["chain_key"].map(lambda k: ec_map.get(k, set()))
    test_default["has_ec_annotation"] = test_default["ec_set"].apply(lambda s: len(s) > 0)
    test_default["OOD_NewFunction"] = test_default["ec_set"].apply(
        lambda s: len(s) > 0 and all(ec not in trainval_ec for ec in s)
    )
    test_default["binding_site_ratio"] = test_default.apply(
        lambda r: binding_site_ratio(r["label"], len(r["aa_seq"])), axis=1
    )
    test_default["OOD_LongTail"] = test_default["binding_site_ratio"] < 0.01

    logger.info(
        f"Default test with EC annotation: {test_default['has_ec_annotation'].sum()}; "
        f"OOD_NewFunction: {test_default['OOD_NewFunction'].sum()}; "
        f"OOD_LongTail: {test_default['OOD_LongTail'].sum()}"
    )

    return test_default[["unique_id", "has_ec_annotation", "OOD_NewFunction", "binding_site_ratio", "OOD_LongTail"]].copy()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out_dir", type=str, default=str(DATASET_DIR))
    parser.add_argument("--recompute_full_ood", action="store_true",
                        help="Regenerate all OOD columns via 11_build_ood_splits.py (slow).")
    args = parser.parse_args()
    out_dir = Path(args.out_dir)

    logger.info("Loading current final splits ...")
    train = pd.read_csv(DATASET_DIR / "ligand_binding_site_train.csv", keep_default_na=False, na_values=[''])
    val = pd.read_csv(DATASET_DIR / "ligand_binding_site_val.csv", keep_default_na=False, na_values=[''])
    test = pd.read_csv(DATASET_DIR / "ligand_binding_site_test.csv", keep_default_na=False, na_values=[''])

    # Boolean flags are stored as strings in the CSVs; normalize to real bools.
    for b in ("Default", "OOD_ExtremeShort", "OOD_ExtremeLong"):
        if b in test.columns and test[b].dtype == object:
            test[b] = test[b].map(lambda v: str(v) == "True")

    # Drop rows with missing ligand_smiles
    train_clean = train[train['ligand_smiles'].ne('') & train['ligand_smiles'].notna()].copy()
    val_clean = val[val['ligand_smiles'].ne('') & val['ligand_smiles'].notna()].copy()
    test_clean = test[test['ligand_smiles'].ne('') & test['ligand_smiles'].notna()].copy()

    logger.info(f"After dropping missing SMILES: train={len(train_clean)}, val={len(val_clean)}, test={len(test_clean)}")

    test_default = test_clean[test_clean['Default']].copy()
    test_extreme = test_clean[~test_clean['Default']].copy()
    logger.info(f"Test Default={len(test_default)}, ExtremeShort={test_extreme['OOD_ExtremeShort'].sum()}, ExtremeLong={test_extreme['OOD_ExtremeLong'].sum()}")

    if args.recompute_full_ood:
        ood_df = recompute_full_ood(train_clean, val_clean, test_default)
        test_default = test_default.drop(columns=[c for c in ood_df.columns if c != "unique_id"], errors="ignore")
        test_default = test_default.merge(ood_df, on="unique_id", how="left")
    else:
        # Refresh only the new OOD columns; all other OOD columns are carried over.
        new_ood_df = add_newfunction_longtail_ood(train_clean, val_clean, test_default)
        test_default = test_default.drop(
            columns=[c for c in new_ood_df.columns if c != "unique_id"],
            errors="ignore"
        )
        test_default = test_default.merge(new_ood_df, on="unique_id", how="left")

    # Extreme rows: other OOD columns False / 0.0
    ood_cols = ["is_time_cutoff"]
    ood_cols += [f"seq_Redundancy_{thr}" for thr in SEQ_THRESHOLDS]
    ood_cols += [f"OOD_NewFold_{thr:.1f}" for thr in TM_THRESHOLDS]
    ood_cols += ["OOD_Orphan", "idr_ratio", "OOD_IDR"]
    ood_cols += ["has_ec_annotation", "OOD_NewFunction", "binding_site_ratio", "OOD_LongTail"]
    ood_cols = [c for c in ood_cols if c in test_default.columns]

    for c in ood_cols:
        if c == "is_time_cutoff":
            continue
        if c in ("idr_ratio", "binding_site_ratio"):
            test_extreme[c] = 0.0
        else:
            test_extreme[c] = False

    test_final = pd.concat([test_default, test_extreme], ignore_index=True)

    bool_ood_cols = [c for c in ood_cols if c not in ("idr_ratio", "binding_site_ratio")] + ["Default", "OOD_ExtremeShort", "OOD_ExtremeLong"]
    for c in bool_ood_cols:
        test_final[c] = test_final[c].fillna(False).astype(bool)

    test_final_cols = CORE_COLS + ["is_time_cutoff", "Default", "OOD_ExtremeShort", "OOD_ExtremeLong"]
    test_final_cols += [f"seq_Redundancy_{thr}" for thr in SEQ_THRESHOLDS]
    test_final_cols += [f"OOD_NewFold_{thr:.1f}" for thr in TM_THRESHOLDS]
    test_final_cols += ["OOD_Orphan", "idr_ratio", "OOD_IDR"]
    test_final_cols += ["has_ec_annotation", "OOD_NewFunction", "binding_site_ratio", "OOD_LongTail"]
    test_final_cols = [c for c in test_final_cols if c in test_final.columns]

    train_final = train_clean[CORE_COLS]
    val_final = val_clean[CORE_COLS]
    test_final = test_final[test_final_cols]

    # Backup current final files
    backup_dir = TASK_DIR / "output" / "intermediate_splits"
    backup_dir.mkdir(parents=True, exist_ok=True)
    for fname in ["ligand_binding_site_train.csv", "ligand_binding_site_val.csv", "ligand_binding_site_test.csv"]:
        src = DATASET_DIR / fname
        if src.exists():
            shutil.move(str(src), str(backup_dir / fname))

    train_final.to_csv(out_dir / "ligand_binding_site_train.csv", index=False)
    val_final.to_csv(out_dir / "ligand_binding_site_val.csv", index=False)
    test_final.to_csv(out_dir / "ligand_binding_site_test.csv", index=False)

    logger.info(f"Saved final cleaned splits: train={len(train_final)}, val={len(val_final)}, test={len(test_final)}")
    logger.info(f"  test Default={test_final['Default'].sum()}, ExtremeShort={test_final['OOD_ExtremeShort'].sum()}, ExtremeLong={test_final['OOD_ExtremeLong'].sum()}")


if __name__ == "__main__":
    main()
