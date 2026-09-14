#!/usr/bin/env python3
"""
准备标准格式的 splits/ 与 pdbs/
1. 合并 test.csv 与 ood_extremelength.csv
2. 为 test 添加 Default / OOD_ExtremeLong / OOD_ExtremeShort 列
3. 整理结构文件到 pdbs/ 目录，统一使用 .cif 格式（若无则回退 .pdb）
4. 从结构文件中提取 aa_seq 替换 uniprot_seq_cut
5. 生成统一 schema 列：unique_id, aa_seq, struct_file, label
   （历史注记：初版产出为 unique_id, file_type, aa_seq, labels，后经全项目
   统一 schema 迁移为当前四列；本脚本现已直接产出统一 schema）
"""
import argparse
import os
import re
import shutil
import json
import warnings
from pathlib import Path
from multiprocessing import Pool, cpu_count

import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
SPLITS_DIR = ROOT / "splits"
PDBS_DIR = ROOT / "pdbs"
OUTPUT_DIR = ROOT / "output"
PDBS_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

META_PATH = DATA_DIR / "metadata.csv"

PDB_DIR = DATA_DIR / "dataset" / "DATASET" / "EnzyBase12k_pbs"
CIF_DIR = DATA_DIR / "dataset_mmcifs" / "EnzyBase12k_mmcifs"

AA3_TO_1 = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
    "MSE": "M", "SEC": "U", "PYL": "O",
}


def clean_seq(seq):
    return re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", seq.upper())


def extract_seq_from_cif_fast(cif_path):
    """快速从 mmCIF 文件中提取 CA 残基序列"""
    residues = []
    seen = set()
    try:
        with open(cif_path) as f:
            for line in f:
                if not line.startswith("ATOM"):
                    continue
                parts = line.split()
                if len(parts) < 9:
                    continue
                atom_id = parts[3]  # label_atom_id
                if atom_id != "CA":
                    continue
                resname = parts[5]  # label_comp_id
                seq_id = parts[8]   # label_seq_id
                key = seq_id
                if key in seen:
                    continue
                seen.add(key)
                if resname in AA3_TO_1:
                    residues.append(AA3_TO_1[resname])
    except Exception:
        return None
    return "".join(residues) if residues else None


def extract_seq_from_pdb_fast(pdb_path):
    """快速从 PDB 文件中提取 CA 残基序列"""
    residues = []
    seen = set()
    try:
        with open(pdb_path) as f:
            for line in f:
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    continue
                if len(line) < 54:
                    continue
                atom_name = line[12:16].strip()
                if atom_name != "CA":
                    continue
                resname = line[17:20].strip()
                seq_id = line[22:26].strip()
                key = seq_id
                if key in seen:
                    continue
                seen.add(key)
                if resname in AA3_TO_1:
                    residues.append(AA3_TO_1[resname])
    except Exception:
        return None
    return "".join(residues) if residues else None


def find_structure_file(uniprot_id):
    """优先返回 .cif 文件路径，否则返回 .pdb"""
    cif_path = CIF_DIR / f"{uniprot_id}.cif"
    pdb_path = PDB_DIR / f"{uniprot_id}.pdb"
    if cif_path.exists():
        return cif_path, "cif"
    elif pdb_path.exists():
        return pdb_path, "pdb"
    return None, None


def process_one_sample(args):
    """处理单个样本：复制结构文件并提取序列"""
    row_dict, is_test, pdbs_dir = args
    uniprot_id = row_dict["uniprot_id"]
    struct_path, file_type = find_structure_file(uniprot_id)

    if struct_path is None:
        return None, f"file_not_found:{uniprot_id}"

    # 提取序列
    if file_type == "cif":
        struct_seq = extract_seq_from_cif_fast(struct_path)
    else:
        struct_seq = extract_seq_from_pdb_fast(struct_path)

    if struct_seq is None or len(clean_seq(struct_seq)) == 0:
        return None, f"seq_extract_failed:{uniprot_id}"

    # 复制结构文件到 pdbs/
    target_name = f"{uniprot_id}.{file_type}"
    target_path = pdbs_dir / target_name
    if not target_path.exists():
        shutil.copy2(struct_path, target_path)

    record = {
        "unique_id": uniprot_id,
        "aa_seq": struct_seq,
        "struct_file": target_name,
        "label": str(row_dict["ph_optimum"]),
    }

    if is_test:
        record["Default"] = row_dict["Default"]
        record["OOD_ExtremeLong"] = row_dict.get("OOD_ExtremeLong", False)
        record["OOD_ExtremeShort"] = row_dict.get("OOD_ExtremeShort", False)

    return record, "ok"


def prepare_split_parallel(df, split_name, pdbs_dir, is_test=False, n_workers=16):
    """并行准备单个划分文件"""
    print(f"[INFO] Preparing {split_name} with {n_workers} workers...")
    tasks = [(row.to_dict(), is_test, pdbs_dir) for _, row in df.iterrows()]

    results = []
    failed = []
    with Pool(n_workers) as pool:
        for i, (record, status) in enumerate(pool.imap_unordered(process_one_sample, tasks), 1):
            if record is not None:
                results.append(record)
            else:
                failed.append(status)
            if i % 1000 == 0:
                print(f"[INFO] {split_name}: processed {i}/{len(tasks)}")

    if failed:
        print(f"[WARN] {split_name}: {len(failed)} samples failed")
        for f in failed[:10]:
            print(f"  {f}")

    return pd.DataFrame(results)


def main():
    parser = argparse.ArgumentParser(description="准备统一 schema 的 splits/ 与 pdbs/")
    parser.add_argument("--input-dir", type=Path, default=SPLITS_DIR,
                        help="02 初分产物（train/val/test/ood_extremelength.csv）所在目录")
    parser.add_argument("--splits-dir", type=Path, default=SPLITS_DIR,
                        help="最终 splits 输出目录")
    parser.add_argument("--pdbs-dir", type=Path, default=PDBS_DIR,
                        help="结构文件输出目录")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR,
                        help="split_mapping.json 输出目录")
    args = parser.parse_args()

    input_dir = args.input_dir
    splits_dir = args.splits_dir
    pdbs_dir = args.pdbs_dir
    splits_dir.mkdir(parents=True, exist_ok=True)
    pdbs_dir.mkdir(parents=True, exist_ok=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    train_in = input_dir / "train.csv"
    val_in = input_dir / "val.csv"
    test_in = input_dir / "test.csv"
    extreme_in = input_dir / "ood_extremelength.csv"

    print("[INFO] Loading splits...")
    train_df = pd.read_csv(train_in)
    val_df = pd.read_csv(val_in)
    test_df = pd.read_csv(test_in)
    extreme_df = pd.read_csv(extreme_in)

    print(f"[INFO] Train: {len(train_df)}, Val: {len(val_df)}, Test: {len(test_df)}, Extreme: {len(extreme_df)}")

    # 1. 合并 test 与 extreme
    test_df["Default"] = True
    test_df["OOD_ExtremeLong"] = False
    test_df["OOD_ExtremeShort"] = False

    extreme_df["Default"] = False
    extreme_df["OOD_ExtremeLong"] = extreme_df["seq_length"] > 1000
    extreme_df["OOD_ExtremeShort"] = extreme_df["seq_length"] < 60

    combined_test = pd.concat([test_df, extreme_df], ignore_index=True)
    print(f"[INFO] Combined test: {len(combined_test)}")
    print(f"[INFO] Default: {combined_test['Default'].sum()}, ExtremeLong: {combined_test['OOD_ExtremeLong'].sum()}, ExtremeShort: {combined_test['OOD_ExtremeShort'].sum()}")

    # 2. 并行准备 train / val / test
    train_out = prepare_split_parallel(train_df, "train", pdbs_dir, is_test=False)
    val_out = prepare_split_parallel(val_df, "val", pdbs_dir, is_test=False)
    test_out = prepare_split_parallel(combined_test, "test", pdbs_dir, is_test=True)

    print(f"[INFO] Train prepared: {len(train_out)}")
    print(f"[INFO] Val prepared: {len(val_out)}")
    print(f"[INFO] Test prepared: {len(test_out)}")

    # 3. 保存为最终标准文件名
    train_out.to_csv(splits_dir / "optimal_ph_prediction_train.csv", index=False)
    val_out.to_csv(splits_dir / "optimal_ph_prediction_val.csv", index=False)
    test_out.to_csv(splits_dir / "optimal_ph_prediction_test.csv", index=False)

    # 仅在正常流水线运行（输入即默认 splits/ 目录）时删除初分中间文件
    if input_dir.resolve() == SPLITS_DIR.resolve():
        extreme_in.unlink()
        train_in.unlink()
        val_in.unlink()
        test_in.unlink()

    print(f"[INFO] Saved final splits to {splits_dir}")
    print(f"[INFO] Copied {len(list(pdbs_dir.glob('*')))} structure files to {pdbs_dir}")

    # 4. 保存元数据映射
    mapping = {
        "n_train": len(train_out),
        "n_val": len(val_out),
        "n_test": len(test_out),
        "test_default_ids": test_out[test_out["Default"] == True]["unique_id"].tolist(),
    }
    with open(args.output_dir / "split_mapping.json", "w") as f:
        json.dump(mapping, f, indent=2)

    print("[INFO] Done")


if __name__ == "__main__":
    main()
