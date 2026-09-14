#!/usr/bin/env python3
"""
解压并校验 EnzyBase12k 结构文件
- 解压 DATASET.zip (pdb/pqr) 与 Dataset_mmcifs.zip (cif)
- 统计文件数量与格式
- 将结构文件名与 metadata.csv 中的 uniprot_id / pdb_id_final 进行匹配
- 输出缺失结构文件的清单
"""
import os
import re
import zipfile
import json
from pathlib import Path
from collections import Counter

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
OUTPUT_DIR = ROOT / "output"
LOG_DIR = ROOT / "logs"
for d in [OUTPUT_DIR, LOG_DIR]:
    d.mkdir(parents=True, exist_ok=True)

META_PATH = DATA_DIR / "metadata.csv"
DATASET_ZIP = DATA_DIR / "DATASET.zip"
MMCIF_ZIP = DATA_DIR / "Dataset_mmcifs.zip"

DATASET_EXTRACT_DIR = DATA_DIR / "dataset"
MMCIF_EXTRACT_DIR = DATA_DIR / "dataset_mmcifs"


def unzip_file(zip_path, extract_to):
    print(f"[INFO] Unzipping {zip_path} -> {extract_to}")
    with zipfile.ZipFile(zip_path, 'r') as z:
        z.extractall(extract_to)
    print(f"[INFO] Unzipped {zip_path}")


def collect_files(directory, extensions):
    """递归收集指定扩展名的文件"""
    files = []
    for ext in extensions:
        files.extend(directory.rglob(f"*{ext}"))
    return files


def extract_id_from_filename(filename):
    """从文件名提取 ID（假设为 uniprot_id 或 pdb_id）"""
    base = filename.stem
    # 移除常见后缀
    base = re.sub(r'_(?:AF|pdb|model|ranked)_?\d*$', '', base, flags=re.IGNORECASE)
    base = re.sub(r'_[A-Za-z]$', '', base)  # 移除链标识
    return base.upper()


def main():
    # 1. 加载元数据
    print(f"[INFO] Loading metadata: {META_PATH}")
    df = pd.read_csv(META_PATH)
    n_meta = len(df)
    print(f"[INFO] Metadata rows: {n_meta}")

    # 2. 解压
    if not DATASET_EXTRACT_DIR.exists() or not any(DATASET_EXTRACT_DIR.iterdir()):
        if DATASET_ZIP.exists():
            unzip_file(DATASET_ZIP, DATASET_EXTRACT_DIR)
        else:
            print(f"[WARN] {DATASET_ZIP} not found, skipping.")
    else:
        print(f"[INFO] {DATASET_EXTRACT_DIR} already exists, skipping unzip.")

    if not MMCIF_EXTRACT_DIR.exists() or not any(MMCIF_EXTRACT_DIR.iterdir()):
        if MMCIF_ZIP.exists():
            unzip_file(MMCIF_ZIP, MMCIF_EXTRACT_DIR)
        else:
            print(f"[WARN] {MMCIF_ZIP} not found, skipping.")
    else:
        print(f"[INFO] {MMCIF_EXTRACT_DIR} already exists, skipping unzip.")

    # 3. 收集结构文件
    pdb_files = collect_files(DATASET_EXTRACT_DIR, [".pdb"]) if DATASET_EXTRACT_DIR.exists() else []
    pqr_files = collect_files(DATASET_EXTRACT_DIR, [".pqr"]) if DATASET_EXTRACT_DIR.exists() else []
    cif_files = collect_files(MMCIF_EXTRACT_DIR, [".cif"]) if MMCIF_EXTRACT_DIR.exists() else []

    print(f"[INFO] Found {len(pdb_files)} .pdb files")
    print(f"[INFO] Found {len(pqr_files)} .pqr files")
    print(f"[INFO] Found {len(cif_files)} .cif files")

    # 4. 文件名统计
    file_stats = {
        "n_metadata_rows": n_meta,
        "n_pdb_files": len(pdb_files),
        "n_pqr_files": len(pqr_files),
        "n_cif_files": len(cif_files),
        "total_structure_files": len(pdb_files) + len(pqr_files) + len(cif_files),
    }

    # 5. 与 metadata 匹配
    # 构建 metadata ID 集合
    uniprot_set = set(df["uniprot_id"].dropna().astype(str).str.upper())
    pdb_set = set(df["pdb_id_final"].dropna().astype(str).str.upper())

    # 提取文件名 ID
    pdb_ids = [extract_id_from_filename(f) for f in pdb_files]
    pqr_ids = [extract_id_from_filename(f) for f in pqr_files]
    cif_ids = [extract_id_from_filename(f) for f in cif_files]

    pdb_uniprot_matched = sum(1 for i in pdb_ids if i in uniprot_set)
    pdb_pdb_matched = sum(1 for i in pdb_ids if i in pdb_set)
    pqr_uniprot_matched = sum(1 for i in pqr_ids if i in uniprot_set)
    pqr_pdb_matched = sum(1 for i in pqr_ids if i in pdb_set)
    cif_uniprot_matched = sum(1 for i in cif_ids if i in uniprot_set)
    cif_pdb_matched = sum(1 for i in cif_ids if i in pdb_set)

    match_stats = {
        "pdb_files_matching_uniprot": pdb_uniprot_matched,
        "pdb_files_matching_pdb_id": pdb_pdb_matched,
        "pdb_files_unmatched": len(pdb_files) - max(pdb_uniprot_matched, pdb_pdb_matched),
        "pqr_files_matching_uniprot": pqr_uniprot_matched,
        "pqr_files_matching_pdb_id": pqr_pdb_matched,
        "pqr_files_unmatched": len(pqr_files) - max(pqr_uniprot_matched, pqr_pdb_matched),
        "cif_files_matching_uniprot": cif_uniprot_matched,
        "cif_files_matching_pdb_id": cif_pdb_matched,
        "cif_files_unmatched": len(cif_files) - max(cif_uniprot_matched, cif_pdb_matched),
    }

    # 6. 缺失文件分析
    # 以 uniprot_id 为主键，检查哪些 metadata 行缺少对应的结构文件
    all_structure_ids = set(pdb_ids + pqr_ids + cif_ids)
    missing_uniprots = sorted(uniprot_set - all_structure_ids)

    missing_by_source = {}
    for source, group in df.groupby("structure_source"):
        source_uniprots = set(group["uniprot_id"].dropna().astype(str).str.upper())
        missing_by_source[source] = {
            "total": len(source_uniprots),
            "missing": len(source_uniprots - all_structure_ids),
            "missing_ratio": round(len(source_uniprots - all_structure_ids) / len(source_uniprots), 4),
        }

    missing_stats = {
        "total_metadata_uniprots": len(uniprot_set),
        "uniprots_with_any_structure_file": len(uniprot_set & all_structure_ids),
        "uniprots_missing_structure_file": len(missing_uniprots),
        "missing_uniprots_sample": missing_uniprots[:50],
        "missing_by_structure_source": missing_by_source,
    }

    # 7. 扩展名分布
    ext_counter = Counter()
    for f in pdb_files + pqr_files + cif_files:
        ext_counter[f.suffix.lower()] += 1

    # 8. 保存报告
    report = {
        "file_stats": file_stats,
        "match_stats": match_stats,
        "missing_stats": missing_stats,
        "extension_counts": dict(ext_counter),
    }

    report_path = OUTPUT_DIR / "structure_verification.json"
    with open(report_path, "w") as f:
        json.dump(report, f, indent=2)
    print(f"[INFO] Verification report saved to {report_path}")

    # 9. 简要输出
    print("\n===== Verification Summary =====")
    for k, v in file_stats.items():
        print(f"  {k}: {v}")
    print("\n-- Match stats --")
    for k, v in match_stats.items():
        print(f"  {k}: {v}")
    print("\n-- Missing stats --")
    for k, v in missing_stats.items():
        if k != "missing_uniprots_sample":
            print(f"  {k}: {v}")
    print(f"  missing_uniprots_sample: {missing_uniprots[:20]}")


if __name__ == "__main__":
    main()
