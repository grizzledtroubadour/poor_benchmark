#!/usr/bin/env python3
"""
脚本：06_build_dataset_csv.py
用途：由 05 的序列产物生成样本级数据集 output/pdbbind_v2020_dataset.csv
      （Stage 2 划分与 OOD 流程的输入）。

      本脚本替代已归档的 08_extract_protein_chains.py 的链级结构路线：
      不再逐链写出独立 PDB 文件（整蛋白结构由 07_copy_whole_protein_pdbs.py 落盘），
      只保留序列提取成功的样本并补充结构文件路径列。
      chain_extraction_success / chain_extraction_error 两列为兼容历史下游
      （08_cluster_time_split.py 等）保留，取值与序列提取成功标记一致。

使用：
    cd ligand_binding_affinity/
    python scripts/06_build_dataset_csv.py

输入：
    output/pdbbind_v2020_metadata_seq.csv

输出：
    output/pdbbind_v2020_dataset.csv
    output/protein_chain_extraction_stats.json
"""

import argparse
import json
from pathlib import Path

import pandas as pd

# 输出列顺序与历史 dataset.csv 保持一致（见 docstring）
FINAL_COLUMNS = [
    "pdb_id", "resolution", "release_year", "affinity_type", "affinity_operator",
    "affinity_value", "pK", "ligand_name", "reference", "year_range",
    "protein_file", "ligand_mol2_file", "ligand_sdf_file", "pocket_file",
    "incomplete_ligand", "covalent_complex", "different_protein", "different_ligand",
    "uncommon_element", "disulfide_bond", "comment", "raw",
    "ligand_parse_success", "ligand_parse_error", "ligand_smiles", "ligand_ecfp4",
    "ligand_num_heavy_atoms",
    "protein_sequence_extract_success", "protein_sequence_error",
    "protein_sequence", "protein_length",
    "chain_extraction_success", "chain_extraction_error",
    "protein_chains", "num_protein_chains", "protein_pdb_file",
]


def main():
    parser = argparse.ArgumentParser(description="Build sample-level dataset CSV from sequence metadata")
    parser.add_argument(
        "--input",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "output" / "pdbbind_v2020_metadata_seq.csv",
        help="Input metadata CSV with per-chain protein sequences (output of 05)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "output" / "pdbbind_v2020_dataset.csv",
        help="Output sample-level dataset CSV",
    )
    parser.add_argument(
        "--stats",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "output" / "protein_chain_extraction_stats.json",
        help="Output stats JSON",
    )
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    total = len(df)
    print(f"[输入] 清洗并解析配体/序列后的样本数: {total}")

    # 仅保留序列提取成功的样本作为最终数据集
    final_df = df[df["protein_sequence_extract_success"] == True].copy()
    failure_df = df[df["protein_sequence_extract_success"] != True].copy()

    # 兼容列：历史链级路线的成功标记（取值与序列提取一致）
    final_df["chain_extraction_success"] = final_df["protein_sequence_extract_success"]
    final_df["chain_extraction_error"] = final_df["protein_sequence_error"]

    # 整蛋白结构文件路径（由 07 落盘）
    final_df["protein_pdb_file"] = final_df["pdb_id"].apply(lambda x: f"pdbs/{x}.pdb")

    final_df = final_df[FINAL_COLUMNS]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    final_df.to_csv(args.output, index=False)
    if not failure_df.empty:
        failure_path = args.output.parent / "protein_chain_extraction_failures.csv"
        failure_df.to_csv(failure_path, index=False)
        print(f"[输出] 提取失败样本 CSV: {failure_path}")

    stats = {
        "total_input": total,
        "chain_extraction_success": int(len(final_df)),
        "chain_extraction_failure": int(len(failure_df)),
        "retention_rate": round(len(final_df) / total, 4) if total > 0 else 0.0,
        "total_chain_files": 0,  # 链级 PDB 路线已废弃，整蛋白结构见 07
        "num_chains_distribution": {
            "min": int(final_df["num_protein_chains"].min()) if len(final_df) else 0,
            "max": int(final_df["num_protein_chains"].max()) if len(final_df) else 0,
            "mean": round(float(final_df["num_protein_chains"].mean()), 2) if len(final_df) else 0.0,
        },
        "note": "链级 PDB 路线已废弃（原 08_extract_protein_chains.py 归档）；"
                "chain_extraction_* 列为兼容保留，整蛋白结构由 07 落盘。",
    }

    with open(args.stats, "w", encoding="utf-8") as f:
        json.dump(stats, f, indent=2)

    print(f"[输出] 保留样本: {len(final_df)}")
    print(f"[输出] 剔除样本: {len(failure_df)}")
    print(f"[输出] 每条样本平均链数: {stats['num_chains_distribution']['mean']}")
    print(f"[输出] 最终数据集 CSV: {args.output}")
    print(f"[输出] 统计 JSON: {args.stats}")


if __name__ == "__main__":
    main()
