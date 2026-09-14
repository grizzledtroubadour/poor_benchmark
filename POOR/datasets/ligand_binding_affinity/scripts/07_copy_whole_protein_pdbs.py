#!/usr/bin/env python3
"""
脚本：07_copy_whole_protein_pdbs.py
用途：删除 pdbs/ 下的链级结构，改为存放整条蛋白结构文件。
      文件名仅使用 {pdb_id}.pdb，不含有配体名；文件中过滤掉非蛋白残基的 HETATM。

使用：
    cd ligand_binding_affinity/
    python scripts/07_copy_whole_protein_pdbs.py

输入：
    output/pdbbind_v2020_dataset.csv

输出：
    pdbs/{pdb_id}.pdb
    output/pdbbind_v2020_dataset.csv（更新蛋白结构文件路径列）
"""

import argparse
import json
from pathlib import Path
from collections import Counter

import pandas as pd

# 标准氨基酸与常见修饰氨基酸
AA_3TO1 = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
}

MODIFIED_AA = {
    "MSE": "M", "SEP": "S", "TPO": "T", "PTR": "Y", "KCX": "K",
    "CSO": "C", "CSD": "C", "CME": "C", "LYZ": "K", "MLY": "K",
    "M3L": "K", "HYP": "P", "DAL": "A", "DVA": "V", "DLE": "L",
    "DIL": "I", "DTY": "Y", "DTR": "W", "DSE": "S", "DTH": "T",
    "DASN": "N", "DGLN": "Q", "DARG": "R", "DLY": "K", "DPR": "P",
    "DPH": "F", "DHI": "H", "DME": "M", "DAS": "D", "DGL": "E",
    "DCY": "C", "DGP": "G",
}

PROTEIN_RESIDUES = set(list(AA_3TO1.keys()) + list(MODIFIED_AA.keys()))


def copy_protein_without_ligand(src_path: Path, dst_path: Path):
    """
    复制蛋白 PDB，排除非蛋白残基的 HETATM（如水、配体小分子）。
    遇到第一个 ENDMDL 后停止，跳过 NMR 后续模型。
    """
    with open(src_path, "r", encoding="utf-8") as fin, \
         open(dst_path, "w", encoding="utf-8") as fout:
        for line in fin:
            if line.startswith("ENDMDL"):
                fout.write("END\n")
                break

            if line.startswith("ATOM  "):
                fout.write(line)
                continue

            if line.startswith("HETATM"):
                if len(line) >= 20:
                    res_name = line[17:20].strip()
                    if res_name in PROTEIN_RESIDUES:
                        fout.write(line)
                continue

            # 其他行（HEADER, COMPND, REMARK, SEQRES, SSBOND, TER, END 等）原样保留
            fout.write(line)

            if line.startswith("END"):
                break


def main():
    parser = argparse.ArgumentParser(description="Copy whole protein PDBs without ligand")
    parser.add_argument(
        "--input",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "output" / "pdbbind_v2020_dataset.csv",
        help="Input dataset CSV",
    )
    parser.add_argument(
        "--pdbs-dir",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "pdbs",
        help="Output directory for whole-protein PDB files",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "data",
        help="Root data directory",
    )
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    total = len(df)
    print(f"[输入] 样本数: {total}")

    args.pdbs_dir.mkdir(parents=True, exist_ok=True)

    # 清空已有的 .pdb 文件
    removed = 0
    for old_pdb in args.pdbs_dir.glob("*.pdb"):
        old_pdb.unlink()
        removed += 1
    print(f"[清理] 删除旧 PDB 文件: {removed} 个")

    success = 0
    failed = []
    project_root = Path(__file__).resolve().parent.parent

    for _, row in df.iterrows():
        pdb_id = row["pdb_id"]
        protein_rel = row.get("protein_file", "")
        src_path = args.data_dir / protein_rel if pd.notna(protein_rel) and protein_rel else None

        if src_path is None or not src_path.exists():
            failed.append((pdb_id, "protein file missing"))
            continue

        dst_path = args.pdbs_dir / f"{pdb_id}.pdb"
        try:
            copy_protein_without_ligand(src_path, dst_path)
            success += 1
        except Exception as e:
            failed.append((pdb_id, str(e)))

    # 更新数据集：去掉链级文件路径列，添加整蛋白结构文件路径
    if "protein_chain_files" in df.columns:
        df = df.drop(columns=["protein_chain_files"])
    df["protein_pdb_file"] = df["pdb_id"].apply(lambda x: f"pdbs/{x}.pdb")

    args.input.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.input, index=False)

    stats = {
        "total_input": total,
        "protein_pdb_written": success,
        "failed": len(failed),
        "old_pdb_removed": removed,
        "failure_reasons": dict(Counter([e for _, e in failed])),
    }
    stats_path = args.input.parent / "whole_protein_pdb_stats.json"
    with open(stats_path, "w", encoding="utf-8") as f:
        json.dump(stats, f, indent=2)

    print(f"[输出] 成功写入整蛋白 PDB: {success}")
    print(f"[输出] 失败: {len(failed)}")
    print(f"[输出] 更新数据集 CSV: {args.input}")
    print(f"[输出] 统计 JSON: {stats_path}")


if __name__ == "__main__":
    main()
