#!/usr/bin/env python3
"""
脚本：04_parse_ligands.py
用途：解析 PDBbind v2020.R1 清洗后数据中的配体分子，生成 canonical SMILES 与 ECFP4 指纹。

使用：
    cd binding_affinity/
    python scripts/04_parse_ligands.py

输入：
    output/pdbbind_v2020_metadata_cleaned.csv

输出：
    output/pdbbind_v2020_metadata_ligand.csv
    output/pdbbind_v2020_metadata_ligand_stats.json
    output/ligand_parse_failures.csv
"""

import argparse
import csv
import json
from pathlib import Path
from collections import Counter

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

# 抑制 RDKit 大量 warning
RDLogger.DisableLog("rdApp.*")


def parse_ligand_from_sdf(path: Path):
    """从 SDF 文件解析配体，返回 (mol, error_msg)。"""
    if not path.exists():
        return None, "SDF file not found"
    try:
        suppl = Chem.SDMolSupplier(str(path), sanitize=True, removeHs=False)
        mols = [m for m in suppl if m is not None]
        if not mols:
            return None, "No valid molecule in SDF"
        mol = mols[0]
        if mol.GetNumHeavyAtoms() == 0:
            return None, "Molecule has no heavy atoms"
        return mol, None
    except Exception as e:
        return None, f"SDF parse error: {e}"


def parse_ligand_from_mol2(path: Path):
    """从 MOL2 文件解析配体，返回 (mol, error_msg)。"""
    if not path.exists():
        return None, "MOL2 file not found"
    try:
        mol = Chem.MolFromMol2File(str(path), removeHs=False)
        if mol is None:
            return None, "MOL2 parse failed"
        if mol.GetNumHeavyAtoms() == 0:
            return None, "Molecule has no heavy atoms"
        return mol, None
    except Exception as e:
        return None, f"MOL2 parse error: {e}"


def get_canonical_smiles(mol):
    """获取 canonical SMILES。"""
    try:
        return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    except Exception:
        return None


def get_ecfp4(mol, radius=2, n_bits=2048):
    """计算 ECFP4（Morgan）指纹，返回为 1 位索引列表。"""
    try:
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=n_bits)
        on_bits = sorted([int(i) for i in fp.GetOnBits()])
        return on_bits
    except Exception:
        return None


def parse_ligand(row, data_dir: Path):
    """尝试解析单条记录的配体。"""
    pdb_id = row["pdb_id"]
    sdf_rel = row.get("ligand_sdf_file", "")
    mol2_rel = row.get("ligand_mol2_file", "")

    sdf_path = data_dir / sdf_rel if pd.notna(sdf_rel) and sdf_rel else None
    mol2_path = data_dir / mol2_rel if pd.notna(mol2_rel) and mol2_rel else None

    # 优先尝试 SDF
    mol, error = None, None
    if sdf_path and sdf_path.exists():
        mol, error = parse_ligand_from_sdf(sdf_path)

    # SDF 失败则尝试 MOL2
    if mol is None and mol2_path and mol2_path.exists():
        mol, error = parse_ligand_from_mol2(mol2_path)
        if mol is not None:
            error = None

    # 两者都失败
    if mol is None:
        return {
            "pdb_id": pdb_id,
            "ligand_parse_success": False,
            "ligand_parse_error": error or "No ligand file available",
            "ligand_smiles": None,
            "ligand_ecfp4": None,
            "ligand_num_heavy_atoms": None,
        }

    smiles = get_canonical_smiles(mol)
    ecfp4 = get_ecfp4(mol)
    num_heavy = mol.GetNumHeavyAtoms()

    return {
        "pdb_id": pdb_id,
        "ligand_parse_success": True,
        "ligand_parse_error": None,
        "ligand_smiles": smiles,
        "ligand_ecfp4": json.dumps(ecfp4) if ecfp4 is not None else None,
        "ligand_num_heavy_atoms": num_heavy,
    }


def main():
    parser = argparse.ArgumentParser(description="Parse ligands for PDBbind v2020.R1")
    parser.add_argument(
        "--input",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "output" / "pdbbind_v2020_metadata_cleaned.csv",
        help="Input cleaned metadata CSV",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "output" / "pdbbind_v2020_metadata_ligand.csv",
        help="Output CSV with ligand features",
    )
    parser.add_argument(
        "--stats",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "output" / "pdbbind_v2020_metadata_ligand_stats.json",
        help="Output stats JSON",
    )
    parser.add_argument(
        "--failures",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "output" / "ligand_parse_failures.csv",
        help="Output failures CSV",
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
    print(f"[输入] 清洗后样本数: {total}")

    results = []
    for _, row in df.iterrows():
        res = parse_ligand(row, args.data_dir)
        results.append(res)

    ligand_df = pd.DataFrame(results)

    # 合并到原表
    merged = df.merge(
        ligand_df,
        on="pdb_id",
        how="left",
    )

    # 过滤解析失败的样本
    success_df = merged[merged["ligand_parse_success"] == True].copy()
    failure_df = merged[merged["ligand_parse_success"] != True].copy()

    # 输出
    args.output.parent.mkdir(parents=True, exist_ok=True)
    success_df.to_csv(args.output, index=False)
    failure_df.to_csv(args.failures, index=False)

    stats = {
        "total_input": total,
        "ligand_parse_success": int(len(success_df)),
        "ligand_parse_failure": int(len(failure_df)),
        "retention_rate": round(len(success_df) / total, 4) if total > 0 else 0.0,
        "failure_reasons": dict(Counter(failure_df["ligand_parse_error"].astype(str))),
        "heavy_atoms_stats": success_df["ligand_num_heavy_atoms"].describe().to_dict() if len(success_df) > 0 else {},
    }

    with open(args.stats, "w", encoding="utf-8") as f:
        json.dump(stats, f, indent=2)

    print(f"[输出] 配体解析成功: {len(success_df)}")
    print(f"[输出] 配体解析失败: {len(failure_df)}")
    print(f"[输出] 保留比例: {stats['retention_rate']*100:.2f}%")
    print(f"[输出] 成功样本 CSV: {args.output}")
    print(f"[输出] 失败样本 CSV: {args.failures}")
    print(f"[输出] 统计 JSON: {args.stats}")
    print()
    print("=== 失败原因分布 ===")
    for reason, count in stats["failure_reasons"].items():
        print(f"  {reason}: {count}")


if __name__ == "__main__":
    main()
