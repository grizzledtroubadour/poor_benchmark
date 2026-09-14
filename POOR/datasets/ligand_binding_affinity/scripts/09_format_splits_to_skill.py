#!/usr/bin/env python3
"""
脚本：09_format_splits_to_skill.py
用途：将数据集按 POOD-Benchmark 统一 schema 输出为划分文件：
      - 文件命名：{task_name}_train.csv / {task_name}_val.csv / {task_name}_test.csv
      - 核心列：unique_id, aa_seq, struct_file, ligand_smiles, ligand_ecfp4, label
      - test 追加 OOD 列：Default, OOD_ExtremeShort, OOD_ExtremeLong,
        seq_Redundancy_90..30, OOD_Orphan
      - label（pK）强制保存为字符串，统一 6 位小数

注意（superseded）：本脚本内的 mmseqs2 序列同源性段为历史首遍实现，
统一口径与重放工具以 .skills/homology-ood-annotation 的
scripts/seq_homology_ood.py 为准（支持 easy-search 与 --m8-cache 重放）。

使用：
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
    """写出 FASTA；mmseqs2 只接受标准氨基酸，因此去除链分隔符 '|'。"""
    with open(path, "w", encoding="utf-8") as f:
        for _, row in df.iterrows():
            uid = row[id_col]
            seq = row[seq_col]
            if pd.isna(seq):
                seq = ""
            f.write(f">{uid}\n{seq.replace('|', '')}\n")


def parse_max_identity(out_tsv):
    """解析 convertalis 结果（query,target,pident,evalue），返回 {query: max_identity}。"""
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
    """运行 mmseqs2 search 并返回最佳 hit 的 identity 字典 {query: max_identity}。

    superseded：统一口径以 homology-ood-annotation skill 的 seq_homology_ood.py 为准。
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
    """从完整元数据构造统一 schema 基础 DataFrame。"""
    out = pd.DataFrame()
    out["unique_id"] = meta_df["pdb_id"].astype(str)
    out["aa_seq"] = meta_df["protein_sequence"]
    out["struct_file"] = meta_df["pdb_id"].astype(str) + ".pdb"
    out["ligand_smiles"] = meta_df["ligand_smiles"]
    out["ligand_ecfp4"] = meta_df["ligand_ecfp4"]
    # label 强制为字符串，统一保留 6 位小数
    out["label"] = meta_df["pK"].apply(lambda x: f"{x:.6f}" if pd.notna(x) else "")
    # 保留后续计算 OOD 所需的辅助列
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
        help="可选：已有的 convertalis 结果 tsv（query,target,pident,evalue），"
             "提供后跳过 mmseqs2 重算，直接解析缓存",
    )
    args = parser.parse_args()

    # 读取完整元数据
    meta_df = pd.read_csv(args.dataset)
    base_df = build_base_df(meta_df)

    # 读取当前划分文件中的 unique_id 列表，确定每个样本属于哪个集合
    def load_split_ids(name):
        path = args.splits_dir / f"{TASK_NAME}_{name}.csv"
        if not path.exists():
            # 兼容旧命名（首次迁移时使用）
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

    # 计算 test 的序列冗余 OOD（superseded：统一口径以 homology-ood-annotation skill 为准）
    print("[OOD] 使用 mmseqs2 计算 test 与 train+val 的序列同源性...")
    train_val_df = pd.concat([train_base, val_base], ignore_index=True)
    query_fasta = args.output_dir / "test_for_mmseqs.fasta"
    target_fasta = args.output_dir / "train_val_for_mmseqs.fasta"
    result_tsv = args.output_dir / "mmseqs_test_vs_trainval.tsv"
    mmseqs_tmp = args.output_dir / "mmseqs_tmp"

    write_fasta(train_val_df, target_fasta)
    write_fasta(test_base, query_fasta)
    if args.mmseqs_cache is not None:
        print(f"[OOD] 使用缓存的 mmseqs2 结果: {args.mmseqs_cache}")
        max_identity = parse_max_identity(args.mmseqs_cache)
    else:
        max_identity = run_mmseqs2(query_fasta, target_fasta, result_tsv, mmseqs_tmp)

    # OOD 列
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

    # Default：常规测试集样本，本任务不对 test 做提前去除，因此全部设为 True
    test_base["Default"] = True

    # 删除辅助列
    train_base = train_base.drop(columns=["_protein_length"])
    val_base = val_base.drop(columns=["_protein_length"])
    test_base = test_base.drop(columns=["_protein_length"])

    # 调整列顺序
    base_cols = ["unique_id", "aa_seq", "struct_file", "ligand_smiles", "ligand_ecfp4", "label"]
    ood_cols_ordered = ["Default", "OOD_ExtremeShort", "OOD_ExtremeLong"] + \
                       [f"seq_Redundancy_{thr}" for thr in SEQ_THRESHOLDS] + \
                       ["OOD_Orphan"]
    test_base = test_base[base_cols + ood_cols_ordered]

    # 保存
    args.splits_dir.mkdir(parents=True, exist_ok=True)
    train_out_path = args.splits_dir / f"{TASK_NAME}_train.csv"
    val_out_path = args.splits_dir / f"{TASK_NAME}_val.csv"
    test_out_path = args.splits_dir / f"{TASK_NAME}_test.csv"

    train_base.to_csv(train_out_path, index=False)
    val_base.to_csv(val_out_path, index=False)
    test_base.to_csv(test_out_path, index=False)

    # 清理旧命名文件（如果存在）
    for old in ["train.csv", "val.csv", "test.csv"]:
        old_path = args.splits_dir / old
        if old_path.exists():
            old_path.unlink()

    # 清理 mmseqs 临时文件
    if mmseqs_tmp.exists():
        shutil.rmtree(mmseqs_tmp)

    print("[输出] 已生成：")
    print(f"  {train_out_path}: {len(train_base)} 行")
    print(f"  {val_out_path}: {len(val_base)} 行")
    print(f"  {test_out_path}: {len(test_base)} 行")
    print(f"\n[Test OOD 统计]")
    print(f"  Default=True: {test_base['Default'].sum()}")
    print(f"  OOD_ExtremeShort=True: {test_base['OOD_ExtremeShort'].sum()}")
    print(f"  OOD_ExtremeLong=True: {test_base['OOD_ExtremeLong'].sum()}")
    print(f"  OOD_Orphan=True: {test_base['OOD_Orphan'].sum()}")
    for thr in SEQ_THRESHOLDS:
        col = f"seq_Redundancy_{thr}"
        print(f"  {col}=True: {test_base[col].sum()}")


if __name__ == "__main__":
    main()
