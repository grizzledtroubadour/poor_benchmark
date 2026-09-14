#!/usr/bin/env python3
"""
脚本：10_compute_structure_idr_ood.py
用途：为测试集计算结构相似性 OOD（Foldseek-Multimer 复合物级对称 TM-score）
      和内在无序区域 OOD（基于 metapredict，阈值 0.1）。

输出列（在 09 产物基础上追加/重排）：
  - TM-score_0.9, ..., TM-score_0.3
  - idr_ratio
  - OOD_IDR
  （Default / OOD_ExtremeShort / OOD_ExtremeLong / seq_Redundancy_* / OOD_Orphan 由 09 生成）

结构 OOD 口径（2026-08-29 起）：foldseek `easy-multimersearch` 复合物比对
（Foldseek-Multimer, Kim et al., Nat Methods 2025），取每对复合物的
complexqtmscore/complexttmscore 对称均值的最大值；复合物级比对最贴合
"复合物-配体亲和力"任务语义。旧口径（单链 `search -a` + qtmscore/ttmscore
对称均值）已于 2026-08-29 替换（备份 output/backup_before_multimer_tm_unify/，
对比见 analysis/output/foldseek_easysearch_audit/lba/multimer/）。
IDR 以 .skills/ood-annotation-toolkit 的 idr_ood.py 为准（本任务阈值 0.1）。

使用：
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
    """运行 Foldseek-Multimer easy-multimersearch（复合物级比对），返回结果 tsv。"""
    out_dir.mkdir(parents=True, exist_ok=True)
    result_tsv = out_dir / "foldseek_multimer_test_vs_trainval.tsv"
    tmp_dir = out_dir / "tmp_multimer"
    q_dir = out_dir / "query_structures"
    t_dir = out_dir / "target_structures"
    # multimersearch 需要结构目录，用符号链接避免拷贝
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
    """foldseek createdb 会把多链 PDB 拆成 '<file>_<chain>' 条目，单链文件保留
    裸文件名；据此把 m8 中的 query/target 名归一回复合物文件 stem。"""
    if chain_name in known_stems:
        return chain_name
    return chain_name.rsplit("_", 1)[0]


def compute_max_tmscore(result_tsv: Path, test_stems: set, ref_stems: set):
    """聚合每个 test 复合物的对称复合物 TM-score（complexqtm/complexttm 平均）最大值。"""
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
        help="Directory containing whole-protein PDB files (struct_file 所在目录)",
    )
    args = parser.parse_args()

    train_df = pd.read_csv(args.splits_dir / f"{TASK_NAME}_train.csv")
    val_df = pd.read_csv(args.splits_dir / f"{TASK_NAME}_val.csv")
    test_df = pd.read_csv(args.splits_dir / f"{TASK_NAME}_test.csv")

    print(f"[输入] train: {len(train_df)}, val: {len(val_df)}, test: {len(test_df)}")

    # 为 Foldseek 创建文件列表（struct_file 为文件名，拼上 pdbs 目录）
    fs_dir = args.output_dir / "foldseek_ood"
    fs_dir.mkdir(parents=True, exist_ok=True)
    trainval_list = fs_dir / "trainval_files.tsv"
    test_list = fs_dir / "test_files.tsv"

    trainval_files = [args.pdbs_dir / f for f in pd.concat([train_df, val_df])["struct_file"]]
    test_files = [args.pdbs_dir / f for f in test_df["struct_file"]]
    write_file_list(trainval_files, trainval_list)
    write_file_list(test_files, test_list)

    # 运行 Foldseek-Multimer
    print("[Foldseek-Multimer] 计算 test 与 train+val 的复合物结构相似性...")
    result_tsv = run_foldseek(trainval_list, test_list, fs_dir)
    test_stems = {Path(f).stem for f in test_df["struct_file"]}
    ref_stems = {Path(f).stem for f in pd.concat([train_df, val_df])["struct_file"]}
    max_tm = compute_max_tmscore(result_tsv, test_stems, ref_stems)

    # 添加 TM-score OOD 列
    for thr in TM_THRESHOLDS:
        col = f"TM-score_{thr:.1f}"
        test_df[col] = test_df["unique_id"].apply(
            lambda uid: max_tm.get(uid, 0.0) < thr
        )

    # 计算 IDR（batch）
    print("[metapredict] 计算 test 序列的内在无序区域比例...")
    # metapredict 只接受标准氨基酸序列，去除链分隔符 '|'
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
            print(f"  已完成 {min(i + batch_size, len(seqs))}/{len(seqs)}")

    test_df["idr_ratio"] = idr_ratios
    test_df["OOD_IDR"] = test_df["idr_ratio"] > IDR_THRESHOLD

    # Default：常规测试集样本，全部设为 True
    test_df["Default"] = True

    # 极端长度（基于氨基酸长度，不含链分隔符 '|'）
    test_lengths = test_df["aa_seq"].str.replace("|", "").str.len().values
    test_df["OOD_ExtremeShort"] = test_lengths < 60
    test_df["OOD_ExtremeLong"] = test_lengths > 1000

    # 序列冗余 OOD（复用 09_format_splits_to_skill.py 结果，若不存在则计算）
    if "OOD_Orphan" not in test_df.columns:
        # 若尚未计算，此处留空；正常情况下 split 文件已包含
        raise RuntimeError("请先运行 09_format_splits_to_skill.py 生成序列同源性 OOD 列")

    # 调整列顺序（已知列按统一 schema 排列，后续步骤追加的其他列原样保留，
    # 以便在当前完整 splits 上重跑时不丢失后续 OOD 列）
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

    # 保存
    test_out_path = args.splits_dir / f"{TASK_NAME}_test.csv"
    test_df.to_csv(test_out_path, index=False)

    # 更新统计 JSON
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

    # 清理 Foldseek 临时文件
    if fs_dir.exists():
        shutil.rmtree(fs_dir)

    print(f"\n[输出] 已更新 {test_out_path}")
    print(f"[OOD] OOD_IDR=True (idr>{IDR_THRESHOLD}): {test_df['OOD_IDR'].sum()}")
    for thr in TM_THRESHOLDS:
        col = f"TM-score_{thr:.1f}"
        print(f"[OOD] {col}=True: {test_df[col].sum()}")
    print(f"[OOD] Default=True: {test_df['Default'].sum()}")


if __name__ == "__main__":
    main()
