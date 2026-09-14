#!/usr/bin/env python3
"""
为 optimal_ph_prediction_test.csv 计算 OOD 标记
- seq_Redundancy_90 ~ seq_Redundancy_30（MMseqs2 easy-search）
- TM-score_0.9 ~ TM-score_0.3（Foldseek easy-search）
- OOD_Orphan（MMseqs2 easy-search，无显著 hit）
- OOD_IDR（metapredict v3 批量预测，IDR_ratio > 0.1，仅 Default 样本；
  本任务历史口径为 0.1，与其他任务的 0.3 不同，口径对照见
  .skills/ood-annotation-toolkit/SKILL.md 的 idr_ood.py 口径表；
  该逻辑已并入本脚本，原独立重算脚本 compute_ood_idr.py 已归档）

Default=False（极端长度）行的上述全部 OOD 列在保存前直接置 NaN
（全项目统一惯例：NaN 表示"未计算"，与 False 区分；仅 OOD_ExtremeShort /
OOD_ExtremeLong 原因列保留取值）。finalize_ood_columns.py --mode nan 仅作
兜底校验，不再承担 NaN 化职责。
"""
import os
import re
import json
import subprocess
import tempfile
from pathlib import Path
from collections import defaultdict

import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SPLITS_DIR = ROOT / "splits"
PDBS_DIR = ROOT / "pdbs"
OUTPUT_DIR = ROOT / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

TRAIN_PATH = SPLITS_DIR / "optimal_ph_prediction_train.csv"
VAL_PATH = SPLITS_DIR / "optimal_ph_prediction_val.csv"
TEST_PATH = SPLITS_DIR / "optimal_ph_prediction_test.csv"

SEQ_REDUNDANCY_THRESHOLDS = [90, 80, 70, 60, 50, 40, 30]
TM_SCORE_THRESHOLDS = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3]


def clean_seq(seq):
    return re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", seq.upper())


def write_fasta(df, fasta_path):
    with open(fasta_path, "w") as f:
        for _, row in df.iterrows():
            seq = clean_seq(row["aa_seq"])
            f.write(f">{row['unique_id']}\n{seq}\n")


def run_mmseqs_search(query_fasta, target_fasta, output_tsv, tmp_dir, sensitivity=7):
    """运行 MMseqs2 easy-search，返回结果 tsv"""
    tmp_dir = Path(tmp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "mmseqs", "easy-search",
        str(query_fasta), str(target_fasta), str(output_tsv), str(tmp_dir),
        "-s", str(sensitivity),
        "--format-output", "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits",
    ]
    print(f"[INFO] Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def run_foldseek_search(query_dir, target_dir, output_tsv, tmp_dir):
    """运行 Foldseek easy-search，返回结果 tsv"""
    tmp_dir = Path(tmp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "foldseek", "easy-search",
        str(query_dir), str(target_dir), str(output_tsv), str(tmp_dir),
        "--format-output", "query,target,alntmscore",
    ]
    print(f"[INFO] Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def parse_mmseqs_max_identity(tsv_path, test_ids):
    """解析 MMseqs2 结果，返回每个 query 的最大 pident"""
    max_identity = {uid: 0.0 for uid in test_ids}
    if not Path(tsv_path).exists():
        return max_identity
    with open(tsv_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            query, pident = parts[0], float(parts[2])
            if query in max_identity:
                max_identity[query] = max(max_identity[query], pident)
    return max_identity


def parse_foldseek_max_tmscore(tsv_path, test_ids):
    """解析 Foldseek 结果，返回每个 query 的最大 alntmscore"""
    max_tmscore = {uid: 0.0 for uid in test_ids}
    if not Path(tsv_path).exists():
        return max_tmscore
    with open(tsv_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            query, tmscore = parts[0], float(parts[2])
            if query in max_tmscore:
                max_tmscore[query] = max(max_tmscore[query], tmscore)
    return max_tmscore


def compute_seq_redundancy(test_df, train_val_fasta, tmp_dir):
    """计算 seq_Redundancy 标记"""
    print("[INFO] Computing seq_Redundancy...")
    test_fasta = tmp_dir / "test.fa"
    result_tsv = tmp_dir / "mmseqs_search.tsv"
    write_fasta(test_df, test_fasta)
    run_mmseqs_search(test_fasta, train_val_fasta, result_tsv, tmp_dir / "mmseqs_tmp")

    max_identity = parse_mmseqs_max_identity(result_tsv, test_df["unique_id"].tolist())

    for thr in SEQ_REDUNDANCY_THRESHOLDS:
        col = f"seq_Redundancy_{thr}"
        # 2026-08-29 起统一为严格小于（与全项目 < 口径一致；此前 <= 曾使
        # pident 恰等于阈值的 35 行多标 True，已重算修正）
        test_df[col] = test_df["unique_id"].map(lambda uid: max_identity.get(uid, 0.0) < thr)

    return test_df


def compute_tm_score(test_df, train_val_dir, tmp_dir):
    """计算 TM-score 标记"""
    print("[INFO] Computing TM-score...")
    # 创建临时 query 目录，包含所有 test 结构文件
    query_dir = tmp_dir / "test_structures"
    query_dir.mkdir(parents=True, exist_ok=True)
    for uid in test_df["unique_id"]:
        src = PDBS_DIR / f"{uid}.cif"
        if not src.exists():
            src = PDBS_DIR / f"{uid}.pdb"
        if src.exists():
            dst = query_dir / f"{uid}{src.suffix}"
            if not dst.exists():
                os.link(src, dst)

    result_tsv = tmp_dir / "foldseek_search.tsv"
    run_foldseek_search(query_dir, train_val_dir, result_tsv, tmp_dir / "foldseek_tmp")

    max_tmscore = parse_foldseek_max_tmscore(result_tsv, test_df["unique_id"].tolist())

    for thr in TM_SCORE_THRESHOLDS:
        col = f"TM-score_{thr:.1f}".replace(".0", "")
        # 注意列名格式统一为 TM-score_0.9
        col = f"TM-score_{thr}"
        # 2026-08-29 起统一为严格小于（实测无恰好等于阈值的样本，数值不变）
        test_df[col] = test_df["unique_id"].map(lambda uid: max_tmscore.get(uid, 0.0) < thr)

    return test_df


def compute_orphan(test_df, train_val_fasta, tmp_dir):
    """计算 OOD_Orphan：test vs train+val 无显著 hit"""
    print("[INFO] Computing OOD_Orphan...")
    test_fasta = tmp_dir / "test_orphan.fa"
    result_tsv = tmp_dir / "mmseqs_orphan.tsv"
    write_fasta(test_df, test_fasta)
    run_mmseqs_search(test_fasta, train_val_fasta, result_tsv, tmp_dir / "mmseqs_orphan_tmp", sensitivity=7)

    has_hit = set()
    if Path(result_tsv).exists():
        with open(result_tsv) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 1:
                    has_hit.add(parts[0])

    test_df["OOD_Orphan"] = ~test_df["unique_id"].isin(has_hit)
    return test_df


def compute_idr(test_df, idr_threshold=0.1):
    """计算 OOD_IDR：metapredict v3 批量预测，仅 Default 样本（非 Default → False）"""
    print("[INFO] Computing OOD_IDR...")
    try:
        from metapredict import predict_disorder_fasta
    except ImportError as e:
        print(f"[ERROR] metapredict not available: {e}, skipping OOD_IDR")
        test_df["OOD_IDR"] = False
        return test_df

    # 写入 FASTA（保留为 output/ 下的中间产物）
    fasta_path = OUTPUT_DIR / "test_sequences_for_idr.fasta"
    print(f"[INFO] Writing FASTA: {fasta_path}")
    with open(fasta_path, "w") as f:
        for _, row in test_df.iterrows():
            seq = clean_seq(row["aa_seq"])
            f.write(f">{row['unique_id']}\n{seq}\n")

    # 批量预测
    print("[INFO] Computing IDR ratios with metapredict v3 (batch)...")
    disorder_dict = predict_disorder_fasta(str(fasta_path), version="V3")

    idr_ratios = {}
    for uid, item in disorder_dict.items():
        # predict_disorder_fasta 返回 [sequence, scores]
        if isinstance(item, (list, tuple)) and len(item) == 2:
            scores = np.array(item[1])
        else:
            scores = np.array(item)
        idr_len = int((scores > 0.5).sum())
        idr_ratio = idr_len / len(scores) if len(scores) > 0 else 0.0
        idr_ratios[uid] = idr_ratio

    test_df["OOD_IDR"] = test_df.apply(
        lambda row: idr_ratios.get(row["unique_id"], 0.0) > idr_threshold if row["Default"] else False,
        axis=1,
    )
    return test_df


def main():
    print("[INFO] Loading splits...")
    train_df = pd.read_csv(TRAIN_PATH)
    val_df = pd.read_csv(VAL_PATH)
    test_df = pd.read_csv(TEST_PATH)

    print(f"[INFO] Train: {len(train_df)}, Val: {len(val_df)}, Test: {len(test_df)}")

    train_val_df = pd.concat([train_df, val_df], ignore_index=True)
    print(f"[INFO] Train+Val: {len(train_val_df)}")

    with tempfile.TemporaryDirectory(prefix="ood_", dir=str(OUTPUT_DIR)) as tmp:
        tmp_dir = Path(tmp)

        # 1. 准备 train+val FASTA
        train_val_fasta = tmp_dir / "train_val.fa"
        write_fasta(train_val_df, train_val_fasta)

        # 2. 准备 train+val 结构目录（用于 Foldseek）
        train_val_dir = tmp_dir / "train_val_structures"
        train_val_dir.mkdir(parents=True, exist_ok=True)
        for uid in train_val_df["unique_id"]:
            src = PDBS_DIR / f"{uid}.cif"
            if not src.exists():
                src = PDBS_DIR / f"{uid}.pdb"
            if src.exists():
                dst = train_val_dir / f"{uid}{src.suffix}"
                if not dst.exists():
                    os.link(src, dst)

        # 3. 计算 OOD 标记
        test_df = compute_seq_redundancy(test_df, train_val_fasta, tmp_dir)
        test_df = compute_tm_score(test_df, train_val_dir, tmp_dir)
        test_df = compute_orphan(test_df, train_val_fasta, tmp_dir)
        test_df = compute_idr(test_df)

    # 4. 确保 OOD 列为布尔值
    ood_cols = [c for c in test_df.columns if c.startswith(("seq_Redundancy_", "TM-score_", "OOD_"))]
    for col in ood_cols:
        test_df[col] = test_df[col].astype(bool)

    # 4b. Default=False（极端长度）行：除 OOD_ExtremeShort/OOD_ExtremeLong 外的
    #     OOD 列直接置 NaN（"未计算"语义，全项目统一惯例，流水线位置写入）
    nd_mask = ~test_df["Default"]
    keep_cols = {"OOD_ExtremeShort", "OOD_ExtremeLong"}
    nan_cols = [c for c in ood_cols if c not in keep_cols]
    for col in nan_cols:
        test_df[col] = test_df[col].astype(object)
        test_df.loc[nd_mask, col] = pd.NA
    print(f"[INFO] Default=False rows -> NaN: {int(nd_mask.sum())} rows x {len(nan_cols)} cols")

    # 5. 保存
    test_df.to_csv(TEST_PATH, index=False)
    print(f"[INFO] Updated {TEST_PATH}")

    # 6. 统计摘要
    summary = {
        "n_test": int(len(test_df)),
        "n_default": int(test_df["Default"].sum()),
        "n_extreme_long": int(test_df["OOD_ExtremeLong"].sum()),
        "n_extreme_short": int(test_df["OOD_ExtremeShort"].sum()),
        "OOD_IDR_threshold": 0.1,
        "OOD_IDR_default_only": True,
        "nondefault_nan_cols": len(nan_cols),
    }

    for col in ood_cols:
        summary[col] = int((test_df[col] == True).sum())  # noqa: E712

    summary_path = OUTPUT_DIR / "ood_markers_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[INFO] Summary saved to {summary_path}")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
