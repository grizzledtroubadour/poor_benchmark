#!/usr/bin/env python3
"""
使用 MMseqs2 对 EnzyBase12k 序列进行聚类去冗余
- 输入：metadata.csv 中的 uniprot_seq_cut
- 输出：聚类结果与统计摘要
"""
import os
import re
import json
import subprocess
import argparse
from pathlib import Path
from collections import Counter

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
OUTPUT_DIR = ROOT / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

META_PATH = DATA_DIR / "metadata.csv"


def clean_seq(seq):
    return re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", seq.upper())


def write_fasta(df, fasta_path):
    print(f"[INFO] Writing FASTA to {fasta_path}")
    with open(fasta_path, "w") as f:
        for _, row in df.iterrows():
            seq = clean_seq(row["uniprot_seq_cut"])
            if len(seq) >= 20:  # MMseqs2 要求最小长度
                f.write(f">{row['uniprot_id']}\n{seq}\n")
    print(f"[INFO] FASTA written")


def run_mmseqs(fasta_path, cluster_out_dir, min_seq_id=0.95, coverage=0.8, cov_mode=0):
    cluster_out_dir = Path(cluster_out_dir)
    cluster_out_dir.mkdir(parents=True, exist_ok=True)

    db_path = cluster_out_dir / "seqDB"
    cluster_db_path = cluster_out_dir / "cluDB"
    cluster_tsv = cluster_out_dir / "mmseqs_clust95.tsv"

    # 创建 MMseqs2 数据库
    print("[INFO] Creating MMseqs2 DB...")
    subprocess.run(
        ["mmseqs", "createdb", str(fasta_path), str(db_path)],
        check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )

    # 聚类
    print(f"[INFO] Clustering with min-seq-id={min_seq_id}, coverage={coverage}, cov-mode={cov_mode}...")
    subprocess.run(
        [
            "mmseqs", "cluster",
            str(db_path), str(cluster_db_path), str(cluster_out_dir / "tmp"),
            "--min-seq-id", str(min_seq_id),
            "-c", str(coverage),
            "--cov-mode", str(cov_mode),
        ],
        check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )

    # 创建 tsv 输出
    print("[INFO] Creating TSV output...")
    subprocess.run(
        ["mmseqs", "createtsv", str(db_path), str(db_path), str(cluster_db_path), str(cluster_tsv)],
        check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )

    print(f"[INFO] Cluster TSV saved to {cluster_tsv}")
    return cluster_tsv


def parse_clusters(cluster_tsv):
    """解析 MMseqs2 tsv 聚类结果：rep -> members"""
    clusters = {}
    with open(cluster_tsv) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            rep, member = parts[0], parts[1]
            clusters.setdefault(rep, []).append(member)
    return clusters


def main(min_seq_id=0.95, coverage=0.8, cov_mode=0):
    print(f"[INFO] Loading metadata: {META_PATH}")
    df = pd.read_csv(META_PATH)
    n_total = len(df)
    print(f"[INFO] Total sequences: {n_total}")

    fasta_path = OUTPUT_DIR / "sequences.fa"
    write_fasta(df, fasta_path)

    cluster_out_dir = OUTPUT_DIR / f"mmseqs_cluster_{int(min_seq_id*100)}"
    cluster_tsv = run_mmseqs(fasta_path, cluster_out_dir, min_seq_id, coverage, cov_mode)

    clusters = parse_clusters(cluster_tsv)
    n_clusters = len(clusters)
    cluster_sizes = [len(members) for members in clusters.values()]

    singletons = sum(1 for s in cluster_sizes if s == 1)
    multi = n_clusters - singletons
    max_cluster_size = max(cluster_sizes)
    removed = n_total - n_clusters

    summary = {
        "parameters": {
            "min_seq_id": min_seq_id,
            "coverage": coverage,
            "cov_mode": cov_mode,
        },
        "input_sequences": n_total,
        "n_clusters": n_clusters,
        "singleton_clusters": singletons,
        "multi_member_clusters": multi,
        "max_cluster_size": max_cluster_size,
        "removed_by_clustering": removed,
        "retention_ratio": round(n_clusters / n_total, 4),
        "reduction_ratio": round(removed / n_total, 4),
        "cluster_size_distribution": dict(Counter(cluster_sizes)),
    }

    summary_path = OUTPUT_DIR / f"mmseqs_cluster_{int(min_seq_id*100)}_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[INFO] Summary saved to {summary_path}")

    print("\n===== Clustering Summary =====")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--min_seq_id", type=float, default=0.95)
    parser.add_argument("--coverage", type=float, default=0.8)
    parser.add_argument("--cov_mode", type=int, default=0)
    args = parser.parse_args()
    main(args.min_seq_id, args.coverage, args.cov_mode)
