#!/usr/bin/env python3
"""
对 chain-centric 数据集进行 95% identity 序列去冗余。

步骤：
1. 从 chain_centric_metadata.csv 读取 22,223 条唯一链序列。
2. 用 mmseqs2 cluster 以 95% identity / 80% coverage 聚类。
3. 每个聚类保留一条代表序列（mmseqs2 自动选择）。

历史版本（dedup_chain_centric_and_link_structures.py）曾从 Pair-Centric 旧线
产物（output/structures/，已归档）为代表链创建结构符号链接到
output/chain_representative_structures/；该环节已被
08_extract_representative_structures.py（直接从 pair 文件提取坐标写入 pdbs/）
完全取代，本脚本已裁掉符号链接死代码，结构来源字段
（has_structure / representative_structure / structure_source_*）
初始化为空，由 08 回填。

输出：
    output/chain_centric_dedup_metadata.csv
    output/chain_centric_dedup_summary.json
"""

import argparse
import csv
import hashlib
import json
import logging
import os
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

Path("logs").mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("logs/dedup_chain_centric.log"),
        logging.StreamHandler(sys.stdout),
    ],
)
logger = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = ROOT / "output"


def seq_hash(seq):
    return hashlib.md5(seq.encode("utf-8")).hexdigest()


def load_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def write_fasta(records, out_path):
    """records: list of dict with 'seq_hash' and 'sequence'."""
    with open(out_path, "w") as f:
        for r in records:
            f.write(f">{r['seq_hash']}\n{r['sequence']}\n")


def run_mmseqs_cluster(fasta_path, tmp_dir, num_workers=32):
    """运行 mmseqs2 cluster 并返回 tsv 聚类结果路径。"""
    os.makedirs(tmp_dir, exist_ok=True)
    db = os.path.join(tmp_dir, "seq_db")
    clu = os.path.join(tmp_dir, "clu")
    tsv = os.path.join(tmp_dir, "cluster.tsv")

    subprocess.run(["mmseqs", "createdb", fasta_path, db], check=True, capture_output=True, text=True)
    subprocess.run(
        ["mmseqs", "cluster", db, clu, tmp_dir,
         "--min-seq-id", "0.95", "-c", "0.8", "--cov-mode", "0",
         "--threads", str(num_workers)],
        check=True, capture_output=True, text=True,
    )
    subprocess.run(["mmseqs", "createtsv", db, db, clu, tsv], check=True, capture_output=True, text=True)
    return tsv


def parse_cluster_tsv(tsv_path):
    """解析 mmseqs2 cluster tsv，返回 representative -> [members] 的 dict。"""
    clusters = defaultdict(list)
    with open(tsv_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            rep, member = parts[0], parts[1]
            clusters[rep].append(member)
    return clusters


def main():
    parser = argparse.ArgumentParser(description="Dedup chain-centric dataset (95% identity)")
    parser.add_argument("--chain_metadata", default="output/chain_centric_metadata.csv")
    parser.add_argument("--num_workers", type=int, default=32)
    args = parser.parse_args()

    chain_records = load_csv(ROOT / args.chain_metadata)
    logger.info("Loaded %d unique chain records", len(chain_records))

    # 1. 写 FASTA
    fasta_path = OUTPUT_DIR / "chain_centric_sequences.fasta"
    write_fasta(chain_records, fasta_path)
    logger.info("Wrote FASTA with %d sequences to %s", len(chain_records), fasta_path)

    # 2. mmseqs2 cluster
    with tempfile.TemporaryDirectory(prefix="chain_cluster_", dir=str(OUTPUT_DIR)) as tmp_dir:
        tsv_path = run_mmseqs_cluster(str(fasta_path), tmp_dir, args.num_workers)
        clusters = parse_cluster_tsv(tsv_path)
    logger.info("Clustering complete: %d clusters", len(clusters))

    # 3. 构建序列哈希 -> 完整记录的映射
    hash_to_record = {r["seq_hash"]: r for r in chain_records}

    # 4. 生成去冗余后的元数据（结构来源字段留空，由 08 回填）
    rep_rows = []

    for rep_hash, members in clusters.items():
        rep_record = hash_to_record.get(rep_hash)
        if rep_record is None:
            logger.warning("Representative hash not found in metadata: %s", rep_hash)
            continue

        # 该代表序列的所有成员记录
        member_records = [hash_to_record[m] for m in members if m in hash_to_record]
        total_occurrences = sum(int(r["num_occurrences"]) for r in member_records)
        total_interface = sum(int(r["num_interface_residues"]) for r in member_records)

        rep_rows.append({
            "seq_hash": rep_hash,
            "sequence": rep_record["sequence"],
            "length": rep_record["length"],
            "cluster_size": len(members),
            "num_occurrences": total_occurrences,
            "num_interface_residues": int(rep_record["num_interface_residues"]),
            "interface_ratio": float(rep_record["interface_ratio"]),
            "as_receptor": int(rep_record["as_receptor"]),
            "as_ligand": int(rep_record["as_ligand"]),
            "member_seq_hashes": ";".join(members[:20]) + (";..." if len(members) > 20 else ""),
            "has_structure": "False",
            "representative_structure": "",
            "structure_source_unique_id": "",
            "structure_source_role": "",
        })

    # 按聚类大小降序排列
    rep_rows.sort(key=lambda r: (-r["cluster_size"], -r["num_occurrences"]))

    out_path = OUTPUT_DIR / "chain_centric_dedup_metadata.csv"
    fieldnames = [
        "seq_hash", "sequence", "length", "cluster_size", "num_occurrences",
        "num_interface_residues", "interface_ratio",
        "as_receptor", "as_ligand",
        "member_seq_hashes", "has_structure", "representative_structure",
        "structure_source_unique_id", "structure_source_role",
    ]
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rep_rows)
    logger.info("Saved dedup metadata: %d representative chains to %s", len(rep_rows), out_path)

    # 5. 统计摘要（结构关联统计由 08 提取代表结构后回填，本阶段不含）
    summary = {
        "original_unique_chains": len(chain_records),
        "representative_chains_after_dedup": len(rep_rows),
        "reduction_ratio": 1.0 - len(rep_rows) / len(chain_records),
        "clusters": len(clusters),
        "mean_cluster_size": sum(r["cluster_size"] for r in rep_rows) / len(rep_rows) if rep_rows else 0,
        "max_cluster_size": max(r["cluster_size"] for r in rep_rows) if rep_rows else 0,
        "mean_interface_ratio": sum(r["interface_ratio"] for r in rep_rows) / len(rep_rows) if rep_rows else 0,
    }
    summary_path = OUTPUT_DIR / "chain_centric_dedup_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    logger.info("Summary: %s", json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
