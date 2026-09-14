#!/usr/bin/env python3
"""
划分 Chain-Centric 数据集：
1. 将链长不在 [60, 1000] 区间内的链作为独立 ExtremeLength 测试集。
2. 剩余链按 72:8:20 随机划分为 Train / Val / Test。

输出：
    output/chain_centric_splits/chain_centric_train.csv
    output/chain_centric_splits/chain_centric_val.csv
    output/chain_centric_splits/chain_centric_test.csv
    output/chain_centric_splits/chain_centric_test_extremelength.csv
    output/chain_centric_split_summary.json
"""

import argparse
import csv
import json
import logging
import os
import random
import sys
from pathlib import Path

Path("logs").mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("logs/split_chain_centric.log"),
        logging.StreamHandler(sys.stdout),
    ],
)
logger = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = ROOT / "output"
SPLIT_DIR = ROOT / "output" / "chain_centric_splits"

OUTPUT_FIELDS = [
    "seq_hash",
    "sequence",
    "length",
    "num_occurrences",
    "num_interface_residues",
    "interface_ratio",
    "as_receptor",
    "as_ligand",
    "representative_structure",
]


def load_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def save_split(name, rows):
    path = SPLIT_DIR / f"chain_centric_{name}.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=OUTPUT_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    logger.info("Saved %s with %d rows to %s", name, len(rows), path)


def main():
    parser = argparse.ArgumentParser(description="Split chain-centric dataset")
    parser.add_argument("--metadata", default="output/chain_centric_dedup_metadata.csv")
    parser.add_argument("--min_len", type=int, default=60)
    parser.add_argument("--max_len", type=int, default=1000)
    parser.add_argument("--ratios", nargs=3, type=float, default=[0.72, 0.08, 0.20])
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    if abs(sum(args.ratios) - 1.0) > 1e-6:
        raise ValueError("ratios 之和必须等于 1.0")

    os.makedirs("logs", exist_ok=True)

    metadata = load_csv(ROOT / args.metadata)
    logger.info("Loaded %d representative chain records", len(metadata))

    # 1. 分离 ExtremeLength
    normal_rows = []
    extreme_rows = []
    for row in metadata:
        length = int(row["length"])
        if length < args.min_len or length > args.max_len:
            extreme_rows.append(row)
        else:
            normal_rows.append(row)

    logger.info(
        "ExtremeLength (%d < len or len > %d): %d; Normal: %d",
        args.min_len, args.max_len, len(extreme_rows), len(normal_rows),
    )

    # 2. 对 normal 链随机 72:8:20 划分
    random.seed(args.seed)
    shuffled = normal_rows.copy()
    random.shuffle(shuffled)

    n_total = len(shuffled)
    n_train = int(n_total * args.ratios[0])
    n_val = int(n_total * args.ratios[1])
    train_rows = shuffled[:n_train]
    val_rows = shuffled[n_train : n_train + n_val]
    test_rows = shuffled[n_train + n_val :]

    # 3. 格式化输出
    def format_row(row):
        return {
            "seq_hash": row["seq_hash"],
            "sequence": row["sequence"],
            "length": row["length"],
            "num_occurrences": row["num_occurrences"],
            "num_interface_residues": row["num_interface_residues"],
            "interface_ratio": row["interface_ratio"],
            "as_receptor": row["as_receptor"],
            "as_ligand": row["as_ligand"],
            "representative_structure": row.get("representative_structure", ""),
        }

    save_split("train", [format_row(r) for r in train_rows])
    save_split("val", [format_row(r) for r in val_rows])
    save_split("test", [format_row(r) for r in test_rows])
    save_split("test_extremelength", [format_row(r) for r in extreme_rows])

    # 4. 统计摘要
    def stats(rows):
        lengths = [int(r["length"]) for r in rows]
        ratios = [float(r["interface_ratio"]) for r in rows]
        return {
            "count": len(rows),
            "mean_length": sum(lengths) / len(lengths) if lengths else 0,
            "min_length": min(lengths) if lengths else 0,
            "max_length": max(lengths) if lengths else 0,
            "mean_interface_ratio": sum(ratios) / len(ratios) if ratios else 0,
        }

    summary = {
        "total_representative_chains": len(metadata),
        "min_len_threshold": args.min_len,
        "max_len_threshold": args.max_len,
        "seed": args.seed,
        "ratios": args.ratios,
        "train": stats(train_rows),
        "val": stats(val_rows),
        "test": stats(test_rows),
        "test_extremelength": stats(extreme_rows),
    }

    summary_path = OUTPUT_DIR / "chain_centric_split_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    logger.info("Summary saved to %s", summary_path)
    logger.info("Summary: %s", json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
