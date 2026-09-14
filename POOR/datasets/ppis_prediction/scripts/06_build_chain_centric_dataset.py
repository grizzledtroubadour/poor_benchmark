#!/usr/bin/env python3
"""
基于原始 42,112 条 pair 构建 chain-centric 数据集。

核心思想：
- 以唯一氨基酸序列作为链的标识。
- 对每条唯一序列，聚合它在所有 pair 中出现时的界面残基（取并集）。
- 得到该序列的"全结合位点"标签。

输出：
    output/chain_centric_metadata.csv
    output/chain_interface_labels/{seq_hash}_labels.npz
    output/chain_centric_summary.json
"""

import argparse
import csv
import hashlib
import json
import logging
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

Path("logs").mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("logs/build_chain_centric_dataset.log"),
        logging.StreamHandler(sys.stdout),
    ],
)
logger = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = ROOT / "output"
INTERFACE_LABEL_DIR = OUTPUT_DIR / "interface_labels"
CHAIN_LABEL_DIR = OUTPUT_DIR / "chain_interface_labels"


def seq_hash(seq):
    """基于序列生成稳定哈希。"""
    return hashlib.md5(seq.encode("utf-8")).hexdigest()


def load_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def load_interface_labels(unique_id):
    """加载 pair 的受体/配体界面标签。"""
    label_path = INTERFACE_LABEL_DIR / f"{unique_id}_labels.npz"
    if not label_path.exists():
        return None, None
    try:
        data = np.load(label_path)
        return data["receptor_labels"], data["ligand_labels"]
    except Exception as e:
        logger.warning("Failed to load labels for %s: %s", unique_id, e)
        return None, None


def build_chain_centric_dataset(metadata):
    """聚合所有链的界面标签。"""
    # chain_data[seq_hash] = {
    #   'seq': sequence,
    #   'interface_positions': set of 0-based positions,
    #   'occurrences': list of (unique_id, role, pdb_id),
    #   'as_receptor': count,
    #   'as_ligand': count,
    # }
    chain_data = defaultdict(lambda: {
        "seq": "",
        "interface_positions": set(),
        "occurrences": [],
        "as_receptor": 0,
        "as_ligand": 0,
    })

    total = len(metadata)
    skipped = 0

    for idx, row in enumerate(metadata, 1):
        uid = row["unique_id"]
        pdb_id = row.get("pdb_id", "")
        r_seq = row.get("receptor_seq", "")
        l_seq = row.get("ligand_seq", "")

        r_labels, l_labels = load_interface_labels(uid)
        if r_labels is None or l_labels is None:
            skipped += 1
            continue

        # 受体链
        if r_seq:
            h = seq_hash(r_seq)
            info = chain_data[h]
            info["seq"] = r_seq
            info["occurrences"].append({"unique_id": uid, "pdb_id": pdb_id, "role": "receptor"})
            info["as_receptor"] += 1
            r_interface = set(np.where(r_labels == 1)[0].tolist())
            info["interface_positions"].update(r_interface)

        # 配体链
        if l_seq:
            h = seq_hash(l_seq)
            info = chain_data[h]
            info["seq"] = l_seq
            info["occurrences"].append({"unique_id": uid, "pdb_id": pdb_id, "role": "ligand"})
            info["as_ligand"] += 1
            l_interface = set(np.where(l_labels == 1)[0].tolist())
            info["interface_positions"].update(l_interface)

        if idx % 5000 == 0:
            logger.info("Processed %d / %d pairs", idx, total)

    logger.info("Processed %d pairs, skipped %d", total, skipped)
    return chain_data


def save_chain_labels(chain_data):
    """保存每条唯一序列的聚合标签为 NPZ。"""
    CHAIN_LABEL_DIR.mkdir(parents=True, exist_ok=True)

    for h, info in chain_data.items():
        seq = info["seq"]
        labels = np.zeros(len(seq), dtype=np.int8)
        if info["interface_positions"]:
            labels[list(info["interface_positions"])] = 1
        out_path = CHAIN_LABEL_DIR / f"{h}_labels.npz"
        np.savez_compressed(out_path, labels=labels, seq=seq)


def save_chain_metadata(chain_data):
    """保存 chain-centric 元数据表。"""
    rows = []
    for h, info in chain_data.items():
        seq = info["seq"]
        interface_pos = info["interface_positions"]
        occurrences = info["occurrences"]
        rows.append({
            "seq_hash": h,
            "sequence": seq,
            "length": len(seq),
            "num_occurrences": len(occurrences),
            "num_interface_residues": len(interface_pos),
            "interface_ratio": len(interface_pos) / len(seq) if seq else 0.0,
            "as_receptor": info["as_receptor"],
            "as_ligand": info["as_ligand"],
            "sample_unique_ids": ";".join(o["unique_id"] for o in occurrences[:10]) + (";..." if len(occurrences) > 10 else ""),
            "pdb_ids": ";".join(sorted(set(o["pdb_id"] for o in occurrences)))[:500],
        })

    # 按出现次数降序，再按序列长度排序
    rows.sort(key=lambda r: (-r["num_occurrences"], -r["length"]))

    out_path = OUTPUT_DIR / "chain_centric_metadata.csv"
    fieldnames = [
        "seq_hash", "sequence", "length", "num_occurrences",
        "num_interface_residues", "interface_ratio",
        "as_receptor", "as_ligand",
        "sample_unique_ids", "pdb_ids",
    ]
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    logger.info("Saved chain-centric metadata: %d rows to %s", len(rows), out_path)
    return rows


def save_summary(chain_data, metadata_rows):
    """保存统计摘要。"""
    total_occurrences = sum(len(info["occurrences"]) for info in chain_data.values())
    multi_occurrence = sum(1 for info in chain_data.values() if len(info["occurrences"]) > 1)
    only_receptor = sum(1 for info in chain_data.values() if info["as_receptor"] > 0 and info["as_ligand"] == 0)
    only_ligand = sum(1 for info in chain_data.values() if info["as_ligand"] > 0 and info["as_receptor"] == 0)
    both_roles = sum(1 for info in chain_data.values() if info["as_receptor"] > 0 and info["as_ligand"] > 0)

    interface_ratios = [r["interface_ratio"] for r in metadata_rows]

    summary = {
        "total_pairs": len(metadata_rows),  # placeholder, corrected below
        "total_pairs_input": sum(1 for r in metadata_rows),  # not accurate
        "unique_chains": len(chain_data),
        "total_chain_occurrences": total_occurrences,
        "multi_occurrence_chains": multi_occurrence,
        "only_receptor": only_receptor,
        "only_ligand": only_ligand,
        "both_roles": both_roles,
        "interface_ratio_mean": float(np.mean(interface_ratios)) if interface_ratios else 0.0,
        "interface_ratio_median": float(np.median(interface_ratios)) if interface_ratios else 0.0,
        "interface_ratio_max": float(max(interface_ratios)) if interface_ratios else 0.0,
    }

    summary_path = OUTPUT_DIR / "chain_centric_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    logger.info("Summary: %s", json.dumps(summary, indent=2))


def main():
    parser = argparse.ArgumentParser(description="Build chain-centric PPIS dataset from 42,112 pairs")
    parser.add_argument("--metadata", default="output/dips_plus_metadata.csv")
    args = parser.parse_args()

    metadata = load_csv(ROOT / args.metadata)
    logger.info("Loaded %d metadata records", len(metadata))

    chain_data = build_chain_centric_dataset(metadata)
    logger.info("Built chain-centric dataset with %d unique chains", len(chain_data))

    save_chain_labels(chain_data)
    logger.info("Saved chain interface labels to %s", CHAIN_LABEL_DIR)

    metadata_rows = save_chain_metadata(chain_data)

    # 修正 summary 中的 total_pairs
    summary = {
        "total_pairs": len(metadata),
        "unique_chains": len(chain_data),
        "total_chain_occurrences": sum(len(info["occurrences"]) for info in chain_data.values()),
        "multi_occurrence_chains": sum(1 for info in chain_data.values() if len(info["occurrences"]) > 1),
        "only_receptor": sum(1 for info in chain_data.values() if info["as_receptor"] > 0 and info["as_ligand"] == 0),
        "only_ligand": sum(1 for info in chain_data.values() if info["as_ligand"] > 0 and info["as_receptor"] == 0),
        "both_roles": sum(1 for info in chain_data.values() if info["as_receptor"] > 0 and info["as_ligand"] > 0),
        "interface_ratio_mean": float(np.mean([r["interface_ratio"] for r in metadata_rows])) if metadata_rows else 0.0,
        "interface_ratio_median": float(np.median([r["interface_ratio"] for r in metadata_rows])) if metadata_rows else 0.0,
        "interface_ratio_max": float(max([r["interface_ratio"] for r in metadata_rows])) if metadata_rows else 0.0,
    }
    summary_path = OUTPUT_DIR / "chain_centric_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    logger.info("Final summary: %s", json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
