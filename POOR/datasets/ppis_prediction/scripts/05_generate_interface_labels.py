#!/usr/bin/env python3
"""
为 DIPS-Plus pair 文件生成残基级界面标签。

默认基于 pair 文件中预采样好的 pos_idx（正例原子对）推导界面残基：
若某残基在 pos_idx 中出现过，则标记为 1，否则为 0。

用法：
    python scripts/generate_interface_labels.py --metadata output/dips_plus_metadata.csv --num_workers 32
"""

import argparse
import csv
import json
import logging
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

Path("logs").mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("logs/generate_interface_labels.log"),
        logging.StreamHandler(sys.stdout),
    ],
)
logger = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = ROOT / "output"


try:
    import numpy as np
except ImportError:
    np = None

try:
    import dill
except ImportError:
    dill = None

try:
    import hickle
except ImportError:
    hickle = None


def load_pair(path, fmt):
    """加载 pair 文件，返回 (df0, df1, sequences, pos_idx, neg_idx)。"""
    if fmt == "dill":
        if dill is None:
            raise ImportError("dill required")
        with open(path, "rb") as f:
            p = dill.load(f)
        return p.df0, p.df1, p.sequences, p.pos_idx, p.neg_idx
    elif fmt in ("hdf5", "h5"):
        if hickle is None:
            raise ImportError("hickle required")
        data = hickle.load(path)
        return data[1], data[2], data[7], data[3], data[4]
    else:
        raise ValueError(f"Unsupported format: {fmt}")


def residue_labels_from_atom_pairs(df, atom_indices, max_residues):
    """根据 atom_indices 在该链中出现的位置，生成残基级标签。"""
    labels = np.zeros(max_residues, dtype=int)
    if len(atom_indices) == 0:
        return labels

    # 将 atom_indices 映射到 residue 编号（按 DataFrame 中的顺序）
    residues = df["residue"].iloc[atom_indices].values
    # 残基在序列中的位置：按首次出现顺序编号 0..N-1
    unique_residues = df["residue"].unique()
    res_to_idx = {res: i for i, res in enumerate(unique_residues)}

    for res in residues:
        idx = res_to_idx.get(res)
        if idx is not None and 0 <= idx < max_residues:
            labels[idx] = 1
    return labels


def process_one(row):
    file_path = ROOT / row["file_path"]
    fmt = row["file_format"]
    result = {
        "unique_id": row["unique_id"],
        "pdb_id": row["pdb_id"],
    }
    try:
        df0, df1, sequences, pos_idx, neg_idx = load_pair(file_path, fmt)

        r_seq = sequences.get("r_b", "")
        l_seq = sequences.get("l_b", "")
        r_len = len(r_seq) if isinstance(r_seq, str) else len(r_seq)
        l_len = len(l_seq) if isinstance(l_seq, str) else len(l_seq)

        # 从 pos_idx 推导界面残基
        r_labels = residue_labels_from_atom_pairs(df0, pos_idx[:, 0], r_len)
        l_labels = residue_labels_from_atom_pairs(df1, pos_idx[:, 1], l_len)

        result["receptor_len"] = int(r_len)
        result["ligand_len"] = int(l_len)
        result["receptor_interface_count"] = int(r_labels.sum())
        result["ligand_interface_count"] = int(l_labels.sum())
        result["receptor_interface_ratio"] = round(float(r_labels.mean()), 4)
        result["ligand_interface_ratio"] = round(float(l_labels.mean()), 4)
        result["status"] = "success"

        # 保存独立标签文件
        label_dir = OUTPUT_DIR / "interface_labels"
        label_dir.mkdir(parents=True, exist_ok=True)
        np.savez(
            label_dir / f"{row['unique_id']}_labels.npz",
            receptor_labels=r_labels,
            ligand_labels=l_labels,
        )
    except Exception as e:
        result["status"] = "failed"
        result["error"] = str(e)
    return result


def main():
    parser = argparse.ArgumentParser(description="Generate interface labels for DIPS-Plus")
    parser.add_argument("--metadata", default="output/dips_plus_metadata.csv", help="元数据 CSV")
    parser.add_argument("--num_workers", type=int, default=32, help="并发 worker 数")
    args = parser.parse_args()

    metadata_path = ROOT / args.metadata
    if not metadata_path.exists():
        logger.error("Metadata file not found: %s", metadata_path)
        sys.exit(1)

    records = []
    with open(metadata_path, newline="") as f:
        reader = csv.DictReader(f)
        records = [r for r in reader if r.get("status") in (None, "", "success")]
    logger.info("Loaded %d records", len(records))

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    results = []
    with ProcessPoolExecutor(max_workers=args.num_workers) as executor:
        futures = {executor.submit(process_one, r): r for r in records}
        for i, future in enumerate(as_completed(futures)):
            if i % 1000 == 0:
                logger.info("Processed %d / %d", i, len(records))
            try:
                results.append(future.result())
            except Exception as e:
                logger.error("Worker error: %s", e)

    # 保存标签分布
    if results:
        out_path = OUTPUT_DIR / "interface_label_distribution.csv"
        with open(out_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)
        logger.info("Label distribution saved to %s", out_path)

    status_counts = {}
    total_pos = 0
    total_len = 0
    for r in results:
        status_counts[r["status"]] = status_counts.get(r["status"], 0) + 1
        if r["status"] == "success":
            total_pos += int(r.get("receptor_interface_count", 0)) + int(r.get("ligand_interface_count", 0))
            total_len += int(r.get("receptor_len", 0)) + int(r.get("ligand_len", 0))

    summary = {
        "total": len(results),
        "status": status_counts,
        "overall_interface_ratio": round(total_pos / total_len, 4) if total_len else 0,
    }
    summary_path = OUTPUT_DIR / "interface_label_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    logger.info("Summary: %s", summary)


if __name__ == "__main__":
    main()
