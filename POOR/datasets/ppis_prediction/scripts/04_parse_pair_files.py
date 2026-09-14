#!/usr/bin/env python3
"""
解析 DIPS-Plus pair 文件，提取序列、长度、坐标可用性、标签等元数据。

用法：
    python scripts/parse_pair_files.py --index output/dips_plus_pair_index.csv --num_workers 16
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
        logging.FileHandler("logs/parse_pair_files.log"),
        logging.StreamHandler(sys.stdout),
    ],
)
logger = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = ROOT / "output"

# 可选依赖，运行前安装：pip install -r requirements.txt
try:
    import dill
except ImportError:
    dill = None

try:
    import hickle
except ImportError:
    hickle = None


def load_pair(path, fmt):
    """加载 pair 文件，返回 (df0, df1, sequences_dict, pos_idx, neg_idx) 五元组。"""
    if fmt == "dill":
        if dill is None:
            raise ImportError("dill package is required for .dill files")
        with open(path, "rb") as f:
            p = dill.load(f)
        return p.df0, p.df1, p.sequences, p.pos_idx, p.neg_idx
    elif fmt in ("hdf5", "h5"):
        if hickle is None:
            raise ImportError("hickle package is required for .hdf5 files")
        data = hickle.load(path)
        # data 是 hickled list：[complex_name, df0, df1, pos_idx, neg_idx, srcs, id, sequences]
        if not isinstance(data, (list, tuple)) or len(data) < 8:
            raise ValueError(f"Unexpected HDF5 structure: {type(data)}")
        return data[1], data[2], data[7], data[3], data[4]
    else:
        raise ValueError(f"Unsupported format: {fmt}")


def safe_seq(seq):
    """将序列转为字符串。"""
    if isinstance(seq, bytes):
        return seq.decode("utf-8", errors="ignore")
    return str(seq)


def parse_one(row):
    """解析单个 pair 文件。"""
    file_path = ROOT / row["file_path"]
    fmt = row["file_format"]
    result = {
        "unique_id": row["unique_id"],
        "pdb_id": row["pdb_id"],
        "model_id": row.get("model_id", ""),
        "pair_index": row.get("pair_index", ""),
        "file_path": row["file_path"],
        "file_format": fmt,
    }
    try:
        df0, df1, sequences, pos_idx, neg_idx = load_pair(file_path, fmt)

        # 序列
        result["receptor_seq"] = safe_seq(sequences.get("r_b", ""))
        result["ligand_seq"] = safe_seq(sequences.get("l_b", ""))
        result["receptor_len"] = len(result["receptor_seq"])
        result["ligand_len"] = len(result["ligand_seq"])

        # 链信息
        result["receptor_chain"] = ",".join(sorted(set(df0["chain"].dropna().astype(str).unique()))) if hasattr(df0, "columns") else ""
        result["ligand_chain"] = ",".join(sorted(set(df1["chain"].dropna().astype(str).unique()))) if hasattr(df1, "columns") else ""

        # 坐标可用性
        has_coords0 = {"x", "y", "z"}.issubset(set(df0.columns)) and len(df0) > 0
        has_coords1 = {"x", "y", "z"}.issubset(set(df1.columns)) and len(df1) > 0
        result["has_receptor_coords"] = has_coords0
        result["has_ligand_coords"] = has_coords1
        result["receptor_atoms"] = len(df0)
        result["ligand_atoms"] = len(df1)
        result["receptor_residues"] = df0["residue"].nunique() if has_coords0 and "residue" in df0.columns else 0
        result["ligand_residues"] = df1["residue"].nunique() if has_coords1 and "residue" in df1.columns else 0

        # 界面样本数
        result["pos_pairs"] = len(pos_idx) if hasattr(pos_idx, "__len__") else 0
        result["neg_pairs"] = len(neg_idx) if hasattr(neg_idx, "__len__") else 0

        # IDR 信息
        result["has_idr_annotations"] = bool(sequences.get("l_b_idr_annotations") and sequences.get("r_b_idr_annotations"))
        result["has_idr_propensities"] = bool(sequences.get("l_b_idr_propensities") and sequences.get("r_b_idr_propensities"))

        result["status"] = "success"
    except Exception as e:
        result["status"] = "failed"
        result["error"] = str(e)
    return result


def main():
    parser = argparse.ArgumentParser(description="Parse DIPS-Plus pair files")
    parser.add_argument("--index", default="output/dips_plus_pair_index.csv", help="pair 索引 CSV")
    parser.add_argument("--num_workers", type=int, default=16, help="并发 worker 数")
    parser.add_argument("--limit", type=int, default=0, help="仅处理前 N 条（测试用）")
    args = parser.parse_args()

    index_path = ROOT / args.index
    if not index_path.exists():
        logger.error("Index file not found: %s", index_path)
        sys.exit(1)

    records = []
    with open(index_path, newline="") as f:
        reader = csv.DictReader(f)
        records = list(reader)
    if args.limit > 0:
        records = records[:args.limit]
    logger.info("Loaded %d records from %s", len(records), index_path)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    results = []
    with ProcessPoolExecutor(max_workers=args.num_workers) as executor:
        futures = {executor.submit(parse_one, r): r for r in records}
        for i, future in enumerate(as_completed(futures)):
            if i % 1000 == 0:
                logger.info("Parsed %d / %d", i, len(records))
            try:
                results.append(future.result())
            except Exception as e:
                logger.error("Worker error: %s", e)

    # 写入元数据表
    if results:
        out_path = OUTPUT_DIR / "dips_plus_metadata.csv"
        with open(out_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)
        logger.info("Metadata saved to %s", out_path)

    # 统计
    status_counts = {}
    for r in results:
        status_counts[r["status"]] = status_counts.get(r["status"], 0) + 1
    summary = {
        "total": len(results),
        "status": status_counts,
    }
    summary_path = OUTPUT_DIR / "parse_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    logger.info("Parse summary: %s", summary)


if __name__ == "__main__":
    main()
