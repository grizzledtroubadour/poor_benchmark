#!/usr/bin/env python3
"""
为 Chain-Centric 测试集标记 OOD 属性。

使用已有的 OOD 定义和阈值：
- OOD_LowHomology_30/40/.../90
- OOD_Orphan
- OOD_NewFold_0.3/0.4/.../0.9
- OOD_ExtremeLength
- OOD_IDR

参考集：Chain-Centric Train + Val。

> SUPERSEDED 说明（2026-08 统一口径后）：
> - 序列同源段（compute_homology_ood，mmseqs2 search）与结构相似段
>   （compute_newfold_ood，foldseek）已被 `.skills/homology-ood-annotation/`
>   （seq_homology_ood.py / struct_homology_ood.py，easy-search 口径）取代，
>   splits 中的现值以 skill 重算为准；本脚本这两段仅作历史任务实现保留，
>   重跑口径不可混用。
> - 列名后经项目统一 schema 迁移改名：OOD_LowHomology_* -> seq_Redundancy_*、
>   OOD_NewFold_* -> TM-score_*、OOD_ExtremeLength -> OOD_ExtremeLong/Short。
> - OOD_IDR 段保留有效：使用 DIPS-Plus 自带 flDPnn 界面残基无序比例，
>   与 ood-annotation-toolkit 的 metapredict 口径不同，未被取代。

输出：
    output/chain_centric_splits/chain_centric_test_ood.csv
    output/ood_chain_centric_summary.json
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
import warnings
from collections import defaultdict
from pathlib import Path

import hickle
import numpy as np

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning)

Path("logs").mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("logs/mark_ood_chain_centric.log"),
        logging.StreamHandler(sys.stdout),
    ],
)
logger = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = ROOT / "output"
SPLIT_DIR = ROOT / "output" / "chain_centric_splits"
STRUCT_DIR = ROOT / "pdbs"

HM_THRESHOLDS = [90, 80, 70, 60, 50, 40, 30]
NF_THRESHOLDS = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3]

RESIDUE_3TO1 = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
}


def load_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def seq_hash(seq):
    return hashlib.md5(seq.encode("utf-8")).hexdigest()


def write_fasta(records, out_path, seq_field="sequence", id_field="seq_hash"):
    with open(out_path, "w") as f:
        for r in records:
            seq = r.get(seq_field, "")
            if seq:
                f.write(f">{r[id_field]}\n{seq}\n")


def run_mmseqs_search(query_fasta, target_fasta, tmp_dir, num_workers=32):
    os.makedirs(tmp_dir, exist_ok=True)
    query_db = os.path.join(tmp_dir, "query_db")
    target_db = os.path.join(tmp_dir, "target_db")
    aln = os.path.join(tmp_dir, "aln")
    m8 = os.path.join(tmp_dir, "result.m8")

    subprocess.run(["mmseqs", "createdb", query_fasta, query_db], check=True, capture_output=True, text=True)
    subprocess.run(["mmseqs", "createdb", target_fasta, target_db], check=True, capture_output=True, text=True)
    subprocess.run(
        ["mmseqs", "search", query_db, target_db, aln, tmp_dir, "--threads", str(num_workers)],
        check=True, capture_output=True, text=True,
    )
    subprocess.run(["mmseqs", "convertalis", query_db, target_db, aln, m8], check=True, capture_output=True, text=True)
    return m8


def parse_m8(m8_path):
    max_pident = defaultdict(float)
    has_hit = defaultdict(bool)
    with open(m8_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 11:
                continue
            query = parts[0]
            try:
                pident_pct = float(parts[2]) * 100.0
            except ValueError:
                continue
            has_hit[query] = True
            if pident_pct > max_pident[query]:
                max_pident[query] = pident_pct
    return max_pident, has_hit


def compute_homology_ood(test, trainval, tmp_root, num_workers=32):
    query_fasta = os.path.join(tmp_root, "query.fasta")
    target_fasta = os.path.join(tmp_root, "target.fasta")
    write_fasta(test, query_fasta)
    write_fasta(trainval, target_fasta)

    m8 = run_mmseqs_search(query_fasta, target_fasta, os.path.join(tmp_root, "search"), num_workers)
    max_pident, has_hit = parse_m8(m8)

    results = {}
    for r in test:
        h = r["seq_hash"]
        pident = max_pident.get(h, 0.0)
        flags = {f"OOD_LowHomology_{thr}": pident < thr for thr in HM_THRESHOLDS}
        results[h] = {
            "max_pident": pident,
            "OOD_Orphan": not has_hit.get(h, False),
            **flags,
        }
    return results


def link_structures(records, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    for r in records:
        h = r["seq_hash"]
        src = Path(r["representative_structure"])
        dst = out_dir / f"{h}.pdb"
        if dst.exists() or dst.is_symlink():
            dst.unlink()
        if src.exists():
            os.symlink(os.path.abspath(str(src)), dst)


def run_foldseek_search(query_dir, target_dir, tmp_dir, num_workers=32):
    os.makedirs(tmp_dir, exist_ok=True)
    query_db = os.path.join(tmp_dir, "query_db")
    target_db = os.path.join(tmp_dir, "target_db")
    aln = os.path.join(tmp_dir, "aln")
    m8 = os.path.join(tmp_dir, "result.m8")

    subprocess.run(["foldseek", "createdb", query_dir, query_db], check=True, capture_output=True, text=True)
    subprocess.run(["foldseek", "createdb", target_dir, target_db], check=True, capture_output=True, text=True)
    subprocess.run(
        ["foldseek", "search", query_db, target_db, aln, tmp_dir, "--threads", str(num_workers), "-a", "1"],
        check=True, capture_output=True, text=True,
    )
    subprocess.run(
        ["foldseek", "convertalis", query_db, target_db, aln, m8,
         "--format-output", "query,target,alntmscore,qtmscore,ttmscore,evalue"],
        check=True, capture_output=True, text=True,
    )
    return m8


def parse_foldseek_max_tmscore(m8_path):
    max_tm = defaultdict(float)
    with open(m8_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue
            query, target = parts[0], parts[1]
            if query == target:
                continue
            try:
                tm = float(parts[2])
            except ValueError:
                continue
            if tm > max_tm[query]:
                max_tm[query] = tm
    return max_tm


def compute_newfold_ood(test, trainval, tmp_root, num_workers=32):
    query_dir = os.path.join(tmp_root, "query")
    target_dir = os.path.join(tmp_root, "target")
    link_structures(test, Path(query_dir))
    link_structures(trainval, Path(target_dir))

    m8 = run_foldseek_search(query_dir, target_dir, os.path.join(tmp_root, "foldseek"), num_workers)
    max_tm = parse_foldseek_max_tmscore(m8)

    results = {}
    for r in test:
        h = r["seq_hash"]
        tm = max_tm.get(h, 0.0)
        flags = {f"OOD_NewFold_{thr}": tm < thr for thr in NF_THRESHOLDS}
        results[h] = {"max_tmscore": tm, **flags}
    return results


def load_chain_labels(seq_hash):
    label_path = OUTPUT_DIR / "chain_interface_labels" / f"{seq_hash}_labels.npz"
    if not label_path.exists():
        return None
    try:
        data = np.load(label_path)
        return data["labels"]
    except Exception:
        return None


def build_pair_path_map():
    """构建 unique_id -> file_path 的映射，避免每次扫描整个 metadata。"""
    metadata_path = OUTPUT_DIR / "dips_plus_metadata.csv"
    mapping = {}
    with open(metadata_path) as f:
        for row in csv.DictReader(f):
            mapping[row["unique_id"]] = ROOT / row["file_path"]
    return mapping


def load_pair_idr_propensities(unique_id, role, path_map=None):
    """从原始 pair 文件加载指定链的 IDR propensities。"""
    if not unique_id:
        return None
    file_path = path_map.get(unique_id) if path_map else None
    if not file_path:
        return None
    try:
        data = hickle.load(str(file_path))
        if isinstance(data, (list, tuple)) and len(data) >= 8:
            seq_info = data[7]
            key = f"{role[0]}_b_idr_propensities"
            return seq_info.get(key)
    except Exception as e:
        logger.warning("Failed to load IDR for %s: %s", unique_id, e)
    return None


def compute_idr_ood(test_records, metadata_by_hash, path_map):
    results = {}
    for r in test_records:
        h = r["seq_hash"]
        labels = load_chain_labels(h)
        meta = metadata_by_hash.get(h, {})
        unique_id = meta.get("structure_source_unique_id", "")
        role = meta.get("structure_source_role", "")

        idr_props = load_pair_idr_propensities(unique_id, role, path_map)
        if idr_props is None or labels is None or len(idr_props) != len(labels):
            results[h] = {"idr_ratio": 0.0, "OOD_IDR": False}
            continue

        interface_indices = np.where(labels == 1)[0]
        if len(interface_indices) == 0:
            idr_ratio = 0.0
        else:
            idr_arr = np.array(idr_props)
            idr_ratio = float(np.mean(idr_arr[interface_indices] > 0.5))

        results[h] = {"idr_ratio": idr_ratio, "OOD_IDR": idr_ratio > 0.3}
    return results


def compute_extreme_length(rows, min_len=60, max_len=1000):
    results = {}
    for r in rows:
        length = int(r["length"])
        results[r["seq_hash"]] = {"OOD_ExtremeLength": length < min_len or length > max_len}
    return results


def main():
    parser = argparse.ArgumentParser(description="Mark OOD for chain-centric test set")
    parser.add_argument("--train", default=str(SPLIT_DIR / "chain_centric_train.csv"))
    parser.add_argument("--val", default=str(SPLIT_DIR / "chain_centric_val.csv"))
    parser.add_argument("--test", default=str(SPLIT_DIR / "chain_centric_test.csv"))
    parser.add_argument("--metadata", default="output/chain_centric_dedup_metadata.csv")
    parser.add_argument("--num_workers", type=int, default=32)
    args = parser.parse_args()

    train = load_csv(args.train)
    val = load_csv(args.val)
    test = load_csv(args.test)
    trainval = train + val

    logger.info("Loaded Train: %d, Val: %d, Test: %d", len(train), len(val), len(test))

    # 加载完整 metadata 以获取 source pair 信息
    metadata = load_csv(ROOT / args.metadata)
    metadata_by_hash = {r["seq_hash"]: r for r in metadata}

    # 1. LowHomology + Orphan
    with tempfile.TemporaryDirectory(prefix="ood_chain_hm_", dir=str(OUTPUT_DIR)) as tmp_root:
        homology = compute_homology_ood(test, trainval, tmp_root, args.num_workers)
    logger.info("Homology OOD computed for %d test chains", len(homology))

    # 2. NewFold
    with tempfile.TemporaryDirectory(prefix="ood_chain_nf_", dir=str(OUTPUT_DIR)) as tmp_root:
        newfold = compute_newfold_ood(test, trainval, tmp_root, args.num_workers)
    logger.info("NewFold OOD computed for %d test chains", len(newfold))

    # 3. IDR
    logger.info("Building pair path map ...")
    path_map = build_pair_path_map()
    logger.info("Pair path map ready: %d entries", len(path_map))
    idr = compute_idr_ood(test, metadata_by_hash, path_map)
    logger.info("IDR OOD computed for %d test chains", len(idr))

    # 4. ExtremeLength
    extreme = compute_extreme_length(test)

    # 5. 写入 CSV
    base_fields = list(test[0].keys()) if test else []
    ood_fields = []
    for thr in HM_THRESHOLDS:
        ood_fields.append(f"OOD_LowHomology_{thr}")
    ood_fields.append("OOD_Orphan")
    for thr in NF_THRESHOLDS:
        ood_fields.append(f"OOD_NewFold_{thr}")
    ood_fields.extend(["OOD_ExtremeLength", "OOD_IDR", "max_pident", "max_tmscore", "idr_ratio"])

    out_path = SPLIT_DIR / "chain_centric_test_ood.csv"
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=base_fields + ood_fields)
        writer.writeheader()
        for r in test:
            h = r["seq_hash"]
            hm = homology.get(h, {})
            nf = newfold.get(h, {})
            idr_info = idr.get(h, {})
            ex = extreme.get(h, {})

            out = dict(r)
            for thr in HM_THRESHOLDS:
                out[f"OOD_LowHomology_{thr}"] = str(hm.get(f"OOD_LowHomology_{thr}", False))
            out["OOD_Orphan"] = str(hm.get("OOD_Orphan", False))
            for thr in NF_THRESHOLDS:
                out[f"OOD_NewFold_{thr}"] = str(nf.get(f"OOD_NewFold_{thr}", False))
            out["OOD_ExtremeLength"] = str(ex.get("OOD_ExtremeLength", False))
            out["OOD_IDR"] = str(idr_info.get("OOD_IDR", False))
            out["max_pident"] = f"{hm.get('max_pident', 0.0):.2f}"
            out["max_tmscore"] = f"{nf.get('max_tmscore', 0.0):.3f}"
            out["idr_ratio"] = f"{idr_info.get('idr_ratio', 0.0):.3f}"
            writer.writerow(out)
    logger.info("OOD test set saved to %s", out_path)

    # 6. 统计
    summary = {"total_test": len(test)}
    for thr in HM_THRESHOLDS:
        col = f"OOD_LowHomology_{thr}"
        summary[col] = sum(1 for r in test if homology.get(r["seq_hash"], {}).get(col, False))
    summary["OOD_Orphan"] = sum(1 for r in test if homology.get(r["seq_hash"], {}).get("OOD_Orphan", False))
    for thr in NF_THRESHOLDS:
        col = f"OOD_NewFold_{thr}"
        summary[col] = sum(1 for r in test if newfold.get(r["seq_hash"], {}).get(col, False))
    summary["OOD_ExtremeLength"] = sum(1 for r in test if extreme.get(r["seq_hash"], {}).get("OOD_ExtremeLength", False))
    summary["OOD_IDR"] = sum(1 for r in test if idr.get(r["seq_hash"], {}).get("OOD_IDR", False))

    summary_path = OUTPUT_DIR / "ood_chain_centric_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    logger.info("Summary: %s", json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
