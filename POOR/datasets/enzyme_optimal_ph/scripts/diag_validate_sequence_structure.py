#!/usr/bin/env python3
"""
序列-结构一致性验证
- 从 .cif 或 .pdb 结构文件中提取氨基酸序列
- 与 metadata.csv 中的 uniprot_seq_cut 做全局比对（edlib NW）
- 计算 identity 与 coverage
"""
import os
import re
import json
import argparse
from pathlib import Path
from collections import defaultdict
from multiprocessing import Pool, cpu_count

import pandas as pd
import numpy as np
import edlib

from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import warnings
warnings.filterwarnings("ignore", category=PDBConstructionWarning)

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
OUTPUT_DIR = ROOT / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

META_PATH = DATA_DIR / "metadata.csv"
PDB_DIR = DATA_DIR / "dataset" / "DATASET" / "EnzyBase12k_pbs"
CIF_DIR = DATA_DIR / "dataset_mmcifs" / "EnzyBase12k_mmcifs"

AA3_TO_1 = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
    "MSE": "M", "SEC": "U", "PYL": "O",
}

pdb_parser = PDBParser(QUIET=True)
cif_parser = MMCIFParser(QUIET=True)


def extract_seq_from_structure(struct_path):
    """从单个结构文件中提取第一条链的氨基酸序列"""
    ext = struct_path.suffix.lower()
    try:
        if ext == ".cif":
            structure = cif_parser.get_structure("id", str(struct_path))
        elif ext == ".pdb":
            structure = pdb_parser.get_structure("id", str(struct_path))
        else:
            return None, "unsupported_extension"
    except Exception as e:
        return None, f"parse_error:{e}"

    residues = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0].strip():  # 跳过异质残基（HETATM）
                    continue
                resname = residue.resname.upper()
                if resname in AA3_TO_1:
                    residues.append(AA3_TO_1[resname])
            # 只取第一条链
            break
        break

    if not residues:
        return None, "no_residues"

    seq = "".join(residues)
    return seq, "ok"


def clean_seq(seq):
    """只保留标准 20 种氨基酸"""
    return re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", seq.upper())


def align_sequences(struct_seq, meta_seq):
    """使用 edlib 全局比对，返回 identity 和 coverage"""
    struct_seq = clean_seq(struct_seq)
    meta_seq = clean_seq(meta_seq)

    if not struct_seq or not meta_seq:
        return 0.0, 0.0

    result = edlib.align(struct_seq, meta_seq, mode="NW", task="path")
    edit_distance = result["editDistance"]
    alignment = edlib.getNiceAlignment(result, struct_seq, meta_seq)

    query_aligned = alignment["query_aligned"]
    target_aligned = alignment["target_aligned"]

    # identity: 匹配位置数 / 比对长度（含 gap）
    matches = sum(1 for q, t in zip(query_aligned, target_aligned) if q == t and q != "-")
    align_len = len(query_aligned)
    identity = matches / align_len if align_len > 0 else 0.0

    # coverage: 结构序列长度 / 元数据序列长度
    coverage = len(struct_seq) / len(meta_seq) if len(meta_seq) > 0 else 0.0

    return round(identity, 4), round(coverage, 4)


def process_one(args):
    uniprot_id, meta_seq = args
    # 优先尝试 cif，再尝试 pdb
    cif_path = CIF_DIR / f"{uniprot_id}.cif"
    pdb_path = PDB_DIR / f"{uniprot_id}.pdb"

    struct_path = None
    file_type = None
    if cif_path.exists():
        struct_path = cif_path
        file_type = "cif"
    elif pdb_path.exists():
        struct_path = pdb_path
        file_type = "pdb"
    else:
        return {
            "uniprot_id": uniprot_id,
            "file_type": None,
            "struct_seq_len": None,
            "meta_seq_len": len(clean_seq(meta_seq)),
            "identity": None,
            "coverage": None,
            "status": "file_not_found",
        }

    struct_seq, status = extract_seq_from_structure(struct_path)
    if struct_seq is None:
        return {
            "uniprot_id": uniprot_id,
            "file_type": file_type,
            "struct_seq_len": None,
            "meta_seq_len": len(clean_seq(meta_seq)),
            "identity": None,
            "coverage": None,
            "status": status,
        }

    identity, coverage = align_sequences(struct_seq, meta_seq)
    return {
        "uniprot_id": uniprot_id,
        "file_type": file_type,
        "struct_seq_len": len(clean_seq(struct_seq)),
        "meta_seq_len": len(clean_seq(meta_seq)),
        "identity": identity,
        "coverage": coverage,
        "status": "ok",
    }


def main(n_workers=None):
    print(f"[INFO] Loading metadata: {META_PATH}")
    df = pd.read_csv(META_PATH)
    print(f"[INFO] Metadata rows: {len(df)}")

    tasks = [(row["uniprot_id"], row["uniprot_seq_cut"]) for _, row in df.iterrows()]

    n_workers = n_workers or min(cpu_count(), 16)
    print(f"[INFO] Using {n_workers} workers")

    results = []
    with Pool(n_workers) as pool:
        for i, res in enumerate(pool.imap_unordered(process_one, tasks), 1):
            results.append(res)
            if i % 1000 == 0:
                print(f"[INFO] Processed {i}/{len(tasks)}")

    results_df = pd.DataFrame(results)
    out_csv = OUTPUT_DIR / "sequence_structure_validation.csv"
    results_df.to_csv(out_csv, index=False)
    print(f"[INFO] Results saved to {out_csv}")

    # Summary stats
    ok_df = results_df[results_df["status"] == "ok"]
    summary = {
        "n_total": int(len(results_df)),
        "n_ok": int(len(ok_df)),
        "n_file_not_found": int((results_df["status"] == "file_not_found").sum()),
        "n_parse_error": int(results_df["status"].str.startswith("parse_error").sum()),
        "n_no_residues": int((results_df["status"] == "no_residues").sum()),
        "identity": {
            "mean": round(float(ok_df["identity"].mean()), 4),
            "median": round(float(ok_df["identity"].median()), 4),
            "min": round(float(ok_df["identity"].min()), 4),
            "max": round(float(ok_df["identity"].max()), 4),
            "q1": round(float(ok_df["identity"].quantile(0.25)), 4),
            "q3": round(float(ok_df["identity"].quantile(0.75)), 4),
            "gte_0.99": int((ok_df["identity"] >= 0.99).sum()),
            "gte_0.95": int((ok_df["identity"] >= 0.95).sum()),
            "gte_0.90": int((ok_df["identity"] >= 0.90).sum()),
            "lt_0.90": int((ok_df["identity"] < 0.90).sum()),
            "lt_0.50": int((ok_df["identity"] < 0.50).sum()),
        },
        "coverage": {
            "mean": round(float(ok_df["coverage"].mean()), 4),
            "median": round(float(ok_df["coverage"].median()), 4),
            "min": round(float(ok_df["coverage"].min()), 4),
            "max": round(float(ok_df["coverage"].max()), 4),
            "gte_0.95": int((ok_df["coverage"] >= 0.95).sum()),
            "gte_0.90": int((ok_df["coverage"] >= 0.90).sum()),
            "lt_0.90": int((ok_df["coverage"] < 0.90).sum()),
        },
    }

    summary_path = OUTPUT_DIR / "sequence_structure_validation_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[INFO] Summary saved to {summary_path}")

    print("\n===== Validation Summary =====")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--n_workers", type=int, default=None)
    args = parser.parse_args()
    main(args.n_workers)
