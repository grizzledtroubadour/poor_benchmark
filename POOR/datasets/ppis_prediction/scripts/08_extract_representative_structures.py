#!/usr/bin/env python3
"""
从原始 pair 文件中为每条代表链提取一个代表性结构。

确保 8,850 条去冗余后的代表链都有对应的 PDB 结构文件。

步骤：
1. 从 chain_centric_dedup_metadata.csv 读取 8,850 条代表链的 seq_hash 与序列。
2. 重新扫描 42,112 条 pair 元数据，为每条唯一序列建立原子数最多的最佳出现
   （seq_hash -> (unique_id, role)）。
3. 对每个代表链，加载对应 pair 文件的 df0/df1，**按序列匹配**选择与目标序列
   对应的 df（DIPS-Plus pair 文件中 df0/df1 与 receptor/ligand 可能互换，
   不能按 role 直取——该逻辑原属一次性补丁 fix_representative_structures.py，
   已融入本脚本；其 resname 列格式 bug 不带入，见 df_to_pdb 注释），
   提取坐标写入单链 PDB（统一链 ID "A"）。

输出：
    pdbs/{seq_hash}.pdb
    output/chain_centric_dedup_metadata.csv（回填结构来源字段）
"""

import argparse
import csv
import hashlib
import logging
import sys
from collections import defaultdict
from multiprocessing import Pool
from functools import partial
from pathlib import Path

import hickle
import numpy as np
import pandas as pd
from tqdm import tqdm

Path("logs").mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("logs/extract_representative_structures.log"),
        logging.StreamHandler(sys.stdout),
    ],
)
logger = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = ROOT / "output"
REP_STRUCT_DIR = ROOT / "pdbs"


def seq_hash(seq):
    return hashlib.md5(seq.encode("utf-8")).hexdigest()


def load_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def load_pair(file_path):
    """加载 HDF5/dill pair 文件，返回 (df0, df1)。"""
    data = hickle.load(str(file_path))
    if isinstance(data, (list, tuple)) and len(data) >= 3:
        return data[1], data[2]
    return data.df0, data.df1


def df_to_pdb(df, out_path, chain_id="A"):
    """将原子 DataFrame 写入 PDB 文件。"""
    records = []
    atom_idx = 1
    for _, row in df.iterrows():
        atom_name = str(row.get("atom_name", "X")).strip()
        resname = str(row.get("resname", "UNK")).strip()
        residue = row.get("residue", 1)
        try:
            resseq = int(residue)
        except (ValueError, TypeError):
            resseq = 1
        x = float(row.get("x", 0.0))
        y = float(row.get("y", 0.0))
        z = float(row.get("z", 0.0))
        element = str(row.get("element", atom_name[0] if atom_name else "X")).strip()

        atom_name_fmt = f"{atom_name:>3} " if len(atom_name) <= 3 else f"{atom_name:>4}"
        # NOTE: resname must be exactly 3 chars; a trailing space here used to
        # shift chain/resseq one column right of the PDB standard, making the
        # files unparseable by Bio.PDB.PDBParser (fixed post-hoc by the
        # archived 18_fix_pdb_column_layout.py, now in output/scripts_archive/).
        resname_fmt = f"{resname:>3}"

        if resseq > 9999:
            resseq = 9999

        # resseq 之后必须留 4 列（iCode col 27 + 空白 cols 28-30），使 x 落在
        # 标准 cols 31-38；少一列会把坐标整体左移（负坐标百位会丢失负号）。
        line = (
            f"ATOM  {atom_idx:5d} {atom_name_fmt} {resname_fmt} {chain_id}{resseq:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}  \n"
        )
        records.append(line)
        atom_idx += 1
        if atom_idx > 99999:
            atom_idx = 1

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        f.writelines(records)
        f.write("END\n")


def extract_sequence_from_df(df):
    """从原子 DataFrame 中提取氨基酸序列（CA 原子、按 residue 去重排序）。"""
    RESIDUE_3TO1 = {
        "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
        "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
        "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
        "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
    }
    ca_df = df[df["atom_name"].str.strip() == "CA"]
    if len(ca_df) == 0:
        return ""
    ca_df = ca_df.drop_duplicates(subset=["residue"]).sort_values("residue")
    return "".join(RESIDUE_3TO1.get(str(r).strip(), "X") for r in ca_df["resname"])


def select_chain_df(df0, df1, target_seq):
    """按序列匹配选择与目标序列对应的 df。

    DIPS-Plus pair 文件中 df0/df1 与元数据 receptor/ligand 的对应关系可能
    互换，必须按序列内容选择：完全匹配优先，其次唯一长度匹配，最后按序列
    相似度。返回 (df, matched_role)。
    """
    seq0 = extract_sequence_from_df(df0)
    seq1 = extract_sequence_from_df(df1)

    if seq0 == target_seq:
        return df0, "df0"
    if seq1 == target_seq:
        return df1, "df1"
    if len(seq0) == len(target_seq) and len(seq1) != len(target_seq):
        return df0, "df0_by_length"
    if len(seq1) == len(target_seq) and len(seq0) != len(target_seq):
        return df1, "df1_by_length"

    def similarity(a, b):
        if len(a) == 0 or len(b) == 0:
            return 0.0
        min_len = min(len(a), len(b))
        matches = sum(1 for i in range(min_len) if a[i] == b[i])
        return matches / max(len(a), len(b))

    sim0 = similarity(seq0, target_seq)
    sim1 = similarity(seq1, target_seq)
    if sim0 >= sim1:
        return df0, f"df0_by_sim_{sim0:.2f}"
    return df1, f"df1_by_sim_{sim1:.2f}"


def find_best_occurrences(metadata):
    """为每个唯一序列找到坐标最完整的出现（原子数最多）。"""
    best_occurrence = {}
    for row in metadata:
        uid = row["unique_id"]
        file_path = ROOT / row["file_path"]
        r_seq = row.get("receptor_seq", "")
        l_seq = row.get("ligand_seq", "")

        for seq, role, atom_key in [
            (r_seq, "receptor", "receptor_atoms"),
            (l_seq, "ligand", "ligand_atoms"),
        ]:
            if not seq:
                continue
            h = seq_hash(seq)
            try:
                atom_count = int(row.get(atom_key, 0))
            except (ValueError, TypeError):
                atom_count = 0

            occ = {
                "unique_id": uid,
                "file_path": str(file_path),
                "role": role,
                "seq": seq,
                "atom_count": atom_count,
            }

            if h not in best_occurrence or atom_count > best_occurrence[h]["atom_count"]:
                best_occurrence[h] = occ
    return best_occurrence


def extract_structure_for_rep(rep_info, occurrence_map):
    """为代表链提取结构，返回 (rep_hash, status, occurrence_info)。"""
    rep_hash = rep_info["seq_hash"]
    target_seq = rep_info.get("sequence", "")
    out_path = REP_STRUCT_DIR / f"{rep_hash}.pdb"

    occ = occurrence_map.get(rep_hash)
    if occ is None:
        return rep_hash, "no_occurrence", None

    file_path = Path(occ["file_path"])
    if not file_path.exists():
        return rep_hash, "missing_pair_file", occ

    try:
        df0, df1 = load_pair(file_path)
    except Exception as e:
        return rep_hash, f"load_error: {e}", occ

    # 按序列匹配选 df（df0/df1 与 receptor/ligand 可能互换）
    if target_seq:
        df, matched_role = select_chain_df(df0, df1, target_seq)
    else:  # 元数据缺序列时回退到按 role 直取
        df = df0 if occ["role"] == "receptor" else df1
        matched_role = "by_role_fallback"

    try:
        # 如果是已存在的符号链接或文件，先删除
        if out_path.is_symlink() or out_path.exists():
            out_path.unlink()
        df_to_pdb(df, out_path, chain_id="A")
        written_seq = extract_sequence_from_df(df)
        if target_seq and len(written_seq) != len(target_seq):
            return rep_hash, f"length_mismatch: written={len(written_seq)} target={len(target_seq)}", occ
        return rep_hash, "ok", {**occ, "matched_role": matched_role}
    except Exception as e:
        return rep_hash, f"write_error: {e}", occ


def main():
    parser = argparse.ArgumentParser(description="Extract representative structures from all pair files")
    parser.add_argument("--rep_metadata", default="output/chain_centric_dedup_metadata.csv")
    parser.add_argument("--metadata", default="output/dips_plus_metadata.csv")
    parser.add_argument("--num_workers", type=int, default=16)
    args = parser.parse_args()

    rep_records = load_csv(ROOT / args.rep_metadata)
    logger.info("Loaded %d representative chain records", len(rep_records))

    metadata = load_csv(ROOT / args.metadata)
    logger.info("Loaded %d pair metadata records", len(metadata))

    logger.info("Finding best occurrence (most atoms) for each unique sequence...")
    best_occurrence = find_best_occurrences(metadata)
    logger.info("Found best occurrences for %d unique sequences", len(best_occurrence))

    REP_STRUCT_DIR.mkdir(parents=True, exist_ok=True)

    logger.info("Extracting representative structures...")
    if args.num_workers > 1:
        worker = partial(extract_structure_for_rep, occurrence_map=best_occurrence)
        with Pool(args.num_workers) as pool:
            results = list(tqdm(
                pool.imap(worker, rep_records),
                total=len(rep_records),
                desc="Extracting rep structures",
            ))
    else:
        results = []
        for rep_info in tqdm(rep_records, desc="Extracting rep structures"):
            results.append(extract_structure_for_rep(rep_info, occurrence_map=best_occurrence))

    ok = sum(1 for _, status, _ in results if status == "ok")
    failed = [(h, s) for h, s, _ in results if s != "ok"]
    logger.info("Extraction complete: %d ok, %d failed", ok, len(failed))
    if failed:
        for h, s in failed[:20]:
            logger.warning("Failed %s: %s", h, s)

    # 更新 dedup metadata 中的结构路径
    update_metadata_paths(results)


def update_metadata_paths(results):
    """更新 chain_centric_dedup_metadata.csv 中的结构路径为实际 PDB 文件。"""
    in_path = OUTPUT_DIR / "chain_centric_dedup_metadata.csv"
    rows = load_csv(in_path)

    # 构建 hash -> 最新结构来源信息
    hash_to_occ = {}
    for rep_hash, status, occ in results:
        if status == "ok" and occ is not None:
            hash_to_occ[rep_hash] = occ

    for r in rows:
        h = r["seq_hash"]
        pdb_path = REP_STRUCT_DIR / f"{h}.pdb"
        if pdb_path.exists():
            r["has_structure"] = "True"
            r["representative_structure"] = str(pdb_path)
            occ = hash_to_occ.get(h)
            if occ:
                r["structure_source_unique_id"] = occ["unique_id"]
                r["structure_source_role"] = f"{occ['role']}->{occ.get('matched_role', '')}"
        else:
            r["has_structure"] = "False"
            r["representative_structure"] = ""
            r["structure_source_unique_id"] = ""
            r["structure_source_role"] = ""

    fieldnames = list(rows[0].keys())
    with open(in_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    logger.info("Updated %s with actual structure paths", in_path)


if __name__ == "__main__":
    main()
