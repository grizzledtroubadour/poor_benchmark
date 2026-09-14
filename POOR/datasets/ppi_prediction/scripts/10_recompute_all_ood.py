#!/usr/bin/env python3
"""
Recompute all OOD columns for the PPI test set after XTAL rebuild.

This script replaces the old scripts 22/18/19 and computes:
  1. seq_Redundancy_90..30 + OOD_Orphan (mmseqs2 easy-search -s 7.5, test vs train+val)
  2. TM-score_0.9..0.3 (Foldseek easy-search, test vs train+val)
  3. OOD_ExtremeShort / OOD_ExtremeLong (length-based)
  4. OOD_IDR (metapredict v3)
  5. OOD_NewEC_L4 / OOD_LongTail_EC_L4 (UniProt + SIFTS EC)
  6. Default column (False for ExtremeLong)
  7. Default=False rows: all OOD/seq/TM columns set to NaN except OOD_ExtremeLong

Also fixes:
  - --max-seqs 1 → 10 for mmseqs2 (get true max identity)
  - Column names aligned with current test CSV
  - Foldseek cache regenerated (test pool changed)
  - Path uses dynamic ROOT instead of hardcoded

NOTE (2026-08 统一整理):
  - §1 mmseqs2 / §2 Foldseek：[superseded] 正式重算以
    `.skills/homology-ood-annotation/scripts/seq_homology_ood.py` /
    `struct_homology_ood.py` 为准（easy-search 口径，缓存可 --m8-cache 重放）。
  - §4 IDR：[superseded] 统一工具为 `.skills/ood-annotation-toolkit/scripts/idr_ood.py`
    （本任务口径 --mode region-ratio --clean map --threshold 0.3）。
  - §5 EC：[superseded] NewEC/LongTailEC 统一为 all-not-in 口径，由
    `.skills/ec-function-ood-annotation/scripts/add_unified_newec.py` /
    `add_unified_longtail_ec.py` 落盘；本段产出的 OOD_LongTail_EC_L4 旧列已删除。
  - §3 长度列 / §6 Default / §7 Default=False 置 NaN 为本脚本保留的独有逻辑；
    NaN/InD 收尾现以 `.skills/ood-annotation-toolkit/scripts/finalize_ood_columns.py` 为准。

Inputs:
  - splits/ppi_{train,val,test}.csv
  - data/uniprot_sequences_all.fasta
  - output/structure_selection.csv
  - output/uniprot_ec_mapping.csv
  - data/sifts/sifts_chain_ec.tsv.gz
  - pdbs/ (Foldseek structure database)

Outputs:
  - Updated splits/ppi_test.csv with all OOD columns
  - output/ood_foldseek_tmscore.csv (regenerated cache)
  - output/idr_predictions.csv (regenerated cache)
  - output/ood_mmseqs_test_vs_trainval.m8 (mmseqs2 result cache)
"""

import gzip
import os
import re
import subprocess
import sys
import tempfile
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SPLITS_DIR = ROOT / "splits"
OUTPUT_DIR = ROOT / "output"
DATA_DIR = ROOT / "data"
PDB_DIR = ROOT / "pdbs"

FASTA_FILE = DATA_DIR / "uniprot_sequences_all.fasta"
STRUCTURE_SELECTION = OUTPUT_DIR / "structure_selection.csv"
EC_CSV = OUTPUT_DIR / "uniprot_ec_mapping.csv"
SIFTS_EC_TSV = ROOT.parent.parent / "data" / "sifts" / "sifts_chain_ec.tsv.gz"
FOLDSEEK_CACHE = OUTPUT_DIR / "ood_foldseek_tmscore.csv"
IDR_CACHE = OUTPUT_DIR / "idr_predictions.csv"
MMSEQS_CACHE = OUTPUT_DIR / "ood_mmseqs_test_vs_trainval.m8"

LOW_HOMOLOGY_THRESHOLDS = [90, 80, 70, 60, 50, 40, 30]
NEWFOLD_THRESHOLDS = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3]
EXTREME_SHORT = 60
EXTREME_LONG = 1000
IDR_DISORDER_THRESHOLD = 0.5
IDR_RATIO_THRESHOLD = 0.3
LONGTAIL_EC_THRESHOLD = 10

# EC level: use L4 (4-level EC, e.g. "1.1.1.1")
EC_LEVEL = 4


# ---------------------------------------------------------------------------
# 1. mmseqs2: seq_Redundancy + OOD_Orphan
# [superseded] 正式重算以 .skills/homology-ood-annotation/scripts/seq_homology_ood.py 为准
# ---------------------------------------------------------------------------

def write_unique_fasta(df, out_fasta):
    records = {}
    for side in ["1", "2"]:
        id_col = f"protein_{('A' if side=='1' else 'B')}_id"
        seq_col = f"aa_seq{side}"
        for pid, seq in zip(df[id_col], df[seq_col]):
            if pd.notna(pid) and pd.notna(seq):
                records[str(pid)] = str(seq)
    with open(out_fasta, "w") as f:
        for pid, seq in records.items():
            f.write(f">{pid}\n{seq}\n")
    return records


def run_mmseqs_search(query_fasta, target_fasta, output_m8, tmp_dir):
    # 2026-08-27 起统一为 easy-search -s 7.5（默认 e-value <= 1e-3 显著性门槛），
    # 与其他任务及 homology-ood-annotation skill 一致；旧式宽松参数
    # (--min-seq-id 0.0 -c 0.0 --cov-mode 0 -e 100) 会把 9-25 aa 弱局部比对计入
    # max pident，口径与其他任务不可比（统一重算备份
    # output/backup_before_easysearch_unify/）。旧缓存（旧参数产物）需删除后重跑。
    subprocess.run(
        ["mmseqs", "easy-search", str(query_fasta), str(target_fasta), str(output_m8),
         str(Path(tmp_dir) / "mmseqs_tmp"),
         "-s", "7.5",
         "--format-output", "query,target,pident,evalue"],
        check=True, capture_output=True)


def parse_m8(m8_file, proteins):
    max_identity = {p: 0.0 for p in proteins}
    with open(m8_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            query, pident = parts[0], float(parts[2])
            if pident > max_identity.get(query, 0.0):
                max_identity[query] = pident
    return max_identity


def compute_seq_redundancy(train, val, test):
    """Returns dict: protein_id -> max_identity vs train+val."""
    test_proteins = set()
    for side in ["A", "B"]:
        test_proteins.update(test[f"protein_{side}_id"].dropna().astype(str).unique())

    if MMSEQS_CACHE.exists():
        print("  Using cached mmseqs2 results")
        max_identity = {p: 0.0 for p in test_proteins}
        with open(MMSEQS_CACHE) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 3:
                    continue
                query, pident = parts[0], float(parts[2])
                if pident > max_identity.get(query, 0.0):
                    max_identity[query] = pident
        return max_identity

    train_val = pd.concat([train, val], ignore_index=True)
    with tempfile.TemporaryDirectory(prefix="ppi_ood_mmseqs_") as tmp_dir:
        qf = str(Path(tmp_dir) / "query.fasta")
        tf = str(Path(tmp_dir) / "target.fasta")
        m8 = str(MMSEQS_CACHE)
        write_unique_fasta(test, qf)
        write_unique_fasta(train_val, tf)
        print("  Running mmseqs2 search (test vs train+val) ...")
        run_mmseqs_search(qf, tf, m8, tmp_dir)
    return parse_m8(m8, test_proteins)


def add_seq_redundancy_cols(test, max_identity):
    test = test.copy()
    test["A_max_id"] = test["protein_A_id"].astype(str).map(lambda x: max_identity.get(x, 0.0))
    test["B_max_id"] = test["protein_B_id"].astype(str).map(lambda x: max_identity.get(x, 0.0))
    for thr in LOW_HOMOLOGY_THRESHOLDS:
        # seq_Redundancy_XX = True means max_identity < XX (i.e. low homology / OOD)
        test[f"seq_Redundancy_{thr}"] = (
            (test["A_max_id"] < thr) | (test["B_max_id"] < thr)
        )
    test["OOD_Orphan"] = (test["A_max_id"] == 0) | (test["B_max_id"] == 0)
    test = test.drop(columns=["A_max_id", "B_max_id"])
    return test


# ---------------------------------------------------------------------------
# 2. Foldseek: TM-score
# [superseded] 正式重算以 .skills/homology-ood-annotation/scripts/struct_homology_ood.py 为准
# ---------------------------------------------------------------------------

def build_protein_to_structure_map():
    sel = pd.read_csv(STRUCTURE_SELECTION)
    protein_to_file = {}
    for _, row in sel.iterrows():
        files = str(row["selected_files"]).split(";")
        pa, pb = str(row["protein_A_id"]), str(row["protein_B_id"])
        for pid, fn in [(pa, files[0] if len(files) > 0 else None),
                        (pb, files[1] if len(files) > 1 else None)]:
            if fn and fn != "nan":
                pdb_path = PDB_DIR / fn.strip()
                if pdb_path.exists() and pid not in protein_to_file:
                    protein_to_file[pid] = str(pdb_path)
    return protein_to_file


def compute_foldseek_tmscore(train, val, test):
    """Returns dict: protein_id -> max_tmscore vs train+val structures."""
    test_proteins = set()
    for side in ["A", "B"]:
        test_proteins.update(test[f"protein_{side}_id"].dropna().astype(str).unique())

    if FOLDSEEK_CACHE.exists():
        cache_df = pd.read_csv(FOLDSEEK_CACHE)
        cached = dict(zip(cache_df["protein_id"].astype(str),
                          cache_df["max_tmscore"].astype(float)))
        # Check if all test proteins are in cache
        missing = test_proteins - set(cached.keys())
        if not missing:
            print("  Using cached Foldseek TM-scores")
            return cached
        print(f"  Foldseek cache missing {len(missing)} proteins, regenerating...")

    print("  Building protein→structure map...")
    protein_to_file = build_protein_to_structure_map()

    train_val = pd.concat([train, val], ignore_index=True)
    train_val_proteins = set()
    for side in ["A", "B"]:
        train_val_proteins.update(
            train_val[f"protein_{side}_id"].dropna().astype(str).unique())

    test_with_struct = sorted([p for p in test_proteins if p in protein_to_file])
    tv_with_struct = sorted([p for p in train_val_proteins if p in protein_to_file])
    print(f"  Test proteins with structure: {len(test_with_struct)}/{len(test_proteins)}")
    print(f"  Train+Val proteins with structure: {len(tv_with_struct)}/{len(train_val_proteins)}")

    max_tmscore = {p: 0.0 for p in test_proteins}
    if not test_with_struct or not tv_with_struct:
        print("  Not enough structures, all TM-scores set to 0")
        return max_tmscore

    with tempfile.TemporaryDirectory(prefix="ppi_foldseek_") as tmp_dir:
        query_dir = os.path.join(tmp_dir, "query_dir")
        target_dir = os.path.join(tmp_dir, "target_dir")
        os.makedirs(query_dir, exist_ok=True)
        os.makedirs(target_dir, exist_ok=True)
        for p in test_with_struct:
            os.symlink(protein_to_file[p], os.path.join(query_dir, f"{p}.pdb"))
        for p in tv_with_struct:
            os.symlink(protein_to_file[p], os.path.join(target_dir, f"{p}.pdb"))

        output_tsv = os.path.join(tmp_dir, "foldseek_result.tsv")
        fs_tmp = os.path.join(tmp_dir, "fs_tmp")
        os.makedirs(fs_tmp, exist_ok=True)

        print("  Running Foldseek easy-search...")
        # 2026-08-27 起统一为 foldseek easy-search（默认 e-value <= 10 门槛，
        # 与其他任务及 homology-ood-annotation skill 一致）；旧式
        # `search -a -e inf` 保留全部弱 hit，会抬高 max alntmscore、压低低阈值
        # OOD（统一重算备份 output/backup_before_foldseek_easysearch_unify/）。
        subprocess.run(["foldseek", "easy-search", query_dir, target_dir,
                         output_tsv, fs_tmp,
                         "--format-output", "query,target,alntmscore"],
                       check=True, capture_output=True)

        if os.path.exists(output_tsv) and os.path.getsize(output_tsv) > 0:
            with open(output_tsv) as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) < 3:
                        continue
                    query_pid = Path(parts[0]).stem
                    tm = float(parts[2])
                    if tm > max_tmscore.get(query_pid, 0.0):
                        max_tmscore[query_pid] = tm

    # Save cache
    cache_df = pd.DataFrame([
        {"protein_id": pid, "max_tmscore": score}
        for pid, score in max_tmscore.items()
    ])
    cache_df.to_csv(FOLDSEEK_CACHE, index=False)
    print(f"  Cached {len(max_tmscore)} TM-scores to {FOLDSEEK_CACHE}")
    return max_tmscore


def add_tmscore_cols(test, max_tmscore):
    test = test.copy()
    test["A_tm"] = test["protein_A_id"].astype(str).map(lambda x: max_tmscore.get(x, 0.0))
    test["B_tm"] = test["protein_B_id"].astype(str).map(lambda x: max_tmscore.get(x, 0.0))
    for thr in NEWFOLD_THRESHOLDS:
        test[f"TM-score_{thr}"] = (test["A_tm"] < thr) | (test["B_tm"] < thr)
    test = test.drop(columns=["A_tm", "B_tm"])
    return test


# ---------------------------------------------------------------------------
# 3. Extreme length
# ---------------------------------------------------------------------------

def add_extreme_length_cols(test):
    test = test.copy()
    len1 = test["aa_seq1"].str.len()
    len2 = test["aa_seq2"].str.len()
    test["OOD_ExtremeShort"] = (len1 < EXTREME_SHORT) | (len2 < EXTREME_SHORT)
    test["OOD_ExtremeLong"] = (len1 > EXTREME_LONG) | (len2 > EXTREME_LONG)
    return test


# ---------------------------------------------------------------------------
# 4. IDR (metapredict v3)
# [superseded] 统一工具为 .skills/ood-annotation-toolkit/scripts/idr_ood.py
# （--mode region-ratio --clean map --threshold 0.3）
# ---------------------------------------------------------------------------

def clean_sequence(seq):
    seq = seq.upper()
    seq = seq.replace("B", "N").replace("Z", "Q").replace("J", "L")
    seq = seq.replace("U", "C").replace("O", "K")
    seq = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "A", seq)
    return seq


def find_continuous_regions(binary_mask):
    regions = []
    current = 0
    for flag in binary_mask:
        if flag:
            current += 1
        else:
            if current > 0:
                regions.append(current)
                current = 0
    if current > 0:
        regions.append(current)
    return regions


def compute_idr_flag(scores, seq_len):
    binary = np.asarray(scores) > IDR_DISORDER_THRESHOLD
    regions = find_continuous_regions(binary)
    if not regions:
        return 0
    max_ratio = max(regions) / seq_len
    return 1 if max_ratio > IDR_RATIO_THRESHOLD else 0


def compute_idr(test):
    try:
        import metapredict
    except ImportError:
        print("  metapredict not available, skipping IDR")
        return test

    print(f"  metapredict version: {metapredict.__version__}")

    all_seqs = set(test["aa_seq1"].dropna()) | set(test["aa_seq2"].dropna())
    clean_to_original = {}
    original_to_clean = {}
    for s in all_seqs:
        c = clean_sequence(s)
        original_to_clean[s] = c
        clean_to_original.setdefault(c, s)
    unique_clean = sorted(set(original_to_clean.values()))

    # Load cache
    clean_predictions = {}
    if IDR_CACHE.exists():
        cache_df = pd.read_csv(IDR_CACHE)
        for _, row in cache_df.iterrows():
            clean_predictions[row["sequence"]] = row["is_idr"]
        print(f"  Loaded {len(clean_predictions)} cached IDR predictions")

    remaining = [s for s in unique_clean if s not in clean_predictions]
    print(f"  Remaining to predict: {len(remaining)}")

    if remaining:
        fasta_dir = OUTPUT_DIR / "idr_fastas"
        fasta_dir.mkdir(exist_ok=True)
        BATCH = 1000
        seq_to_id = {s: f"seq_{i:07d}" for i, s in enumerate(remaining)}
        for batch_idx in range(0, len(remaining), BATCH):
            batch = remaining[batch_idx:batch_idx + BATCH]
            fasta_path = fasta_dir / f"batch_{batch_idx//BATCH:05d}.fasta"
            with open(fasta_path, "w") as fh:
                for s in batch:
                    fh.write(f">{seq_to_id[s]}\n{s}\n")
            result = metapredict.predict_disorder_fasta(
                str(fasta_path), show_progress_bar=False,
                invalid_sequence_action="convert")
            for _, (seq, scores) in result.items():
                clean_predictions[seq] = compute_idr_flag(scores, len(seq))
            # Save cache
            df_cache = pd.DataFrame(
                [(s, f) for s, f in clean_predictions.items()],
                columns=["sequence", "is_idr"])
            df_cache.to_csv(IDR_CACHE, index=False)
            print(f"    Batch {batch_idx//BATCH+1}: {len(batch)} seqs, "
                  f"total cached: {len(clean_predictions)}/{len(unique_clean)}")

    original_predictions = {}
    for orig, clean in original_to_clean.items():
        original_predictions[orig] = clean_predictions.get(clean, 0)

    test = test.copy()
    test["A_idr"] = test["aa_seq1"].map(original_predictions).fillna(0).astype(int)
    test["B_idr"] = test["aa_seq2"].map(original_predictions).fillna(0).astype(int)
    test["OOD_IDR"] = (test["A_idr"] == 1) | (test["B_idr"] == 1)
    test = test.drop(columns=["A_idr", "B_idr"])
    return test


# ---------------------------------------------------------------------------
# 5. EC OOD
# [superseded] 统一为 all-not-in 口径（.skills/ec-function-ood-annotation/scripts/
# add_unified_newec.py / add_unified_longtail_ec.py）；OOD_LongTail_EC_L4 旧列已删除
# ---------------------------------------------------------------------------

def load_uniprot_ec():
    df = pd.read_csv(EC_CSV)
    df["ec_numbers"] = df["ec_numbers"].fillna("")
    ec_map = {}
    for _, row in df.iterrows():
        ec_str = str(row["ec_numbers"]).strip()
        ec_map[row["uniprot_id"]] = {x.strip() for x in ec_str.split(";") if x.strip()} if ec_str else set()
    return ec_map


def load_sifts_ec(relevant_pdbs):
    sifts_ec = defaultdict(set)
    with gzip.open(SIFTS_EC_TSV, "rt") as f:
        next(f)
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                pdb_id = parts[0].lower()
                chain = parts[1]
                ec = parts[3]
                if pdb_id in relevant_pdbs:
                    sifts_ec[(pdb_id, chain)].add(ec)
    return sifts_ec


def truncate_ec(ec_set, level):
    """Truncate EC numbers to given level (e.g. 4 → '1.1.1.1')."""
    result = set()
    for ec in ec_set:
        parts = ec.split(".")
        if len(parts) >= level:
            result.append(".".join(parts[:level])) if hasattr(result, 'append') else result.add(".".join(parts[:level]))
        else:
            result.add(ec)  # keep partial EC as-is
    return result


def get_row_ecs(row, uniprot_ec, sifts_ec):
    """Get EC numbers for a row. Uses unique_id for BIO/XTAL, protein_id for random."""
    uid = str(row["unique_id"])
    pair_type = row.get("pair_type", "bio")
    ecs = set()

    if pair_type in ("bio", "xtal"):
        # Parse PDB/chain from unique_id
        for part in uid.split("--"):
            sub = part.split("_")
            if len(sub) < 4:
                continue
            pdb_id = sub[0].lower()
            chain_letter = sub[-2][0]
            uniprot_id = sub[-1]
            ecs.update(uniprot_ec.get(uniprot_id, set()))
            ecs.update(sifts_ec.get((pdb_id, chain_letter), set()))
    else:
        # random negative: only UniProt EC
        for side in ["A", "B"]:
            pid = row.get(f"protein_{side}_id")
            if pd.notna(pid):
                ecs.update(uniprot_ec.get(str(pid), set()))
    return ecs


def compute_ec_ood(train, val, test):
    uniprot_ec = load_uniprot_ec()
    print(f"  UniProt EC: {len(uniprot_ec)} IDs, "
          f"{sum(1 for v in uniprot_ec.values() if v)} with EC")

    # Collect PDB IDs for SIFTS
    relevant_pdbs = set()
    for df in [train, val, test]:
        for uid in df["unique_id"]:
            if str(uid).startswith("neg_"):
                continue
            for part in str(uid).split("--"):
                sub = part.split("_")
                if len(sub) >= 4:
                    relevant_pdbs.add(sub[0].lower())
    sifts_ec = load_sifts_ec(relevant_pdbs)
    print(f"  SIFTS EC: {len(sifts_ec)} PDB/chain pairs")

    # Build train+val EC frequency (L4)
    trainval_ec_counts = Counter()
    for df in [train, val]:
        for _, row in df.iterrows():
            ecs = get_row_ecs(row, uniprot_ec, sifts_ec)
            for ec in ecs:
                parts = ec.split(".")
                if len(parts) >= EC_LEVEL:
                    trainval_ec_counts[".".join(parts[:EC_LEVEL])] += 1

    trainval_ec_set = set(trainval_ec_counts.keys())
    longtail_ecs = {ec for ec, cnt in trainval_ec_counts.items() if cnt < LONGTAIL_EC_THRESHOLD}
    print(f"  Train+Val unique EC L4: {len(trainval_ec_set)}, long-tail: {len(longtail_ecs)}")

    new_ec_flags = []
    longtail_flags = []
    for _, row in test.iterrows():
        test_ecs = get_row_ecs(row, uniprot_ec, sifts_ec)
        test_ecs_l4 = set()
        for ec in test_ecs:
            parts = ec.split(".")
            if len(parts) >= EC_LEVEL:
                test_ecs_l4.add(".".join(parts[:EC_LEVEL]))
            else:
                test_ecs_l4.add(ec)
        new_ec_flags.append(bool(test_ecs_l4 and (test_ecs_l4 - trainval_ec_set)))
        longtail_flags.append(bool(test_ecs_l4 & longtail_ecs))

    test = test.copy()
    test["OOD_NewEC_L4"] = new_ec_flags
    test["OOD_LongTail_EC_L4"] = longtail_flags
    return test


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=== Loading splits ===")
    train = pd.read_csv(SPLITS_DIR / "ppi_train.csv")
    val = pd.read_csv(SPLITS_DIR / "ppi_val.csv")
    test = pd.read_csv(SPLITS_DIR / "ppi_test.csv")
    print(f"  train={len(train)}, val={len(val)}, test={len(test)}")

    # 1. seq_Redundancy + OOD_Orphan
    print("\n=== 1. mmseqs2: seq_Redundancy + OOD_Orphan ===")
    max_identity = compute_seq_redundancy(train, val, test)
    test = add_seq_redundancy_cols(test, max_identity)
    for thr in [90, 50, 30]:
        print(f"  seq_Redundancy_{thr}: {test[f'seq_Redundancy_{thr}'].sum()}")
    print(f"  OOD_Orphan: {test['OOD_Orphan'].sum()}")

    # 2. TM-score (Foldseek)
    print("\n=== 2. Foldseek: TM-score ===")
    # Delete old cache to force regeneration
    if FOLDSEEK_CACHE.exists():
        FOLDSEEK_CACHE.unlink()
    max_tmscore = compute_foldseek_tmscore(train, val, test)
    test = add_tmscore_cols(test, max_tmscore)
    for thr in [0.9, 0.5, 0.3]:
        print(f"  TM-score_{thr}: {test[f'TM-score_{thr}'].sum()}")

    # 3. Extreme length
    print("\n=== 3. Extreme length ===")
    test = add_extreme_length_cols(test)
    print(f"  OOD_ExtremeShort: {test['OOD_ExtremeShort'].sum()}")
    print(f"  OOD_ExtremeLong: {test['OOD_ExtremeLong'].sum()}")

    # Update Default: False for ExtremeLong
    if "Default" in test.columns:
        test["Default"] = ~test["OOD_ExtremeLong"]
    else:
        test["Default"] = ~test["OOD_ExtremeLong"]
    print(f"  Default=True: {test['Default'].sum()}, Default=False: {(~test['Default']).sum()}")

    # 4. IDR
    print("\n=== 4. IDR (metapredict v3) ===")
    # Delete old cache
    if IDR_CACHE.exists():
        IDR_CACHE.unlink()
    test = compute_idr(test)
    print(f"  OOD_IDR: {test['OOD_IDR'].sum()}")

    # 5. EC OOD
    print("\n=== 5. EC OOD ===")
    test = compute_ec_ood(train, val, test)
    print(f"  OOD_NewEC_L4: {test['OOD_NewEC_L4'].sum()}")
    print(f"  OOD_LongTail_EC_L4: {test['OOD_LongTail_EC_L4'].sum()}")

    # 6. Default=False（预挑出极长行）收尾：除 OOD_ExtremeLong（预挑出原因列）
    #    外全部 OOD/同源/结构列置 NaN（含 OOD_ExtremeShort），表示"未计算"。
    #    2026-08-27 起与其他任务 Default=False 行惯例统一（此前 seq/TM 列对
    #    全行计算、OOD_ExtremeShort 保留，备份 output/backup_before_nondefault_full_nan/）。
    print("\n=== 6. Default=False 行置 NaN（仅保留 OOD_ExtremeLong）===")
    nd = ~test["Default"]
    nan_cols = ([c for c in test.columns
                 if c.startswith("OOD_") and c != "OOD_ExtremeLong"]
                + [f"seq_Redundancy_{t}" for t in LOW_HOMOLOGY_THRESHOLDS]
                + [f"TM-score_{t}" for t in NEWFOLD_THRESHOLDS])
    for c in nan_cols:
        if c in test.columns:
            ser = test[c].astype(object)
            ser.loc[nd] = np.nan
            test[c] = ser
    print(f"  {int(nd.sum())} 行 x {len(nan_cols)} 列置 NaN")

    # --- Reorder columns ---
    base_cols = ["unique_id", "protein_A_id", "protein_B_id",
                 "aa_seq1", "aa_seq2", "struct_file1", "struct_file2",
                 "label", "Default", "pair_type"]
    ood_cols = [
        "OOD_IDR", "OOD_NewEC_L4", "OOD_LongTail_EC_L4",
        "OOD_ExtremeShort", "OOD_ExtremeLong",
    ] + [f"seq_Redundancy_{t}" for t in LOW_HOMOLOGY_THRESHOLDS] \
      + ["OOD_Orphan"] \
      + [f"TM-score_{t}" for t in NEWFOLD_THRESHOLDS]

    final_cols = base_cols + [c for c in ood_cols if c in test.columns]
    # Add any remaining columns
    for c in test.columns:
        if c not in final_cols:
            final_cols.append(c)
    test = test[final_cols]

    test.to_csv(SPLITS_DIR / "ppi_test.csv", index=False)
    test_path = SPLITS_DIR / "ppi_test.csv"
    print(f"\n=== Saved {test_path}: {len(test)} rows, {len(test.columns)} cols ===")
    print(f"Columns: {list(test.columns)}")

    # Summary
    print("\n=== OOD Summary ===")
    summary = {"test_total": len(test)}
    for c in ood_cols:
        if c in test.columns:
            summary[c] = int(test[c].sum())
    summary["Default_True"] = int(test["Default"].sum())
    for k, v in summary.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
