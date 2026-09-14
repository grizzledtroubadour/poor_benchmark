#!/usr/bin/env python3
"""Build final train/val/test splits with OOD labels for ligand-binding-site.

Inputs:
  - time-based train/val/test CSVs
  - extracted chain PDBs in ligand_binding_site/pdbs/

Outputs:
  - splits/ligand_binding_site_{train,val,test}.csv
  - intermediate caches in output/ood/

NOTE (2026-08 统一整理): 列名已适配统一 schema（labels->label, pdb_file->struct_file,
file_type 列取消）；缺失的 pdb_id/chain_id/protein_length 会从 unique_id/aa_seq 派生。
序列同源 OOD 的正式重算口径以 `.skills/homology-ood-annotation/scripts/seq_homology_ood.py`
为准（本脚本内嵌 mmseqs 段仅作历史流程留存，视为 superseded）；EC/holdout 等后续列
由对应 skill 脚本落盘（见 DATA_PROCESS.md 阶段四）。
"""
import argparse
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import ast
import gzip
import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s')
logger = logging.getLogger(__name__)

TASK_DIR = Path(__file__).resolve().parent.parent  # datasets/ligand_binding_site
ROOT = TASK_DIR.parent.parent                      # 项目根（POOR/）
DATASET_DIR = TASK_DIR / "splits"
PDB_DIR = TASK_DIR / "pdbs"
OOD_DIR = TASK_DIR / "output" / "ood"
SIFTS_EC_FILE = ROOT / "data" / "sifts" / "sifts_chain_ec.tsv.gz"

SEQ_THRESHOLDS = [90, 80, 70, 60, 50, 40, 30]
TM_THRESHOLDS = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3]


def chain_key(row):
    return f"{row['pdb_id'].upper()}_{row['chain_id']}"


def load_ec_annotations(ec_file):
    """Load SIFTS chain-level EC annotations into a chain_key -> set(EC_NUMBER) map."""
    if not ec_file.exists():
        logger.warning(f"EC annotation file not found: {ec_file}")
        return {}
    logger.info(f"Loading EC annotations from {ec_file}")
    df = pd.read_csv(ec_file, sep="\t", comment="#", keep_default_na=False, na_values=[""])
    df["PDB"] = df["PDB"].str.upper()
    df["chain_key"] = df["PDB"] + "_" + df["CHAIN"]
    # Some EC numbers contain '-' placeholders (e.g. 1.1.1.-); keep them as-is.
    ec_map = (
        df.groupby("chain_key")["EC_NUMBER"]
        .apply(lambda x: set(v for v in x if v and v.strip()))
        .to_dict()
    )
    logger.info(f"Loaded EC annotations for {len(ec_map)} unique chains")
    return ec_map


def binding_site_ratio(labels_str, seq_len):
    """Parse labels string and return fraction of binding residues."""
    try:
        labels = ast.literal_eval(str(labels_str))
    except Exception:
        return 0.0
    if not isinstance(labels, (list, tuple)) or len(labels) == 0:
        return 0.0
    if seq_len <= 0:
        return 0.0
    return float(sum(1 for v in labels if v)) / seq_len


def write_unique_fasta(df, out_fasta, key_col="chain_key", seq_col="aa_seq"):
    """Write one FASTA entry per unique key, using the first sequence seen."""
    records = {}
    for k, seq in zip(df[key_col], df[seq_col]):
        k = str(k)
        if pd.notna(seq) and k not in records:
            records[k] = str(seq)
    with open(out_fasta, "w") as f:
        for k, seq in records.items():
            f.write(f">{k}\n{seq}\n")
    return records


def run_mmseqs_search(query_fasta, target_fasta, output_m8, tmp_dir):
    # 2026-08-27 起统一为 easy-search -s 7 口径（默认 e-value <= 1e-3 显著性门槛），
    # 与 func_prediction 等任务及 homology-ood-annotation skill 一致；
    # 旧式 `search --min-seq-id 0.0 -c 0.0 --cov-mode 0 -s 7.5 -e 100` 会把 9-25 aa 弱
    # 局部比对计入 max pident 并破坏 Orphan⊆seq30 嵌套（修复记录见归档的
    # output/scripts_archive/24_fix_seq_homology_easysearch.py）。
    logger.info("Running mmseqs2 easy-search -s 7 ...")
    subprocess.run(
        [
            "mmseqs", "easy-search", str(query_fasta), str(target_fasta), str(output_m8), str(tmp_dir),
            "-s", "7",
            "--format-output", "query,target,pident,evalue",
        ],
        check=True,
        capture_output=True,
    )


def parse_mmseqs_m8(m8_file, query_keys):
    max_identity = {k: 0.0 for k in query_keys}
    best_evalue = {k: np.inf for k in query_keys}
    if not m8_file.exists():
        return max_identity, best_evalue
    with open(m8_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            q, t, pident, evalue = parts[0], parts[1], float(parts[2]), float(parts[3])
            if q == t:
                continue
            if pident > max_identity.get(q, 0.0):
                max_identity[q] = pident
            if evalue < best_evalue.get(q, np.inf):
                best_evalue[q] = evalue
    return max_identity, best_evalue


def symlink_unique_structures(df, link_dir, pdb_dir):
    """Create one symlink per chain_key pointing to the first available PDB file."""
    link_dir.mkdir(parents=True, exist_ok=True)
    seen = set()
    links = 0
    missing = 0
    for _, row in df.iterrows():
        key = row["chain_key"]
        if key in seen:
            continue
        seen.add(key)
        src = pdb_dir / f"{row['unique_id']}.pdb"
        if not src.exists():
            missing += 1
            continue
        dst = link_dir / f"{key}.pdb"
        if dst.exists() or dst.is_symlink():
            os.remove(dst)
        os.symlink(src.resolve(), dst)
        links += 1
    logger.info(f"Created {links} structure symlinks in {link_dir} (missing {missing})")
    return links


def run_foldseek_search(query_dir, target_dir, output_m8, tmp_dir):
    # 2026-08-28 起统一为 easy-search 口径（默认 e-value <= 10、sensitivity 9.5），
    # 与 ppi 等任务及 homology-ood-annotation skill 一致；旧式
    # `search -a -e inf --max-seqs 1000` + convertalis 的缓存含 21% 非法 >1 值
    # 且整体系统性偏高，曾严重低估 TM-score_* 各档计数（对比见
    # output/foldseek_easysearch_compare/compare_report.txt）。
    logger.info("Running foldseek easy-search ...")
    subprocess.run(
        ["foldseek", "easy-search", str(query_dir), str(target_dir), str(output_m8), str(tmp_dir),
         "--format-output", "query,target,alntmscore"],
        check=True,
    )


def parse_foldseek_m8(m8_file, query_keys):
    max_tm = {k: 0.0 for k in query_keys}
    if not m8_file.exists():
        return max_tm
    with open(m8_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            q, t, score = parts[0], parts[1], float(parts[2])
            # easy-search 直接输出 chain_key；旧式 convertalis 可能带 .pdb 后缀
            q = Path(q).stem
            t = Path(t).stem
            if q == t:
                continue
            if score > max_tm.get(q, 0.0):
                max_tm[q] = score
    return max_tm


def clean_sequence(seq):
    seq = str(seq).upper()
    seq = seq.replace("B", "N").replace("Z", "Q").replace("J", "L")
    seq = seq.replace("U", "C").replace("O", "K")
    seq = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "A", seq)
    return seq


def continuous_regions(mask):
    lengths = []
    cur = 0
    for flag in mask:
        if flag:
            cur += 1
        else:
            if cur:
                lengths.append(cur)
                cur = 0
    if cur:
        lengths.append(cur)
    return lengths


def compute_idr(scores, seq_len):
    arr = np.asarray(scores)
    binary = arr > 0.5
    total_ratio = float(binary.sum()) / seq_len
    regions = continuous_regions(binary)
    max_cont_ratio = max(regions) / seq_len if regions else 0.0
    is_idr = max_cont_ratio > 0.3
    return total_ratio, is_idr


def compute_idr_for_df(df, cache_file, batch_size=1000):
    try:
        import metapredict
    except ImportError as exc:
        logger.error(f"metapredict not installed: {exc}")
        raise

    cache = {}
    if cache_file.exists():
        cache_df = pd.read_csv(cache_file)
        cache = dict(zip(cache_df["sequence"], cache_df["idr_ratio"]))
        logger.info(f"Loaded {len(cache)} cached IDR predictions")

    unique_seqs = df["aa_seq"].dropna().unique()
    clean_map = {s: clean_sequence(s) for s in unique_seqs}
    clean_unique = sorted(set(clean_map.values()))
    remaining = [s for s in clean_unique if s not in cache]
    logger.info(f"IDR: {len(clean_unique)} unique cleaned sequences, {len(remaining)} to predict")

    fasta_dir = OOD_DIR / "idr_fastas"
    fasta_dir.mkdir(parents=True, exist_ok=True)

    for i in range(0, len(remaining), batch_size):
        batch = remaining[i : i + batch_size]
        fasta_path = fasta_dir / f"batch_{i:06d}.fasta"
        with open(fasta_path, "w") as fh:
            for j, s in enumerate(batch):
                fh.write(f">seq_{i+j}\n{s}\n")
        logger.info(f"IDR batch {i//batch_size + 1}/{(len(remaining)+batch_size-1)//batch_size}: {len(batch)} seqs")
        result = metapredict.predict_disorder_fasta(
            str(fasta_path),
            show_progress_bar=False,
            invalid_sequence_action="convert",
        )
        for _, (seq, scores) in result.items():
            total_ratio, _ = compute_idr(scores, len(seq))
            cache[seq] = total_ratio
        pd.DataFrame({"sequence": list(cache.keys()), "idr_ratio": list(cache.values())}).to_csv(cache_file, index=False)

    # Map back to original sequences.
    seq_to_total = {orig: cache.get(clean_map[orig], 0.0) for orig in unique_seqs}
    df["idr_ratio"] = df["aa_seq"].map(seq_to_total)
    df["OOD_IDR"] = df["idr_ratio"] > 0.3
    return df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--train", type=str, default=str(DATASET_DIR / "ligand_binding_site_train.csv"))
    parser.add_argument("--val", type=str, default=str(DATASET_DIR / "ligand_binding_site_val.csv"))
    parser.add_argument("--test", type=str, default=str(DATASET_DIR / "ligand_binding_site_test.csv"))
    parser.add_argument("--pdb_dir", type=str, default=str(PDB_DIR))
    parser.add_argument("--out_dir", type=str, default=str(DATASET_DIR))
    parser.add_argument("--skip_structures", action="store_true", help="Skip Foldseek structural OOD (use if pdbs not ready)")
    args = parser.parse_args()

    OOD_DIR.mkdir(parents=True, exist_ok=True)
    pdb_dir = Path(args.pdb_dir)
    out_dir = Path(args.out_dir)

    train = pd.read_csv(args.train, low_memory=False, dtype={"label": str}, keep_default_na=False, na_values=[''])
    val = pd.read_csv(args.val, low_memory=False, dtype={"label": str}, keep_default_na=False, na_values=[''])
    test = pd.read_csv(args.test, low_memory=False, dtype={"label": str}, keep_default_na=False, na_values=[''])

    for df in [train, val, test]:
        # Derive metadata columns when absent (current splits keep only core columns).
        if "pdb_id" not in df.columns:
            df["pdb_id"] = df["unique_id"].str.split("_").str[0]
        if "chain_id" not in df.columns:
            df["chain_id"] = df["unique_id"].str.split("_").str[1]
        if "protein_length" not in df.columns:
            df["protein_length"] = df["aa_seq"].str.len()
        # Preserve 'NA' chain IDs which pandas otherwise interprets as NaN.
        df["chain_id"] = df["chain_id"].fillna("NA").astype(str)
        df["chain_key"] = df.apply(chain_key, axis=1)
        df["struct_file"] = df["unique_id"].apply(lambda uid: f"{uid}.pdb")

    trainval = pd.concat([train, val], ignore_index=True)
    test_keys = set(test["chain_key"].unique())
    trainval_keys = set(trainval["chain_key"].unique())
    logger.info(f"Train rows: {len(train)}, Val rows: {len(val)}, Test rows: {len(test)}")
    logger.info(f"Test unique chains: {len(test_keys)}, Train+Val unique chains: {len(trainval_keys)}")

    # ------------------------------------------------------------------
    # Sequence homology (mmseqs2)
    # [superseded] 正式重算请用 .skills/homology-ood-annotation/scripts/seq_homology_ood.py；
    # 本段仅保留为历史流程记录（已对齐 easy-search -s 7 口径）。
    # ------------------------------------------------------------------
    seq_cache = OOD_DIR / "mmseqs_max_identity.csv"
    if seq_cache.exists():
        logger.info(f"Loading cached sequence homology from {seq_cache}")
        seq_df = pd.read_csv(seq_cache)
        max_identity = dict(zip(seq_df["chain_key"], seq_df["max_identity"]))
        best_evalue = dict(zip(seq_df["chain_key"], seq_df["best_evalue"]))
    else:
        with tempfile.TemporaryDirectory(prefix="lbs_ood_seq_", dir=OOD_DIR) as tmp:
            tmp = Path(tmp)
            query_fa = tmp / "query.fasta"
            target_fa = tmp / "target.fasta"
            out_m8 = tmp / "result.m8"
            write_unique_fasta(test, query_fa)
            write_unique_fasta(trainval, target_fa)
            run_mmseqs_search(query_fa, target_fa, out_m8, tmp)
            max_identity, best_evalue = parse_mmseqs_m8(out_m8, test_keys)
        seq_df = pd.DataFrame({
            "chain_key": list(max_identity.keys()),
            "max_identity": [max_identity[k] for k in max_identity.keys()],
            "best_evalue": [best_evalue[k] for k in max_identity.keys()],
        })
        seq_df.to_csv(seq_cache, index=False)
        logger.info(f"Saved sequence homology cache: {seq_cache}")

    test["max_seq_identity"] = test["chain_key"].map(max_identity).fillna(0.0)
    test["best_evalue"] = test["chain_key"].map(best_evalue).fillna(np.inf)
    for thr in SEQ_THRESHOLDS:
        test[f"seq_Redundancy_{thr}"] = test["max_seq_identity"] < thr
    test["OOD_Orphan"] = test["best_evalue"] > 1e-3

    # ------------------------------------------------------------------
    # Structural homology (Foldseek)
    # ------------------------------------------------------------------
    if args.skip_structures:
        logger.warning("Skipping Foldseek structural OOD; all TM-score / NewFold columns set to True")
        for thr in TM_THRESHOLDS:
            test[f"TM-score_{thr:.1f}"] = True
            test[f"OOD_NewFold_{thr:.1f}"] = True
    else:
        tm_cache = OOD_DIR / "foldseek_max_tmscore.csv"
        if tm_cache.exists():
            logger.info(f"Loading cached structural homology from {tm_cache}")
            tm_df = pd.read_csv(tm_cache)
            max_tm = dict(zip(tm_df["chain_key"], tm_df["max_tmscore"]))
        else:
            with tempfile.TemporaryDirectory(prefix="lbs_ood_struct_", dir=OOD_DIR) as tmp:
                tmp = Path(tmp)
                query_dir = tmp / "query_structures"
                target_dir = tmp / "target_structures"
                symlink_unique_structures(test, query_dir, pdb_dir)
                symlink_unique_structures(trainval, target_dir, pdb_dir)
                out_m8 = tmp / "foldseek_result.m8"
                run_foldseek_search(query_dir, target_dir, out_m8, tmp)
                max_tm = parse_foldseek_m8(out_m8, test_keys)
            tm_df = pd.DataFrame({
                "chain_key": list(max_tm.keys()),
                "max_tmscore": [max_tm[k] for k in max_tm.keys()],
            })
            tm_df.to_csv(tm_cache, index=False)
            logger.info(f"Saved structural homology cache: {tm_cache}")

        test["max_tmscore"] = test["chain_key"].map(max_tm).fillna(0.0)
        for thr in TM_THRESHOLDS:
            test[f"OOD_NewFold_{thr:.1f}"] = test["max_tmscore"] < thr

    # ------------------------------------------------------------------
    # IDR (metapredict)
    # ------------------------------------------------------------------
    idr_cache = OOD_DIR / "idr_predictions.csv"
    # Compute IDR on unique test chain sequences.
    test_unique_seq = test.drop_duplicates(subset=["chain_key"]).copy()
    test_unique_seq = compute_idr_for_df(test_unique_seq, idr_cache, batch_size=1000)
    idr_map = dict(zip(test_unique_seq["chain_key"], test_unique_seq["idr_ratio"]))
    idr_flag_map = dict(zip(test_unique_seq["chain_key"], test_unique_seq["OOD_IDR"]))
    test["idr_ratio"] = test["chain_key"].map(idr_map).fillna(0.0)
    test["OOD_IDR"] = test["chain_key"].map(idr_flag_map).fillna(False)

    # ------------------------------------------------------------------
    # EC-based New Function OOD
    # ------------------------------------------------------------------
    ec_map = load_ec_annotations(SIFTS_EC_FILE)
    trainval["ec_set"] = trainval["chain_key"].map(lambda k: ec_map.get(k, set()))
    test["ec_set"] = test["chain_key"].map(lambda k: ec_map.get(k, set()))
    test["has_ec_annotation"] = test["ec_set"].apply(lambda s: len(s) > 0)
    trainval_ec = set()
    for ecs in trainval["ec_set"]:
        trainval_ec.update(ecs)
    logger.info(f"Train+Val unique EC numbers: {len(trainval_ec)}")
    logger.info(f"Test proteins with EC annotation: {test['has_ec_annotation'].sum()} / {len(test)}")
    test["OOD_NewFunction"] = test["ec_set"].apply(
        lambda s: len(s) > 0 and all(ec not in trainval_ec for ec in s)
    )

    # ------------------------------------------------------------------
    # Long-tail binding site OOD
    # ------------------------------------------------------------------
    test["binding_site_ratio"] = test.apply(
        lambda r: binding_site_ratio(r["label"], len(r["aa_seq"])), axis=1
    )
    test["OOD_LongTail"] = test["binding_site_ratio"] < 0.01

    # ------------------------------------------------------------------
    # Extreme length
    # ------------------------------------------------------------------
    test["OOD_ExtremeLength"] = (test["protein_length"] < 60) | (test["protein_length"] > 1000)
    test["is_time_cutoff"] = True

    # ------------------------------------------------------------------
    # Compose final CSVs
    # ------------------------------------------------------------------
    core_meta_cols = [
        "unique_id", "aa_seq", "label",
        "pdb_id", "chain_id", "ligand_id", "ligand_name", "ligand_chain", "ligand_resnum",
        "num_binding_residues", "num_ligand_heavy_atoms", "protein_length",
        "struct_file",
    ]
    train_out = train[[c for c in core_meta_cols if c in train.columns]]
    val_out = val[[c for c in core_meta_cols if c in val.columns]]

    ood_cols = ["is_time_cutoff"]
    ood_cols += [f"seq_Redundancy_{thr}" for thr in SEQ_THRESHOLDS]
    ood_cols += [f"OOD_NewFold_{thr:.1f}" for thr in TM_THRESHOLDS]
    ood_cols += ["OOD_Orphan", "OOD_ExtremeLength", "idr_ratio", "OOD_IDR"]
    ood_cols += ["has_ec_annotation", "OOD_NewFunction", "binding_site_ratio", "OOD_LongTail"]

    test_out = test[[c for c in core_meta_cols if c in test.columns] + ood_cols]

    train_out.to_csv(out_dir / "ligand_binding_site_train.csv", index=False)
    val_out.to_csv(out_dir / "ligand_binding_site_val.csv", index=False)
    test_out.to_csv(out_dir / "ligand_binding_site_test.csv", index=False)

    logger.info(f"Saved final splits to {out_dir}")

    # Summary
    logger.info("OOD summary:")
    for thr in SEQ_THRESHOLDS:
        col = f"seq_Redundancy_{thr}"
        logger.info(f"  {col}: {test[col].sum()} / {len(test)} ({100*test[col].mean():.2f}%)")
    if not args.skip_structures:
        for thr in TM_THRESHOLDS:
            col = f"OOD_NewFold_{thr:.1f}"
            logger.info(f"  {col}: {test[col].sum()} / {len(test)} ({100*test[col].mean():.2f}%)")
    logger.info(f"  OOD_Orphan: {test['OOD_Orphan'].sum()} / {len(test)} ({100*test['OOD_Orphan'].mean():.2f}%)")
    logger.info(f"  OOD_ExtremeLength: {test['OOD_ExtremeLength'].sum()} / {len(test)} ({100*test['OOD_ExtremeLength'].mean():.2f}%)")
    logger.info(f"  OOD_IDR: {test['OOD_IDR'].sum()} / {len(test)} ({100*test['OOD_IDR'].mean():.2f}%)")
    logger.info(f"  OOD_NewFunction: {test['OOD_NewFunction'].sum()} / {len(test)} ({100*test['OOD_NewFunction'].mean():.2f}%)")
    logger.info(f"  OOD_LongTail: {test['OOD_LongTail'].sum()} / {len(test)} ({100*test['OOD_LongTail'].mean():.2f}%)")


if __name__ == "__main__":
    main()
