#!/usr/bin/env python3
"""Exact post-cleanup of mmseqs 95% ligand-wise dedup using edlib pairwise alignment."""
import argparse
import gzip
import logging
import os
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import edlib
import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s')
logger = logging.getLogger(__name__)

TASK_DIR = Path(__file__).resolve().parent.parent  # datasets/ligand_binding_site
ROOT = TASK_DIR.parent.parent                      # 项目根（POOR/）
OUTPUT_DIR = str(TASK_DIR / "output")
SPLIT_DIR = str(TASK_DIR / "splits")
SIFTS_CSV = str(ROOT / "data" / "sifts" / "annotated_chains.csv")
PDB_DIR = str(ROOT / "data" / "pdb")

EDLIB_IDENTITY_THRESHOLD = 0.95
KMER_K = 3
KMER_JACCARD_CUTOFF = 0.75
MMSSEQS_LARGE_LIGAND_THRESHOLD = 50


def load_resolution(pdb_ids, sifts_csv=SIFTS_CSV, pdb_dir=PDB_DIR):
    """Load resolution for each PDB ID, preferring sifts then parsing files."""
    logger.info("Loading resolution information...")
    try:
        sifts = pd.read_csv(sifts_csv, usecols=['pdb_id', 'resolution'])
        sifts = sifts.dropna(subset=['resolution'])
        sifts_res = sifts.groupby('pdb_id')['resolution'].first().to_dict()
    except Exception as e:
        logger.warning(f"Could not load sifts resolution: {e}")
        sifts_res = {}

    res_map = {}
    missing = []
    for pid in pdb_ids:
        pid_upper = pid.upper()
        if pid_upper in sifts_res:
            res_map[pid] = sifts_res[pid_upper]
        else:
            missing.append(pid)

    if missing:
        logger.info(f"Parsing resolution from files for {len(missing)} PDBs...")
        for i, pid in enumerate(missing):
            r = parse_pdb_resolution(Path(pdb_dir) / f"{pid.upper()}.pdb.gz")
            if r is None:
                r = parse_cif_resolution(Path(pdb_dir) / f"{pid.upper()}.cif.gz")
            if r is not None:
                res_map[pid] = r
            if (i + 1) % 10000 == 0:
                logger.info(f"Parsed {i+1}/{len(missing)} PDBs")

    unresolved = set(pdb_ids) - set(res_map.keys())
    if unresolved:
        logger.warning(f"Could not resolve {len(unresolved)} PDBs; treating as worst resolution")
    logger.info(f"Resolution loaded for {len(res_map)}/{len(pdb_ids)} PDBs")
    return res_map


def parse_pdb_resolution(pdb_file):
    if not pdb_file.exists():
        return None
    try:
        with gzip.open(pdb_file, 'rt', errors='ignore') as f:
            for line in f:
                if line.startswith('REMARK   2 RESOLUTION'):
                    parts = line.split()
                    for p in parts:
                        try:
                            return float(p)
                        except ValueError:
                            continue
    except Exception:
        pass
    return None


def parse_cif_resolution(cif_file):
    if not cif_file.exists():
        return None
    try:
        with gzip.open(cif_file, 'rt', errors='ignore') as f:
            for line in f:
                if '_refine.ls_d_res_high' in line:
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        try:
                            return float(parts[-1])
                        except ValueError:
                            continue
    except Exception:
        pass
    return None


def parse_cigar(cigar):
    """Parse edlib CIGAR into counts per operation."""
    counts = defaultdict(int)
    for n, op in re.findall(r'(\d+)([=XID])', cigar):
        counts[op] += int(n)
    return counts


def kmer_set(seq, k=KMER_K):
    """Return set of k-mers."""
    return set(seq[i:i+k] for i in range(len(seq) - k + 1))


def jaccard(a, b):
    if not a or not b:
        return 0.0
    inter = len(a & b)
    return inter / len(a | b)


def edlib_identity(seq1, seq2):
    """Compute global alignment identity (matches / alignment length)."""
    result = edlib.align(seq1, seq2, task="path")
    counts = parse_cigar(result['cigar'])
    aln_len = sum(counts.values())
    matches = counts['=']
    return matches / aln_len if aln_len > 0 else 0.0


def pick_best(component, quality_df):
    """Pick the best representative from a component by quality metrics.

    Priority:
      1. Most binding residues
      2. Highest resolution (lowest value)
      3. Longest protein chain
      4. Most ligand heavy atoms
    """
    rows = quality_df.loc[component]
    rows_sorted = rows.sort_values(
        by=['num_binding_residues', 'resolution', 'protein_length', 'num_ligand_heavy_atoms'],
        ascending=[False, True, False, False]
    )
    return rows_sorted.index[0]


def exact_pairwise_dedup(items, quality_df=None, threshold=EDLIB_IDENTITY_THRESHOLD):
    """Remove sequences within a ligand that are >= threshold identical.

    items: list of (unique_id, aa_seq)
    quality_df: DataFrame indexed by unique_id with quality columns
    Returns: set of unique_ids to keep.
    """
    n = len(items)
    if n <= 1:
        return {uid for uid, _ in items}

    ids = [uid for uid, _ in items]
    seqs = [seq for _, seq in items]
    kmers = [kmer_set(seq) for seq in seqs]

    # Build adjacency list for pairs above threshold
    adj = defaultdict(set)
    for i in range(n):
        len_i = len(seqs[i])
        for j in range(i + 1, n):
            len_j = len(seqs[j])
            # Length prefilter
            if min(len_i, len_j) / max(len_i, len_j) < threshold - 0.03:
                continue
            # K-mer Jaccard prefilter
            if jaccard(kmers[i], kmers[j]) < KMER_JACCARD_CUTOFF:
                continue
            identity = edlib_identity(seqs[i], seqs[j])
            if identity >= threshold:
                adj[ids[i]].add(ids[j])
                adj[ids[j]].add(ids[i])

    # Find connected components and keep one per component
    visited = set()
    keep = set()
    for uid in ids:
        if uid in visited:
            continue
        component = [uid]
        visited.add(uid)
        stack = [uid]
        while stack:
            cur = stack.pop()
            for nb in adj[cur]:
                if nb not in visited:
                    visited.add(nb)
                    stack.append(nb)
                    component.append(nb)
        if quality_df is not None:
            keep.add(pick_best(component, quality_df))
        else:
            keep.add(component[0])
    return keep


def mmseqs_per_ligand_dedup(items, quality_df=None, threads=8):
    """Run mmseqs easy-cluster on a single ligand's sequences."""
    tmp = Path(tempfile.mkdtemp(prefix="mmseqs_lig_"))
    fasta = tmp / "seqs.fasta"
    with open(fasta, "w") as f:
        for idx, (uid, seq) in enumerate(items):
            f.write(f">{idx}\n{seq}\n")
    prefix = tmp / "clust"
    subprocess.run([
        "mmseqs", "easy-cluster", str(fasta), str(prefix), str(tmp / "tmp"),
        "--min-seq-id", "0.95", "-c", "0.95", "--cov-mode", "1",
        "--threads", str(threads), "-v", "0"
    ], check=True)

    # Read clusters and map integer indices back to unique_ids
    clusters = defaultdict(list)
    with open(f"{prefix}_cluster.tsv") as f:
        for line in f:
            rep, mem = line.strip().split("\t")
            clusters[int(rep)].append(int(mem))

    keep = set()
    uids = [uid for uid, _ in items]
    for rep, idxs in clusters.items():
        if quality_df is not None:
            member_uids = [uids[i] for i in idxs]
            keep.add(pick_best(member_uids, quality_df))
        else:
            keep.add(uids[rep])
    shutil.rmtree(tmp)
    return keep


def main():
    parser = argparse.ArgumentParser(description="Exact cleanup of 95% ligand-wise dedup")
    parser.add_argument("--input_csv", type=str, default=os.path.join(OUTPUT_DIR, "ligand_binding_sites_std_dedup95.csv"))
    parser.add_argument("--output_csv", type=str, default=os.path.join(OUTPUT_DIR, "ligand_binding_sites_std_dedup95.csv"))
    parser.add_argument("--split_dir", type=str, default=SPLIT_DIR)
    parser.add_argument("--threads", type=int, default=32)
    parser.add_argument("--selection", type=str, default="first", choices=["first", "quality"],
                        help="Representative selection strategy: first or quality")
    args = parser.parse_args()

    logger.info(f"Loading {args.input_csv}")
    df = pd.read_csv(args.input_csv)
    logger.info(f"Loaded {len(df)} pairs, {df['ligand_name'].nunique()} ligands")

    # Build quality lookup if needed
    if args.selection == "quality":
        res_map = load_resolution(df['pdb_id'].unique())
        df['_resolution'] = df['pdb_id'].map(res_map).fillna(np.inf)
        quality_df = df.set_index('unique_id')[
            ['num_binding_residues', '_resolution', 'protein_length', 'num_ligand_heavy_atoms']
        ].rename(columns={'_resolution': 'resolution'})
    else:
        quality_df = None

    keep_uids = set()
    ligand_counts = df['ligand_name'].value_counts()

    for idx, (ligand, count) in enumerate(ligand_counts.items(), 1):
        sub = df[df['ligand_name'] == ligand]
        items = list(zip(sub['unique_id'], sub['aa_seq']))
        if count <= MMSSEQS_LARGE_LIGAND_THRESHOLD:
            kept = exact_pairwise_dedup(items, quality_df=quality_df)
        else:
            kept = mmseqs_per_ligand_dedup(items, quality_df=quality_df, threads=max(1, args.threads // 4))
        keep_uids.update(kept)
        if idx % 1000 == 0:
            logger.info(f"Processed {idx}/{len(ligand_counts)} ligands, kept {len(keep_uids)} so far")

    cleaned_df = df[df['unique_id'].isin(keep_uids)].copy()
    logger.info(f"After exact cleanup: {len(cleaned_df)} pairs from {cleaned_df['pdb_id'].nunique()} PDBs, {cleaned_df['ligand_name'].nunique()} ligands")

    std_cols = [
        "unique_id", "file_type", "aa_seq", "labels", "pdb_id", "chain_id",
        "ligand_id", "ligand_name", "ligand_chain", "ligand_resnum",
        "num_binding_residues", "num_ligand_heavy_atoms", "protein_length"
    ]
    cleaned_df[std_cols].to_csv(args.output_csv, index=False)

    # Re-split by pdb_id
    split_dir = Path(args.split_dir)
    split_dir.mkdir(parents=True, exist_ok=True)
    unique_pdbs = cleaned_df['pdb_id'].unique()
    rng = np.random.default_rng(42)
    shuffled = rng.permutation(unique_pdbs)
    n_train = int(len(shuffled) * 0.8)
    n_val = int(len(shuffled) * 0.1)
    train_pdbs = shuffled[:n_train]
    val_pdbs = shuffled[n_train:n_train + n_val]
    test_pdbs = shuffled[n_train + n_val:]

    cleaned_df[cleaned_df['pdb_id'].isin(train_pdbs)][std_cols].to_csv(split_dir / "ligand_binding_site_train_dedup95.csv", index=False)
    cleaned_df[cleaned_df['pdb_id'].isin(val_pdbs)][std_cols].to_csv(split_dir / "ligand_binding_site_val_dedup95.csv", index=False)
    cleaned_df[cleaned_df['pdb_id'].isin(test_pdbs)][std_cols].to_csv(split_dir / "ligand_binding_site_test_dedup95.csv", index=False)

    # Final mmseqs2 convergence check to catch any local-alignment redundancy
    # that edlib global alignment might have missed
    logger.info("Running final mmseqs2 convergence check...")
    script_dir = Path(__file__).parent
    subprocess.run([
        "python", str(script_dir / "03_dedup_by_ligand_95.py"),
        "--input_csv", str(args.output_csv),
        "--output_dir", str(Path(args.output_csv).parent),
        "--split_dir", str(args.split_dir),
        "--threads", str(args.threads),
        "--max_iterations", "10",
        "--selection", args.selection
    ], check=True)

    logger.info("Done")


if __name__ == "__main__":
    main()
