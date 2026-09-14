#!/usr/bin/env python3
"""Deduplicate chain-ligand pairs by ligand at 95% sequence identity using mmseqs."""
import argparse
import gzip
import hashlib
import logging
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

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

TAG_LENGTH = 150
AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"


def ligand_tag(ligand_name, length=TAG_LENGTH):
    """Create a ligand-specific amino-acid tag of `length` residues."""
    ligand_name = str(ligand_name)
    h = hashlib.sha256(ligand_name.encode()).hexdigest()
    mapped = "".join(AA_ALPHABET[int(ch, 16) % len(AA_ALPHABET)] for ch in h)
    tag = (mapped * ((length // len(mapped)) + 1))[:length]
    return tag


def load_resolution(pdb_ids, sifts_csv=SIFTS_CSV, pdb_dir=PDB_DIR):
    """Load resolution for each PDB ID, preferring sifts then parsing files."""
    logger.info("Loading resolution information...")
    # Try sifts first
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
    """Parse resolution from PDB REMARK 2."""
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
    """Parse resolution from mmCIF _refine.ls_d_res_high."""
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


def write_fasta(df, fasta_path):
    """Write sequences to FASTA with integer headers to avoid truncation issues."""
    with open(fasta_path, "w") as f:
        for idx, row in df.iterrows():
            seq = row['_prefixed_seq']
            f.write(f">{idx}\n{seq}\n")
    logger.info(f"Wrote {len(df)} sequences to {fasta_path}")


def run_mmseqs(fasta_path, prefix, tmp_dir, threads=32):
    """Run mmseqs easy-cluster at 95% seq id / 95% coverage."""
    cmd = [
        "mmseqs", "easy-cluster",
        str(fasta_path),
        str(prefix),
        str(tmp_dir),
        "--min-seq-id", "0.95",
        "-c", "0.95",
        "--cov-mode", "1",
        "--threads", str(threads),
        "--cluster-mode", "1",
        "--split-memory-limit", "200G",
        "-v", "1",
    ]
    logger.info(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    logger.info("mmseqs clustering completed")


def parse_clusters(cluster_tsv):
    """Parse mmseqs cluster TSV: representative \t member."""
    clusters = {}
    with open(cluster_tsv) as f:
        for line in f:
            rep, member = line.strip().split("\t")
            clusters.setdefault(rep, []).append(member)
    logger.info(f"Parsed {len(clusters)} clusters")
    return clusters


def select_first_representatives(clusters, df):
    """For each cluster, keep one sequence per ligand (first occurrence)."""
    selected = set()
    for rep, members in clusters.items():
        seen_ligands = set()
        for m in members:
            idx = int(m)
            ligand = df.iloc[idx]['ligand_name']
            if ligand not in seen_ligands:
                selected.add(idx)
                seen_ligands.add(ligand)
    logger.info(f"Selected {len(selected)} representatives (first-occurrence)")
    return selected


def select_quality_representatives(clusters, df):
    """For each cluster, keep one sequence per ligand based on quality metrics.

    Priority:
      1. Most binding residues
      2. Highest resolution (lowest value)
      3. Longest protein chain
      4. Most ligand heavy atoms
    """
    # Compute quality ranks on current df (lower rank = better)
    quality_keys = df[['num_binding_residues', '_resolution', 'protein_length', 'num_ligand_heavy_atoms']].copy()
    quality_keys['num_binding_residues'] = -quality_keys['num_binding_residues']
    quality_keys['protein_length'] = -quality_keys['protein_length']
    quality_keys['num_ligand_heavy_atoms'] = -quality_keys['num_ligand_heavy_atoms']
    quality_ranks = quality_keys.apply(tuple, axis=1).rank(method='min').values

    selected = set()
    for rep, members in clusters.items():
        # Group members by ligand
        ligand_members = {}
        for m in members:
            idx = int(m)
            ligand = df.iloc[idx]['ligand_name']
            ligand_members.setdefault(ligand, []).append(idx)

        for ligand, idxs in ligand_members.items():
            best_idx = min(idxs, key=lambda i: quality_ranks[i])
            selected.add(best_idx)
    logger.info(f"Selected {len(selected)} representatives (quality-based)")
    return selected


def dedup_one_round(df, tmp_root, iteration, threads, selection_func):
    """Run one round of mmseqs clustering and select representatives."""
    iter_dir = tmp_root / f"iter_{iteration:02d}"
    iter_dir.mkdir(exist_ok=True)
    fasta_path = iter_dir / "all_seqs.fasta"
    prefix = iter_dir / "cluster"
    # Reset integer index for this round
    df = df.reset_index(drop=True)
    write_fasta(df, fasta_path)
    run_mmseqs(fasta_path, prefix, iter_dir / "tmp", threads=threads)
    cluster_tsv = f"{prefix}_cluster.tsv"
    clusters = parse_clusters(cluster_tsv)
    selected = selection_func(clusters, df)
    dedup_df = df.iloc[list(selected)].copy()
    return dedup_df


def main():
    parser = argparse.ArgumentParser(description="Deduplicate by ligand at 95% sequence identity")
    parser.add_argument("--input_csv", type=str, default=os.path.join(OUTPUT_DIR, "ligand_binding_sites_std.csv"))
    parser.add_argument("--output_dir", type=str, default=OUTPUT_DIR)
    parser.add_argument("--split_dir", type=str, default=SPLIT_DIR)
    parser.add_argument("--train_ratio", type=float, default=0.8)
    parser.add_argument("--val_ratio", type=float, default=0.1)
    parser.add_argument("--test_ratio", type=float, default=0.1)
    parser.add_argument("--random_seed", type=int, default=42)
    parser.add_argument("--threads", type=int, default=32)
    parser.add_argument("--max_iterations", type=int, default=10)
    parser.add_argument("--selection", type=str, default="first", choices=["first", "quality"],
                        help="Representative selection strategy: first or quality")
    args = parser.parse_args()

    assert abs(args.train_ratio + args.val_ratio + args.test_ratio - 1.0) < 1e-6

    output_dir = Path(args.output_dir)
    split_dir = Path(args.split_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    split_dir.mkdir(parents=True, exist_ok=True)

    # Load full dataset
    logger.info(f"Loading {args.input_csv}")
    df = pd.read_csv(args.input_csv)
    logger.info(f"Loaded {len(df)} chain-ligand pairs, {df['ligand_name'].nunique()} ligands")

    # Prepare helper columns
    df['_prefixed_seq'] = df['ligand_name'].apply(ligand_tag) + df['aa_seq']

    # Build quality lookup if needed
    if args.selection == "quality":
        res_map = load_resolution(df['pdb_id'].unique())
        df['_resolution'] = df['pdb_id'].map(res_map).fillna(np.inf)
        selection_func = select_quality_representatives
    else:
        selection_func = select_first_representatives

    # Iterative deduplication
    current = df.copy()
    tmp_root = Path(tempfile.mkdtemp(prefix="mmseqs_dedup_", dir=output_dir))
    for iteration in range(1, args.max_iterations + 1):
        logger.info(f"=== Iteration {iteration}: {len(current)} sequences ===")
        next_df = dedup_one_round(current, tmp_root, iteration, args.threads, selection_func)
        if len(next_df) == len(current):
            logger.info(f"Converged at iteration {iteration}")
            break
        logger.info(f"Removed {len(current) - len(next_df)} redundant sequences")
        current = next_df
    else:
        logger.warning(f"Did not converge within {args.max_iterations} iterations")

    drop_cols = ['_prefixed_seq']
    if '_resolution' in current.columns:
        drop_cols.append('_resolution')
    dedup_df = current.drop(columns=drop_cols)
    logger.info(f"After dedup: {len(dedup_df)} pairs from {dedup_df['pdb_id'].nunique()} PDBs, {dedup_df['ligand_name'].nunique()} ligands")

    # Report reduction per ligand
    logger.info("Top ligands before/after dedup:")
    before_counts = df['ligand_name'].value_counts().head(10)
    after_counts = dedup_df['ligand_name'].value_counts()
    for ligand in before_counts.index:
        logger.info(f"  {ligand}: {before_counts[ligand]} -> {after_counts.get(ligand, 0)}")

    # Save deduplicated std CSV
    std_cols = [
        "unique_id", "file_type", "aa_seq", "labels", "pdb_id", "chain_id",
        "ligand_id", "ligand_name", "ligand_chain", "ligand_resnum",
        "num_binding_residues", "num_ligand_heavy_atoms", "protein_length"
    ]
    dedup_df[std_cols].to_csv(output_dir / "ligand_binding_sites_std_dedup95.csv", index=False)

    # Re-split by pdb_id
    unique_pdbs = dedup_df['pdb_id'].unique()
    rng = np.random.default_rng(args.random_seed)
    shuffled = rng.permutation(unique_pdbs)
    n_train = int(len(shuffled) * args.train_ratio)
    n_val = int(len(shuffled) * args.val_ratio)
    train_pdbs = shuffled[:n_train]
    val_pdbs = shuffled[n_train:n_train + n_val]
    test_pdbs = shuffled[n_train + n_val:]

    train_df = dedup_df[dedup_df['pdb_id'].isin(train_pdbs)][std_cols]
    val_df = dedup_df[dedup_df['pdb_id'].isin(val_pdbs)][std_cols]
    test_df = dedup_df[dedup_df['pdb_id'].isin(test_pdbs)][std_cols]

    train_df.to_csv(split_dir / "ligand_binding_site_train_dedup95.csv", index=False)
    val_df.to_csv(split_dir / "ligand_binding_site_val_dedup95.csv", index=False)
    test_df.to_csv(split_dir / "ligand_binding_site_test_dedup95.csv", index=False)

    logger.info(f"Splits: train={len(train_df)}, val={len(val_df)}, test={len(test_df)}")

    # Cleanup tmp
    shutil.rmtree(tmp_root)
    logger.info(f"Cleaned up {tmp_root}")

    logger.info("Done")


if __name__ == "__main__":
    main()
