#!/usr/bin/env python3
"""Create time-based Train/Val/Test splits (7:1:2) for the dedup95 dataset."""
import argparse
import gzip
import logging
import os
import re
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s')
logger = logging.getLogger(__name__)

TASK_DIR = Path(__file__).resolve().parent.parent  # datasets/ligand_binding_site
ROOT = TASK_DIR.parent.parent                      # 项目根（POOR/）
OUTPUT_DIR = str(TASK_DIR / "output")
SPLIT_DIR = str(TASK_DIR / "splits")
PDB_DIR = str(ROOT / "data" / "pdb")

MONTH_MAP = {
    'JAN': 1, 'FEB': 2, 'MAR': 3, 'APR': 4, 'MAY': 5, 'JUN': 6,
    'JUL': 7, 'AUG': 8, 'SEP': 9, 'OCT': 10, 'NOV': 11, 'DEC': 12
}


def parse_pdb_date(pdb_file):
    """Parse deposition date from PDB format HEADER line."""
    try:
        with gzip.open(pdb_file, 'rt', errors='ignore') as f:
            for line in f:
                if line.startswith('HEADER'):
                    date_str = line[50:59].strip()
                    if date_str:
                        day, mon, yy = date_str.split('-')
                        year = 1900 + int(yy) if int(yy) >= 50 else 2000 + int(yy)
                        return f"{year:04d}-{MONTH_MAP[mon]:02d}-{int(day):02d}"
    except Exception as e:
        logger.debug(f"Failed to parse PDB date from {pdb_file}: {e}")
    return None


def parse_cif_date(cif_file):
    """Parse initial deposition date from mmCIF _pdbx_database_status."""
    try:
        with gzip.open(cif_file, 'rt', errors='ignore') as f:
            for line in f:
                if '_pdbx_database_status.recvd_initial_deposition_date' in line:
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        return parts[-1]
    except Exception as e:
        logger.debug(f"Failed to parse mmCIF date from {cif_file}: {e}")
    return None


def get_pdb_date(pdb_id, pdb_dir):
    """Get deposition date for a PDB ID from local files."""
    pdb_id = pdb_id.upper()
    pdb_file = Path(pdb_dir) / f"{pdb_id}.pdb.gz"
    cif_file = Path(pdb_dir) / f"{pdb_id}.cif.gz"

    if pdb_file.exists():
        d = parse_pdb_date(pdb_file)
        if d:
            return d
    if cif_file.exists():
        d = parse_cif_date(cif_file)
        if d:
            return d
    return None


def main():
    parser = argparse.ArgumentParser(description="Create time-based 7:1:2 splits")
    parser.add_argument("--input_csv", type=str, default=os.path.join(OUTPUT_DIR, "ligand_binding_sites_std_dedup95.csv"))
    parser.add_argument("--pdb_dir", type=str, default=PDB_DIR)
    parser.add_argument("--split_dir", type=str, default=SPLIT_DIR)
    parser.add_argument("--train_ratio", type=float, default=0.7)
    parser.add_argument("--val_ratio", type=float, default=0.1)
    parser.add_argument("--test_ratio", type=float, default=0.2)
    args = parser.parse_args()

    assert abs(args.train_ratio + args.val_ratio + args.test_ratio - 1.0) < 1e-6

    logger.info(f"Loading {args.input_csv}")
    df = pd.read_csv(args.input_csv)
    logger.info(f"Loaded {len(df)} pairs from {df['pdb_id'].nunique()} PDBs")

    # Get unique PDB IDs and their dates
    pdb_ids = df['pdb_id'].unique()
    logger.info(f"Parsing deposition dates for {len(pdb_ids)} PDBs...")

    pdb_dates = {}
    missing = []
    for i, pid in enumerate(pdb_ids):
        d = get_pdb_date(pid, args.pdb_dir)
        if d:
            pdb_dates[pid] = d
        else:
            missing.append(pid)
        if (i + 1) % 10000 == 0:
            logger.info(f"Parsed {i+1}/{len(pdb_ids)} PDBs, missing {len(missing)}")

    if missing:
        logger.warning(f"Could not parse dates for {len(missing)} PDBs: {missing[:10]}...")
        logger.info(f"These PDBs will be assigned to the test set")
    else:
        missing = []

    # Sort PDBs with known dates by deposition date
    sorted_pdbs = sorted(pdb_dates.keys(), key=lambda p: pdb_dates[p])
    n_known = len(sorted_pdbs)
    n_train = int(n_known * args.train_ratio)
    n_val = int(n_known * args.val_ratio)
    # Known-date test portion + all missing-date PDBs go to test
    known_test_pdbs = sorted_pdbs[n_train + n_val:]

    train_pdbs = sorted_pdbs[:n_train]
    val_pdbs = sorted_pdbs[n_train:n_train + n_val]
    test_pdbs = known_test_pdbs + missing

    logger.info(f"Split by PDB count: train={len(train_pdbs)}, val={len(val_pdbs)}, test={len(test_pdbs)} (including {len(missing)} missing-date PDBs in test)")

    # Date ranges
    logger.info(f"Train date range: {pdb_dates[train_pdbs[0]]} ~ {pdb_dates[train_pdbs[-1]]}")
    logger.info(f"Val   date range: {pdb_dates[val_pdbs[0]]} ~ {pdb_dates[val_pdbs[-1]]}")
    known_test_dates = [pdb_dates[p] for p in test_pdbs if p in pdb_dates]
    if known_test_dates:
        logger.info(f"Test  date range (known dates): {min(known_test_dates)} ~ {max(known_test_dates)}, plus {len(missing)} PDBs without dates")

    # Assign splits
    pdb_to_split = {}
    for p in train_pdbs:
        pdb_to_split[p] = 'train'
    for p in val_pdbs:
        pdb_to_split[p] = 'val'
    for p in test_pdbs:
        pdb_to_split[p] = 'test'

    df['_split'] = df['pdb_id'].map(pdb_to_split)

    std_cols = [
        "unique_id", "file_type", "aa_seq", "labels", "pdb_id", "chain_id",
        "ligand_id", "ligand_name", "ligand_chain", "ligand_resnum",
        "num_binding_residues", "num_ligand_heavy_atoms", "protein_length"
    ]

    split_dir = Path(args.split_dir)
    split_dir.mkdir(parents=True, exist_ok=True)

    for split_name in ['train', 'val', 'test']:
        sub = df[df['_split'] == split_name][std_cols]
        out_path = split_dir / f"ligand_binding_site_{split_name}_dedup95_time.csv"
        sub.to_csv(out_path, index=False)
        logger.info(f"{split_name}: {len(sub)} pairs from {sub['pdb_id'].nunique()} PDBs -> {out_path}")

    # Also save the split mapping
    mapping_df = pd.DataFrame([
        {'pdb_id': p, 'deposition_date': pdb_dates.get(p, 'NA'), 'split': s}
        for p, s in pdb_to_split.items()
    ])
    mapping_path = split_dir / "ligand_binding_site_dedup95_time_split_mapping.csv"
    mapping_df.to_csv(mapping_path, index=False)
    logger.info(f"Saved split mapping to {mapping_path}")

    # Leakage check
    train_pdb = set(train_pdbs)
    val_pdb = set(val_pdbs)
    test_pdb = set(test_pdbs)
    logger.info(f"PDB overlap train/val: {len(train_pdb & val_pdb)}, train/test: {len(train_pdb & test_pdb)}, val/test: {len(val_pdb & test_pdb)}")

    logger.info("Done")


if __name__ == "__main__":
    main()
