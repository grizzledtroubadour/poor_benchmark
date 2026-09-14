#!/usr/bin/env python3
"""Aggregate all binding_sites_batch_*.csv files into final datasets."""
import os
import glob
import pandas as pd
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s')
logger = logging.getLogger(__name__)

TASK_DIR = Path(__file__).resolve().parent.parent  # datasets/ligand_binding_site
OUTPUT_DIR = str(TASK_DIR / "output")

def main():
    batch_files = sorted(glob.glob(os.path.join(OUTPUT_DIR, "binding_sites_batch_*.csv")))
    logger.info(f"Found {len(batch_files)} batch files")
    
    if not batch_files:
        logger.error("No batch files found")
        return
    
    # Read and concatenate all batch CSVs
    dfs = []
    for f in batch_files:
        df = pd.read_csv(f)
        dfs.append(df)
    
    combined = pd.concat(dfs, ignore_index=True)
    logger.info(f"Combined rows: {len(combined)}")
    
    # Deduplicate by unique_id
    before_dedup = len(combined)
    combined = combined.drop_duplicates(subset=["unique_id"], keep="first")
    after_dedup = len(combined)
    logger.info(f"Deduplicated: {before_dedup} -> {after_dedup} rows")
    
    # Write raw format (same columns as batch files)
    raw_path = os.path.join(OUTPUT_DIR, "ligand_binding_sites_raw.csv")
    combined.to_csv(raw_path, index=False)
    logger.info(f"Wrote {raw_path}: {len(combined)} rows")
    
    # Write std format (different column order)
    std_cols = [
        "unique_id", "file_type", "aa_seq", "labels", "pdb_id", "chain_id",
        "ligand_id", "ligand_name", "ligand_chain", "ligand_resnum",
        "num_binding_residues", "num_ligand_heavy_atoms", "protein_length"
    ]
    std_path = os.path.join(OUTPUT_DIR, "ligand_binding_sites_std.csv")
    combined[std_cols].to_csv(std_path, index=False)
    logger.info(f"Wrote {std_path}: {len(combined)} rows")
    
    # Summary stats
    logger.info(f"Unique PDB files: {combined['pdb_id'].nunique()}")
    logger.info(f"Unique ligands: {combined['ligand_name'].nunique()}")
    logger.info(f"Total chain-ligand pairs: {len(combined)}")
    logger.info("Top 10 ligands:")
    for ligand, count in combined['ligand_name'].value_counts().head(10).items():
        logger.info(f"  {ligand}: {count}")

if __name__ == "__main__":
    main()
