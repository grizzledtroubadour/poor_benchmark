#!/usr/bin/env python3
"""Filter test/val samples by ligand Tanimoto similarity to training set ligands."""
import argparse
import gzip
import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s')
logger = logging.getLogger(__name__)

TASK_DIR = Path(__file__).resolve().parent.parent  # datasets/ligand_binding_site
ROOT = TASK_DIR.parent.parent                      # 项目根（POOR/）
OUTPUT_DIR = str(TASK_DIR / "output")
SPLIT_DIR = str(TASK_DIR / "splits")
CCD_PATH = str(ROOT / "data" / "pdb" / "components.cif.gz")


def parse_ccd_smiles(cif_path):
    """Extract ligand_id -> SMILES mapping from PDB Chemical Component Dictionary."""
    logger.info(f"Parsing SMILES from {cif_path}")
    smiles_map = {}
    in_descriptor_loop = False
    cols = []
    with gzip.open(cif_path, 'rt', errors='ignore') as f:
        for line in f:
            line = line.strip()
            if line.startswith('data_'):
                in_descriptor_loop = False
                cols = []
            elif line.startswith('loop_'):
                in_descriptor_loop = False
                cols = []
            elif line.startswith('_pdbx_chem_comp_descriptor.'):
                if not in_descriptor_loop:
                    in_descriptor_loop = True
                    cols = []
                cols.append(line.split('.')[1])
            elif in_descriptor_loop and line and not line.startswith('_') and not line.startswith('#') and not line.startswith('data_'):
                parts = line.split()
                if len(parts) >= len(cols):
                    row = dict(zip(cols, parts))
                    comp_id = row.get('comp_id')
                    dtype = row.get('type', '')
                    descriptor = row.get('descriptor', '')
                    if 'SMILES' in dtype and descriptor and descriptor != '?' and comp_id:
                        # CIF descriptors are often quoted; strip surrounding quotes
                        descriptor = descriptor.strip('"').strip("'")
                        if comp_id not in smiles_map:
                            smiles_map[comp_id] = descriptor
                else:
                    in_descriptor_loop = False
            elif in_descriptor_loop and (line.startswith('#') or line.startswith('data_')):
                in_descriptor_loop = False
    logger.info(f"Parsed {len(smiles_map)} ligand SMILES")
    return smiles_map


def compute_fingerprint(smiles, radius=2, n_bits=2048):
    """Compute ECFP4 Morgan fingerprint from SMILES."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=n_bits)
        return fp
    except Exception:
        return None


def main():
    parser = argparse.ArgumentParser(description="Filter test/val by ligand similarity to train")
    parser.add_argument("--train_csv", type=str, default=os.path.join(SPLIT_DIR, "ligand_binding_site_train_dedup95_time.csv"))
    parser.add_argument("--val_csv", type=str, default=os.path.join(SPLIT_DIR, "ligand_binding_site_val_dedup95_time.csv"))
    parser.add_argument("--test_csv", type=str, default=os.path.join(SPLIT_DIR, "ligand_binding_site_test_dedup95_time.csv"))
    parser.add_argument("--output_dir", type=str, default=SPLIT_DIR)
    parser.add_argument("--output_suffix", type=str, default="dedup95_time_ligsim60",
                        help="Filename suffix for outputs")
    parser.add_argument("--ccd_path", type=str, default=CCD_PATH)
    parser.add_argument("--threshold", type=float, default=0.6)
    parser.add_argument("--radius", type=int, default=2)
    parser.add_argument("--n_bits", type=int, default=2048)
    args = parser.parse_args()

    smiles_map = parse_ccd_smiles(args.ccd_path)

    logger.info("Loading splits")
    train_df = pd.read_csv(args.train_csv)
    val_df = pd.read_csv(args.val_csv)
    test_df = pd.read_csv(args.test_csv)

    # Compute fingerprints for all unique train ligands
    train_ligands = train_df['ligand_name'].unique()
    logger.info(f"Computing fingerprints for {len(train_ligands)} train ligands")
    train_fps = {}
    failed_train = []
    for lig in train_ligands:
        smi = smiles_map.get(lig)
        if not smi:
            failed_train.append(lig)
            continue
        fp = compute_fingerprint(smi, radius=args.radius, n_bits=args.n_bits)
        if fp is not None:
            train_fps[lig] = fp
        else:
            failed_train.append(lig)

    if failed_train:
        logger.warning(f"Could not compute fingerprints for {len(failed_train)} train ligands: {failed_train[:10]}...")

    if not train_fps:
        raise ValueError("No valid train fingerprints; cannot filter")

    train_fp_list = list(train_fps.values())
    logger.info(f"Valid train fingerprints: {len(train_fp_list)}")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for split_name, split_df in [('val', val_df), ('test', test_df)]:
        logger.info(f"Processing {split_name}: {len(split_df)} pairs")
        max_sims = []
        keep_mask = []
        missing_smiles = 0
        failed_fp = 0

        for _, row in split_df.iterrows():
            lig = row['ligand_name']
            smi = smiles_map.get(lig)
            if not smi:
                max_sims.append(-1.0)
                keep_mask.append(False)
                missing_smiles += 1
                continue
            fp = compute_fingerprint(smi, radius=args.radius, n_bits=args.n_bits)
            if fp is None:
                max_sims.append(-1.0)
                keep_mask.append(False)
                failed_fp += 1
                continue
            # Compute max Tanimoto to train fingerprints
            max_sim = DataStructs.BulkTanimotoSimilarity(fp, train_fp_list)
            max_sim = max(max_sim) if max_sim else 0.0
            max_sims.append(max_sim)
            keep_mask.append(max_sim > args.threshold)

        split_df = split_df.copy()
        split_df['max_ligand_tanimoto_to_train'] = max_sims
        filtered_df = split_df[keep_mask].drop(columns=['max_ligand_tanimoto_to_train'])

        logger.info(f"{split_name}: kept {len(filtered_df)}/{len(split_df)} pairs "
                    f"(missing SMILES: {missing_smiles}, failed FP: {failed_fp}, "
                    f"filtered by threshold: {len(split_df) - len(filtered_df) - missing_smiles - failed_fp})")

        # Save filtered split
        out_path = output_dir / f"ligand_binding_site_{split_name}_{args.output_suffix}.csv"
        filtered_df.to_csv(out_path, index=False)
        logger.info(f"Saved {out_path}")

        # Save similarity report
        report_path = output_dir / f"ligand_binding_site_{split_name}_{args.output_suffix}_report.csv"
        split_df[['unique_id', 'pdb_id', 'chain_id', 'ligand_name', 'max_ligand_tanimoto_to_train']].to_csv(report_path, index=False)
        logger.info(f"Saved similarity report {report_path}")

    logger.info("Done")


if __name__ == "__main__":
    main()
