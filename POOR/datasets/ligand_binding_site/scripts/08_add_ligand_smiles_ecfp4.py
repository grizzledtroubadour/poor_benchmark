#!/usr/bin/env python3
"""Add canonical SMILES and ECFP4 fingerprint columns directly to split CSVs."""
import argparse
import gzip
import logging
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s')
logger = logging.getLogger(__name__)

TASK_DIR = Path(__file__).resolve().parent.parent  # datasets/ligand_binding_site
ROOT = TASK_DIR.parent.parent                      # 项目根（POOR/）
SPLIT_DIR = TASK_DIR / "splits"
CCD_PATH = ROOT / "data" / "pdb" / "components.cif.gz"


def parse_ccd_smiles(cif_path):
    """Extract comp_id -> list of SMILES/InChI candidates from PDB CCD."""
    logger.info(f"Parsing descriptors from {cif_path}")
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
                    if descriptor and descriptor not in ('?', '.') and comp_id:
                        descriptor = descriptor.strip('"').strip("'")
                        if 'SMILES' in dtype or dtype == 'InChI':
                            smiles_map.setdefault(comp_id, []).append(descriptor)
                else:
                    in_descriptor_loop = False
            elif in_descriptor_loop and (line.startswith('#') or line.startswith('data_')):
                in_descriptor_loop = False
    logger.info(f"Parsed descriptors for {len(smiles_map)} ligands")
    return smiles_map


def _mol_to_fp(mol, radius=2, n_bits=2048):
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=n_bits)
    bitstring = fp.ToBitString()
    # store as a compact sparse list of on-bit indices
    indices = [str(i) for i, ch in enumerate(bitstring) if ch == '1']
    return '[' + ','.join(indices) + ']'


def _try_smiles(smiles, radius=2, n_bits=2048):
    if not smiles:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            canon = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
            return canon, _mol_to_fp(mol, radius=radius, n_bits=n_bits)
    except Exception:
        pass

    # Fallback: some aromatic/charged CCD SMILES fail kekulization.
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is not None:
            Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
            canon = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
            return canon, _mol_to_fp(mol, radius=radius, n_bits=n_bits)
    except Exception:
        pass
    return None


def _try_inchi(inchi, radius=2, n_bits=2048):
    if not inchi or not inchi.startswith('InChI='):
        return None
    try:
        mol = Chem.MolFromInchi(inchi, treatWarningAsError=False)
        if mol is not None:
            Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
            canon = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
            return canon, _mol_to_fp(mol, radius=radius, n_bits=n_bits)
    except Exception:
        pass
    return None


def canonical_smiles_and_ecfp4(candidates, radius=2, n_bits=2048):
    """Return (canonical_smiles, ecfp4_sparse_list) or ('', '[]') on failure."""
    if not candidates:
        return '', '[]'
    if isinstance(candidates, str):
        candidates = [candidates]
    for cand in candidates:
        res = _try_smiles(cand, radius=radius, n_bits=n_bits)
        if res is not None:
            return res
    # If no SMILES works, try any InChI candidate.
    for cand in candidates:
        res = _try_inchi(cand, radius=radius, n_bits=n_bits)
        if res is not None:
            return res
    return '', '[]'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--split_dir", type=str, default=str(SPLIT_DIR))
    parser.add_argument("--ccd_path", type=str, default=str(CCD_PATH))
    parser.add_argument("--radius", type=int, default=2)
    parser.add_argument("--n_bits", type=int, default=2048)
    args = parser.parse_args()

    split_dir = Path(args.split_dir)
    smiles_map = parse_ccd_smiles(args.ccd_path)

    for split in ['train', 'val', 'test']:
        path = split_dir / f"ligand_binding_site_{split}.csv"
        logger.info(f"Processing {path}")
        df = pd.read_csv(path, keep_default_na=False, na_values=[''])

        # ligand_name is the 3rd component of unique_id
        df['ligand_name'] = df['unique_id'].str.split('_').str[2]

        unique_ligands = df['ligand_name'].unique()
        logger.info(f"  {split}: {len(df)} rows, {len(unique_ligands)} unique ligands")

        cache = {}
        failed = []
        for lig in unique_ligands:
            cands = smiles_map.get(lig, [])
            if not cands:
                failed.append(lig)
                cache[lig] = ('', '[]')
                continue
            canon, fp = canonical_smiles_and_ecfp4(cands, radius=args.radius, n_bits=args.n_bits)
            if not canon:
                failed.append(lig)
                cache[lig] = ('', '[]')
                continue
            cache[lig] = (canon, fp)

        if failed:
            logger.warning(f"  Failed to parse SMILES/FP for {len(failed)} ligands: {failed[:10]}...")

        df['ligand_smiles'] = df['ligand_name'].map(lambda x: cache[x][0])
        df['ligand_ecfp4'] = df['ligand_name'].map(lambda x: cache[x][1])

        # Drop temporary ligand_name column (user requested no extra metadata columns)
        df = df.drop(columns=['ligand_name'])

        # Place ligand info right after unique_id for readability
        cols = list(df.columns)
        core = ['unique_id']
        ligand_cols = ['ligand_smiles', 'ligand_ecfp4']
        other = [c for c in cols if c not in core + ligand_cols]
        df = df[core + ligand_cols + other]

        df.to_csv(path, index=False)
        missing = (df['ligand_smiles'] == '').sum()
        logger.info(f"  Saved {path}; missing SMILES: {missing}")


if __name__ == "__main__":
    main()
