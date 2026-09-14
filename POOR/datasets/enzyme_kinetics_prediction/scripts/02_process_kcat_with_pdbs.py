#!/usr/bin/env python3
"""
Process kcat dataset:
1. Filter out entries without local PDB files
2. Copy corresponding PDB files to enzyme_kinetics_prediction/pdbs/
3. Check PDB sequence consistency with dataset sequences
4. Analyze substrate type distribution
"""
import pandas as pd
import numpy as np
import gzip
import shutil
import json
from pathlib import Path
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio.Data.IUPACData import protein_letters_3to1

# Paths
TASK_DIR = Path(__file__).resolve().parent.parent  # datasets/enzyme_kinetics_prediction
PROJECT_ROOT = TASK_DIR.parent.parent              # POOD repo root

KCAT_CSV = TASK_DIR / "data/catpred-db-github/datasets/processed/kcat_max_wt_singleSeqs_wpdbs.csv"
AF_DIR = PROJECT_ROOT / "data/alphafold_structures"
PDB_OUT_DIR = TASK_DIR / "pdbs"
OUT_DIR = TASK_DIR / "data"
OUT_DIR.mkdir(parents=True, exist_ok=True)
PDB_OUT_DIR.mkdir(parents=True, exist_ok=True)

# Standard amino acid 3-to-1 mapping
AA_MAP = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
}


def parse_pdb_sequence(pdb_file):
    """Parse sequence from PDB ATOM records."""
    if str(pdb_file).endswith('.gz'):
        fopen = gzip.open
        mode = 'rt'
    else:
        fopen = open
        mode = 'r'
    
    residues = []
    seen = set()
    with fopen(pdb_file, mode) as f:
        for line in f:
            if line.startswith('ATOM'):
                res_name = line[17:20].strip()
                chain = line[21].strip()
                res_seq = int(line[22:26])
                key = (chain, res_seq)
                if key not in seen and res_name in AA_MAP:
                    seen.add(key)
                    residues.append(AA_MAP[res_name])
    return ''.join(residues)


def check_one_uniprot(args):
    """Check sequence consistency for one UniProt."""
    uniprot, dataset_seq, pdb_file = args
    try:
        pdb_seq = parse_pdb_sequence(pdb_file)
        if len(dataset_seq) == 0:
            return uniprot, 'empty_dataset_seq', 0, 0, 0
        
        # Check exact match
        exact_match = (pdb_seq == dataset_seq)
        
        # Check length difference
        len_diff = len(pdb_seq) - len(dataset_seq)
        
        # Check identity if lengths differ
        if len_diff == 0:
            identity = sum(a == b for a, b in zip(pdb_seq, dataset_seq)) / len(dataset_seq)
        else:
            # Substring check
            if dataset_seq in pdb_seq:
                identity = len(dataset_seq) / len(pdb_seq)
            elif pdb_seq in dataset_seq:
                identity = len(pdb_seq) / len(dataset_seq)
            else:
                # Use local alignment-ish: count matching positions in overlap
                min_len = min(len(pdb_seq), len(dataset_seq))
                max_len = max(len(pdb_seq), len(dataset_seq))
                matches = sum(a == b for a, b in zip(pdb_seq[:min_len], dataset_seq[:min_len]))
                identity = matches / max_len
        
        status = 'exact_match' if exact_match else ('length_diff' if len_diff != 0 else 'mismatch')
        return uniprot, status, len(pdb_seq), len(dataset_seq), identity
    except Exception as e:
        return uniprot, f'error:{e}', 0, 0, 0


def main():
    print("=" * 80)
    print("Step 1: Load kcat dataset and identify UniProts with local PDB files")
    print("=" * 80)
    
    df = pd.read_csv(KCAT_CSV)
    print(f"Original kcat dataset: {len(df):,} rows, {df['uniprot'].nunique():,} unique UniProts")
    
    # Find local PDB files
    local_pdbs = {}
    for f in AF_DIR.glob('AF-*.pdb.gz'):
        # AF-<uniprot>-F1-model_v6.pdb.gz
        parts = f.name.split('-')
        if len(parts) >= 4:
            uniprot = parts[1]
            local_pdbs[uniprot] = f
    
    print(f"Local PDB files: {len(local_pdbs):,}")
    
    # Filter dataset
    has_pdb = df['uniprot'].isin(local_pdbs)
    df_filtered = df[has_pdb].copy()
    missing_uniprots = set(df['uniprot'].unique()) - set(local_pdbs.keys())
    
    print(f"Rows with PDB: {len(df_filtered):,} / {len(df):,}")
    print(f"Unique UniProts with PDB: {df_filtered['uniprot'].nunique():,} / {df['uniprot'].nunique():,}")
    print(f"Missing UniProts: {len(missing_uniprots)}")
    
    # Save filtered dataset
    filtered_csv = OUT_DIR / 'kcat_filtered_with_pdb.csv'
    df_filtered.to_csv(filtered_csv, index=False)
    print(f"Saved filtered dataset to: {filtered_csv}")
    
    # Save missing uniprots
    missing_file = OUT_DIR / 'kcat_uniprots_removed_no_pdb.txt'
    with open(missing_file, 'w') as f:
        for u in sorted(missing_uniprots):
            f.write(u + '\n')
    print(f"Saved removed UniProts to: {missing_file}")
    
    print("\n" + "=" * 80)
    print("Step 2: Copy PDB files to enzyme_kinetics_prediction/pdbs/")
    print("=" * 80)
    
    needed_uniprots = df_filtered['uniprot'].unique()
    copied = 0
    for uniprot in needed_uniprots:
        src = local_pdbs[uniprot]
        dst = PDB_OUT_DIR / src.name
        if not dst.exists():
            shutil.copy2(src, dst)
            copied += 1
    
    print(f"Copied {copied:,} PDB files to {PDB_OUT_DIR}")
    print(f"Total PDB files in {PDB_OUT_DIR}: {len(list(PDB_OUT_DIR.glob('*.pdb.gz'))):,}")
    
    print("\n" + "=" * 80)
    print("Step 3: Check PDB sequence consistency with dataset sequences")
    print("=" * 80)
    
    # Get representative sequence per UniProt
    uniprot_to_seq = df_filtered.groupby('uniprot')['sequence'].first().to_dict()
    
    args_list = [(u, uniprot_to_seq[u], local_pdbs[u]) for u in needed_uniprots]
    
    results = []
    with ProcessPoolExecutor(max_workers=16) as executor:
        for result in executor.map(check_one_uniprot, args_list):
            results.append(result)
    
    consistency_df = pd.DataFrame(results, columns=['uniprot', 'status', 'pdb_seq_len', 'dataset_seq_len', 'identity'])
    consistency_csv = OUT_DIR / 'kcat_pdb_sequence_consistency.csv'
    consistency_df.to_csv(consistency_csv, index=False)
    
    print(f"Saved consistency check to: {consistency_csv}")
    print(f"\nConsistency summary:")
    print(consistency_df['status'].value_counts())
    print(f"\nIdentity stats:")
    print(consistency_df['identity'].describe())
    print(f"Median identity: {consistency_df['identity'].median():.4f}")
    print(f"Proteins with identity < 0.95: {(consistency_df['identity'] < 0.95).sum()}")
    print(f"Proteins with identity < 0.90: {(consistency_df['identity'] < 0.90).sum()}")
    print(f"Proteins with identity < 0.80: {(consistency_df['identity'] < 0.80).sum()}")
    
    print("\n" + "=" * 80)
    print("Step 4: Analyze substrate type distribution")
    print("=" * 80)
    
    # Basic substrate stats
    n_unique_reactions = df_filtered['reaction_smiles'].nunique()
    n_non_null_reactions = df_filtered['reaction_smiles'].notna().sum()
    print(f"Total rows: {len(df_filtered):,}")
    print(f"Rows with reaction_smiles: {n_non_null_reactions:,} ({n_non_null_reactions/len(df_filtered)*100:.1f}%)")
    print(f"Unique reaction_smiles: {n_unique_reactions:,}")
    
    # Most common reactions
    print(f"\nTop 20 most common reactions:")
    reaction_counts = df_filtered['reaction_smiles'].value_counts().head(20)
    for i, (rxn, count) in enumerate(reaction_counts.items(), 1):
        print(f"  {i:2d}. count={count:5d}  rxn={rxn[:100]}{'...' if len(rxn) > 100 else ''}")
    
    # EC distribution
    print(f"\nTop 20 EC codes:")
    ec_counts = df_filtered['ec'].value_counts().head(20)
    for i, (ec, count) in enumerate(ec_counts.items(), 1):
        print(f"  {i:2d}. EC={ec:<15s} count={count:5d}")
    
    # EC class distribution (first digit)
    df_filtered['ec_class'] = df_filtered['ec'].astype(str).str.split('.').str[0]
    print(f"\nEC class distribution:")
    print(df_filtered['ec_class'].value_counts().sort_index())
    
    # Save substrate stats
    substrate_stats = {
        'total_rows': len(df_filtered),
        'rows_with_reaction_smiles': int(n_non_null_reactions),
        'unique_reaction_smiles': int(n_unique_reactions),
        'top_reactions': reaction_counts.to_dict(),
        'top_ec_codes': ec_counts.to_dict(),
        'ec_class_distribution': df_filtered['ec_class'].value_counts().sort_index().to_dict(),
    }
    with open(OUT_DIR / 'kcat_substrate_stats.json', 'w') as f:
        json.dump(substrate_stats, f, indent=2)
    print(f"\nSaved substrate stats to: {OUT_DIR / 'kcat_substrate_stats.json'}")
    
    print("\n" + "=" * 80)
    print("All steps completed!")
    print("=" * 80)


if __name__ == '__main__':
    main()
