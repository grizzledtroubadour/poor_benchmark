#!/usr/bin/env python3
"""
Filter kcat dataset to keep only entries with valid reactant SMILES,
compute ECFP4 fingerprints for all reactants (combined via bitwise OR),
and exclude rows where no fingerprint can be generated.
"""
import pandas as pd
import numpy as np
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger

# Suppress RDKit warnings
RDLogger.DisableLog('rdApp.*')

# Paths
TASK_DIR = Path(__file__).resolve().parent.parent  # datasets/enzyme_kinetics_prediction
DATA_DIR = TASK_DIR / "data"
OUT_DIR = TASK_DIR / "output"
OUT_DIR.mkdir(parents=True, exist_ok=True)

INPUT_CSV = DATA_DIR / "kcat_filtered_with_pdb.csv"

# ECFP4 parameters
FP_RADIUS = 2  # ECFP4
FP_NBITS = 2048


def compute_combined_ecfp4(reactant_smiles):
    """
    Compute combined ECFP4 fingerprint for all reactant molecules.
    Returns (fingerprint_array, status, list_of_mol_info).
    """
    if pd.isna(reactant_smiles) or str(reactant_smiles).strip() == '':
        return None, 'empty_reactant', []
    
    reactant_smiles = str(reactant_smiles).strip()
    
    # Split reactants by '.'
    mol_smiles_list = [s.strip() for s in reactant_smiles.split('.') if s.strip()]
    
    if len(mol_smiles_list) == 0:
        return None, 'empty_reactant', []
    
    combined_fp = None
    mol_info = []
    valid_mols = 0
    
    for smi in mol_smiles_list:
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                mol_info.append({'smiles': smi, 'valid': False, 'reason': 'parse_failed'})
                continue
            
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=FP_RADIUS, nBits=FP_NBITS)
            fp_array = np.array(fp)
            
            if combined_fp is None:
                combined_fp = fp_array.copy()
            else:
                combined_fp = np.maximum(combined_fp, fp_array)  # bitwise OR
            
            valid_mols += 1
            mol_info.append({'smiles': smi, 'valid': True, 'num_bits': int(fp_array.sum())})
        except Exception as e:
            mol_info.append({'smiles': smi, 'valid': False, 'reason': str(e)})
    
    if combined_fp is None:
        return None, 'no_valid_fingerprint', mol_info
    
    return combined_fp, 'success', mol_info


def main():
    print("=" * 80)
    print("Step 1: Load filtered kcat dataset")
    print("=" * 80)
    
    df = pd.read_csv(INPUT_CSV)
    print(f"Original rows: {len(df):,}")
    
    # Step 1: Filter rows with non-empty reactant_smiles
    has_reactant = df['reactant_smiles'].notna() & (df['reactant_smiles'].astype(str).str.len() > 0)
    df = df[has_reactant].copy()
    print(f"After excluding empty reactant_smiles: {len(df):,}")
    
    print("\n" + "=" * 80)
    print("Step 2: Compute ECFP4 fingerprints")
    print("=" * 80)
    
    fingerprints = []
    statuses = []
    fp_bit_counts = []
    valid_mol_counts = []
    total_mol_counts = []
    
    for idx, reactant_smiles in enumerate(df['reactant_smiles']):
        fp, status, mol_info = compute_combined_ecfp4(reactant_smiles)
        fingerprints.append(fp)
        statuses.append(status)
        fp_bit_counts.append(int(fp.sum()) if fp is not None else 0)
        valid_mol_counts.append(sum(1 for m in mol_info if m['valid']))
        total_mol_counts.append(len(mol_info))
        
        if (idx + 1) % 2000 == 0:
            print(f"  Processed {idx + 1}/{len(df)} rows...")
    
    df['fp_status'] = statuses
    df['fp_bit_count'] = fp_bit_counts
    df['fp_valid_mols'] = valid_mol_counts
    df['fp_total_mols'] = total_mol_counts
    
    print(f"\nFingerprint status distribution:")
    print(pd.Series(statuses).value_counts())
    
    # Filter rows with successful fingerprints
    df_with_fp = df[df['fp_status'] == 'success'].copy()
    print(f"\nAfter excluding rows without valid fingerprints: {len(df_with_fp):,}")
    
    # Add fingerprint columns efficiently
    fp_matrix = np.vstack([fp for fp in fingerprints if fp is not None])
    
    # Create fingerprint DataFrame separately to avoid fragmentation
    fp_cols = [f'ecfp4_bit_{i}' for i in range(FP_NBITS)]
    fp_df = pd.DataFrame(fp_matrix, columns=fp_cols, index=df_with_fp.index)
    df_with_fp = pd.concat([df_with_fp, fp_df], axis=1)
    
    # Save fingerprint as string for reference
    df_with_fp['ecfp4'] = [''.join(map(str, row)) for row in fp_matrix]
    
    # Save dataset
    output_csv = DATA_DIR / 'kcat_with_reactant_ecfp4.csv'
    df_with_fp.to_csv(output_csv, index=False)
    print(f"\nSaved dataset with ECFP4 to: {output_csv}")
    
    # Save summary
    summary = {
        'original_rows': int(len(pd.read_csv(INPUT_CSV))),
        'after_reactant_filter': int(len(df)),
        'after_fp_filter': int(len(df_with_fp)),
        'fp_radius': FP_RADIUS,
        'fp_nbits': FP_NBITS,
        'status_distribution': pd.Series(statuses).value_counts().to_dict(),
        'mean_bit_count': float(np.mean(fp_bit_counts)),
        'median_bit_count': float(np.median(fp_bit_counts)),
    }
    import json
    with open(OUT_DIR / 'kcat_ecfp4_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"Saved summary to: {OUT_DIR / 'kcat_ecfp4_summary.json'}")
    
    print("\n" + "=" * 80)
    print("All steps completed!")
    print("=" * 80)


if __name__ == '__main__':
    main()
