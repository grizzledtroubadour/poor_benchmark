#!/usr/bin/env python3
"""
Step 3b: Extract chain-level CIF files for CIF-format extreme-length chains.

Uses direct text-based extraction (no BioPython parsing) for speed.
Filters _atom_site lines by chain_id (auth_asym_id).
"""

import gzip
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parent.parent  # func_prediction/
POOR_ROOT = BASE.parent.parent  # POOR/
PDB_RAW_DIR = POOR_ROOT / 'data' / 'pdb'
PDBS_DIR = BASE / 'pdbs'


def extract_chain_cif_text(cif_path, target_chain_id, output_path):
    """Extract a single chain from a mmCIF.gz file using text-based filtering."""
    atom_site_cols = {}
    in_atom_site = False
    atom_count = 0
    chain_col = None
    comp_col = None

    try:
        with gzip.open(cif_path, 'rt') as f_in:
            with open(output_path, 'w') as f_out:
                for line in f_in:
                    line_stripped = line.rstrip('\n')

                    # Track _atom_site column definitions
                    if line_stripped.startswith('_atom_site.'):
                        if not in_atom_site:
                            in_atom_site = True
                        field = line_stripped.split('.', 1)[1].strip()
                        atom_site_cols[field] = len(atom_site_cols)
                        if field == 'auth_asym_id':
                            chain_col = len(atom_site_cols) - 1
                        elif field == 'label_comp_id':
                            comp_col = len(atom_site_cols) - 1
                        # Write the column definition to output
                        f_out.write(line)
                        continue

                    # End of atom_site block
                    if line_stripped.startswith('#') or line_stripped.startswith('loop_'):
                        in_atom_site = False
                        # Write non-atom_site lines to output for metadata
                        f_out.write(line)
                        continue

                    # Skip non-atom_site content
                    if not in_atom_site:
                        if not line_stripped.startswith('_'):
                            f_out.write(line)
                        continue

                    # Skip atom_site data lines that don't start with ATOM/HETATM
                    parts = line_stripped.split()
                    if len(parts) < 5:
                        continue
                    if parts[0] not in ('ATOM', 'HETATM'):
                        continue

                    # Filter by chain_id
                    if chain_col is not None and parts[chain_col] != target_chain_id:
                        continue

                    f_out.write(line)
                    atom_count += 1

        if atom_count == 0:
            os.unlink(output_path)
            return False
        return True
    except Exception as e:
        if output_path.exists():
            os.unlink(output_path)
        print(f"Error: {e}")
        return False


def process_chain(args):
    pdb_id, chain_id = args
    chain_key = f'{pdb_id}_{chain_id}'
    output_path = PDBS_DIR / f'{chain_key}.cif'

    if output_path.exists():
        return chain_key, True

    cif_path = PDB_RAW_DIR / f'{pdb_id}.cif.gz'
    if not cif_path.exists():
        cif_path = PDB_RAW_DIR / f'{pdb_id}.cif'
    if not cif_path.exists():
        return chain_key, False

    success = extract_chain_cif_text(cif_path, chain_id, output_path)
    return chain_key, success


def main():
    print("=" * 60)
    print("Step 3b: Extract CIF files (text-based, fast)")
    print("=" * 60)

    # Load extreme chains
    ext_df = pd.read_csv(BASE / 'output' / 'extreme_length_chains.csv')

    # Get chains that were actually added to test sets
    needed = set()
    for task in ['ec', 'go_mf', 'go_cc', 'go_bp']:
        test_df = pd.read_csv(BASE / 'splits' / f'{task}_test.csv')
        new_keys = test_df[test_df['OOD_ExtremeShort'] | test_df['OOD_ExtremeLong']]['unique_id']
        needed.update(new_keys)

    # Filter to CIF chains not yet extracted
    cif_chains = ext_df[(ext_df['is_cif'] == True) & (ext_df['chain_key'].isin(needed))]
    cif_to_extract = []
    for _, row in cif_chains.iterrows():
        fp = PDBS_DIR / f'{row["chain_key"]}.cif'
        if not fp.exists():
            cif_to_extract.append((row['pdb_id'], row['chain_id']))

    print(f"CIF chains to extract: {len(cif_to_extract)}")

    # Update struct_file in test CSVs for CIF chains
    cif_keys = set(cif_chains['chain_key'])

    success = 0
    failed = 0
    with ProcessPoolExecutor(max_workers=16) as executor:
        futures = {executor.submit(process_chain, t): t for t in cif_to_extract}
        for i, future in enumerate(as_completed(futures)):
            if i % 20 == 0:
                print(f"  Processed {i}/{len(cif_to_extract)}...")
            chain_key, ok = future.result()
            if ok:
                success += 1
            else:
                failed += 1
                print(f"  FAILED: {chain_key}")

    print(f"\nResults: {success} extracted, {failed} failed")
    print(f"Total PDB files in pdbs/: {len(os.listdir(PDBS_DIR))}")

    # Update struct_file in test CSVs for CIF chains
    for task in ['ec', 'go_mf', 'go_cc', 'go_bp']:
        test_path = BASE / 'splits' / f'{task}_test.csv'
        df = pd.read_csv(test_path)
        mask = df['unique_id'].isin(cif_keys)
        df.loc[mask, 'struct_file'] = df.loc[mask, 'unique_id'] + '.cif'
        df.to_csv(test_path, index=False)
        updated = mask.sum()
        if updated > 0:
            print(f"  {task}: updated struct_file for {updated} CIF chains")


if __name__ == '__main__':
    main()
