#!/usr/bin/env python3
"""
Step 3: Extract chain-level PDB files for extreme-length chains.

For each chain in the extreme_length_chains.csv that was added to test sets,
extract the chain from the raw PDB file and save as a chain-level PDB file
in the pdbs/ directory.

Uses full-chain extraction (C-grade): extract all ATOM/HETATM records
for the target chain_id from the raw PDB file.
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


def extract_chain_pdb(pdb_raw_path, chain_id, output_path):
    """Extract a single chain from a raw PDB file."""
    opener = gzip.open if pdb_raw_path.suffix == '.gz' else open
    lines_written = 0
    try:
        with opener(pdb_raw_path, 'rt') as f_in:
            with open(output_path, 'w') as f_out:
                for line in f_in:
                    # Multi-model (NMR) entries: keep model 1 only.
                    if line.startswith('ENDMDL'):
                        break
                    if line.startswith('ATOM  ') or line.startswith('HETATM'):
                        chain = line[21:22].strip()
                        if chain == chain_id:
                            f_out.write(line)
                            lines_written += 1
                    elif line.startswith('TER') and lines_written > 0:
                        chain = line[21:22].strip()
                        if chain == chain_id:
                            f_out.write(line)
                            break
        if lines_written == 0:
            os.unlink(output_path)
            return False
        return True
    except Exception as e:
        if os.path.exists(output_path):
            os.unlink(output_path)
        print(f"Error extracting {pdb_raw_path.name} chain {chain_id}: {e}")
        return False


def process_chain(args):
    """Process a single chain extraction."""
    pdb_id, chain_id, is_cif = args
    chain_key = f'{pdb_id}_{chain_id}'
    ext = '.cif' if is_cif else '.pdb'
    output_path = PDBS_DIR / f'{chain_key}{ext}'

    if output_path.exists():
        return chain_key, True  # Already exists

    # Find raw PDB file
    pdb_path = PDB_RAW_DIR / f'{pdb_id}.pdb.gz'
    if not pdb_path.exists():
        pdb_path = PDB_RAW_DIR / f'{pdb_id}.cif.gz'
    if not pdb_path.exists():
        pdb_path = PDB_RAW_DIR / f'{pdb_id}.pdb'
    if not pdb_path.exists():
        pdb_path = PDB_RAW_DIR / f'{pdb_id}.cif'
    if not pdb_path.exists():
        return chain_key, False

    # For CIF files, extract the chain
    # (mmCIF chain extraction is more complex; for now, we use the PDB format)
    success = extract_chain_pdb(pdb_path, chain_id, output_path)
    return chain_key, success


def main():
    print("=" * 60)
    print("Step 3: Extract PDB files for extreme-length chains")
    print("=" * 60)

    # Load extreme chains
    ext_df = pd.read_csv(BASE / 'output' / 'extreme_length_chains.csv')

    # Get chains that were actually added to test sets
    needed = set()
    for task in ['ec', 'go_mf', 'go_cc', 'go_bp']:
        test_df = pd.read_csv(BASE / 'splits' / f'{task}_test.csv')
        new_keys = test_df[test_df['OOD_ExtremeShort'] | test_df['OOD_ExtremeLong']]['unique_id']
        needed.update(new_keys)

    print(f"Chains needing PDB files: {len(needed)}")

    # Prepare tasks
    tasks = []
    for _, row in ext_df.iterrows():
        ck = row['chain_key']
        if ck not in needed:
            continue
        tasks.append((row['pdb_id'], row['chain_id'], row.get('is_cif', False)))

    print(f"Processing {len(tasks)} extractions...")

    success = 0
    failed = 0
    with ProcessPoolExecutor(max_workers=16) as executor:
        futures = {executor.submit(process_chain, t): t for t in tasks}
        for i, future in enumerate(as_completed(futures)):
            if i % 200 == 0:
                print(f"  Processed {i}/{len(tasks)}...")
            chain_key, ok = future.result()
            if ok:
                success += 1
            else:
                failed += 1
                print(f"  FAILED: {chain_key}")

    print(f"\nResults: {success} extracted, {failed} failed")
    print(f"Total PDB files in pdbs/: {len(os.listdir(PDBS_DIR))}")


if __name__ == '__main__':
    main()
