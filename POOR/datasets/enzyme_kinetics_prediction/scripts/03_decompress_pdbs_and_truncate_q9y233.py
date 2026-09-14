#!/usr/bin/env python3
"""
1. Decompress all .pdb.gz files in enzyme_kinetics_prediction/pdbs/
2. Process Q9Y233 structure to keep only residues matching dataset sequence
"""
import gzip
import shutil
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

TASK_DIR = Path(__file__).resolve().parent.parent  # datasets/enzyme_kinetics_prediction

PDBS_DIR = TASK_DIR / "pdbs"
DATA_DIR = TASK_DIR / "data"


def decompress_one(gz_file):
    """Decompress a single .pdb.gz file and remove the archive."""
    pdb_file = gz_file.with_suffix('')  # removes .gz
    try:
        with gzip.open(gz_file, 'rb') as f_in, open(pdb_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        gz_file.unlink()
        return gz_file.name, "success"
    except Exception as e:
        return gz_file.name, f"error:{e}"


def parse_pdb_atoms(pdb_file):
    """Parse ATOM/HETATM records from PDB, return list of dicts."""
    atoms = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom = {
                    'record': line[:6].strip(),
                    'atom_serial': int(line[6:11]),
                    'atom_name': line[12:16].strip(),
                    'alt_loc': line[16].strip(),
                    'res_name': line[17:20].strip(),
                    'chain_id': line[21].strip(),
                    'res_seq': int(line[22:26]),
                    'insertion': line[26].strip(),
                    'x': float(line[30:38]),
                    'y': float(line[38:46]),
                    'z': float(line[46:54]),
                    'occupancy': line[54:60].strip(),
                    'temp_factor': line[60:66].strip(),
                    'element': line[76:78].strip(),
                    'charge': line[78:80].strip(),
                    'line': line,
                }
                atoms.append(atom)
    return atoms


def get_pdb_sequence_with_residues(atoms):
    """Extract sequence and residue ranges from PDB atoms."""
    AA_MAP = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
    }
    residues = []
    seen = set()
    for atom in atoms:
        if atom['res_name'] not in AA_MAP:
            continue
        key = (atom['chain_id'], atom['res_seq'], atom['insertion'])
        if key not in seen:
            seen.add(key)
            residues.append({
                'chain_id': atom['chain_id'],
                'res_seq': atom['res_seq'],
                'insertion': atom['insertion'],
                'res_name': atom['res_name'],
                'aa': AA_MAP[atom['res_name']],
            })
    seq = ''.join(r['aa'] for r in residues)
    return seq, residues


def find_best_window_position(pdb_seq, target_seq):
    """Find the best matching window of target_seq in pdb_seq."""
    if len(target_seq) > len(pdb_seq):
        return None, 0
    best_pos = -1
    best_matches = -1
    for i in range(len(pdb_seq) - len(target_seq) + 1):
        matches = sum(a == b for a, b in zip(pdb_seq[i:i+len(target_seq)], target_seq))
        if matches > best_matches:
            best_matches = matches
            best_pos = i
    return best_pos, best_matches


def extract_residues_to_pdb(input_pdb, output_pdb, start_res_idx, end_res_idx, residues):
    """Extract PDB lines for residues from start_res_idx to end_res_idx (inclusive)."""
    target_residues = set()
    for i in range(start_res_idx, end_res_idx + 1):
        r = residues[i]
        target_residues.add((r['chain_id'], r['res_seq'], r['insertion']))
    
    # Read original PDB and filter
    with open(input_pdb, 'r') as f_in, open(output_pdb, 'w') as f_out:
        for line in f_in:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                chain_id = line[21].strip()
                res_seq = int(line[22:26])
                insertion = line[26].strip()
                if (chain_id, res_seq, insertion) in target_residues:
                    f_out.write(line)
            elif line.startswith('TER') or line.startswith('END'):
                # Write END at the end only
                pass
        f_out.write('END\n')


def main():
    print("=" * 80)
    print("Step 1: Decompress all .pdb.gz files")
    print("=" * 80)
    
    gz_files = sorted(PDBS_DIR.glob('*.pdb.gz'))
    print(f"Found {len(gz_files)} .pdb.gz files")
    
    success = 0
    failed = 0
    with ProcessPoolExecutor(max_workers=16) as executor:
        for name, status in executor.map(decompress_one, gz_files):
            if status == "success":
                success += 1
            else:
                failed += 1
                print(f"Failed {name}: {status}")
    
    print(f"Decompressed: {success}, Failed: {failed}")
    print(f"Total .pdb files now: {len(list(PDBS_DIR.glob('*.pdb')))}")
    print(f"Remaining .pdb.gz files: {len(list(PDBS_DIR.glob('*.pdb.gz')))}")
    
    print("\n" + "=" * 80)
    print("Step 2: Process Q9Y233 structure")
    print("=" * 80)
    
    # Load dataset sequence for Q9Y233
    df = pd.read_csv(DATA_DIR / 'kcat_filtered_with_pdb.csv')
    q9y_rows = df[df['uniprot'] == 'Q9Y233']
    dataset_seq = q9y_rows['sequence'].iloc[0]
    print(f"Q9Y233 dataset sequence length: {len(dataset_seq)}")
    print(f"Q9Y233 dataset rows: {len(q9y_rows)}")
    
    # Parse PDB
    q9y_pdb = PDBS_DIR / 'AF-Q9Y233-F1-model_v6.pdb'
    atoms = parse_pdb_atoms(q9y_pdb)
    pdb_seq, residues = get_pdb_sequence_with_residues(atoms)
    print(f"Q9Y233 PDB sequence length: {len(pdb_seq)}")
    
    # Find best matching window of dataset sequence in PDB sequence
    pos, matches = find_best_window_position(pdb_seq, dataset_seq)
    if pos is None:
        print("ERROR: Dataset sequence not found in PDB sequence!")
        return
    
    identity = matches / len(dataset_seq)
    print(f"Best match: position {pos}, matches {matches}/{len(dataset_seq)} ({identity*100:.1f}%)")
    
    start_res_idx = pos
    end_res_idx = pos + len(dataset_seq) - 1
    start_res = residues[start_res_idx]
    end_res = residues[end_res_idx]
    
    print(f"Dataset sequence found at PDB residue index {start_res_idx} to {end_res_idx}")
    print(f"Corresponding PDB residue range: chain {start_res['chain_id']}, "
          f"residue {start_res['res_seq']}{start_res['insertion']} to "
          f"{end_res['res_seq']}{end_res['insertion']}")
    
    # Extract matching region
    output_pdb = PDBS_DIR / 'AF-Q9Y233-F1-model_v6_truncated.pdb'
    extract_residues_to_pdb(q9y_pdb, output_pdb, start_res_idx, end_res_idx, residues)
    print(f"Extracted structure saved to: {output_pdb}")
    
    # Verify extracted sequence
    extracted_atoms = parse_pdb_atoms(output_pdb)
    extracted_seq, _ = get_pdb_sequence_with_residues(extracted_atoms)
    print(f"Extracted PDB sequence length: {len(extracted_seq)}")
    print(f"Extracted sequence matches dataset: {extracted_seq == dataset_seq}")
    
    # Save mapping info
    mapping_info = {
        'uniprot': 'Q9Y233',
        'dataset_seq_len': len(dataset_seq),
        'pdb_seq_len': len(pdb_seq),
        'dataset_start_idx_in_pdb': start_res_idx,
        'dataset_end_idx_in_pdb': end_res_idx,
        'pdb_start_residue': {
            'chain': start_res['chain_id'],
            'res_seq': start_res['res_seq'],
            'insertion': start_res['insertion'],
            'res_name': start_res['res_name'],
        },
        'pdb_end_residue': {
            'chain': end_res['chain_id'],
            'res_seq': end_res['res_seq'],
            'insertion': end_res['insertion'],
            'res_name': end_res['res_name'],
        },
        'extracted_pdb': str(output_pdb),
    }
    import json
    mapping_json = DATA_DIR / 'q9y233_truncation_mapping.json'
    with open(mapping_json, 'w') as f:
        json.dump(mapping_info, f, indent=2)
    print(f"Saved mapping info to: {mapping_json}")
    
    print("\n" + "=" * 80)
    print("All steps completed!")
    print("=" * 80)


if __name__ == '__main__':
    import pandas as pd
    main()
