"""
Extract amino acid sequences from PDB/CIF files and count actual residues.
Filter out chains with < 60 residues.

从任务根目录（datasets/func_prediction/）运行；产物写入 output/。
"""
import os
import sys
import json
from collections import OrderedDict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

PDB_DIR = 'pdbs'
OUTPUT_FASTA = 'output/filtered_sequences.fa'
OUTPUT_JSON = 'output/chain_residue_counts.json'
MIN_RESIDUES = 60

os.makedirs('output', exist_ok=True)

# 3-letter to 1-letter amino acid code
AA_MAP = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
    'MSE': 'M',  # Selenomethionine
    'SEC': 'C',  # Selenocysteine
    'PYL': 'K',  # Pyrrolysine
    'UNK': 'X',  # Unknown
}


def extract_from_pdb(filepath):
    """Extract sequence and count residues from a PDB file."""
    residues = OrderedDict()  # (res_seq, icode) -> res_name
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if not (line.startswith('ATOM  ') or line.startswith('HETATM')):
                    continue
                res_name = line[17:20].strip()
                res_seq = line[22:26].strip()
                icode = line[26:27].strip()
                key = (res_seq, icode)
                if key not in residues:
                    residues[key] = res_name
    except Exception as e:
        return None, str(e)
    
    seq = ''.join(AA_MAP.get(r, '') for r in residues.values())
    return {'count': len(residues), 'seq': seq, 'seq_len': len(seq)}


def extract_from_cif(filepath):
    """Extract sequence and count residues from a mmCIF file."""
    residues = OrderedDict()
    in_atom_site = False
    
    # Column indices
    group_col = None
    comp_col = None
    seq_col = None
    icode_col = None
    
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.rstrip('\n')
                
                if line.startswith('_atom_site.'):
                    in_atom_site = True
                    if 'group_PDB' in line:
                        group_col = len([l for l in line[:line.index('group_PDB')].split() if l.startswith('_atom_site.')])
                    elif 'label_comp_id' in line and 'auth' not in line.lower():
                        comp_col = len([l for l in line[:line.index('label_comp_id')].split() if l.startswith('_atom_site.')])
                    elif 'auth_seq_id' in line:
                        seq_col = len([l for l in line[:line.index('auth_seq_id')].split() if l.startswith('_atom_site.')])
                    elif 'pdbx_PDB_ins_code' in line:
                        icode_col = len([l for l in line[:line.index('pdbx_PDB_ins_code')].split() if l.startswith('_atom_site.')])
                    continue
                
                if not in_atom_site:
                    continue
                if line.startswith('#') or line.startswith('loop_'):
                    in_atom_site = False
                    continue
                if not line.strip() or line.startswith('_'):
                    continue
                
                parts = line.split()
                if len(parts) < 5:
                    continue
                
                # Simple heuristic: look for ATOM/HETATM in first few columns
                record_type = parts[0] if parts[0] in ('ATOM', 'HETATM') else None
                if not record_type:
                    continue
                
                # Find res_name - it's typically around column 5 in mmCIF
                # Use a simpler approach: look for standard amino acid names
                res_name = None
                for p in parts[1:8]:
                    if p.upper() in AA_MAP:
                        res_name = p.upper()
                        break
                
                if not res_name:
                    continue
                
                # Find sequence id - usually a numeric or mixed field
                seq_id = None
                for p in parts[2:10]:
                    if p.replace('-', '').replace('?', '').isdigit() or (p.startswith('-') and p[1:].isdigit()):
                        seq_id = p
                        break
                
                if not seq_id:
                    continue
                
                # Find insertion code - usually after seq_id
                icode = ''
                for i, p in enumerate(parts):
                    if p == seq_id and i + 1 < len(parts):
                        next_p = parts[i + 1]
                        if len(next_p) == 1 and next_p.isalpha():
                            icode = next_p
                            break
                
                key = (seq_id, icode)
                if key not in residues:
                    residues[key] = res_name
                    
    except Exception as e:
        return None, str(e)
    
    seq = ''.join(AA_MAP.get(r, '') for r in residues.values())
    return {'count': len(residues), 'seq': seq, 'seq_len': len(seq)}


def process_file(filename):
    """Process a single PDB/CIF file."""
    filepath = os.path.join(PDB_DIR, filename)
    chain_key = filename.replace('.pdb', '').replace('.cif', '')
    
    if filename.endswith('.pdb'):
        result = extract_from_pdb(filepath)
    elif filename.endswith('.cif'):
        result = extract_from_cif(filepath)
    else:
        return None
    
    if result is None or isinstance(result, tuple):
        return {'chain_key': chain_key, 'error': result[1] if isinstance(result, tuple) else 'unknown'}
    
    result['chain_key'] = chain_key
    result['filename'] = filename
    return result


def main():
    files = [f for f in os.listdir(PDB_DIR) if f.endswith(('.pdb', '.cif'))]
    print(f"Total files: {len(files)}")
    
    results = []
    errors = []
    
    with ProcessPoolExecutor(max_workers=32) as executor:
        futures = {executor.submit(process_file, f): f for f in files}
        for i, future in enumerate(as_completed(futures)):
            if i % 5000 == 0:
                print(f"  Processed {i}/{len(files)}...")
            result = future.result()
            if result is None:
                continue
            if 'error' in result:
                errors.append(result)
            else:
                results.append(result)
    
    print(f"\nSuccessful: {len(results)}, Errors: {len(errors)}")
    if errors:
        print(f"First 5 errors: {errors[:5]}")
    
    # Count residues and filter
    counts = {}
    filtered_seqs = []
    short_chains = []
    
    for r in results:
        ck = r['chain_key']
        count = r['count']
        seq = r['seq']
        counts[ck] = {'count': count, 'seq_len': len(seq), 'filename': r['filename']}
        
        if count >= MIN_RESIDUES and len(seq) >= MIN_RESIDUES:
            filtered_seqs.append((ck, seq))
        else:
            short_chains.append((ck, count, len(seq)))
    
    print(f"\nChains with >= {MIN_RESIDUES} residues: {len(filtered_seqs)}")
    print(f"Chains with < {MIN_RESIDUES} residues: {len(short_chains)}")
    
    # Save counts
    with open(OUTPUT_JSON, 'w') as f:
        json.dump(counts, f, indent=2)
    
    # Save filtered sequences
    with open(OUTPUT_FASTA, 'w') as f:
        for ck, seq in filtered_seqs:
            f.write(f">{ck}\n{seq}\n")
    
    print(f"\nSaved to {OUTPUT_FASTA} and {OUTPUT_JSON}")
    
    # Show short chain stats
    if short_chains:
        counts_only = [c for _, c, _ in short_chains]
        seq_counts = [s for _, _, s in short_chains]
        print(f"\nShort chains - residue count: min={min(counts_only)}, max={max(counts_only)}, mean={sum(counts_only)/len(counts_only):.1f}")
        print(f"Short chains - seq length: min={min(seq_counts)}, max={max(seq_counts)}, mean={sum(seq_counts)/len(seq_counts):.1f}")


if __name__ == '__main__':
    main()
