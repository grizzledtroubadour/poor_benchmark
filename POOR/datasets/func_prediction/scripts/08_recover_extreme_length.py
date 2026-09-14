#!/usr/bin/env python3
"""
Step 1: Recover extreme-length chains that were filtered out during stage 2.3.1.

These chains have actual residue counts < 60 (ExtremeShort) or > 1000 (ExtremeLong),
and were removed from the final nonredundant set. This script recovers them from
the original SIFTS annotated_chains.csv + raw PDB files.

Output: output/extreme_length_chains.csv
"""

import gzip
import json
import os
import re
from collections import OrderedDict
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parent.parent  # func_prediction/
POOR_ROOT = BASE.parent.parent  # POOR/
DATA_SIFTS = POOR_ROOT / 'data' / 'sifts'
OUTPUT_DIR = BASE / 'output'
PDB_RAW_DIR = POOR_ROOT / 'data' / 'pdb'

AA_MAP = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
    'MSE': 'M', 'SEC': 'C', 'PYL': 'K',
}

MONTH_MAP = {
    'JAN': 1, 'FEB': 2, 'MAR': 3, 'APR': 4, 'MAY': 5, 'JUN': 6,
    'JUL': 7, 'AUG': 8, 'SEP': 9, 'OCT': 10, 'NOV': 11, 'DEC': 12,
}


def extract_chain_seq_from_pdb(filepath, target_chain_id):
    """Extract amino acid sequence and residue count from a specific chain in a PDB.gz file."""
    residues = OrderedDict()
    opener = gzip.open if filepath.suffix == '.gz' else open
    try:
        with opener(filepath, 'rt') as f:
            for line in f:
                if not (line.startswith('ATOM  ') or line.startswith('HETATM')):
                    continue
                chain_id = line[21:22].strip()
                if chain_id != target_chain_id:
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
    return {'count': len(residues), 'seq': seq, 'seq_len': len(seq)}, None


def extract_chain_seq_from_cif(filepath, target_chain_id):
    """Extract amino acid sequence from a specific chain in a mmCIF.gz file."""
    residues = OrderedDict()
    in_atom_site = False
    col_map = {}
    opener = gzip.open if filepath.suffix == '.gz' else open
    try:
        with opener(filepath, 'rt') as f:
            for line in f:
                line = line.rstrip('\n')
                if line.startswith('_atom_site.'):
                    if not in_atom_site:
                        in_atom_site = True
                    field = line.split('.', 1)[1].strip()
                    col_map[field] = len(col_map)
                    continue
                if not in_atom_site:
                    continue
                if line.startswith('#') or line.startswith('loop_'):
                    in_atom_site = False
                    col_map = {}
                    continue
                if not line.strip() or line.startswith('_'):
                    continue
                parts = line.split()
                if len(parts) < 5:
                    continue
                group = parts[0] if parts[0] in ('ATOM', 'HETATM') else None
                if not group:
                    continue
                # auth_atom_site uses auth_asym_id for chain, label uses label_asym_id
                # Try auth_asym_id first (matches PDB chain ID)
                chain_col = col_map.get('auth_asym_id')
                comp_col = col_map.get('label_comp_id')
                seq_col = col_map.get('auth_seq_id')
                icode_col = col_map.get('pdbx_PDB_ins_code')
                if chain_col is not None and parts[chain_col] != target_chain_id:
                    continue
                if comp_col is None or seq_col is None:
                    continue
                res_name = parts[comp_col].upper() if comp_col < len(parts) else ''
                if res_name not in AA_MAP:
                    continue
                seq_id = parts[seq_col] if seq_col < len(parts) else ''
                icode = parts[icode_col].strip() if icode_col and icode_col < len(parts) else ''
                key = (seq_id, icode)
                if key not in residues:
                    residues[key] = res_name
    except Exception as e:
        return None, str(e)
    seq = ''.join(AA_MAP.get(r, '') for r in residues.values())
    return {'count': len(residues), 'seq': seq, 'seq_len': len(seq)}, None


def parse_deposition_date(filepath):
    """Parse deposition date from PDB HEADER line."""
    opener = gzip.open if filepath.suffix == '.gz' else open
    try:
        with opener(filepath, 'rt') as f:
            for line in f:
                if line.startswith('HEADER'):
                    # Format: HEADER    COMPOUND    DD-MMM-YY   IDCODE
                    # Date is at fixed position: columns 50-59
                    date_str = line[50:59].strip()
                    if date_str:
                        parts = date_str.split('-')
                        if len(parts) == 3:
                            day, month, year = parts
                            month = MONTH_MAP.get(month.upper())
                            if month:
                                if len(year) == 2:
                                    year = '20' + year if int(year) < 50 else '19' + year
                                return f'{year}-{month:02d}-{int(day):02d}'
                    return None
    except Exception:
        pass
    return None


def process_chain(args):
    """Process a single chain candidate."""
    pdb_id, chain_id, uniprot, resolution, mapped_length, has_ec, has_go = args
    chain_key = f'{pdb_id}_{chain_id}'

    # Find PDB file
    pdb_path = PDB_RAW_DIR / f'{pdb_id}.pdb.gz'
    is_cif = False
    if not pdb_path.exists():
        pdb_path = PDB_RAW_DIR / f'{pdb_id}.cif.gz'
        is_cif = True
    if not pdb_path.exists():
        pdb_path = PDB_RAW_DIR / f'{pdb_id}.pdb'
        if not pdb_path.exists():
            pdb_path = PDB_RAW_DIR / f'{pdb_id}.cif'
            is_cif = True
    if not pdb_path.exists():
        return {'chain_key': chain_key, 'error': f'PDB file not found for {pdb_id}'}

    # Extract sequence
    if pdb_path.suffix == '.cif' or (pdb_path.suffix == '.gz' and pdb_path.stem.endswith('.cif')):
        result, err = extract_chain_seq_from_cif(pdb_path, chain_id)
    else:
        result, err = extract_chain_seq_from_pdb(pdb_path, chain_id)

    if result is None:
        return {'chain_key': chain_key, 'error': err}
    if result['seq_len'] == 0:
        return {'chain_key': chain_key, 'error': 'empty sequence'}

    actual_len = result['seq_len']

    # Length filter
    if actual_len > 2000:
        return {'chain_key': chain_key, 'error': f'actual_len={actual_len} > 2000, discarded'}
    if 60 <= actual_len <= 1000:
        return {'chain_key': chain_key, 'error': f'actual_len={actual_len} in normal range, skip'}

    extreme_type = 'ExtremeShort' if actual_len < 60 else 'ExtremeLong'

    # Parse deposition date
    dep_date = parse_deposition_date(pdb_path)

    return {
        'chain_key': chain_key,
        'uniprot': uniprot,
        'pdb_id': pdb_id,
        'chain_id': chain_id,
        'aa_seq': result['seq'],
        'actual_len': actual_len,
        'deposition_date': dep_date,
        'extreme_type': extreme_type,
        'mapped_length': mapped_length,
        'resolution': resolution,
        'has_ec': has_ec,
        'has_go': has_go,
        'pdb_file': pdb_path.name,
        'is_cif': is_cif,
    }


def main():
    print("=" * 60)
    print("Step 1: Recover extreme-length chains")
    print("=" * 60)

    # Load final nonredundant set
    final_df = pd.read_csv(OUTPUT_DIR / 'final_nonredundant_representatives.csv')
    final_uniprots = set(final_df['uniprot'])
    print(f"Final nonredundant UniProts: {len(final_uniprots)}")

    # Load annotated chains
    ac = pd.read_csv(DATA_SIFTS / 'annotated_chains.csv', low_memory=False)
    print(f"Annotated chains: {len(ac)} chains, {ac['uniprot'].nunique()} UniProts")

    # Filter: has EC or GO, UniProt NOT in final
    candidates = ac[
        (ac['has_ec'] | ac['has_go']) &
        (~ac['uniprot'].isin(final_uniprots))
    ].copy()
    print(f"Candidates (has EC/GO, not in final): {len(candidates)} chains, {candidates['uniprot'].nunique()} UniProts")

    # UniProt dedup: keep best chain per UniProt
    # Sort: has_ec desc, resolution asc, mapped_length desc, pdb_id+chain_id asc
    candidates['_has_ec_sort'] = candidates['has_ec'].astype(int)
    candidates = candidates.sort_values(
        ['_has_ec_sort', 'resolution', 'mapped_length', 'pdb_id', 'chain_id'],
        ascending=[False, True, False, True, True]
    ).drop_duplicates(subset='uniprot', keep='first')
    print(f"After UniProt dedup: {len(candidates)} chains")

    # Prepare tasks for parallel processing
    tasks = []
    for _, row in candidates.iterrows():
        tasks.append((
            row['pdb_id'], row['chain_id'], row['uniprot'],
            row['resolution'], row['mapped_length'],
            row['has_ec'], row['has_go']
        ))

    print(f"\nExtracting sequences from PDB files ({len(tasks)} chains)...")
    results = []
    errors = []
    with ProcessPoolExecutor(max_workers=16) as executor:
        futures = {executor.submit(process_chain, t): t for t in tasks}
        for i, future in enumerate(as_completed(futures)):
            if i % 500 == 0:
                print(f"  Processed {i}/{len(tasks)}...")
            result = future.result()
            if 'error' in result:
                errors.append(result)
            else:
                results.append(result)

    print(f"\nResults: {len(results)} extreme-length chains found")
    print(f"Errors/skipped: {len(errors)}")

    # Separate by type
    short = [r for r in results if r['extreme_type'] == 'ExtremeShort']
    long = [r for r in results if r['extreme_type'] == 'ExtremeLong']
    print(f"  ExtremeShort (< 60 residues): {len(short)}")
    print(f"  ExtremeLong (1000 < residues <= 2000): {len(long)}")

    # Load EC/GO annotations
    print("\nLoading EC annotations...")
    with open(DATA_SIFTS / 'sifts_ec_annotations.json') as f:
        ec_data = json.load(f)

    print("Loading GO annotations...")
    with open(DATA_SIFTS / 'sifts_go_annotations.json') as f:
        go_data = json.load(f)

    print("Loading GO namespace mapping...")
    with open(OUTPUT_DIR / 'go_id_to_namespace_full.json') as f:
        go2ns = json.load(f)

    print("Loading UniProt metadata (for reviewed status)...")
    with open(OUTPUT_DIR / 'uniprot_metadata.json') as f:
        uniprot_meta = json.load(f)

    # Add EC/GO labels to results
    EXP_EVIDENCE = {'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'TAS', 'IC', 'HTP', 'HDA'}
    for r in results:
        ck = r['chain_key']
        uni = r['uniprot']

        # EC labels
        ec_labels = []
        if ck in ec_data:
            for ann in ec_data[ck].get('annotations', []):
                ec = ann.get('ec_number', '')
                if ec:
                    # Normalize to 4-part format
                    parts = ec.split('.')
                    while len(parts) < 4:
                        parts.append('-')
                    ec_labels.append('.'.join(parts))
        r['ec_labels'] = ';'.join(ec_labels) if ec_labels else ''

        # GO labels by namespace
        go_mf, go_cc, go_bp = [], [], []
        has_exp_go = False
        if ck in go_data:
            for ann in go_data[ck].get('annotations', []):
                go_id = ann.get('go_id', '')
                evidence = ann.get('evidence', '')
                if not go_id:
                    continue
                ns = go2ns.get(go_id, '')
                if evidence in EXP_EVIDENCE:
                    has_exp_go = True
                if ns == 'molecular_function':
                    if go_id not in go_mf:
                        go_mf.append(go_id)
                elif ns == 'cellular_component':
                    if go_id not in go_cc:
                        go_cc.append(go_id)
                elif ns == 'biological_process':
                    if go_id not in go_bp:
                        go_bp.append(go_id)
        r['go_mf_labels'] = ';'.join(go_mf) if go_mf else ''
        r['go_cc_labels'] = ';'.join(go_cc) if go_cc else ''
        r['go_bp_labels'] = ';'.join(go_bp) if go_bp else ''
        r['has_exp_go'] = has_exp_go

        # Reviewed status
        meta = uniprot_meta.get(uni, {})
        r['is_reviewed'] = meta.get('reviewed', False)

    # Save to CSV
    df = pd.DataFrame(results)
    output_file = OUTPUT_DIR / 'extreme_length_chains.csv'
    df.to_csv(output_file, index=False)
    print(f"\nSaved to {output_file}")
    print(f"Total extreme-length chains: {len(df)}")

    # Summary
    print(f"\n{'='*60}")
    print("Summary")
    print(f"{'='*60}")
    print(f"ExtremeShort: {len(df[df['extreme_type']=='ExtremeShort'])}")
    print(f"ExtremeLong:  {len(df[df['extreme_type']=='ExtremeLong'])}")
    print(f"  With EC:    {len(df[df['ec_labels']!=''])}")
    print(f"  With GO-MF: {len(df[df['go_mf_labels']!=''])}")
    print(f"  With GO-CC: {len(df[df['go_cc_labels']!=''])}")
    print(f"  With GO-BP: {len(df[df['go_bp_labels']!=''])}")
    print(f"  Reviewed:   {len(df[df['is_reviewed']==True])}")
    print(f"  Has exp GO: {len(df[df['has_exp_go']==True])}")


if __name__ == '__main__':
    main()
