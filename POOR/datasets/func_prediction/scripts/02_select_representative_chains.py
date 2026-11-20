#!/usr/bin/env python3
"""UniProt-level representative-chain selection (reimplements DATA_PROCESS.md 2.2).

Group by UniProt ID and retain the best PDB chain per group using these priorities:
  1. Extraction tier: A > B (tier C is not included as a fallback).
  2. Source quality: Swiss-Prot > TrEMBL (constant within each group, so not discriminative).
  3. EC availability: annotated > unannotated.
  4. Number of experimental GO terms: descending.
  5. Total GO terms: descending.
  6. Resolution: ascending; missing last.
  7. R-free: ascending; missing last.
  8. Mapped length: descending.
  9. PDB ID + chain ID: ascending lexicographic deterministic tiebreak.

Build the candidate pool from current outputs:
  - data/sifts/annotated_chains.csv: initial metadata (mapped_length in [60,1000] and EC/GO).
  - data/sifts/sifts_uniprot_mapping.json: SIFTS segments determine tier A/B; valid pdb_beg/pdb_end implies A.
  - data/sifts/sifts_go_annotations.json: experimental GO counts using ijson and ten evidence codes.
  - output/uniprot_metadata.json: reviewed status (Swiss-Prot > TrEMBL).

Output: output/representative_chains.csv, consumed by 03_extract_representative_structures.py;
the historical run selected 55,124 representative chains.

Known deviation (see DATA_PROCESS.md notes):
  - Historically, B-tier candidates first passed inferred-boundary identity >=0.90 (419,370/439,427 passed).
    data/intermediate/b_boundary_inference.jsonl was not migrated. This script includes all B-tier candidates by default.
    After rerunning inference, supply --b-keep-list with passing chain_keys (one per line or a CSV chain_key column)
    to reproduce the historical filter.
"""

import argparse
import csv
import json
from pathlib import Path

import ijson
import pandas as pd

TASK_DIR = Path(__file__).resolve().parent.parent  # datasets/func_prediction/
POOR_ROOT = TASK_DIR.parent.parent                 # POOR/
SIFTS_DIR = POOR_ROOT / 'data' / 'sifts'
OUT_DIR = TASK_DIR / 'output'

EXP_EVIDENCE = {'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'TAS', 'IC', 'HTP', 'HDA'}


def load_filtered_chains():
    """Initial filter: mapped_length in [60,1000] and EC or GO annotations (document section 1.5)."""
    df = pd.read_csv(SIFTS_DIR / 'annotated_chains.csv', low_memory=False)
    df = df[(df['has_ec'] | df['has_go']) &
            (df['mapped_length'] >= 60) & (df['mapped_length'] <= 1000)].copy()
    # Exclude missing chain_id/pdb_id rows, whose structures cannot be located.
    df = df[df['pdb_id'].notna() & df['chain_id'].notna()].copy()
    df['chain_key'] = df['pdb_id'] + '_' + df['chain_id']
    return df


def classify_extract_mode(chain_keys):
    """Tier A requires a SIFTS segment with valid pdb_beg and pdb_end; otherwise tier B."""
    with open(SIFTS_DIR / 'sifts_uniprot_mapping.json') as f:
        mapping = json.load(f)
    modes = {}
    for ck in chain_keys:
        info = mapping.get(ck)
        mode = 'B'
        if info:
            for seg in info.get('segments', []):
                if seg.get('pdb_beg') and seg.get('pdb_end'):
                    mode = 'A'
                    break
        modes[ck] = mode
    return modes


def count_go_evidence(chain_keys):
    """Stream SIFTS GO annotations; count experimental and total GO terms per chain."""
    go_count, exp_go_count = {}, {}
    with open(SIFTS_DIR / 'sifts_go_annotations.json', 'rb') as f:
        for chain_key, record in ijson.kvitems(f, ''):
            if chain_key not in chain_keys:
                continue
            seen, seen_exp = set(), set()
            for ann in record.get('annotations', []):
                go_id = ann.get('go_id', '')
                if not go_id:
                    continue
                seen.add(go_id)
                if ann.get('evidence', '') in EXP_EVIDENCE:
                    seen_exp.add(go_id)
            go_count[chain_key] = len(seen)
            exp_go_count[chain_key] = len(seen_exp)
    return go_count, exp_go_count


def load_b_keep_list(path):
    """Read passing B-tier chain_keys, one per line or from a CSV chain_key column."""
    with open(path) as f:
        first = f.readline().strip()
    if first == 'chain_key' or 'chain_key' in first.split(','):
        return set(pd.read_csv(path)['chain_key'])
    with open(path) as f:
        return {line.strip() for line in f if line.strip()}


def main():
    ap = argparse.ArgumentParser(description='UniProt-level representative chain selection')
    ap.add_argument('--out', type=Path, default=OUT_DIR / 'representative_chains.csv',
                    help='Output CSV (default: output/representative_chains.csv)')
    ap.add_argument('--b-keep-list', type=Path, default=None,
                    help='Optional chain_keys passing B-tier boundary identity >=0.90, to reproduce historical filtering')
    args = ap.parse_args()

    print('[1/4] Loading filtered candidate chains...')
    df = load_filtered_chains()
    print(f'  candidates: {len(df):,} chains, {df["uniprot"].nunique():,} UniProts')

    print('[2/4] Classifying extract mode (A/B) from SIFTS segments...')
    chain_keys = set(df['chain_key'])
    df['extract_mode'] = df['chain_key'].map(classify_extract_mode(chain_keys))
    if args.b_keep_list:
        keep = load_b_keep_list(args.b_keep_list)
        before = len(df)
        df = df[(df['extract_mode'] == 'A') | (df['chain_key'].isin(keep))].copy()
        print(f'  B-keep filter: {before:,} -> {len(df):,} chains')
    print(f'  A: {(df["extract_mode"] == "A").sum():,}, B: {(df["extract_mode"] == "B").sum():,}')

    print('[3/4] Counting GO evidence (streaming) + reviewed status...')
    go_count, exp_go_count = count_go_evidence(chain_keys)
    df['go_count'] = df['chain_key'].map(go_count).fillna(0).astype(int)
    df['exp_go_count'] = df['chain_key'].map(exp_go_count).fillna(0).astype(int)
    with open(OUT_DIR / 'uniprot_metadata.json') as f:
        uniprot_meta = json.load(f)
    df['is_reviewed'] = df['uniprot'].map(
        lambda u: bool(uniprot_meta.get(u, {}).get('reviewed', False)))

    print('[4/4] Selecting one representative chain per UniProt...')
    df['_mode_pri'] = df['extract_mode'].map({'A': 0, 'B': 1})
    df['_has_ec'] = df['has_ec'].astype(int)
    df['_reviewed'] = df['is_reviewed'].astype(int)
    df = df.sort_values(
        ['uniprot', '_mode_pri', '_reviewed', '_has_ec', 'exp_go_count', 'go_count',
         'resolution', 'r_free', 'mapped_length', 'pdb_id', 'chain_id'],
        ascending=[True, True, False, False, False, False,
                   True, True, False, True, True],
        na_position='last',
        kind='mergesort',  # stable sort for deterministic tiebreak
    )
    rep = df.drop_duplicates(subset='uniprot', keep='first').copy()

    out_cols = ['uniprot', 'pdb_id', 'chain_id', 'chain_key', 'extract_mode',
                'is_reviewed', 'ec_count', 'exp_go_count', 'go_count',
                'resolution', 'r_free', 'mapped_length', 'method']
    rep = rep[out_cols].sort_values('uniprot').reset_index(drop=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    rep.to_csv(args.out, index=False)

    print(f'\nRepresentative chains: {len(rep):,} '
          f'(A: {(rep["extract_mode"] == "A").sum():,}, B: {(rep["extract_mode"] == "B").sum():,})')
    print(f'  Swiss-Prot: {rep["is_reviewed"].sum():,}, with EC: {(rep["ec_count"] > 0).sum():,}, '
          f'with exp GO: {(rep["exp_go_count"] > 0).sum():,}')
    print(f'Saved to {args.out}')


if __name__ == '__main__':
    main()
