#!/usr/bin/env python3
"""Add a struct_label column to a residue-level task CSV.

struct_label[i] is the residue-level label for the i-th residue of the PDB
structure file (residues enumerated in ATOM-record file order, keyed by
(chain, resseq, icode) at first occurrence). The PDB residue sequence is
globally aligned to the CSV sequence column; each structure residue inherits
the label of the aligned sequence position. Structure residues that do not
align to any sequence position (e.g. expression tags present only in the
structure) are labeled -1 (ignore index).

Usage:
    python add_struct_label.py --csv splits/ssp_train.csv --pdb-dir pdbs \
        [--seq-col aa_seq] [--label-col label] [--struct-col struct_file] \
        [--out-col struct_label] [--workers 32] [--dry-run]
"""
import argparse
import os
import sys
from multiprocessing import Pool

import pandas as pd
from Bio.Align import PairwiseAligner

RES3 = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G',
    'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V',
    'TRP': 'W', 'TYR': 'Y', 'SEC': 'U', 'PYL': 'O',
    # common modified amino acids (kept in sync with
    # ligand_binding_site/scripts/extract_binding_sites_core.py)
    'MSE': 'M', 'SEP': 'S', 'TPO': 'T', 'PTR': 'Y', 'KCX': 'K', 'MLY': 'K',
    'M3L': 'K', 'MLZ': 'K', 'HLY': 'K', 'FME': 'M', 'CME': 'C', 'CSD': 'C',
    'CSW': 'C', 'OCY': 'C', 'CSO': 'C', 'HYP': 'P', 'NLN': 'L', 'NLE': 'L',
    'ALC': 'A', 'SAC': 'S', 'PCA': 'E', 'CXM': 'M', 'MVA': 'V', 'ORN': 'R',
    'PFF': 'F', 'TPQ': 'Y', 'TYS': 'Y', 'SME': 'M', 'SNC': 'C', 'YCM': 'C',
    'MEN': 'N', 'DM0': 'K', 'DMH': 'K', 'LLP': 'K',
    # additional D-/modified residues (from lbs 13_verify_seq_structure_consistency.py)
    'AIB': 'A', 'NRQ': 'N', 'DAL': 'A', 'DLE': 'L', 'DVA': 'V', 'DSE': 'S',
    'DTH': 'T', 'DSN': 'N', 'DGL': 'E', 'DLY': 'K', 'DAR': 'R', 'DGN': 'Q',
    'DTR': 'W', 'DHI': 'H', 'DIL': 'I', 'DPH': 'F', 'DPR': 'P', 'DTY': 'Y',
    'DCY': 'C', 'DAS': 'D',
}

IGNORE = '-1'


def pdb_residue_seq(path):
    """Return residue one-letter sequence of a PDB file in ATOM file order.

    Residues are keyed by (chain, resseq, icode) at first occurrence so that
    residues lacking a CA atom are still included.
    """
    seen = set()
    seq = []
    with open(path) as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            key = (line[21], line[22:26], line[26])
            if key in seen:
                continue
            seen.add(key)
            seq.append(RES3.get(line[17:20].strip(), 'X'))
    return ''.join(seq)


def _make_aligner():
    a = PairwiseAligner()
    a.mode = 'global'
    a.match_score = 2.0
    a.mismatch_score = -1.0
    a.open_gap_score = -2.0
    a.extend_gap_score = -0.5
    return a


_ALIGNER = None


def _get_aligner():
    global _ALIGNER
    if _ALIGNER is None:
        _ALIGNER = _make_aligner()
    return _ALIGNER


def map_row(args):
    """Return (struct_label_str, n_res, n_ignore, n_nonx_mismatch) for one row."""
    uid, pdb_path, seq, labels = args
    pseq = pdb_residue_seq(pdb_path)
    n_res = len(pseq)
    if pseq == seq:
        return ' '.join(labels), n_res, 0, 0
    if len(pseq) == len(seq):
        # equal length: positional mapping (substitutions stay aligned)
        n_mm = sum(1 for a, b in zip(pseq, seq) if a != 'X' and a != b)
        return ' '.join(labels), n_res, 0, n_mm
    aln = _get_aligner().align(seq, pseq)[0]
    pos_map = {}
    a_blocks, b_blocks = aln.aligned
    for (a_start, a_end), (b_start, b_end) in zip(a_blocks, b_blocks):
        for k in range(int(a_end) - int(a_start)):
            pos_map[int(b_start) + k] = int(a_start) + k
    out = []
    n_ignore = 0
    n_mm = 0
    for i, pc in enumerate(pseq):
        j = pos_map.get(i)
        if j is None:
            out.append(IGNORE)
            n_ignore += 1
        else:
            out.append(labels[j])
            if pc != 'X' and pc != seq[j]:
                n_mm += 1
    return ' '.join(out), n_res, n_ignore, n_mm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', required=True)
    ap.add_argument('--pdb-dir', required=True)
    ap.add_argument('--seq-col', default='aa_seq')
    ap.add_argument('--label-col', default='label')
    ap.add_argument('--struct-col', default='struct_file')
    ap.add_argument('--out-col', default='struct_label')
    ap.add_argument('--workers', type=int, default=32)
    ap.add_argument('--dry-run', action='store_true',
                    help='compute and report stats only, do not write the csv')
    ap.add_argument('--overwrite', action='store_true',
                    help='recompute and replace an existing out column '
                         '(e.g. after pdb files were fixed)')
    args = ap.parse_args()

    df = pd.read_csv(args.csv, dtype=str, keep_default_na=False)
    if args.out_col in df.columns and not args.overwrite:
        sys.exit(f'column {args.out_col!r} already exists in {args.csv} '
                 f'(pass --overwrite to recompute it)')

    jobs = []
    for _, r in df.iterrows():
        pdb_path = os.path.join(args.pdb_dir, r[args.struct_col])
        labels = r[args.label_col].strip().lstrip('[').rstrip(']').split()
        if len(labels) != len(r[args.seq_col]):
            sys.exit(f'{r.name}: len({args.label_col}) != len({args.seq_col}) '
                     f'for {r.get("unique_id", "?")}')
        jobs.append((r.get('unique_id', ''), pdb_path, r[args.seq_col], labels))

    with Pool(args.workers) as pool:
        results = pool.map(map_row, jobs, chunksize=64)

    struct_labels = []
    n_rows_gap = n_rows_mm = n_ignore_tot = 0
    for (sl, n_res, n_ign, n_mm), (uid, _, seq, _) in zip(results, jobs):
        toks = sl.split(' ')
        assert len(toks) == n_res, f'{uid}: struct_label length mismatch'
        struct_labels.append('[ ' + sl + ' ]')
        if n_ign:
            n_rows_gap += 1
            n_ignore_tot += n_ign
        if n_mm:
            n_rows_mm += 1

    print(f'{os.path.basename(args.csv)}: rows={len(df)} '
          f'rows_with_-1={n_rows_gap} total_-1={n_ignore_tot} '
          f'rows_with_nonX_mismatch={n_rows_mm}')

    if not args.dry_run:
        df[args.out_col] = struct_labels
        df.to_csv(args.csv, index=False, lineterminator='\n')
        print(f'  wrote {args.out_col} -> {args.csv}')


if __name__ == '__main__':
    main()
