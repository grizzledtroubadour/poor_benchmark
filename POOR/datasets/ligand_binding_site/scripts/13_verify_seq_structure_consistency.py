#!/usr/bin/env python3
"""Verify that aa_seq in CSV matches the sequence extracted from pdbs/*.pdb files."""
import argparse
import logging
import sys
from collections import OrderedDict
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from extract_binding_sites_core import MODIFIED_AA_1LETTER  # noqa: E402

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s')
logger = logging.getLogger(__name__)

RES3 = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G',
    'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V',
    'TRP': 'W', 'TYR': 'Y', 'SEC': 'U', 'PYL': 'O',
    # common modified amino acids
    'MSE': 'M', 'SEP': 'S', 'TPO': 'T', 'PTR': 'Y', 'KCX': 'K', 'CSD': 'C',
    'CME': 'C', 'CSO': 'C', 'LLP': 'K', 'MLE': 'L', 'FME': 'M', 'HYP': 'P',
    'SAC': 'S', 'AIB': 'A', 'NRQ': 'N', 'DAL': 'A', 'DLE': 'L', 'DVA': 'V',
    'DSE': 'S', 'DTH': 'T', 'DSN': 'N', 'DGL': 'E', 'DLY': 'K', 'DAR': 'R',
    'DGN': 'Q', 'DTR': 'W', 'DHI': 'H', 'DIL': 'I', 'DPH': 'F', 'DPR': 'P',
    'DTY': 'Y', 'DCY': 'C', 'DAS': 'D',
}
RES3.update(MODIFIED_AA_1LETTER)


def seq_from_pdb(path):
    residues = OrderedDict()
    with open(path) as f:
        for line in f:
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                continue
            if len(line) < 27:
                continue
            resname = line[17:20].strip()
            try:
                resseq = int(line[22:26])
            except ValueError:
                continue
            icode = line[26]
            key = (resseq, icode)
            if key not in residues:
                residues[key] = resname
    return ''.join(RES3.get(r, 'X') for r in residues.values())


def main():
    parser = argparse.ArgumentParser()
    task_dir = Path(__file__).resolve().parent.parent  # datasets/ligand_binding_site
    parser.add_argument('--train', type=str, default=str(task_dir / 'splits' / 'ligand_binding_site_train.csv'))
    parser.add_argument('--val', type=str, default=str(task_dir / 'splits' / 'ligand_binding_site_val.csv'))
    parser.add_argument('--test', type=str, default=str(task_dir / 'splits' / 'ligand_binding_site_test.csv'))
    parser.add_argument('--pdb_dir', type=str, default=str(task_dir / 'pdbs'))
    parser.add_argument('--out', type=str, default=str(task_dir / 'output' / 'seq_structure_mismatches.csv'))
    args = parser.parse_args()

    dfs = [pd.read_csv(p, keep_default_na=False, na_values=['']) for p in [args.train, args.val, args.test]]
    df = pd.concat(dfs, ignore_index=True)
    pdb_dir = Path(args.pdb_dir)

    exact = 0
    len_only = 0
    mismatches = []
    empty = 0
    missing_file = 0

    for _, row in df.iterrows():
        uid = row['unique_id']
        aa = row['aa_seq']
        pdb_path = pdb_dir / f"{uid}.pdb"
        if not pdb_path.exists():
            missing_file += 1
            mismatches.append((uid, 'missing_file', len(aa), 0, aa[:50], ''))
            continue
        pdb_seq = seq_from_pdb(pdb_path)
        if len(pdb_seq) == 0:
            empty += 1
            mismatches.append((uid, 'empty_structure', len(aa), 0, aa[:50], ''))
            continue
        if pdb_seq == aa:
            exact += 1
        else:
            if len(pdb_seq) == len(aa):
                category = 'same_length_diff_seq'
            elif abs(len(pdb_seq) - len(aa)) <= 5:
                category = 'small_len_diff'
            else:
                category = 'large_len_diff'
            mismatches.append((uid, category, len(aa), len(pdb_seq), aa[:50], pdb_seq[:50]))
            len_only += 1

    total = len(df)
    logger.info(f'Total rows: {total}')
    logger.info(f'Exact sequence match: {exact} / {total} ({100*exact/total:.2f}%)')
    logger.info(f'Mismatches: {len_only} / {total} ({100*len_only/total:.2f}%)')
    logger.info(f'Empty structures: {empty}')
    logger.info(f'Missing pdb files: {missing_file}')

    # category breakdown
    cat_counts = {}
    for m in mismatches:
        cat_counts[m[1]] = cat_counts.get(m[1], 0) + 1
    logger.info('Mismatch categories: ' + ', '.join(f'{k}={v}' for k, v in cat_counts.items()))

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(mismatches, columns=['unique_id', 'category', 'aa_seq_len', 'pdb_seq_len', 'aa_seq_head', 'pdb_seq_head']).to_csv(out_path, index=False)
    logger.info(f'Saved mismatch details to {out_path}')


if __name__ == '__main__':
    main()
