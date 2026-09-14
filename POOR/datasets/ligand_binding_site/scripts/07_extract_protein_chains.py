#!/usr/bin/env python3
"""Fast extraction of protein chain PDBs for the ligand-binding-site dataset.

This version groups rows by PDB ID, parses each structure file only once,
and writes all requested chain-ligand files in one pass.  It uses a manual
PDB parser (fast) and a lightweight mmCIF parser for the ATOM_SITE loop.
"""
import argparse
import gzip
import logging
import os
import re
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s')
logger = logging.getLogger(__name__)

TASK_DIR = Path(__file__).resolve().parent.parent  # datasets/ligand_binding_site
ROOT = TASK_DIR.parent.parent                      # 项目根（POOR/）
PDB_DIR = ROOT / "data" / "pdb"
OUT_DIR = TASK_DIR / "pdbs"

ATOM_LINE = (
    "ATOM  {serial:5d} {name:4s}{altloc:1s}{resname:3s} {chain:1s}"
    "{resseq:4d}{icode:1s}   {x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{bfactor:6.2f}"
    "          {element:>2s}\n"
)

# Protein residue filters, matching extract_binding_sites_core.py
PROTEIN_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL",
    # modified amino acids
    "MSE", "SEP", "TPO", "PTR", "MEN", "FME", "CME", "CSD",
    "CSW", "OCY", "KCX", "MLY", "M3L", "DM0", "DMH", "HYP",
    "NLN", "ALC", "SAC", "PCA", "MLZ", "HLY", "CSO", "LLP",
    "CXM", "MVA", "NLE", "ORN", "PFF", "TPQ", "TYS", "SME",
    "SNC", "YCM",
}


def format_atom(atom):
    name = atom['name']
    if len(name) < 4:
        name = " " + name
    # Standard PDB chain field is a single character; truncate if necessary.
    chain = atom['chain'][0] if atom['chain'] else 'A'
    return ATOM_LINE.format(
        serial=atom['serial'],
        name=name,
        altloc=atom.get('altloc', ' '),
        resname=atom['resname'],
        chain=chain,
        resseq=atom['resseq'],
        icode=atom.get('icode', ' '),
        x=atom['x'],
        y=atom['y'],
        z=atom['z'],
        occupancy=atom.get('occupancy', 1.0),
        bfactor=atom.get('bfactor', 0.0),
        element=atom.get('element', atom['name'][0]),
    )


def parse_pdb_file(path):
    """Parse ATOM records from a PDB file; return dict chain->list of atom dicts."""
    chain_atoms = defaultdict(list)
    serial = 1
    resseq_re = re.compile(r'(\d+)([A-Za-z]?)')
    with gzip.open(path, 'rt', errors='ignore') as fh:
        for line in fh:
            # Multi-model (NMR) entries: keep model 1 only. Without this,
            # all models were concatenated into one chain, producing
            # duplicated residues and atom-serial overflow (>99,999).
            if line.startswith("ENDMDL"):
                break
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            if len(line) < 54:
                continue
            resname = line[17:20].strip()
            if resname not in PROTEIN_RESIDUES:
                continue
            try:
                resseq = int(line[22:26])
                icode = line[26] if len(line) > 26 else ' '
            except ValueError:
                m = resseq_re.search(line[22:27])
                if not m:
                    continue
                resseq = int(m.group(1))
                icode = m.group(2) if m.group(2) else ' '
            try:
                atom = {
                    'serial': serial,
                    'name': line[12:16].strip(),
                    'altloc': line[16] if len(line) > 16 else ' ',
                    'resname': resname,
                    'chain': line[21] if len(line) > 21 else ' ',
                    'resseq': resseq,
                    'icode': icode,
                    'x': float(line[30:38]),
                    'y': float(line[38:46]),
                    'z': float(line[46:54]),
                    'occupancy': float(line[54:60]) if len(line) >= 60 else 1.0,
                    'bfactor': float(line[60:66]) if len(line) >= 66 else 0.0,
                    'element': line[76:78].strip() if len(line) >= 78 else '',
                }
            except Exception:
                continue
            if not atom['name']:
                continue
            if not atom['element']:
                atom['element'] = atom['name'][0]
            chain_atoms[atom['chain']].append(atom)
            serial += 1
    return chain_atoms


def parse_cif_file(path):
    """Parse ATOM records from an mmCIF ATOM_SITE loop."""
    chain_atoms = defaultdict(list)
    serial = 1
    resseq_re = re.compile(r'(\d+)([A-Za-z]?)')
    with gzip.open(path, 'rt', errors='ignore') as fh:
        lines = fh.readlines()

    idx = 0
    n = len(lines)
    while idx < n:
        line = lines[idx].strip()
        if line.startswith('loop_'):
            # collect header
            idx += 1
            cols = []
            while idx < n:
                l = lines[idx].strip()
                if l.startswith('_atom_site.'):
                    cols.append(l.split('.')[1])
                    idx += 1
                else:
                    break
            if not cols:
                continue
            # find column indices we need
            need = {
                'group_PDB': None,
                'id': None,
                'label_atom_id': None,
                'auth_atom_id': None,
                'label_alt_id': None,
                'label_comp_id': None,
                'auth_comp_id': None,
                'auth_asym_id': None,
                'auth_seq_id': None,
                'pdbx_PDB_ins_code': None,
                'Cartn_x': None,
                'Cartn_y': None,
                'Cartn_z': None,
                'occupancy': None,
                'B_iso_or_equiv': None,
                'type_symbol': None,
                'pdbx_PDB_model_num': None,
            }
            for c in need:
                if c in cols:
                    need[c] = cols.index(c)
            if need['group_PDB'] is None or need['auth_asym_id'] is None:
                # skip data until next loop/#
                while idx < n and not lines[idx].strip().startswith('#') and not lines[idx].strip().startswith('loop_'):
                    idx += 1
                if idx < n and lines[idx].strip().startswith('#'):
                    idx += 1
                continue
            # parse data rows
            first_model = None
            while idx < n:
                l = lines[idx].strip()
                if l.startswith('#') or l.startswith('loop_') or l.startswith('data_'):
                    break
                if not l:
                    idx += 1
                    continue
                parts = l.split()
                if len(parts) < len(cols):
                    idx += 1
                    continue
                # Multi-model entries: keep the first model only
                if need['pdbx_PDB_model_num'] is not None:
                    mnum = parts[need['pdbx_PDB_model_num']]
                    if first_model is None:
                        first_model = mnum
                    elif mnum != first_model:
                        idx += 1
                        continue
                resname = parts[need['auth_comp_id']] if need['auth_comp_id'] is not None else parts[need['label_comp_id']]
                if resname not in PROTEIN_RESIDUES:
                    idx += 1
                    continue
                m = resseq_re.search(parts[need['auth_seq_id']])
                if not m:
                    idx += 1
                    continue
                try:
                    resseq = int(m.group(1))
                    icode = m.group(2)
                    if not icode and need['pdbx_PDB_ins_code'] is not None:
                        raw_icode = parts[need['pdbx_PDB_ins_code']]
                        icode = raw_icode if raw_icode not in ('?', '.') else ' '
                    atom = {
                        'serial': serial,
                        # mmCIF quotes atom names containing ' (e.g. "C2'"); strip quotes
                        'name': (parts[need['auth_atom_id']] if need['auth_atom_id'] is not None else parts[need['label_atom_id']]).strip('"'),
                        'altloc': parts[need['label_alt_id']] if need['label_alt_id'] is not None and parts[need['label_alt_id']] != '.' else ' ',
                        'resname': resname,
                        'chain': parts[need['auth_asym_id']],
                        'resseq': resseq,
                        'icode': icode if icode else ' ',
                        'x': float(parts[need['Cartn_x']]),
                        'y': float(parts[need['Cartn_y']]),
                        'z': float(parts[need['Cartn_z']]),
                        'occupancy': float(parts[need['occupancy']]) if need['occupancy'] is not None else 1.0,
                        'bfactor': float(parts[need['B_iso_or_equiv']]) if need['B_iso_or_equiv'] is not None else 0.0,
                        'element': parts[need['type_symbol']] if need['type_symbol'] is not None else '',
                    }
                except Exception:
                    idx += 1
                    continue
                if atom['altloc'] == '.':
                    atom['altloc'] = ' '
                if not atom['element']:
                    atom['element'] = atom['name'][0]
                chain_atoms[atom['chain']].append(atom)
                serial += 1
                idx += 1
            # skip the terminating #
            if idx < n and lines[idx].strip().startswith('#'):
                idx += 1
        else:
            idx += 1
    return chain_atoms


def locate_source_file(pdb_id, ext):
    """Find source file case-insensitively."""
    for name in (f"{pdb_id}.{ext}.gz", f"{pdb_id.lower()}.{ext}.gz"):
        src = PDB_DIR / name
        if src.exists():
            return src
    return None


def process_pdb(args):
    pdb_id, chain_to_uids, ext, skip_existing = args

    src = locate_source_file(pdb_id, ext)
    if src is None and ext == 'pdb':
        # Fallback to cif if pdb not found.
        src = locate_source_file(pdb_id, 'cif')
    if src is None:
        missing = {chain: len(uids) for chain, uids in chain_to_uids.items()}
        return (pdb_id, 0, missing, f"missing_source:{pdb_id}.{ext}.gz")

    # If all expected output files already exist, skip parsing entirely.
    all_uids = [uid for uids in chain_to_uids.values() for uid in uids]
    if skip_existing and all((OUT_DIR / f"{uid}.pdb").exists() for uid in all_uids):
        return (pdb_id, len(all_uids), {}, None)

    actual_ext = src.suffixes[0].lstrip('.')
    try:
        if actual_ext == 'cif':
            chain_atoms = parse_cif_file(src)
        else:
            chain_atoms = parse_pdb_file(src)
    except Exception as exc:
        missing = {chain: len(uids) for chain, uids in chain_to_uids.items()}
        return (pdb_id, 0, missing, f"parse_error:{exc}")

    written = 0
    missing_chains = []
    for chain_id, uids in chain_to_uids.items():
        atoms = chain_atoms.get(chain_id)
        if not atoms:
            missing_chains.append(chain_id)
            continue
        # Renumber atoms per chain to keep PDB serial field within standard 5 columns.
        lines = [format_atom({**a, 'serial': i + 1}) for i, a in enumerate(atoms)]
        content = "".join(lines)
        for uid in uids:
            out_file = OUT_DIR / f"{uid}.pdb"
            if skip_existing and out_file.exists():
                written += 1
                continue
            out_file.parent.mkdir(parents=True, exist_ok=True)
            with open(out_file, "w") as fh:
                fh.write(content)
            written += 1
    err = None
    if missing_chains:
        err = f"missing_chains:{','.join(missing_chains)}"
    return (pdb_id, written, {}, err)


def main():
    global OUT_DIR
    parser = argparse.ArgumentParser()
    parser.add_argument("--train", type=str, default=str(TASK_DIR / "splits" / "ligand_binding_site_train.csv"))
    parser.add_argument("--val", type=str, default=str(TASK_DIR / "splits" / "ligand_binding_site_val.csv"))
    parser.add_argument("--test", type=str, default=str(TASK_DIR / "splits" / "ligand_binding_site_test.csv"))
    parser.add_argument("--workers", type=int, default=12)
    parser.add_argument("--out_dir", type=str, default=str(OUT_DIR))
    parser.add_argument("--skip_existing", action="store_true", default=True)
    parser.add_argument("--overwrite", action="store_true", help="Force overwrite of existing PDB files")
    args = parser.parse_args()

    OUT_DIR = Path(args.out_dir)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if args.overwrite:
        args.skip_existing = False

    dfs = []
    for path in [args.train, args.val, args.test]:
        dfs.append(pd.read_csv(path, low_memory=False, keep_default_na=False, na_values=['']))
    df = pd.concat(dfs, ignore_index=True)
    logger.info(f"Total rows: {len(df)}")

    # Build pdb -> chain -> list of unique_ids
    pdb_chains = defaultdict(lambda: defaultdict(list))
    pdb_ext = {}
    for _, row in df.iterrows():
        pdb_id = str(row['pdb_id']).upper()
        chain_id = str(row['chain_id'])
        pdb_chains[pdb_id][chain_id].append(row['unique_id'])
        if pdb_id not in pdb_ext:
            # Determine actual source extension by file existence (case-insensitive).
            if locate_source_file(pdb_id, 'pdb'):
                pdb_ext[pdb_id] = 'pdb'
            elif locate_source_file(pdb_id, 'cif'):
                pdb_ext[pdb_id] = 'cif'
            else:
                pdb_ext[pdb_id] = 'pdb'  # will be reported missing later

    tasks = [(pdb_id, dict(chain_to_uids), pdb_ext[pdb_id], args.skip_existing) for pdb_id, chain_to_uids in pdb_chains.items()]
    logger.info(f"Unique PDB files to parse: {len(tasks)}")

    total_written = 0
    failures = []
    missing_total = 0
    with Pool(args.workers) as pool:
        for pdb_id, written, missing_map, err in pool.imap_unordered(process_pdb, tasks, chunksize=1):
            total_written += written
            missing_total += sum(missing_map.values())
            if err:
                failures.append((pdb_id, err))
            if (total_written + len(failures)) % 5000 == 0:
                logger.info(f"Progress: written={total_written}, failures={len(failures)}")

    logger.info(f"Done. written={total_written}, missing={missing_total}, failures={len(failures)}")
    log_file = OUT_DIR.parent / "output" / "structure_extraction_log.csv"
    log_file.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(failures, columns=["pdb_id", "reason"]).to_csv(log_file, index=False)
    logger.info(f"Failure log saved to {log_file}")


if __name__ == "__main__":
    main()
