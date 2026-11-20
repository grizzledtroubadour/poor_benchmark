#!/usr/bin/env python3
"""
Extract all CATH S95 domain structures.
Extract 62,915 representative S95 domains from source PDB/mmCIF files.
Output directory: fold_classification/output/cath_s95_domains/
"""

import gzip
import os
import sys
import json
import time
import pickle
import csv
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count

from Bio.PDB import MMCIFParser
from tqdm import tqdm

# ── Path configuration ──
TASK_DIR = Path(__file__).resolve().parents[1]  # datasets/fold_classification
ROOT = Path(__file__).resolve().parents[3]      # Repository root (POOR)
OUTPUT_DIR = TASK_DIR / "output" / "cath_s95_domains"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

DATA_CATH = ROOT / "data" / "CATHv44" / "cath-classification-data"
DATA_PDB = ROOT / "data" / "pdb"

CATH_DOMAIN_LIST_S95 = DATA_CATH / "cath-domain-list-S95.txt"
CATH_BOUNDARIES = DATA_CATH / "cath-domain-boundaries.txt"
CATH_SEQ_S95 = ROOT / "data" / "CATHv44" / "sequence-data" / "cath-domain-seqs-S95.fa"

CHECKPOINT_FILE = TASK_DIR / "output" / "extract_checkpoint.pkl"
LOG_FILE = TASK_DIR / "output" / "extract_s95.log"

CIF_PARSER = MMCIFParser(QUIET=True)

THREE_TO_ONE = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
    'SEC': 'U', 'PYL': 'O', 'UNK': 'X', 'MSE': 'M',
}


# ── Logging ──

def log(msg):
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{timestamp}] {msg}"
    print(line)
    with open(LOG_FILE, "a") as f:
        f.write(line + "\n")


# ── Data loading ──

def build_pdb_index():
    """Prebuild a {pdb_id_lower: Path} index."""
    index = {}
    for f in DATA_PDB.iterdir():
        if f.is_file() and f.suffix == '.gz':
            pdb_id = f.name.split('.')[0].lower()[:4]
            if pdb_id not in index:
                index[pdb_id] = f
    return index


def parse_cath_boundaries(boundaries_file: Path) -> dict:
    """Parse cath-domain-boundaries.txt (CDF Format 2.0)."""
    domain_map = {}
    with open(boundaries_file) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            tokens = line.split()
            if len(tokens) < 4:
                continue

            chain_name = tokens[0]
            n_domains = int(tokens[1][1:])
            n_fragments = int(tokens[2][1:])
            idx = 3

            pdb_id = chain_name[:4]
            chain_char = chain_name[4] if len(chain_name) > 4 else '0'

            for dom_i in range(1, n_domains + 1):
                n_segments = int(tokens[idx])
                idx += 1
                segments = []
                for _ in range(n_segments):
                    c_start = tokens[idx]
                    start = tokens[idx+1]
                    insert_start = tokens[idx+2]
                    c_end = tokens[idx+3]
                    end = tokens[idx+4]
                    insert_end = tokens[idx+5]

                    try:
                        start_num = int(start)
                    except ValueError:
                        start_num = None
                    try:
                        end_num = int(end)
                    except ValueError:
                        end_num = None

                    segments.append({
                        'chain': c_start if c_start != '0' else chain_char,
                        'start': start_num,
                        'start_icode': insert_start if insert_start != '-' else ' ',
                        'end': end_num,
                        'end_icode': insert_end if insert_end != '-' else ' ',
                    })
                    idx += 6

                domain_id = f"{pdb_id}{chain_char}{dom_i:02d}"
                domain_map[domain_id] = (chain_char, segments)

            for _ in range(n_fragments):
                idx += 6
                if idx < len(tokens) and tokens[idx].startswith('('):
                    idx += 1

    return domain_map


def load_s95_domain_list(list_file: Path) -> list:
    """Load the S95 list."""
    domains = []
    with open(list_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            tokens = line.split()
            if len(tokens) < 12:
                continue
            domains.append({
                'domain_id': tokens[0],
                'pdb_id': tokens[0][:4].lower(),
                'chain': tokens[0][4],
                'dom_num': tokens[0][5:7],
                'C': int(tokens[1]),
                'A': int(tokens[2]),
                'T': int(tokens[3]),
                'H': int(tokens[4]),
                'S35': int(tokens[5]),
                'S60': int(tokens[6]),
                'S95': int(tokens[7]),
                'S100': int(tokens[8]),
                'count': int(tokens[9]),
                'length': int(tokens[10]),
                'resolution': float(tokens[11]),
            })
    return domains


def load_fasta_sequences(fasta_file: Path) -> dict:
    """Load FASTA sequences."""
    seqs = {}
    current_id = None
    current_seq = []
    with open(fasta_file) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if current_id is not None:
                    seqs[current_id] = ''.join(current_seq)
                header = line[1:].strip()
                last_part = header.split('|')[-1].strip()
                current_id = last_part.split('/')[0].strip()
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_id is not None:
            seqs[current_id] = ''.join(current_seq)
    return seqs


# ── Domain extraction core ──

def parse_resnum_from_pdb_line(line: str) -> tuple:
    """Parse residue number and insertion code from a PDB ATOM record."""
    resnum_str = line[22:27]
    try:
        resnum = int(resnum_str)
    except ValueError:
        resnum = None
    icode = line[26] if len(line) > 26 else ' '
    return resnum, icode


def res_in_segments(resnum, icode, segments):
    """Check whether a residue belongs to any segment range."""
    for seg in segments:
        start = seg['start']
        end = seg['end']
        start_icode = seg.get('start_icode', ' ')
        end_icode = seg.get('end_icode', ' ')

        if start is not None and end is not None:
            if start < resnum < end:
                return True
            if resnum == start and icode == start_icode:
                return True
            if resnum == end and icode == end_icode:
                return True
        elif start is not None:
            if resnum > start:
                return True
            if resnum == start and icode == start_icode:
                return True
        elif end is not None:
            if resnum < end:
                return True
            if resnum == end and icode == end_icode:
                return True
    return False


def extract_from_pdb_gz(pdb_path: Path, chain_target: str, segments: list, whole_chain: bool) -> list:
    """Extract ATOM records from gzipped PDB."""
    import subprocess
    lines_out = []
    target_chains = [chain_target]
    if chain_target == '0':
        target_chains = [' ', 'A']

    proc = subprocess.Popen(['zcat', str(pdb_path)], stdout=subprocess.PIPE, text=True, errors='replace')
    try:
        for line in proc.stdout:
            # Multi-model (NMR) entries: keep model 1 only. Without this,
            # all models were concatenated, duplicating residues
            # (fixed post-hoc by output/scripts_archive/fix_multimodel_pdbs.py).
            if line.startswith("ENDMDL"):
                break
            if not line.startswith("ATOM"):
                continue
            chain = line[21]
            if chain not in target_chains:
                continue
            if whole_chain:
                lines_out.append(line)
                continue
            resnum, icode = parse_resnum_from_pdb_line(line)
            if resnum is None:
                continue
            if res_in_segments(resnum, icode, segments):
                lines_out.append(line)
    finally:
        proc.stdout.close()
        proc.wait()

    return lines_out


def extract_from_cif_gz(cif_path: Path, chain_target: str, segments: list, whole_chain: bool) -> list:
    """Extract from gzipped mmCIF and emit PDB-format ATOM records."""
    try:
        with gzip.open(cif_path, 'rt') as f:
            structure = CIF_PARSER.get_structure('tmp', f)
    except Exception:
        return []

    target_chains = [chain_target]
    if chain_target == '0':
        target_chains = [' ', 'A']

    lines_out = []
    atom_serial = 1

    for model in structure:
        for chain in model:
            if chain.id not in target_chains:
                continue
            for residue in chain:
                res_id = residue.get_id()
                res_seq = res_id[1]
                res_icode = res_id[2] if res_id[2] not in ('', ' ') else ' '

                if not whole_chain and not res_in_segments(res_seq, res_icode, segments):
                    continue

                res_name = residue.resname
                for atom in residue:
                    x, y, z = atom.coord
                    occ = atom.occupancy if atom.occupancy is not None else 1.0
                    bfac = atom.bfactor if atom.bfactor is not None else 0.0
                    elem = atom.element if atom.element else atom.name[0]
                    aname = atom.name.strip()
                    if len(aname) <= 3:
                        aname_fmt = f" {aname:>3s}"
                    else:
                        aname_fmt = f"{aname:<4s}"
                    line = (
                        f"ATOM  {atom_serial:5d} {aname_fmt}{' ':1s}{res_name:3s} {chain.id:1s}"
                        f"{res_seq:4d}{res_icode:1s}   {x:8.3f}{y:8.3f}{z:8.3f}"
                        f"{occ:6.2f}{bfac:6.2f}          {elem:>2s}  \n"
                    )
                    lines_out.append(line)
                    atom_serial += 1

    return lines_out


def pdb_lines_to_sequence(atom_lines: list) -> str:
    """Extract amino-acid sequence from PDB ATOM records."""
    seen = set()
    seq = []
    for line in atom_lines:
        if not line.startswith("ATOM"):
            continue
        res_name = line[17:20].strip()
        chain = line[21]
        res_num_str = line[22:27]
        res_id = (chain, res_num_str)
        if res_id not in seen:
            seen.add(res_id)
            seq.append(THREE_TO_ONE.get(res_name, 'X'))
    return ''.join(seq)


# ── Summary generation, incorporating former fix_summary.py logic ──

SUMMARY_FIELDS = [
    'domain_id', 'pdb_id', 'chain', 'dom_num', 'type',
    'C', 'A', 'T', 'H', 'S35', 'S60', 'S95', 'S100',
    'length_expected', 'n_residues', 'n_atoms', 'resolution',
]


def count_residues_and_atoms(pdb_path: Path) -> tuple:
    """Read a PDB and count residues and atoms."""
    seen_residues = set()
    n_atoms = 0
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM"):
                n_atoms += 1
                res_name = line[17:20].strip()
                chain = line[21]
                res_num = line[22:27]
                res_id = (chain, res_num, res_name)
                seen_residues.add(res_id)
    return len(seen_residues), n_atoms


def regenerate_summary():
    """Scan extracted files and rebuild the complete output/cath_s95_summary.csv.

    Historical bug: earlier code appended CSV checkpoints inside the parallel extraction loop;
    processed %500 ==0 rarely triggered, saving only some batches.
    Now rebuild from all output files after extraction, or when all domains hit checkpoints,
    producing a complete summary directly.
    """
    domains_dir = OUTPUT_DIR
    if not domains_dir.is_dir() or not any(domains_dir.glob('*.pdb')):
        # Extracted files have moved to the final pdbs/ structure library.
        domains_dir = TASK_DIR / "pdbs"

    meta = {d['domain_id']: d for d in load_s95_domain_list(CATH_DOMAIN_LIST_S95)}
    pdb_files = sorted(domains_dir.glob('*.pdb'))
    print(f"\n[*] Rebuilding summary: scanning {len(pdb_files)} PDB files ({domains_dir})")

    summary_csv = TASK_DIR / "output" / "cath_s95_summary.csv"
    n_written = 0
    with open(summary_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        for pdb_path in pdb_files:
            domain_id = pdb_path.stem
            m = meta.get(domain_id)
            if m is None:
                continue
            n_res, n_atoms = count_residues_and_atoms(pdb_path)
            writer.writerow({
                'domain_id': domain_id,
                'pdb_id': m['pdb_id'],
                'chain': m['chain'],
                'dom_num': m['dom_num'],
                'type': 'whole' if m['dom_num'] == '00' else 'chopped',
                'C': m['C'], 'A': m['A'], 'T': m['T'], 'H': m['H'],
                'S35': m['S35'], 'S60': m['S60'], 'S95': m['S95'], 'S100': m['S100'],
                'length_expected': m['length'],
                'n_residues': n_res,
                'n_atoms': n_atoms,
                'resolution': m['resolution'],
            })
            n_written += 1
    print(f"    Summary saved: {summary_csv} ({n_written} rows)")


# ── Worker functions ──

def process_batch(batch_data):
    """Process one domain batch and return results."""
    domain_batch, boundaries, pdb_index = batch_data
    results = []

    for d in domain_batch:
        domain_id = d['domain_id']
        pdb_id = d['pdb_id']
        chain = d['chain']
        dom_num = d['dom_num']

        result = {
            'domain_id': domain_id,
            'pdb_id': pdb_id,
            'chain': chain,
            'dom_num': dom_num,
            'type': 'whole' if dom_num == '00' else 'chopped',
            'length_expected': d['length'],
            'resolution': d['resolution'],
            'C': d['C'], 'A': d['A'], 'T': d['T'], 'H': d['H'],
            'S35': d['S35'], 'S60': d['S60'], 'S95': d['S95'], 'S100': d['S100'],
            'status': 'pending',
            'n_atoms': 0,
            'n_residues': 0,
            'seq_extracted': '',
            'error': '',
        }

        pdb_file = pdb_index.get(pdb_id)
        if pdb_file is None:
            result['status'] = 'missing_pdb'
            result['error'] = f'PDB file not found for {pdb_id}'
            results.append(result)
            continue

        is_pdb = '.pdb' in pdb_file.name
        whole_chain = (dom_num == '00')

        if whole_chain:
            if is_pdb:
                atom_lines = extract_from_pdb_gz(pdb_file, chain, [], True)
            else:
                atom_lines = extract_from_cif_gz(pdb_file, chain, [], True)
        else:
            if domain_id in boundaries:
                _, segments = boundaries[domain_id]
                if is_pdb:
                    atom_lines = extract_from_pdb_gz(pdb_file, chain, segments, False)
                else:
                    atom_lines = extract_from_cif_gz(pdb_file, chain, segments, False)
            else:
                result['status'] = 'no_boundary'
                result['error'] = f'No PDB boundary for {domain_id}'
                results.append(result)
                continue

        if not atom_lines:
            result['status'] = 'empty'
            result['error'] = 'No atoms extracted'
            results.append(result)
            continue

        out_pdb = OUTPUT_DIR / f"{domain_id}.pdb"
        try:
            with open(out_pdb, 'w') as f:
                f.writelines(atom_lines)
        except Exception as e:
            result['status'] = 'write_error'
            result['error'] = str(e)
            results.append(result)
            continue

        result['n_atoms'] = len(atom_lines)
        result['seq_extracted'] = pdb_lines_to_sequence(atom_lines)
        result['n_residues'] = len(result['seq_extracted'])
        result['status'] = 'success'
        results.append(result)

    return results


# ── Main workflow ──

def main():
    print("=" * 70)
    print("Full CATH S95 domain extraction")
    print("=" * 70)

    # Load data.
    print("\n[1/5] Loading S95 domain list...")
    s95_domains = load_s95_domain_list(CATH_DOMAIN_LIST_S95)
    print(f"      Total S95 domains: {len(s95_domains)}")

    print("\n[2/5] Loading domain boundaries...")
    boundaries = parse_cath_boundaries(CATH_BOUNDARIES)
    print(f"      Domains with PDB boundaries: {len(boundaries)}")

    print("\n[3/5] Building PDB file index...")
    pdb_index = build_pdb_index()
    print(f"      Local PDB files: {len(pdb_index)}")

    # Check resume state.
    completed_ids = set()
    if CHECKPOINT_FILE.exists():
        print(f"\n[*] Checkpoint found; loading completed IDs...")
        try:
            with open(CHECKPOINT_FILE, 'rb') as f:
                completed_ids = pickle.load(f)
            print(f"    Completed: {len(completed_ids)}")
        except Exception:
            completed_ids = set()

    # Select remaining domains.
    todo_domains = [d for d in s95_domains if d['domain_id'] not in completed_ids]
    print(f"\n[4/5] Pending domains: {len(todo_domains)} / {len(s95_domains)}")

    if len(todo_domains) == 0:
        print("\nAll domains have been processed!")
        regenerate_summary()
        return

    # Divide into batches.
    n_workers = max(1, cpu_count() - 1)
    batch_size = max(1, len(todo_domains) // (n_workers * 4))
    batches = []
    for i in range(0, len(todo_domains), batch_size):
        batches.append((todo_domains[i:i+batch_size], boundaries, pdb_index))

    print(f"\n[5/5] Starting parallel extraction: {n_workers} workers, {len(batches)} batches")

    all_results = []
    processed = len(completed_ids)

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {executor.submit(process_batch, batch): i for i, batch in enumerate(batches)}

        with tqdm(total=len(todo_domains), desc="Extraction", initial=0) as pbar:
            for future in as_completed(futures):
                batch_results = future.result()
                all_results.extend(batch_results)
                processed += len(batch_results)
                pbar.update(len(batch_results))

                # Save progress.
                if processed % 500 == 0:
                    new_completed = completed_ids | {r['domain_id'] for r in all_results if r['status'] == 'success'}
                    with open(CHECKPOINT_FILE, 'wb') as f:
                        pickle.dump(new_completed, f)

    # Final save
    final_completed = completed_ids | {r['domain_id'] for r in all_results if r['status'] == 'success'}
    with open(CHECKPOINT_FILE, 'wb') as f:
        pickle.dump(final_completed, f)

    # Statistics
    print("\n" + "=" * 70)
    print("Extraction complete")
    print("=" * 70)

    status_counts = defaultdict(int)
    for r in all_results:
        status_counts[r['status']] += 1

    print(f"\nStatus distribution:")
    for status, count in sorted(status_counts.items(), key=lambda x: -x[1]):
        print(f"  {status}: {count}")

    success_count = status_counts.get('success', 0)
    print(f"\nSuccessfully extracted: {len(final_completed)} / {len(s95_domains)} ({len(final_completed)/len(s95_domains)*100:.1f}%)")
    print(f"Output directory: {OUTPUT_DIR}")

    # Rebuild the entire summary from output files after extraction, replacing broken incremental writes.
    regenerate_summary()


if __name__ == "__main__":
    main()
