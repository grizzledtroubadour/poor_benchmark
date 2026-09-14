import gzip
import json
import os
import re
import sys
import time
import pandas as pd
from collections import defaultdict
from Bio.PDB import MMCIFParser, MMCIFIO, PDBParser, PDBIO, Select

# ==================== Config ====================
# 脚本内路径为运行目录相对路径，需在含 data/pdb/（项目级原始结构）与 pdbs/ 的工作根下执行
SOURCE_DIR = 'data/pdb/'
OUTPUT_DIR = 'pdbs/'
CHECKPOINT_FILE = 'output/extract_checkpoint.json'
LOG_FILE = 'logs/extract_representative.log'

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs('output', exist_ok=True)
os.makedirs('logs', exist_ok=True)

# ==================== Residue Range Helper ====================
def parse_resid(resid_str):
    """Parse residue identifier like '21A' or '-1' into (resseq_int, icode_str)."""
    resid_str = str(resid_str).strip()
    m = re.match(r'(-?\d+)(\D?)', resid_str)
    if m:
        return int(m.group(1)), m.group(2)
    return None, None

def resid_in_range(resseq, icode, beg_str, end_str):
    """Check if (resseq, icode) is within [beg_str, end_str] inclusive."""
    beg_seq, beg_icode = parse_resid(beg_str)
    end_seq, end_icode = parse_resid(end_str)
    if beg_seq is None or end_seq is None:
        return True
    def make_tuple(seq, icode):
        return (seq, icode.ljust(1) if icode else ' ')
    return make_tuple(beg_seq, beg_icode) <= make_tuple(resseq, icode) <= make_tuple(end_seq, end_icode)

# ==================== Extraction Classes ====================
class ResidueRangeSelect(Select):
    """BioPython Select for filtering by residue ranges."""
    def __init__(self, chain_id, ranges):
        self.target_chain = chain_id
        self.ranges = ranges
    def accept_chain(self, chain):
        return chain.id == self.target_chain
    def accept_residue(self, residue):
        for beg_str, end_str in self.ranges:
            if resid_in_range(residue.id[1], residue.id[2], beg_str, end_str):
                return True
        return False


def extract_from_pdb_text(pdb_path, tasks, output_dir):
    """Extract chains from PDB format (.pdb.gz) using text processing. Output .pdb."""
    if pdb_path.endswith('.gz'):
        with gzip.open(pdb_path, 'rt') as f:
            lines = f.readlines()
    else:
        with open(pdb_path, 'r') as f:
            lines = f.readlines()
    
    results = []
    for task in tasks:
        chain_id = task['chain_id']
        ranges = task['ranges']
        output_path = os.path.join(output_dir, f"{task['chain_key']}.pdb")
        
        out_lines = []
        for line in lines:
            # Multi-model (NMR) entries: keep model 1 only（根因修复，截断于首个 ENDMDL）
            if line.startswith('ENDMDL'):
                break
            if line.startswith('ATOM  ') or line.startswith('HETATM'):
                # Standard PDB: chain at col 22 (0-indexed 21), resseq at 23-26, icode at 27
                file_chain = line[21:22].strip()
                if file_chain != chain_id:
                    continue
                try:
                    resseq = int(line[22:26].strip())
                    icode = line[26:27] if len(line) > 26 else ' '
                except (ValueError, IndexError):
                    continue
                if any(resid_in_range(resseq, icode, b, e) for b, e in ranges):
                    out_lines.append(line)
            elif line.startswith('TER'):
                if out_lines:
                    out_lines.append(line)
        
        if out_lines:
            with open(output_path, 'w') as f:
                f.writelines(out_lines)
            results.append((task['chain_key'], True, len(out_lines)))
        else:
            results.append((task['chain_key'], False, 0))
    return results


def extract_from_cif(cif_path, tasks, output_dir):
    """Extract chains from mmCIF format (.cif.gz) using BioPython. Output .cif."""
    parser = MMCIFParser(QUIET=True)
    if cif_path.endswith('.gz'):
        with gzip.open(cif_path, 'rt') as f:
            structure = parser.get_structure('struct', f)
    else:
        structure = parser.get_structure('struct', cif_path)
    
    io = MMCIFIO()
    results = []
    for task in tasks:
        chain_id = task['chain_id']
        ranges = task['ranges']
        output_path = os.path.join(output_dir, f"{task['chain_key']}.cif")
        
        select = ResidueRangeSelect(chain_id, ranges)
        io.set_structure(structure)
        try:
            io.save(output_path, select=select)
            if os.path.exists(output_path) and os.path.getsize(output_path) > 0:
                results.append((task['chain_key'], True, os.path.getsize(output_path)))
            else:
                results.append((task['chain_key'], False, 0))
        except Exception as e:
            print(f"ERROR extracting {task['chain_key']}: {e}", file=sys.stderr)
            results.append((task['chain_key'], False, 0))
    return results


def process_pdb(pdb_id, tasks):
    """Process all chains for a single PDB."""
    pdb_path = os.path.join(SOURCE_DIR, f"{pdb_id}.pdb.gz")
    cif_path = os.path.join(SOURCE_DIR, f"{pdb_id}.cif.gz")
    
    if os.path.exists(pdb_path):
        return extract_from_pdb_text(pdb_path, tasks, OUTPUT_DIR)
    elif os.path.exists(cif_path):
        return extract_from_cif(cif_path, tasks, OUTPUT_DIR)
    else:
        print(f"ERROR: No source file for {pdb_id}", file=sys.stderr)
        return [(t['chain_key'], False, 0) for t in tasks]


# ==================== Main ====================
if __name__ == '__main__':
    print("=" * 60)
    print("Representative Structure Extraction")
    print("=" * 60)
    
    # Load representative chains（02_select_representative_chains.py 产物）
    print("\n[1/4] Loading representative chains...")
    rep_df = pd.read_csv('output/representative_chains.csv')
    
    # Load SIFTS mapping for A-class
    print("[2/4] Loading SIFTS mapping...")
    with open('data/sifts/sifts_uniprot_mapping.json', 'r') as f:
        sifts = json.load(f)
    
    # Load B-boundary inference for B-class
    # 注意：b_boundary_inference.jsonl 为 B 档边界推断中间产物，未随仓库迁移；
    # 重跑需先按 DATA_PROCESS.md 2.1 的子序列比对法重新生成（放到 output/ 下）。
    print("[3/4] Loading B-boundary inference...")
    b_boundaries = {}
    with open('output/b_boundary_inference.jsonl', 'r') as f:
        for line in f:
            d = json.loads(line.strip())
            key = f"{d['pdb_id']}_{d['chain_id']}"
            b_boundaries[key] = d
    
    # Build extraction tasks
    print("[4/4] Building extraction parameters...")
    extract_tasks = []
    for _, row in rep_df.iterrows():
        ck = row['chain_key']
        mode = row['extract_mode']
        
        if mode == 'A':
            info = sifts.get(ck, {})
            segments = info.get('segments', [])
            ranges = [(s.get('pdb_beg'), s.get('pdb_end')) 
                      for s in segments if s.get('pdb_beg') and s.get('pdb_end')]
            if not ranges:
                print(f"WARNING: A-class {ck} no valid segments", file=sys.stderr)
                continue
        else:
            info = b_boundaries.get(ck)
            if not info:
                print(f"WARNING: B-class {ck} no boundary", file=sys.stderr)
                continue
            ranges = [(str(info.get('inferred_pdb_beg')), str(info.get('inferred_pdb_end')))]
        
        extract_tasks.append({
            'chain_key': ck,
            'pdb_id': row['pdb_id'],
            'chain_id': row['chain_id'],
            'mode': mode,
            'ranges': ranges
        })
    
    print(f"\nTotal extraction tasks: {len(extract_tasks):,}")
    
    # Group by pdb_id
    pdb_tasks = defaultdict(list)
    for task in extract_tasks:
        pdb_tasks[task['pdb_id'].upper()].append(task)
    print(f"Unique PDBs: {len(pdb_tasks):,}")
    
    # Checkpoint
    completed_pdbs = set()
    if os.path.exists(CHECKPOINT_FILE):
        with open(CHECKPOINT_FILE, 'r') as f:
            completed_pdbs = set(json.load(f))
        print(f"Resuming: {len(completed_pdbs):,} PDBs already done")
    
    # Process
    all_pdb_ids = sorted(pdb_tasks.keys())
    remaining = [p for p in all_pdb_ids if p not in completed_pdbs]
    print(f"Remaining: {len(remaining):,} PDBs\n")
    
    success_count = 0
    fail_count = 0
    start_time = time.time()
    
    with open(LOG_FILE, 'a') as logf:
        for i, pdb_id in enumerate(remaining):
            tasks = pdb_tasks[pdb_id]
            results = process_pdb(pdb_id, tasks)
            
            for chain_key, ok, size in results:
                if ok:
                    success_count += 1
                else:
                    fail_count += 1
                    logf.write(f"FAIL {chain_key}\n")
            
            completed_pdbs.add(pdb_id)
            
            if (i + 1) % 500 == 0:
                with open(CHECKPOINT_FILE, 'w') as f:
                    json.dump(sorted(list(completed_pdbs)), f)
                elapsed = time.time() - start_time
                rate = (i + 1) / elapsed
                print(f"  {i+1:,}/{len(remaining):,} | Success: {success_count:,} | "
                      f"Fail: {fail_count:,} | Rate: {rate:.1f} PDBs/s")
    
    # Final checkpoint
    with open(CHECKPOINT_FILE, 'w') as f:
        json.dump(sorted(list(completed_pdbs)), f)
    
    print(f"\n{'='*60}")
    print("COMPLETE")
    print(f"{'='*60}")
    print(f"Total chains: {len(extract_tasks):,}")
    print(f"Success: {success_count:,}")
    print(f"Fail: {fail_count:,}")
    print(f"Time: {time.time() - start_time:.1f}s")
