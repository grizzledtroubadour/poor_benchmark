#!/usr/bin/env python3
"""
CATH S95 结构域全量提取脚本。
从 PDB/mmCIF 源文件中提取 62,915 个 S95 代表域的结构。
输出目录: fold_classification/output/cath_s95_domains/
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

# ── 路径配置 ──────────────────────────────────────────────────────────
TASK_DIR = Path(__file__).resolve().parents[1]  # datasets/fold_classification
ROOT = Path(__file__).resolve().parents[3]      # 仓库根（POOR）
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


# ── 日志 ──────────────────────────────────────────────────────────────

def log(msg):
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{timestamp}] {msg}"
    print(line)
    with open(LOG_FILE, "a") as f:
        f.write(line + "\n")


# ── 数据加载 ──────────────────────────────────────────────────────────

def build_pdb_index():
    """预建 {pdb_id_lower: Path} 索引"""
    index = {}
    for f in DATA_PDB.iterdir():
        if f.is_file() and f.suffix == '.gz':
            pdb_id = f.name.split('.')[0].lower()[:4]
            if pdb_id not in index:
                index[pdb_id] = f
    return index


def parse_cath_boundaries(boundaries_file: Path) -> dict:
    """解析 cath-domain-boundaries.txt (CDF Format 2.0)"""
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
    """加载 S95 列表"""
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
    """加载 FASTA 序列"""
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


# ── 结构域提取核心 ────────────────────────────────────────────────────

def parse_resnum_from_pdb_line(line: str) -> tuple:
    """从 PDB ATOM 行解析残基编号和 insertion code"""
    resnum_str = line[22:27]
    try:
        resnum = int(resnum_str)
    except ValueError:
        resnum = None
    icode = line[26] if len(line) > 26 else ' '
    return resnum, icode


def res_in_segments(resnum, icode, segments):
    """检查残基是否在任何 segment 范围内"""
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
    """从 gzipped PDB 提取 ATOM 行"""
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
    """从 gzipped mmCIF 提取，输出 PDB 格式 ATOM 行"""
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
    """从 PDB ATOM 行提取氨基酸序列"""
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


# ── 汇总表生成（原 fix_summary.py 逻辑，已融入本脚本）─────────────────────

SUMMARY_FIELDS = [
    'domain_id', 'pdb_id', 'chain', 'dom_num', 'type',
    'C', 'A', 'T', 'H', 'S35', 'S60', 'S95', 'S100',
    'length_expected', 'n_residues', 'n_atoms', 'resolution',
]


def count_residues_and_atoms(pdb_path: Path) -> tuple:
    """读取 PDB 文件，统计残基数和原子数"""
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
    """扫描提取产物，重新生成完整汇总表 output/cath_s95_summary.csv。

    历史缺陷：早期版本在并行提取循环内按 checkpoint 追加写 CSV
    （`processed % 500 == 0` 判定几乎不命中，仅保存了部分批次）。
    现统一在提取结束后（或全部命中 checkpoint 时）从产物文件全量重建，
    直接产出完整汇总表。
    """
    domains_dir = OUTPUT_DIR
    if not domains_dir.is_dir() or not any(domains_dir.glob('*.pdb')):
        # 提取产物已迁移为最终结构库 pdbs/（当前仓库状态）
        domains_dir = TASK_DIR / "pdbs"

    meta = {d['domain_id']: d for d in load_s95_domain_list(CATH_DOMAIN_LIST_S95)}
    pdb_files = sorted(domains_dir.glob('*.pdb'))
    print(f"\n[*] 重新生成汇总表: 扫描 {len(pdb_files)} 个 PDB 文件 ({domains_dir})")

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
    print(f"    汇总表已保存: {summary_csv} ({n_written} 行)")


# ── Worker 函数 ───────────────────────────────────────────────────────

def process_batch(batch_data):
    """处理一批域，返回结果列表"""
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


# ── 主流程 ────────────────────────────────────────────────────────────

def main():
    print("=" * 70)
    print("CATH S95 结构域全量提取")
    print("=" * 70)
    
    # 加载数据
    print("\n[1/5] 加载 S95 域列表...")
    s95_domains = load_s95_domain_list(CATH_DOMAIN_LIST_S95)
    print(f"      S95 域总数: {len(s95_domains)}")
    
    print("\n[2/5] 加载结构域边界...")
    boundaries = parse_cath_boundaries(CATH_BOUNDARIES)
    print(f"      PDB 边界域数: {len(boundaries)}")
    
    print("\n[3/5] 构建 PDB 文件索引...")
    pdb_index = build_pdb_index()
    print(f"      本地 PDB 文件数: {len(pdb_index)}")
    
    # 检查断点续传
    completed_ids = set()
    if CHECKPOINT_FILE.exists():
        print(f"\n[*] 发现检查点文件，加载已完成列表...")
        try:
            with open(CHECKPOINT_FILE, 'rb') as f:
                completed_ids = pickle.load(f)
            print(f"    已完成: {len(completed_ids)}")
        except Exception:
            completed_ids = set()
    
    # 筛选待处理域
    todo_domains = [d for d in s95_domains if d['domain_id'] not in completed_ids]
    print(f"\n[4/5] 待处理域: {len(todo_domains)} / {len(s95_domains)}")

    if len(todo_domains) == 0:
        print("\n所有域已处理完成！")
        regenerate_summary()
        return
    
    # 分 batch
    n_workers = max(1, cpu_count() - 1)
    batch_size = max(1, len(todo_domains) // (n_workers * 4))
    batches = []
    for i in range(0, len(todo_domains), batch_size):
        batches.append((todo_domains[i:i+batch_size], boundaries, pdb_index))
    
    print(f"\n[5/5] 启动并行提取: {n_workers} workers, {len(batches)} batches")

    all_results = []
    processed = len(completed_ids)

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {executor.submit(process_batch, batch): i for i, batch in enumerate(batches)}

        with tqdm(total=len(todo_domains), desc="提取", initial=0) as pbar:
            for future in as_completed(futures):
                batch_results = future.result()
                all_results.extend(batch_results)
                processed += len(batch_results)
                pbar.update(len(batch_results))

                # 保存进度
                if processed % 500 == 0:
                    new_completed = completed_ids | {r['domain_id'] for r in all_results if r['status'] == 'success'}
                    with open(CHECKPOINT_FILE, 'wb') as f:
                        pickle.dump(new_completed, f)
    
    # 最终保存
    final_completed = completed_ids | {r['domain_id'] for r in all_results if r['status'] == 'success'}
    with open(CHECKPOINT_FILE, 'wb') as f:
        pickle.dump(final_completed, f)
    
    # 统计
    print("\n" + "=" * 70)
    print("提取完成")
    print("=" * 70)
    
    status_counts = defaultdict(int)
    for r in all_results:
        status_counts[r['status']] += 1
    
    print(f"\n状态分布:")
    for status, count in sorted(status_counts.items(), key=lambda x: -x[1]):
        print(f"  {status}: {count}")
    
    success_count = status_counts.get('success', 0)
    print(f"\n成功提取: {len(final_completed)} / {len(s95_domains)} ({len(final_completed)/len(s95_domains)*100:.1f}%)")
    print(f"输出目录: {OUTPUT_DIR}")

    # 提取结束后从产物文件全量重建汇总表（替代有缺陷的增量追加写）
    regenerate_summary()


if __name__ == "__main__":
    main()
