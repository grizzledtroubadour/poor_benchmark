#!/usr/bin/env python3
"""
批量解析原始 PDB/mmCIF 文件的 deposition date，映射到 chain_key。
"""
import gzip
import json
from datetime import datetime
from pathlib import Path
from tqdm import tqdm

TASK_DIR = Path(__file__).resolve().parent.parent  # datasets/func_prediction/
POOR_ROOT = TASK_DIR.parent.parent  # POOR/
PDB_DIR = POOR_ROOT / "data/pdb"
OUT_JSON = TASK_DIR / "output/deposition_dates.json"

def parse_pdb_date_from_gz(filepath):
    """从 gzipped PDB 文件的 HEADER 行解析日期。"""
    month_map = {
        'JAN': 1, 'FEB': 2, 'MAR': 3, 'APR': 4, 'MAY': 5, 'JUN': 6,
        'JUL': 7, 'AUG': 8, 'SEP': 9, 'OCT': 10, 'NOV': 11, 'DEC': 12
    }
    with gzip.open(filepath, 'rt') as f:
        for line in f:
            if line.startswith('HEADER'):
                stripped = line.strip()
                parts = stripped.split()
                for part in parts:
                    if len(part) == 9 and part[2] == '-' and part[6] == '-':
                        try:
                            day = int(part[:2])
                            month = month_map.get(part[3:6].upper())
                            year_str = part[7:9]
                            if month:
                                year = int(year_str)
                                if year >= 50:
                                    year += 1900
                                else:
                                    year += 2000
                                return f"{year:04d}-{month:02d}-{day:02d}"
                        except:
                            pass
    return None


def parse_cif_date_from_gz(filepath):
    """从 gzipped mmCIF 文件解析 deposition date。"""
    with gzip.open(filepath, 'rt') as f:
        for line in f:
            line_stripped = line.strip()
            if line_stripped.startswith('_pdbx_database_status.recvd_initial_deposition_date'):
                parts = line_stripped.split()
                if len(parts) >= 2:
                    date_str = parts[1]
                    try:
                        dt = datetime.strptime(date_str, '%Y-%m-%d')
                        return dt.strftime('%Y-%m-%d')
                    except ValueError:
                        try:
                            dt = datetime.strptime(date_str, '%Y/%m/%d')
                            return dt.strftime('%Y-%m-%d')
                        except ValueError:
                            pass
    return None


def main():
    # Load all chain_keys from final representatives
    import csv
    chain_keys = []
    with open(TASK_DIR / 'output/final_nonredundant_representatives.csv', 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            chain_keys.append(row['chain_key'])
    
    dates = {}
    failed = []
    
    for chain_key in tqdm(chain_keys, desc="Parsing deposition dates"):
        pdb_id = chain_key.split('_')[0]
        
        # Try PDB first
        pdb_path = PDB_DIR / f"{pdb_id}.pdb.gz"
        cif_path = PDB_DIR / f"{pdb_id}.cif.gz"
        
        date = None
        if pdb_path.exists():
            date = parse_pdb_date_from_gz(pdb_path)
        elif cif_path.exists():
            date = parse_cif_date_from_gz(cif_path)
        
        if date:
            dates[chain_key] = date
        else:
            failed.append(chain_key)
    
    with open(OUT_JSON, 'w') as f:
        json.dump(dates, f, indent=2)
    
    print(f"Parsed: {len(dates)}/{len(chain_keys)}")
    if failed:
        print(f"Failed ({len(failed)}): {failed[:10]}...")


if __name__ == '__main__':
    main()
