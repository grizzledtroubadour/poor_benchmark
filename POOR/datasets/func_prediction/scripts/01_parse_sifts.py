#!/usr/bin/env python3
"""
Stage 1: SIFTS Processing
Parse SIFTS flat files to generate PDB chain to UniProt mappings with annotations.

输入：项目级共享 SIFTS flat files（POOR/data/sifts/*.tsv.gz，需先自行下载）。
输出：POOR/data/sifts/ 下的 sifts_pdb_chains.csv、sifts_uniprot_mapping.json、
sifts_go_annotations.json、sifts_ec_annotations.json（项目级共享产物）。
"""

import os
import sys
import json
import gzip
import csv
import logging
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Optional
import requests
from tqdm import tqdm
import time

# Paths
TASK_DIR = Path(__file__).resolve().parent.parent  # datasets/func_prediction/
POOR_ROOT = TASK_DIR.parent.parent                 # POOR/
SIFTS_DIR = POOR_ROOT / "data" / "sifts"           # 项目级共享 SIFTS 目录
OUTPUT_DIR = POOR_ROOT / "data" / "sifts"
UNIPROT_DIR = TASK_DIR / "output" / "uniprot"      # UniProt 参考序列与元数据
LOG_DIR = TASK_DIR / "output" / "logs"

for d in [OUTPUT_DIR, UNIPROT_DIR, LOG_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[
        logging.FileHandler(LOG_DIR / "stage1.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)


def stage1_process_sifts(limit: Optional[int] = None, sifts_dir: Optional[Path] = None):
    """
    Process SIFTS flat files to generate:
    - sifts_pdb_chains.csv
    - sifts_uniprot_mapping.json
    - sifts_go_annotations.json
    - sifts_ec_annotations.json
    
    Args:
        limit: Limit number of records for testing
        sifts_dir: Custom SIFTS directory (default: DATA_DIR/sifts)
    
    Returns:
        Tuple of (chains_df, uniprot_mapping, go_annotations, ec_annotations)
    """
    if sifts_dir is None:
        sifts_dir = SIFTS_DIR
    
    logger.info("=== Stage 1: SIFTS Processing ===")
    logger.info(f"SIFTS directory: {sifts_dir}")
    
    # Check required files
    required_files = [
        sifts_dir / "pdb_chain_uniprot.tsv.gz",
        sifts_dir / "sifts_chain_go.tsv.gz",
        sifts_dir / "sifts_chain_ec.tsv.gz",
    ]
    for f in required_files:
        if not f.exists():
            raise FileNotFoundError(f"Required SIFTS file not found: {f}")
    
    # ---- Step 1: Parse uniprot segments observed (residue-level mapping) ----
    logger.info("Step 1: Parsing SIFTS uniprot segments...")
    pdb_chain_uniprot = defaultdict(list)
    
    with gzip.open(sifts_dir / "pdb_chain_uniprot.tsv.gz", 'rt') as f:
        # Skip comment lines starting with #
        while True:
            pos = f.tell()
            line = f.readline()
            if not line.startswith('#'):
                f.seek(pos)
                break
        
        reader = csv.DictReader(f, delimiter='\t')
        for i, row in enumerate(tqdm(reader, desc="Uniprot segments")):
            if limit and i >= limit:
                break
            
            pdb = row['PDB'].upper()
            chain = row['CHAIN']
            uniprot = row['SP_PRIMARY']
            pdb_beg = row.get('PDB_BEG', '')
            pdb_end = row.get('PDB_END', '')
            sp_beg = row.get('SP_BEG', '')
            sp_end = row.get('SP_END', '')
            
            key = (pdb, chain)
            pdb_chain_uniprot[key].append({
                'uniprot': uniprot,
                'pdb_beg': pdb_beg,
                'pdb_end': pdb_end,
                'sp_beg': sp_beg,
                'sp_end': sp_end
            })
    
    logger.info(f"  -> Found {len(pdb_chain_uniprot)} PDB-chain combinations")
    
    # ---- Step 2: Parse GO annotations ----
    logger.info("Step 2: Parsing SIFTS GO annotations...")
    pdb_chain_go = defaultdict(list)
    
    with gzip.open(sifts_dir / "sifts_chain_go.tsv.gz", 'rt') as f:
        while True:
            pos = f.tell()
            line = f.readline()
            if not line.startswith('#'):
                f.seek(pos)
                break
        
        reader = csv.DictReader(f, delimiter='\t')
        for row in tqdm(reader, desc="GO annotations"):
            pdb = row['PDB'].upper()
            chain = row['CHAIN']
            uniprot = row['SP_PRIMARY']
            go_id = row['GO_ID']
            evidence = row['EVIDENCE']
            with_str = row.get('WITH_STRING', '')
            
            key = (pdb, chain)
            # deduplicate by GO_ID
            existing = {g['go_id'] for g in pdb_chain_go[key]}
            if go_id not in existing:
                pdb_chain_go[key].append({
                    'uniprot': uniprot,
                    'go_id': go_id,
                    'evidence': evidence,
                    'with_string': with_str
                })
    
    logger.info(f"  -> Found GO annotations for {len(pdb_chain_go)} PDB-chain combinations")
    
    # ---- Step 3: Parse EC annotations ----
    logger.info("Step 3: Parsing SIFTS EC annotations...")
    pdb_chain_ec = defaultdict(list)
    
    with gzip.open(sifts_dir / "sifts_chain_ec.tsv.gz", 'rt') as f:
        while True:
            pos = f.tell()
            line = f.readline()
            if not line.startswith('#'):
                f.seek(pos)
                break
        
        reader = csv.DictReader(f, delimiter='\t')
        for row in tqdm(reader, desc="EC annotations"):
            pdb = row['PDB'].upper()
            chain = row['CHAIN']
            uniprot = row['ACCESSION']
            ec = row['EC_NUMBER']
            
            key = (pdb, chain)
            existing = {e['ec_number'] for e in pdb_chain_ec[key]}
            if ec not in existing:
                pdb_chain_ec[key].append({
                    'uniprot': uniprot,
                    'ec_number': ec
                })
    
    logger.info(f"  -> Found EC annotations for {len(pdb_chain_ec)} PDB-chain combinations")
    
    # ---- Step 4: Get unique PDB entries for metadata query ----
    unique_pdbs = sorted({k[0] for k in pdb_chain_uniprot.keys()})
    logger.info(f"Step 4: Total unique PDB entries: {len(unique_pdbs)}")
    
    # ---- Step 5: Fetch PDB metadata via RCSB API (batch) ----
    logger.info("Step 5: Fetching PDB metadata from RCSB GraphQL API...")
    pdb_metadata = fetch_pdb_metadata_batch(unique_pdbs)
    
    # ---- Step 6: Build output files ----
    logger.info("Step 6: Building Stage 1 outputs...")
    
    import pandas as pd
    
    # sifts_pdb_chains.csv
    chains_records = []
    for (pdb, chain), mappings in pdb_chain_uniprot.items():
        meta = pdb_metadata.get(pdb, {})
        # Calculate mapped length (sum of all segments)
        mapped_len = 0
        for m in mappings:
            try:
                sp_beg = int(m['sp_beg']) if m['sp_beg'] else 0
                sp_end = int(m['sp_end']) if m['sp_end'] else 0
                mapped_len += max(0, sp_end - sp_beg + 1)
            except ValueError:
                pass
        
        chains_records.append({
            'pdb_id': pdb,
            'chain_id': chain,
            'uniprot': mappings[0]['uniprot'],
            'resolution': meta.get('resolution'),
            'method': meta.get('method'),
            'r_work': meta.get('r_work'),
            'r_free': meta.get('r_free'),
            'mapped_length': mapped_len,
            'num_segments': len(mappings)
        })
    
    df_chains = pd.DataFrame(chains_records)
    df_chains.to_csv(OUTPUT_DIR / "sifts_pdb_chains.csv", index=False)
    logger.info(f"  -> Written {OUTPUT_DIR / 'sifts_pdb_chains.csv'} ({len(df_chains)} chains)")
    
    # sifts_uniprot_mapping.json
    uniprot_mapping = {}
    for (pdb, chain), mappings in pdb_chain_uniprot.items():
        key = f"{pdb}_{chain}"
        uniprot_mapping[key] = {
            'pdb': pdb,
            'chain': chain,
            'uniprot': mappings[0]['uniprot'],
            'segments': mappings
        }
    
    with open(OUTPUT_DIR / "sifts_uniprot_mapping.json", 'w') as f:
        json.dump(uniprot_mapping, f)
    logger.info(f"  -> Written {OUTPUT_DIR / 'sifts_uniprot_mapping.json'}")
    
    # sifts_go_annotations.json
    go_annotations = {}
    for (pdb, chain), annotations in pdb_chain_go.items():
        key = f"{pdb}_{chain}"
        go_annotations[key] = {
            'pdb': pdb,
            'chain': chain,
            'annotations': annotations
        }
    
    with open(OUTPUT_DIR / "sifts_go_annotations.json", 'w') as f:
        json.dump(go_annotations, f)
    logger.info(f"  -> Written {OUTPUT_DIR / 'sifts_go_annotations.json'}")
    
    # sifts_ec_annotations.json
    ec_annotations = {}
    for (pdb, chain), annotations in pdb_chain_ec.items():
        key = f"{pdb}_{chain}"
        ec_annotations[key] = {
            'pdb': pdb,
            'chain': chain,
            'annotations': annotations
        }
    
    with open(OUTPUT_DIR / "sifts_ec_annotations.json", 'w') as f:
        json.dump(ec_annotations, f)
    logger.info(f"  -> Written {OUTPUT_DIR / 'sifts_ec_annotations.json'}")
    
    # Summary stats
    total_pdbs = len(unique_pdbs)
    total_chains = len(df_chains)
    total_go = len(go_annotations)
    total_ec = len(ec_annotations)
    logger.info(f"Stage 1 Complete: {total_pdbs} PDBs, {total_chains} chains, {total_go} with GO, {total_ec} with EC")
    
    return df_chains, uniprot_mapping, go_annotations, ec_annotations


def fetch_pdb_metadata_batch(pdb_ids: List[str], batch_size: int = 50) -> Dict:
    """
    Fetch PDB metadata from RCSB GraphQL API in batches.
    
    Args:
        pdb_ids: List of PDB IDs to fetch
        batch_size: Number of IDs per batch
    
    Returns:
        Dict mapping PDB ID -> metadata dict
    """
    metadata = {}
    url = "https://data.rcsb.org/graphql"
    
    # Filter to only valid PDB entries (4-char codes starting with digit)
    valid_pdbs = [p for p in pdb_ids if len(p) == 4 and p[0].isdigit()]
    logger.info(f"Fetching metadata for {len(valid_pdbs)} valid PDB entries (filtered from {len(pdb_ids)} total)")
    
    for i in tqdm(range(0, len(valid_pdbs), batch_size), desc="PDB metadata batches"):
        batch = valid_pdbs[i:i+batch_size]
        ids_str = ', '.join([f'"{pdb.lower()}"' for pdb in batch])
        
        query = f"""
        {{
          entries(entry_ids: [{ids_str}]) {{
            rcsb_id
            exptl {{
              method
            }}
            refine {{
              ls_d_res_high
              ls_R_factor_R_work
              ls_R_factor_R_free
            }}
            em_3d_reconstruction {{
              resolution
            }}
          }}
        }}
        """
        
        try:
            resp = requests.post(url, json={"query": query}, timeout=120)
            resp.raise_for_status()
            data = resp.json()
            
            if 'errors' in data:
                logger.debug(f"Batch {i//batch_size} errors: {data['errors'][:2]}")
            
            entries = data.get('data', {}).get('entries', [])
            for entry in entries:
                pdb_id = entry['rcsb_id'].upper()
                exptl = entry.get('exptl', [{}])[0] if entry.get('exptl') else {}
                refine = entry.get('refine', [{}])[0] if entry.get('refine') else {}
                em_recon = entry.get('em_3d_reconstruction', [{}])[0] if entry.get('em_3d_reconstruction') else {}
                
                # Resolution: prefer refine (X-ray), fallback to em_3d_reconstruction (EM)
                resolution = refine.get('ls_d_res_high')
                if resolution is None and em_recon:
                    resolution = em_recon.get('resolution')
                
                metadata[pdb_id] = {
                    'method': exptl.get('method'),
                    'resolution': resolution,
                    'r_work': refine.get('ls_R_factor_R_work'),
                    'r_free': refine.get('ls_R_factor_R_free'),
                }
        except Exception as e:
            logger.warning(f"Batch {i//batch_size} failed: {e}")
        
        time.sleep(0.3)  # Polite delay
    
    logger.info(f"  -> Fetched metadata for {len(metadata)} out of {len(valid_pdbs)} valid PDB entries")
    return metadata


# Standalone entry point
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Stage 1: SIFTS Processing")
    parser.add_argument("--limit", type=int, default=None, help="Limit records for testing")
    parser.add_argument("--sifts-dir", type=Path, default=None, help="Custom SIFTS directory")
    args = parser.parse_args()
    
    stage1_process_sifts(limit=args.limit, sifts_dir=args.sifts_dir)