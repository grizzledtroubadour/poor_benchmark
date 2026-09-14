#!/usr/bin/env python3
"""Download missing AlphaFold structures for kcat task via AlphaFold DB API."""
import json
import gzip
import requests
import time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

TASK_DIR = Path(__file__).resolve().parent.parent  # datasets/enzyme_kinetics_prediction
PROJECT_ROOT = TASK_DIR.parent.parent              # POOD repo root

MISSING_FILE = TASK_DIR / "data/kcat_missing_alphafold_uniprots.txt"
OUT_DIR = PROJECT_ROOT / "data/alphafold_structures"
LOG_FILE = TASK_DIR / "data/alphafold_download_progress.log"

OUT_DIR.mkdir(parents=True, exist_ok=True)

API_URL = "https://alphafold.ebi.ac.uk/api/prediction/{uniprot}"

# Read missing uniprots
with open(MISSING_FILE) as f:
    uniprots = [line.strip() for line in f if line.strip()]

print(f"Total missing UniProts to download: {len(uniprots)}")

# Load progress
if LOG_FILE.exists():
    with open(LOG_FILE) as f:
        done = {line.strip().split()[0] for line in f if line.strip() and not line.startswith("#")}
else:
    done = set()

print(f"Already downloaded: {len(done)}")
remaining = [u for u in uniprots if u not in done]
print(f"Remaining: {len(remaining)}")

session = requests.Session()

def download_one(uniprot):
    out_path = OUT_DIR / f"AF-{uniprot}-F1-model_v6.pdb.gz"
    if out_path.exists():
        return uniprot, "exists"
    
    try:
        # Query API
        resp = session.get(API_URL.format(uniprot=uniprot), timeout=30)
        if resp.status_code == 404:
            return uniprot, "api_404"
        resp.raise_for_status()
        data = resp.json()
        if not data:
            return uniprot, "api_empty"
        
        entry = data[0]
        pdb_url = entry.get("pdbUrl")
        if not pdb_url:
            return uniprot, "no_pdb_url"
        
        # Download PDB
        pdb_resp = session.get(pdb_url, timeout=60)
        pdb_resp.raise_for_status()
        
        # Gzip and save
        with gzip.open(out_path, "wt") as f:
            f.write(pdb_resp.text)
        
        return uniprot, "success"
    except Exception as e:
        return uniprot, f"error:{e}"

# Use moderate concurrency to avoid rate limiting
MAX_WORKERS = 8
success = 0
failed = 0
not_found = 0
already = 0

with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
    futures = {executor.submit(download_one, u): u for u in remaining}
    for future in as_completed(futures):
        uniprot, status = future.result()
        with open(LOG_FILE, "a") as f:
            f.write(f"{uniprot}\t{status}\t{time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        
        if status == "success":
            success += 1
        elif status == "exists":
            already += 1
        elif status in ("api_404", "api_empty", "no_pdb_url"):
            not_found += 1
        else:
            failed += 1
        
        total_processed = success + failed + not_found + already
        if total_processed % 50 == 0:
            print(f"Progress: {total_processed}/{len(remaining)} success={success} already={already} not_found={not_found} failed={failed}")

print(f"\nDone. success={success}, already={already}, not_found={not_found}, failed={failed}")
