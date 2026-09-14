#!/usr/bin/env python3
"""
Fetch PDB entry-level metadata directly from RCSB PDB GraphQL API.
No dependency on func_prediction files.

Outputs:
    ss_prediction/output/pdb_metadata_api.csv
    ss_prediction/output/api_checkpoint.json
    ss_prediction/output/fetch_metadata.log

Usage:
    python ss_prediction/scripts/01_fetch_pdb_metadata_api.py
"""

import json
import csv
import time
import logging
from pathlib import Path
from typing import List, Dict, Optional
import requests
from tqdm import tqdm

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
TASK_DIR = Path(__file__).resolve().parent.parent  # datasets/ss_prediction
ROOT = TASK_DIR.parent.parent                      # 仓库根（POOR）
OUTPUT_DIR = TASK_DIR / "output"
OUTPUT_CSV = OUTPUT_DIR / "pdb_metadata_api.csv"
CHECKPOINT_FILE = OUTPUT_DIR / "api_checkpoint.json"
LOG_FILE = OUTPUT_DIR / "fetch_metadata.log"

API_URL = "https://data.rcsb.org/graphql"
BATCH_SIZE = 50
CHECKPOINT_INTERVAL = 500  # save checkpoint every N entries
REQUEST_DELAY = 0.3  # seconds between requests to be polite
MAX_RETRIES = 5

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(LOG_FILE, mode="w"),
        logging.StreamHandler(),
    ],
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# GraphQL query
# ---------------------------------------------------------------------------
GRAPHQL_QUERY = """
query($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    exptl {
      method
    }
    refine {
      ls_d_res_high
      ls_R_factor_R_work
      ls_R_factor_R_free
    }
    em_3d_reconstruction {
      resolution
    }
    rcsb_accession_info {
      deposit_date
    }
  }
}
"""


def parse_api_response(data: dict) -> List[Dict]:
    """Parse GraphQL response into flat records."""
    records = []
    entries = data.get("data", {}).get("entries", [])
    for entry in entries:
        if entry is None:
            continue
        pdb_id = entry.get("rcsb_id", "").upper()

        # Method — take first if multiple experiments
        exptl_list = entry.get("exptl") or []
        method = None
        if isinstance(exptl_list, list) and len(exptl_list) > 0:
            method = exptl_list[0].get("method")
        elif isinstance(exptl_list, dict):
            method = exptl_list.get("method")

        # Resolution: prefer X-ray, fallback to EM
        resolution = None
        refine = entry.get("refine")
        if isinstance(refine, list) and len(refine) > 0:
            resolution = refine[0].get("ls_d_res_high")
        elif isinstance(refine, dict):
            resolution = refine.get("ls_d_res_high")

        if resolution is None:
            em = entry.get("em_3d_reconstruction")
            if isinstance(em, list) and len(em) > 0:
                resolution = em[0].get("resolution")
            elif isinstance(em, dict):
                resolution = em.get("resolution")

        # R-factors
        r_work = None
        r_free = None
        if isinstance(refine, list) and len(refine) > 0:
            r_work = refine[0].get("ls_R_factor_R_work")
            r_free = refine[0].get("ls_R_factor_R_free")
        elif isinstance(refine, dict):
            r_work = refine.get("ls_R_factor_R_work")
            r_free = refine.get("ls_R_factor_R_free")

        # Deposition date
        deposit_date = None
        acc = entry.get("rcsb_accession_info")
        if isinstance(acc, dict):
            deposit_date = acc.get("deposit_date")

        records.append({
            "pdb_id": pdb_id,
            "method": method,
            "resolution": resolution,
            "r_work": r_work,
            "r_free": r_free,
            "deposition_date": deposit_date,
        })
    return records


def fetch_batch(ids: List[str]) -> Optional[List[Dict]]:
    """Query API for a batch of PDB IDs. Returns records or None on failure."""
    payload = {"query": GRAPHQL_QUERY, "variables": {"ids": ids}}
    for attempt in range(MAX_RETRIES):
        try:
            resp = requests.post(
                API_URL,
                json=payload,
                headers={"Content-Type": "application/json"},
                timeout=60,
            )
            if resp.status_code == 429:
                wait = 2 ** attempt
                logger.warning(f"Rate limited, sleeping {wait}s...")
                time.sleep(wait)
                continue
            resp.raise_for_status()
            data = resp.json()
            if "errors" in data:
                logger.warning(f"GraphQL errors: {data['errors']}")
                return None
            return parse_api_response(data)
        except requests.exceptions.RequestException as e:
            logger.warning(f"Request failed (attempt {attempt + 1}/{MAX_RETRIES}): {e}")
            time.sleep(2 ** attempt)
    logger.error(f"Failed to fetch batch after {MAX_RETRIES} retries: {ids[:5]}...")
    return None


def load_checkpoint() -> int:
    if CHECKPOINT_FILE.exists():
        with open(CHECKPOINT_FILE) as f:
            ck = json.load(f)
        return ck.get("processed", 0)
    return 0


def save_checkpoint(processed: int):
    with open(CHECKPOINT_FILE, "w") as f:
        json.dump({"processed": processed}, f)


def main():
    # Gather all local PDB IDs from filenames
    pdb_dir = ROOT / "data" / "pdb"
    all_ids = sorted(
        set(
            f.name.split(".")[0].upper()
            for f in pdb_dir.iterdir()
            if f.suffix == ".gz"
        )
    )
    logger.info(f"Total local PDB entries: {len(all_ids)}")

    # Load existing results if any
    existing_records = []
    if OUTPUT_CSV.exists():
        with open(OUTPUT_CSV, newline="") as f:
            reader = csv.DictReader(f)
            existing_records = list(reader)
        logger.info(f"Loaded {len(existing_records)} existing records from {OUTPUT_CSV}")

    processed = load_checkpoint()
    logger.info(f"Resuming from checkpoint: {processed} already processed")

    records = existing_records.copy()
    failed_batches = []

    # Slice remaining IDs
    remaining_ids = all_ids[processed:]
    total_batches = (len(remaining_ids) + BATCH_SIZE - 1) // BATCH_SIZE

    pbar = tqdm(total=len(remaining_ids), desc="Fetching metadata", unit="entry")
    for i in range(total_batches):
        batch_ids = remaining_ids[i * BATCH_SIZE : (i + 1) * BATCH_SIZE]
        batch_records = fetch_batch(batch_ids)

        if batch_records is not None:
            records.extend(batch_records)
            processed += len(batch_ids)
        else:
            failed_batches.extend(batch_ids)
            processed += len(batch_ids)  # mark as processed to avoid infinite loop

        pbar.update(len(batch_ids))

        # Checkpoint
        if processed % CHECKPOINT_INTERVAL < BATCH_SIZE or i == total_batches - 1:
            save_checkpoint(processed)
            # Write incremental CSV
            with open(OUTPUT_CSV, "w", newline="") as f:
                writer = csv.DictWriter(
                    f,
                    fieldnames=["pdb_id", "method", "resolution", "r_work", "r_free", "deposition_date"],
                )
                writer.writeheader()
                writer.writerows(records)
            logger.info(f"Checkpoint saved: {processed}/{len(all_ids)} entries")

        time.sleep(REQUEST_DELAY)

    pbar.close()

    # Final save
    with open(OUTPUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["pdb_id", "method", "resolution", "r_work", "r_free", "deposition_date"],
        )
        writer.writeheader()
        writer.writerows(records)

    logger.info(f"Done. Total records: {len(records)}")
    if failed_batches:
        failed_path = OUTPUT_DIR / "failed_batches.txt"
        with open(failed_path, "w") as f:
            f.write("\n".join(failed_batches))
        logger.info(f"Failed batches saved to {failed_path}: {len(failed_batches)} entries")

    # Quick stats
    import pandas as pd

    df = pd.DataFrame(records)
    logger.info(f"Non-null resolution: {df['resolution'].notna().sum()}")
    logger.info(f"Non-null R-free: {df['r_free'].notna().sum()}")
    logger.info(f"Method dist:\n{df['method'].value_counts().head(10)}")


if __name__ == "__main__":
    main()
