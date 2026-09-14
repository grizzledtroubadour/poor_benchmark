#!/usr/bin/env python3
"""
Stage 1.5: Entry-level quality filtering from PDB metadata.

Inputs:
    ss_prediction/output/pdb_metadata_api.csv

Outputs:
    ss_prediction/output/stage1_filtered_entries.csv
    ss_prediction/output/stage1_filtering_stats.json

Filter criteria:
    - method == 'X-RAY DIFFRACTION'
    - resolution <= 2.5
    - r_free <= 0.25
"""

import json
import logging
import pandas as pd
from pathlib import Path

TASK_DIR = Path(__file__).resolve().parent.parent  # datasets/ss_prediction
OUTPUT_DIR = TASK_DIR / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger(__name__)


def main():
    meta_path = OUTPUT_DIR / "pdb_metadata_api.csv"
    if not meta_path.exists():
        logger.error(f"Metadata not found: {meta_path}")
        logger.error("Please run 01_fetch_pdb_metadata_api.py first.")
        return

    df = pd.read_csv(meta_path)
    logger.info(f"Loaded {len(df)} entries from {meta_path}")

    # Convert numeric columns
    df["resolution"] = pd.to_numeric(df["resolution"], errors="coerce")
    df["r_free"] = pd.to_numeric(df["r_free"], errors="coerce")

    total = len(df)

    # Step-wise filtering stats
    stats = {"total_entries": total}

    # Filter 1: X-ray only
    df_xray = df[df["method"] == "X-RAY DIFFRACTION"].copy()
    stats["after_xray_filter"] = len(df_xray)

    # Filter 2: resolution <= 2.5
    df_res = df_xray[df_xray["resolution"] <= 2.5].copy()
    stats["after_resolution_filter"] = len(df_res)

    # Filter 3: R-free <= 0.25
    df_final = df_res[df_res["r_free"] <= 0.25].copy()
    stats["after_rfree_filter"] = len(df_final)

    # Distribution stats
    stats["resolution_mean"] = float(df_final["resolution"].mean())
    stats["resolution_median"] = float(df_final["resolution"].median())
    stats["r_free_mean"] = float(df_final["r_free"].mean())
    stats["r_free_median"] = float(df_final["r_free"].median())
    stats["deposit_date_earliest"] = str(df_final["deposition_date"].min())
    stats["deposit_date_latest"] = str(df_final["deposition_date"].max())

    # Save filtered entries
    out_csv = OUTPUT_DIR / "stage1_filtered_entries.csv"
    df_final[["pdb_id", "method", "resolution", "r_work", "r_free", "deposition_date"]].to_csv(
        out_csv, index=False
    )
    logger.info(f"Saved {len(df_final)} filtered entries to {out_csv}")

    # Save stats
    out_stats = OUTPUT_DIR / "stage1_filtering_stats.json"
    with open(out_stats, "w") as f:
        json.dump(stats, f, indent=2)
    logger.info(f"Saved stats to {out_stats}")

    # Print summary
    print("\n=== Filtering Summary ===")
    for k, v in stats.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
