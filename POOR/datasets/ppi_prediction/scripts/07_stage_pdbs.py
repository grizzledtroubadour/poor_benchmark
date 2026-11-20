#!/usr/bin/env python3
"""Stage selected monomer structures from data/pdb_structures_all/ into pdbs/.

Reconstruct a missing step (added 2026-08): pdbs/ was historically populated by copying
split-referenced files from data/pdb_structures_all/ without a standalone script. This formalizes it.

Required files =
            union of struct_file1 / struct_file2 in splits/ppi_{train,val,test}.csv
            (09_rebuild_with_xtal_negatives.py generates these fields from selected structures in
                  output/structure_selection.csv)
      union files referenced by output/splits_pre_xtal_backup/ (historical pre-XTAL splits)
            (included by default; disable with --no-pre-xtal). Existing pdbs/ has 46,402 files:
                  44,178 current references plus 2,224 historical pre-XTAL references.

Idempotent copy rules:
      - Destination exists: skip.
      - Source exists: copy.
      - Source missing: record in missing list. data/pdb_structures_all/ now retains only
            10,690 files; most historical files were moved elsewhere. For fresh reproduction, run
            06_download_structures_and_filter_test.py to redownload all 50,306 selected structures.

Usage:
  python 07_stage_pdbs.py                      # Idempotent staging into pdbs/.
  python 07_stage_pdbs.py --dry-run            # Report expected and missing files without writing.
  python 07_stage_pdbs.py --dest DIR [--limit N]  # Redirect destination for validation.
"""
import argparse
import json
import logging
import shutil
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

TASK_DIR = Path(__file__).resolve().parent.parent  # datasets/ppi_prediction
SPLITS_DIR = TASK_DIR / "splits"
PRE_XTAL_DIR = TASK_DIR / "output" / "splits_pre_xtal_backup"
SRC_DIR = TASK_DIR / "data" / "pdb_structures_all"
DEST_DIR = TASK_DIR / "pdbs"
SUMMARY_JSON = TASK_DIR / "output" / "stage_pdbs_summary.json"

STRUCT_COLS = ["struct_file1", "struct_file2"]


def collect_split_files(csv_path):
    df = pd.read_csv(csv_path, usecols=lambda c: c in STRUCT_COLS)
    files = set()
    for col in df.columns:
        files.update(f for f in df[col].dropna().astype(str) if f and f != "nan")
    return files


def collect_needed_files(include_pre_xtal):
    needed = set()
    for split in ["train", "val", "test"]:
        path = SPLITS_DIR / f"ppi_{split}.csv"
        files = collect_split_files(path)
        logger.info(f"{path.name}: {len(files)} unique structure files")
        needed.update(files)
    n_current = len(needed)
    if include_pre_xtal and PRE_XTAL_DIR.exists():
        for path in sorted(PRE_XTAL_DIR.glob("*.csv")):
            files = collect_split_files(path)
            logger.info(f"pre-XTAL {path.name}: {len(files)} unique structure files")
            needed.update(files)
        logger.info(
            f"needed total: {len(needed)} (current splits {n_current}, "
            f"pre-XTAL extra {len(needed) - n_current})"
        )
    return needed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--src", type=str, default=str(SRC_DIR), help="Source structure directory")
    parser.add_argument("--dest", type=str, default=str(DEST_DIR), help="Destination directory (default pdbs/)")
    parser.add_argument("--no-pre-xtal", action="store_true",
                        help="Exclude files referenced by historical pre-XTAL splits")
    parser.add_argument("--dry-run", action="store_true", help="Report counts only; do not copy")
    parser.add_argument("--limit", type=int, default=0, help="Copy at most N files for validation (0=unlimited)")
    args = parser.parse_args()

    src = Path(args.src)
    dest = Path(args.dest)
    needed = collect_needed_files(include_pre_xtal=not args.no_pre_xtal)

    existing = {f.name for f in dest.glob("*.pdb")} if dest.exists() else set()
    to_copy = sorted(needed - existing)
    missing_src = [f for f in to_copy if not (src / f).exists()]
    copyable = [f for f in to_copy if (src / f).exists()]
    if args.limit:
        copyable = copyable[: args.limit]

    logger.info(
        f"needed={len(needed)}, already in dest={len(needed & existing)}, "
        f"to copy={len(to_copy)} (src available {len(to_copy) - len(missing_src)}, "
        f"src missing {len(missing_src)})"
    )

    n_copied = 0
    if not args.dry_run:
        dest.mkdir(parents=True, exist_ok=True)
        for f in copyable:
            shutil.copy2(src / f, dest / f)
            n_copied += 1
        logger.info(f"copied {n_copied} files -> {dest}")
    else:
        logger.info("dry-run: no files written")

    summary = {
        "needed_files": len(needed),
        "already_in_dest": len(needed & existing),
        "copied": n_copied,
        "src_missing": len(missing_src),
        "src_missing_files": missing_src,
    }
    out = SUMMARY_JSON if dest == DEST_DIR else dest / "stage_pdbs_summary.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as fh:
        json.dump(summary, fh, indent=2)
    logger.info(f"summary -> {out}")


if __name__ == "__main__":
    main()
