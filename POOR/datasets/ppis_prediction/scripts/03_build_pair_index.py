#!/usr/bin/env python3
"""
Build an index of DIPS-Plus pair files.

Scan dill/HDF5 pairs under final/raw/ and extract key metadata;
also try to read official split lists (pairs-postprocessed-*.txt).

Usage:
    python scripts/build_pair_index.py
"""

import csv
import json
import logging
import re
import sys
from pathlib import Path

Path("logs").mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("logs/build_pair_index.log"),
        logging.StreamHandler(sys.stdout),
    ],
)
logger = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data"
OUTPUT_DIR = ROOT / "output"

DIPS_RAW_DIR_CANDIDATES = [
    DATA_DIR / "raw",
    DATA_DIR / "project" / "datasets" / "DIPS" / "final" / "raw",
    DATA_DIR / "DIPS" / "final" / "raw",
    DATA_DIR,
]


def find_raw_dir():
    """Locate final/raw/."""
    for cand in DIPS_RAW_DIR_CANDIDATES:
        if cand.exists() and any(cand.iterdir()):
            logger.info("Located raw pair directory: %s", cand)
            return cand
    raise FileNotFoundError("Could not locate DIPS final/raw directory")


def extract_pdb_info(filename):
    """Infer PDB ID and model information from a pair filename.

    DIPS-Plus v1.3.0 filenames have the form:
        1ad1.pdb1_0.dill
        2ad7.pdb1_2.hdf5
    Here 1ad1 is the PDB ID, pdb1 is the model, and _0/_2 is the pair index.
    """
    stem = Path(filename).stem
    # Match four-character PDB ID + .pdb + model + _ + pair_index.
    m = re.match(r"([0-9][a-z0-9]{3})\.pdb(\d+)_(\d+)", stem, re.IGNORECASE)
    if m:
        return m.group(1).lower(), m.group(2), m.group(3)
    # Fallback: match the PDB ID only.
    m = re.match(r"([0-9][a-z0-9]{3})", stem, re.IGNORECASE)
    if m:
        return m.group(1).lower(), "", ""
    return None, "", ""


def scan_pair_files(raw_dir):
    """Scan pair files, preferring HDF5."""
    records = []
    seen = set()
    # Scan HDF5 first so it wins when dill and HDF5 share the same name.
    for fmt in ("*.hdf5", "*.h5", "*.dill"):
        for path in raw_dir.rglob(fmt):
            rel = path.relative_to(ROOT)
            unique_id = path.stem
            if unique_id in seen:
                continue
            seen.add(unique_id)
            pdb_id, model_id, pair_index = extract_pdb_info(path.name)
            records.append({
                "unique_id": unique_id,
                "pdb_id": pdb_id or "",
                "model_id": model_id,
                "pair_index": pair_index,
                "file_name": path.name,
                "file_path": str(rel),
                "file_format": path.suffix.lstrip("."),
            })
    return records


def load_official_split_lists(raw_dir):
    """Read official split-list files."""
    splits = {}
    split_files = {
        "train": "pairs-postprocessed-train.txt",
        "val": "pairs-postprocessed-val.txt",
        "test": "pairs-postprocessed-test.txt",
        "train_before_filter": "pairs-postprocessed-train-before-structure-based-filtering.txt",
        "val_before_filter": "pairs-postprocessed-val-before-structure-based-filtering.txt",
    }
    for split_name, filename in split_files.items():
        path = raw_dir / filename
        if path.exists():
            logger.info("Found official split list: %s", path)
            splits[split_name] = {line.strip() for line in path.read_text().splitlines() if line.strip()}
    return splits


def assign_splits(records, splits):
    """Merge official split labels into the index.

    Prefer final filtered train/val lists; before-filter lists provide supplementary information.
    """
    split_map = {}
    # Read before-filter lists first so final lists override them.
    order = ["train_before_filter", "val_before_filter", "test", "train", "val"]
    for split_name in order:
        if split_name not in splits:
            continue
        for fn in splits[split_name]:
            base = Path(fn).stem
            name = Path(fn).name
            split_map[base] = split_name
            split_map[name] = split_name
            split_map[fn] = split_name

    for rec in records:
        rec["official_split"] = split_map.get(rec["unique_id"], "")
        if not rec["official_split"]:
            rec["official_split"] = split_map.get(rec["file_name"], "")
        # Support legacy extensionless lists.
        if not rec["official_split"]:
            rec["official_split"] = split_map.get(rec["unique_id"].replace(".dill", "").replace(".hdf5", "").replace(".h5", ""), "")
    return records


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    raw_dir = find_raw_dir()
    records = scan_pair_files(raw_dir)
    logger.info("Scanned %d pair files", len(records))

    if not records:
        logger.warning("No pair files found. Please check extraction path: %s", raw_dir)
        return

    splits = load_official_split_lists(raw_dir)
    records = assign_splits(records, splits)

    # Save the index.
    index_path = OUTPUT_DIR / "dips_plus_pair_index.csv"
    with open(index_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=records[0].keys())
        writer.writeheader()
        writer.writerows(records)
    logger.info("Pair index saved to %s", index_path)

    # Split index
    split_records = [
        {"unique_id": r["unique_id"], "file_path": r["file_path"], "split": r["official_split"]}
        for r in records if r["official_split"]
    ]
    if split_records:
        split_path = OUTPUT_DIR / "dips_plus_split_index.csv"
        with open(split_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=split_records[0].keys())
            writer.writeheader()
            writer.writerows(split_records)
        logger.info("Split index saved to %s", split_path)

    # Summary statistics
    fmt_counts = {}
    split_counts = {}
    for r in records:
        fmt_counts[r["file_format"]] = fmt_counts.get(r["file_format"], 0) + 1
        sp = r["official_split"] or "unassigned"
        split_counts[sp] = split_counts.get(sp, 0) + 1

    summary = {
        "total_unique_complexes": len(records),
        "by_format": fmt_counts,
        "by_official_split": split_counts,
        "raw_dir": str(raw_dir),
    }
    summary_path = OUTPUT_DIR / "dips_plus_index_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    logger.info("Summary: %s", summary)


if __name__ == "__main__":
    main()
