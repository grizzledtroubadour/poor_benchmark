#!/usr/bin/env python3
"""
diag_validate_structures.py (formerly 14_validate_structures.py; renamed as diagnostic in 2026-08).

Check downloaded monomer PDB completeness and parsability; diagnostic only, outside the main pipeline.

Checks:
1. Reject zero-byte files.
2. Parse using Bio.PDB.
3. Require at least one chain and some ATOM residues.
4. Write validation report and failed-file list.

Input:
- data/pdb_structures_all/*.pdb

Output:
- output/structure_validation_summary.json: validation summary
- output/structure_validation_failed.csv: parsing failures
"""

import os
import sys
import json
import glob
import csv
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio import PDB


ROOT = Path(__file__).resolve().parent.parent
PDB_DIR = ROOT / "data" / "pdb_structures_all"
OUTPUT_DIR = ROOT / "output"
OUTPUT_DIR.mkdir(exist_ok=True)

SUMMARY_FILE = OUTPUT_DIR / "structure_validation_summary.json"
FAILED_FILE = OUTPUT_DIR / "structure_validation_failed.csv"


def validate_one(pdb_path: str):
    """Return (pdb_path, status, error_msg, n_chains, n_residues)."""
    rel_path = os.path.relpath(pdb_path, PDB_DIR)
    try:
        size = os.path.getsize(pdb_path)
        if size == 0:
            return rel_path, "empty", "File size is 0", 0, 0

        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("model", pdb_path)

        n_chains = 0
        n_residues = 0
        has_atom = False
        for model in structure:
            for chain in model:
                n_chains += 1
                for residue in chain:
                    if residue.id[0].strip() == "":  # Count ATOM only; skip HETATM/water.
                        n_residues += 1
                        has_atom = True

        if not has_atom or n_residues == 0:
            return rel_path, "no_atoms", "No ATOM residues found", n_chains, n_residues

        return rel_path, "ok", "", n_chains, n_residues

    except Exception as e:
        return rel_path, "parse_error", str(e), 0, 0


def main():
    pdb_files = sorted(glob.glob(str(PDB_DIR / "*.pdb")))
    total = len(pdb_files)
    print(f"Found {total} PDB files under {PDB_DIR}")

    if total == 0:
        print("No PDB files to validate.")
        sys.exit(1)

    failed = []
    summary = {"total": total, "ok": 0, "empty": 0, "no_atoms": 0, "parse_error": 0}

    # Validate in parallel.
    max_workers = min(os.cpu_count() or 4, 16)
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(validate_one, f): f for f in pdb_files}
        for i, future in enumerate(as_completed(futures), 1):
            rel_path, status, error_msg, n_chains, n_residues = future.result()
            summary[status] += 1
            if status != "ok":
                failed.append({
                    "file": rel_path,
                    "status": status,
                    "error": error_msg,
                })
            if i % 5000 == 0 or i == total:
                print(f"  validated {i}/{total} ... ok={summary['ok']} failed={summary['total']-summary['ok']}")

    # Save summary.
    with open(SUMMARY_FILE, "w") as f:
        json.dump(summary, f, indent=2)

    # Save failed-file list.
    with open(FAILED_FILE, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["file", "status", "error"])
        writer.writeheader()
        writer.writerows(failed)

    print("\n=== Validation Summary ===")
    print(json.dumps(summary, indent=2))
    print(f"\nFailed list saved to: {FAILED_FILE}")
    print(f"Summary saved to: {SUMMARY_FILE}")


if __name__ == "__main__":
    main()
