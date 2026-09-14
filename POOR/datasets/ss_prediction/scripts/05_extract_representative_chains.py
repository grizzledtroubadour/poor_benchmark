#!/usr/bin/env python3
"""
Extract single-chain PDB files for representative chains.
Handles both .pdb.gz and .cif.gz source files.
Adds minimal headers required for DSSP v4 compatibility.
"""

import csv
import gzip
import os
import sys
import tempfile
from multiprocessing import Process, Queue
from pathlib import Path

from Bio.PDB import MMCIFParser, PDBIO, PDBParser, Select

# Config
TASK_DIR = Path(__file__).resolve().parent.parent  # datasets/ss_prediction
ROOT = TASK_DIR.parent.parent                      # 仓库根（POOR）
RAW_PDB_DIR = ROOT / "data" / "pdb"
OUTPUT_DIR = TASK_DIR / "pdbs"
REPRESENTATIVES_CSV = TASK_DIR / "output" / "stage4_cluster_representatives.csv"
NUM_WORKERS = min(16, os.cpu_count() or 4)
CHUNK_SIZE = 500  # Each worker processes this many entries before exiting

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


class ChainSelect(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        return residue.id[0] == " "


def fix_pdb_headers(input_path, output_path, pdb_id, chain_id):
    """Add minimal required PDB records for DSSP v4."""
    with open(input_path) as f:
        lines = f.readlines()

    has_header = any(l.startswith("HEADER") for l in lines)
    has_cryst1 = any(l.startswith("CRYST1") for l in lines)
    has_end = any(l.startswith("END") for l in lines)

    with open(output_path, "w") as f:
        if not has_header:
            f.write(
                f"HEADER    PROTEIN                                 01-JAN-00   {pdb_id.upper():>4s}\n"
            )
            f.write("TITLE     SINGLE CHAIN EXTRACT\n")
            f.write("EXPDTA    X-RAY DIFFRACTION\n")
        if not has_cryst1:
            f.write(
                "CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1           1\n"
            )

        for line in lines:
            if (line.startswith("ATOM") or line.startswith("HETATM")) and len(line) > 21 and line[21] == " ":
                line = line[:21] + chain_id + line[22:]
            f.write(line)

        if not has_end:
            f.write("END\n")


def extract_one(row):
    pdb_id = row["pdb_id"].lower()
    chain_id = row["chain_id"]
    file_type = row["file_type"]
    out_file = OUTPUT_DIR / f"{pdb_id.upper()}_{chain_id}.pdb"

    if out_file.exists() and out_file.stat().st_size > 0:
        return (pdb_id, chain_id, "skipped")

    if file_type == "pdb":
        src_path = RAW_PDB_DIR / f"{pdb_id.upper()}.pdb.gz"
        parser = PDBParser(QUIET=True)
    elif file_type == "cif":
        src_path = RAW_PDB_DIR / f"{pdb_id.upper()}.cif.gz"
        parser = MMCIFParser(QUIET=True)
    else:
        return (pdb_id, chain_id, f"unknown_type:{file_type}")

    if not src_path.exists():
        alt = RAW_PDB_DIR / f"{pdb_id}.pdb.gz"
        if alt.exists():
            src_path = alt
            parser = PDBParser(QUIET=True)
        else:
            alt2 = RAW_PDB_DIR / f"{pdb_id}.cif.gz"
            if alt2.exists():
                src_path = alt2
                parser = MMCIFParser(QUIET=True)
            else:
                return (pdb_id, chain_id, f"missing:{src_path}")

    try:
        with gzip.open(src_path, "rt") as handle:
            structure = parser.get_structure(pdb_id.upper(), handle)

        model = structure[0]
        if chain_id not in [c.id for c in model]:
            return (pdb_id, chain_id, "chain_not_found")

        io = PDBIO()
        # Keep model 1 only: passing the whole Structure writes every model of
        # NMR entries (with MODEL/ENDMDL records), which downstream pipelines
        # do not expect (fixed post-hoc by output/scripts_archive/fix_multimodel_pdbs.py).
        io.set_structure(structure[0])

        with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as tmp:
            tmp_path = tmp.name

        io.save(tmp_path, ChainSelect(chain_id))
        fix_pdb_headers(tmp_path, str(out_file), pdb_id.upper(), chain_id)
        os.unlink(tmp_path)

        atom_count = 0
        with open(out_file) as f:
            for line in f:
                if line.startswith("ATOM"):
                    atom_count += 1
        if atom_count == 0:
            out_file.unlink(missing_ok=True)
            return (pdb_id, chain_id, "no_atoms")

        return (pdb_id, chain_id, "success")

    except Exception as e:
        return (pdb_id, chain_id, f"error:{type(e).__name__}:{e}")


def worker_process(rows_chunk, result_queue):
    for row in rows_chunk:
        result_queue.put(extract_one(row))
    result_queue.put(None)  # Sentinel


def main():
    rows = []
    with open(REPRESENTATIVES_CSV, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)

    print(f"Total chains to extract: {len(rows)}")
    print(f"Using {NUM_WORKERS} workers, chunk size {CHUNK_SIZE}")

    # Split into chunks
    chunks = []
    for i in range(0, len(rows), CHUNK_SIZE):
        chunks.append(rows[i : i + CHUNK_SIZE])

    success = 0
    skipped = 0
    failed = 0
    errors = []

    for batch_start in range(0, len(chunks), NUM_WORKERS):
        batch = chunks[batch_start : batch_start + NUM_WORKERS]
        result_queue = Queue()
        processes = []

        for chunk in batch:
            p = Process(target=worker_process, args=(chunk, result_queue))
            p.start()
            processes.append(p)

        finished_workers = 0
        total_expected = sum(len(c) for c in batch)
        processed = 0

        while finished_workers < len(processes):
            result = result_queue.get()
            if result is None:
                finished_workers += 1
                continue

            processed += 1
            pdb_id, chain_id, status = result
            if status == "success":
                success += 1
            elif status == "skipped":
                skipped += 1
            else:
                failed += 1
                errors.append(result)

        for p in processes:
            p.join()

        total_done = success + skipped + failed
        print(
            f"Batch {batch_start+1}-{min(batch_start+NUM_WORKERS, len(chunks))}/{len(chunks)} | "
            f"total_done={total_done} success={success} skipped={skipped} failed={failed}",
            flush=True,
        )

    if errors:
        err_path = OUTPUT_DIR / "../output/extraction_errors.csv"
        err_path.parent.mkdir(parents=True, exist_ok=True)
        with open(err_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["pdb_id", "chain_id", "status"])
            writer.writerows(errors)
        print(f"Wrote {len(errors)} errors to {err_path}")

    print(f"\nDone: success={success}, skipped={skipped}, failed={failed}")


if __name__ == "__main__":
    main()
