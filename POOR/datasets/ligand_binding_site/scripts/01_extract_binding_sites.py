#!/usr/bin/env python3
"""
Extract protein-ligand binding sites from PDB/mmCIF structure files.

See extract_binding_sites_core.py for the core processing logic.
"""

import argparse
import logging
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

from extract_binding_sites_core import process_file_wrapper


def main():
    parser = argparse.ArgumentParser(description="Extract protein-ligand binding sites from PDB/mmCIF files.")
    parser.add_argument("--input_dir", type=str, required=True)
    parser.add_argument("--output_dir", type=str, required=True)
    parser.add_argument("--distance_threshold", type=float, default=4.5)
    parser.add_argument("--min_ligand_atoms", type=int, default=5)
    parser.add_argument("--min_binding_residues", type=int, default=3)
    parser.add_argument("--num_workers", type=int, default=16)
    parser.add_argument("--batch_size", type=int, default=5000)
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    log_path = output_dir / "01_extract_binding_sites.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_path, mode="w"),
            logging.StreamHandler(sys.stdout),
        ],
    )
    logger = logging.getLogger(__name__)

    from extract_binding_sites_core import (
        AA_NAMES, MODIFIED_AA_1LETTER, CRYSTALLIZATION_ADDITIVE_NAMES
    )

    logger.info("Starting binding site extraction")
    logger.info(f"Input directory: {input_dir}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Distance threshold: {args.distance_threshold} Å")
    logger.info(f"Min ligand heavy atoms: {args.min_ligand_atoms}")
    logger.info(f"Min binding residues: > {args.min_binding_residues}")
    logger.info(f"Crystallization additives filtered: {len(CRYSTALLIZATION_ADDITIVE_NAMES)} HET codes")
    logger.info(f"Modified amino acids mapped to standard AAs: {len(MODIFIED_AA_1LETTER)} codes")
    logger.info(f"Number of workers: {args.num_workers}")

    structure_files = sorted(input_dir.glob("*.pdb.gz")) + sorted(input_dir.glob("*.cif.gz"))
    if not structure_files:
        structure_files = sorted(input_dir.glob("*.pdb")) + sorted(input_dir.glob("*.cif"))
    logger.info(f"Found {len(structure_files)} structure files")

    processed_log = output_dir / "processed_files.txt"
    processed_set = set()
    if processed_log.exists():
        with open(processed_log, "r") as f:
            processed_set = set(line.strip() for line in f if line.strip())
        logger.info(f"Resuming: {len(processed_set)} files already processed")

    remaining_files = [f for f in structure_files if str(f) not in processed_set]
    logger.info(f"Remaining files to process: {len(remaining_files)}")

    work_args = [
        (str(f), args.distance_threshold, args.min_ligand_atoms, args.min_binding_residues)
        for f in remaining_files
    ]

    all_records = []
    failed_records = []
    processed_since_last_batch = []
    stats = {
        "total_files": len(structure_files),
        "processed_files": len(processed_set),
        "successful_files": 0,
        "failed_files": 0,
        "total_chain_ligand_pairs": 0,
        "total_protein_chains": 0,
        "total_ligands": 0,
    }

    existing_batches = sorted(output_dir.glob("binding_sites_batch_*.csv"))
    batch_counter = len(existing_batches)
    logger.info(f"Found {len(existing_batches)} existing batch files; new batches start at {batch_counter:05d}")

    def _save_batch():
        nonlocal all_records, processed_since_last_batch, batch_counter
        if not all_records:
            return
        batch_df = pd.DataFrame(all_records)
        batch_path = output_dir / f"binding_sites_batch_{batch_counter:05d}.csv"
        batch_df.to_csv(batch_path, index=False)
        with open(processed_log, "a") as proc_fh:
            for fp in processed_since_last_batch:
                proc_fh.write(f"{fp}\n")
        logger.info(
            f"Saved batch {batch_counter}: {len(batch_df)} records from "
            f"{len(processed_since_last_batch)} files to {batch_path.name}"
        )
        all_records = []
        processed_since_last_batch = []
        batch_counter += 1

    start_time = time.time()

    with open(processed_log, "a") as proc_fh, open(output_dir / "failed_files.txt", "a") as fail_fh:
        with ProcessPoolExecutor(max_workers=args.num_workers) as executor:
            futures = {executor.submit(process_file_wrapper, wa): wa for wa in work_args}

            for future in as_completed(futures):
                filepath, result = future.result()
                stats["processed_files"] += 1
                processed_since_last_batch.append(filepath)

                if result["error"] is not None:
                    failed_records.append({"filepath": filepath, "error": result["error"]})
                    fail_fh.write(f"{filepath}\t{result['error']}\n")
                    fail_fh.flush()
                    stats["failed_files"] += 1
                else:
                    stats["successful_files"] += 1
                    stats["total_chain_ligand_pairs"] += len(result["records"])
                    stats["total_protein_chains"] += result["n_chains"]
                    stats["total_ligands"] += result["n_ligands"]
                    all_records.extend(result["records"])

                if len(all_records) >= args.batch_size:
                    _save_batch()

                if stats["processed_files"] % 1000 == 0:
                    elapsed = time.time() - start_time
                    rate = stats["processed_files"] / elapsed if elapsed > 0 else 0
                    logger.info(
                        f"Progress: {stats['processed_files']}/{stats['total_files']} "
                        f"({100*stats['processed_files']/stats['total_files']:.1f}%) | "
                        f"Successful: {stats['successful_files']} | Failed: {stats['failed_files']} | "
                        f"Pairs: {stats['total_chain_ligand_pairs']} | "
                        f"Rate: {rate:.1f} files/s"
                    )

        if all_records or processed_since_last_batch:
            _save_batch()

    if all_records:
        batch_df = pd.DataFrame(all_records)
        batch_path = output_dir / f"binding_sites_batch_{batch_counter:05d}.csv"
        batch_df.to_csv(batch_path, index=False)
        logger.info(f"Saved final batch {batch_counter}: {len(batch_df)} records to {batch_path.name}")
        batch_counter += 1

    batch_files = sorted(output_dir.glob("binding_sites_batch_*.csv"))
    if batch_files:
        logger.info(f"Aggregating {len(batch_files)} batch files into final outputs...")
        df = pd.concat([pd.read_csv(f) for f in batch_files], ignore_index=True)

        main_output = output_dir / "ligand_binding_sites_raw.csv"
        df.to_csv(main_output, index=False)
        logger.info(f"Saved raw binding sites to {main_output}: {len(df)} rows")

        std_cols = ["unique_id", "file_type", "aa_seq", "labels", "pdb_id", "chain_id",
                    "ligand_id", "ligand_name", "ligand_chain", "ligand_resnum",
                    "num_binding_residues", "num_ligand_heavy_atoms", "protein_length"]
        std_df = df[std_cols].copy()
        std_output = output_dir / "ligand_binding_sites_std.csv"
        std_df.to_csv(std_output, index=False)
        logger.info(f"Saved standard-format binding sites to {std_output}: {len(std_df)} rows")

    if failed_records:
        failed_df = pd.DataFrame(failed_records)
        failed_output = output_dir / "failed_files.csv"
        failed_df.to_csv(failed_output, index=False)
        logger.info(f"Saved {len(failed_df)} failed files to {failed_output}")

    import json
    stats["elapsed_seconds"] = time.time() - start_time
    stats_path = output_dir / "extraction_stats.json"
    with open(stats_path, "w") as f:
        json.dump(stats, f, indent=2)
    logger.info(f"Saved statistics to {stats_path}")
    logger.info("Binding site extraction completed")
    logger.info(f"Final stats: {json.dumps(stats, indent=2)}")


if __name__ == "__main__":
    main()
