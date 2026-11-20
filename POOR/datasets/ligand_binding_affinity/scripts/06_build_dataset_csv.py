#!/usr/bin/env python3
"""
Script: 06_build_dataset_csv.py
Purpose: Build output/pdbbind_v2020_dataset.csv from step-05 sequence outputs
      as input to the Stage 2 splitting and OOD workflow.

      Replaces the archived chain-structure route in 08_extract_protein_chains.py:
      no separate chain PDBs are written (07_copy_whole_protein_pdbs.py writes whole proteins).
      Keep samples with successful sequence extraction and add structure paths.
      chain_extraction_success / chain_extraction_error are retained for legacy downstream
      consumers (including 08_cluster_time_split.py), matching sequence-extraction success.

Usage:
    cd ligand_binding_affinity/
    python scripts/06_build_dataset_csv.py

Input:
    output/pdbbind_v2020_metadata_seq.csv

Output:
    output/pdbbind_v2020_dataset.csv
    output/protein_chain_extraction_stats.json
"""

import argparse
import json
from pathlib import Path

import pandas as pd

# Preserve the historical dataset.csv column order (see docstring).
FINAL_COLUMNS = [
    "pdb_id", "resolution", "release_year", "affinity_type", "affinity_operator",
    "affinity_value", "pK", "ligand_name", "reference", "year_range",
    "protein_file", "ligand_mol2_file", "ligand_sdf_file", "pocket_file",
    "incomplete_ligand", "covalent_complex", "different_protein", "different_ligand",
    "uncommon_element", "disulfide_bond", "comment", "raw",
    "ligand_parse_success", "ligand_parse_error", "ligand_smiles", "ligand_ecfp4",
    "ligand_num_heavy_atoms",
    "protein_sequence_extract_success", "protein_sequence_error",
    "protein_sequence", "protein_length",
    "chain_extraction_success", "chain_extraction_error",
    "protein_chains", "num_protein_chains", "protein_pdb_file",
]


def main():
    parser = argparse.ArgumentParser(description="Build sample-level dataset CSV from sequence metadata")
    parser.add_argument(
        "--input",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "output" / "pdbbind_v2020_metadata_seq.csv",
        help="Input metadata CSV with per-chain protein sequences (output of 05)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "output" / "pdbbind_v2020_dataset.csv",
        help="Output sample-level dataset CSV",
    )
    parser.add_argument(
        "--stats",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "output" / "protein_chain_extraction_stats.json",
        help="Output stats JSON",
    )
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    total = len(df)
    print(f"[Input] Samples after cleaning and ligand/sequence parsing: {total}")

    # Keep successfully extracted sequences in the final dataset.
    final_df = df[df["protein_sequence_extract_success"] == True].copy()
    failure_df = df[df["protein_sequence_extract_success"] != True].copy()

    # Compatibility flag for the historical chain route; matches sequence extraction.
    final_df["chain_extraction_success"] = final_df["protein_sequence_extract_success"]
    final_df["chain_extraction_error"] = final_df["protein_sequence_error"]

    # Whole-protein structure path (written by step 07).
    final_df["protein_pdb_file"] = final_df["pdb_id"].apply(lambda x: f"pdbs/{x}.pdb")

    final_df = final_df[FINAL_COLUMNS]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    final_df.to_csv(args.output, index=False)
    if not failure_df.empty:
        failure_path = args.output.parent / "protein_chain_extraction_failures.csv"
        failure_df.to_csv(failure_path, index=False)
        print(f"[Output] Failed extraction CSV: {failure_path}")

    stats = {
        "total_input": total,
        "chain_extraction_success": int(len(final_df)),
        "chain_extraction_failure": int(len(failure_df)),
        "retention_rate": round(len(final_df) / total, 4) if total > 0 else 0.0,
        "total_chain_files": 0,  # Chain-level PDB export is obsolete; see step 07.
        "num_chains_distribution": {
            "min": int(final_df["num_protein_chains"].min()) if len(final_df) else 0,
            "max": int(final_df["num_protein_chains"].max()) if len(final_df) else 0,
            "mean": round(float(final_df["num_protein_chains"].mean()), 2) if len(final_df) else 0.0,
        },
        "note": "Chain-level PDB export is obsolete (08_extract_protein_chains.py is archived); "
                "chain_extraction_* columns are retained for compatibility; step 07 writes whole proteins.",
    }

    with open(args.stats, "w", encoding="utf-8") as f:
        json.dump(stats, f, indent=2)

    print(f"[Output] Retained samples: {len(final_df)}")
    print(f"[Output] Excluded samples: {len(failure_df)}")
    print(f"[Output] Mean chains per sample: {stats['num_chains_distribution']['mean']}")
    print(f"[Output] Final dataset CSV: {args.output}")
    print(f"[Output] Statistics JSON: {args.stats}")


if __name__ == "__main__":
    main()
