#!/usr/bin/env python3
"""
Script: 05_extract_protein_sequences.py
Purpose: Extract protein sequences by chain from PDBbind v2020.R1 protein PDBs.
      Separate chains with `|` to align with the final whole-protein structures in pdbs/.
      Incorporates chain extraction from the former 17_reextract_sequences_from_pdbs.py.
      Deduplication by (chain, res_seq, res_name) and protein-only filtering ensure agreement
      between extraction from source {pdb_id}_protein.pdb and step-07-filtered pdbs/.

Usage:
    cd ligand_binding_affinity/
    python scripts/05_extract_protein_sequences.py

Input:
    output/pdbbind_v2020_metadata_ligand.csv

Output:
    output/pdbbind_v2020_metadata_seq.csv
    output/pdbbind_v2020_metadata_seq_stats.json
    output/sequence_extraction_failures.csv
"""

import argparse
import json
from pathlib import Path
from collections import Counter, OrderedDict

import pandas as pd

# Standard amino-acid mapping
AA_3TO1 = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
}

# Map common modified residues to standard amino acids.
MODIFIED_AA = {
    "MSE": "M",  # Selenomethionine
    "SEP": "S",  # Phosphoserine
    "TPO": "T",  # Phosphothreonine
    "PTR": "Y",  # Phosphotyrosine
    "KCX": "K",  # Carboxylated lysine
    "CSO": "C",  # Oxidized cysteine
    "CSD": "C",  # Cysteine sulfoxide
    "CME": "C",  # S,S-(2-hydroxyethyl)thiocysteine
    "LYZ": "K",  # Hydroxylysine
    "MLY": "K",  # Methyllysine
    "M3L": "K",  # Trimethyllysine
    "HYP": "P",  # Hydroxyproline
    "DAL": "A", "DVA": "V", "DLE": "L", "DIL": "I", "DTY": "Y",
    "DTR": "W", "DSE": "S", "DTH": "T", "DASN": "N", "DGLN": "Q",
    "DARG": "R", "DLY": "K", "DPR": "P", "DPH": "F", "DHI": "H",
    "DME": "M", "DAS": "D", "DGL": "E", "DCY": "C", "DGP": "G",
}


PROTEIN_RES = set(list(AA_3TO1.keys()) + list(MODIFIED_AA.keys()))


def extract_chain_sequences(pdb_path: Path):
    """
    Extract protein sequences from a PDB by chain.
    Preserve first-appearance chain order and the first occurrence of each residue
    keyed by (chain, res_seq, res_name); duplicates across NMR models collapse naturally.
    Return (chain_ids, chain_sequences), or (None, error_msg) on failure.
    """
    if not pdb_path.exists():
        return None, None, "PDB file not found"

    chains = OrderedDict()  # chain_id -> [aa, ...], preserving first appearance.
    seen = set()

    try:
        with open(pdb_path, "r", encoding="utf-8") as f:
            for line in f:
                if not (line.startswith("ATOM  ") or line.startswith("HETATM")):
                    continue
                if len(line) < 54:
                    continue

                res_name = line[17:20].strip()
                if res_name not in PROTEIN_RES:
                    continue

                chain_id = line[21].strip() if len(line) > 21 else ""
                res_seq_str = line[22:27].strip()  # Includes insertion codes.

                key = (chain_id, res_seq_str, res_name)
                if key in seen:
                    continue
                seen.add(key)

                aa = AA_3TO1.get(res_name) or MODIFIED_AA.get(res_name)
                if aa:
                    if chain_id not in chains:
                        chains[chain_id] = []
                    chains[chain_id].append(aa)
    except Exception as e:
        return None, None, f"PDB read error: {e}"

    if not chains or not any(chains.values()):
        return None, None, "No protein residues found"

    chain_ids = list(chains.keys())
    chain_seqs = ["".join(chains[ch]) for ch in chain_ids]
    return chain_ids, chain_seqs, None


def main():
    parser = argparse.ArgumentParser(description="Extract protein sequences for PDBbind v2020.R1")
    parser.add_argument(
        "--input",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "output" / "pdbbind_v2020_metadata_ligand.csv",
        help="Input metadata CSV with ligand features",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "output" / "pdbbind_v2020_metadata_seq.csv",
        help="Output CSV with protein sequences",
    )
    parser.add_argument(
        "--stats",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "output" / "pdbbind_v2020_metadata_seq_stats.json",
        help="Output stats JSON",
    )
    parser.add_argument(
        "--failures",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "output" / "sequence_extraction_failures.csv",
        help="Output failures CSV",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "data",
        help="Root data directory",
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=1,
        help="Minimum protein sequence length",
    )
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    # Process only samples with successfully parsed ligands.
    df = df[df["ligand_parse_success"] == True].copy()
    total = len(df)
    print(f"[Input] Successfully parsed ligand samples: {total}")

    results = []
    for _, row in df.iterrows():
        pdb_id = row["pdb_id"]
        protein_rel = row.get("protein_file", "")
        protein_path = args.data_dir / protein_rel if pd.notna(protein_rel) and protein_rel else None

        chain_ids, chain_seqs, error = extract_chain_sequences(protein_path)

        if error is None:
            # Separate chains with `|`; omit separators for single-chain proteins.
            seq = "|".join(chain_seqs)
            total_len = sum(len(s) for s in chain_seqs)
            success = total_len >= args.min_length
            results.append({
                "pdb_id": pdb_id,
                "protein_sequence_extract_success": success,
                "protein_sequence_error": None if success else "Sequence too short",
                "protein_sequence": seq if success else None,
                "protein_length": total_len,
                "protein_chains": ";".join(chain_ids),
                "num_protein_chains": len(chain_ids),
            })
        else:
            results.append({
                "pdb_id": pdb_id,
                "protein_sequence_extract_success": False,
                "protein_sequence_error": error,
                "protein_sequence": None,
                "protein_length": None,
                "protein_chains": None,
                "num_protein_chains": None,
            })

    seq_df = pd.DataFrame(results)
    merged = df.merge(seq_df, on="pdb_id", how="left")

    success_df = merged[merged["protein_sequence_extract_success"] == True].copy()
    failure_df = merged[merged["protein_sequence_extract_success"] != True].copy()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    success_df.to_csv(args.output, index=False)
    failure_df.to_csv(args.failures, index=False)

    stats = {
        "total_input": total,
        "sequence_extract_success": int(len(success_df)),
        "sequence_extract_failure": int(len(failure_df)),
        "retention_rate": round(len(success_df) / total, 4) if total > 0 else 0.0,
        "multi_chain_samples": int((success_df["num_protein_chains"] > 1).sum()),
        "failure_reasons": dict(Counter(failure_df["protein_sequence_error"].astype(str))),
        "protein_length_stats": success_df["protein_length"].describe().to_dict() if len(success_df) > 0 else {},
    }

    with open(args.stats, "w", encoding="utf-8") as f:
        json.dump(stats, f, indent=2)

    print(f"[Output] Successful sequence extractions: {len(success_df)}")
    print(f"[Output] Failed sequence extractions: {len(failure_df)}")
    print(f"[Output] Retention rate: {stats['retention_rate']*100:.2f}%")
    print(f"[Output] Successful samples CSV: {args.output}")
    print(f"[Output] Failed samples CSV: {args.failures}")
    print(f"[Output] Statistics JSON: {args.stats}")
    print()
    print("=== Failure reason distribution ===")
    for reason, count in stats["failure_reasons"].items():
        print(f"  {reason}: {count}")


if __name__ == "__main__":
    main()
