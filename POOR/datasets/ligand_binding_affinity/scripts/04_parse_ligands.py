#!/usr/bin/env python3
"""
Script: 04_parse_ligands.py
Purpose: Parse ligands in cleaned PDBbind v2020.R1 data to generate canonical SMILES and ECFP4 fingerprints.

Usage:
    cd binding_affinity/
    python scripts/04_parse_ligands.py

Input:
    output/pdbbind_v2020_metadata_cleaned.csv

Output:
    output/pdbbind_v2020_metadata_ligand.csv
    output/pdbbind_v2020_metadata_ligand_stats.json
    output/ligand_parse_failures.csv
"""

import argparse
import csv
import json
from pathlib import Path
from collections import Counter

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

# Suppress verbose RDKit warnings.
RDLogger.DisableLog("rdApp.*")


def parse_ligand_from_sdf(path: Path):
    """Parse a ligand from SDF and return (mol, error_msg)."""
    if not path.exists():
        return None, "SDF file not found"
    try:
        suppl = Chem.SDMolSupplier(str(path), sanitize=True, removeHs=False)
        mols = [m for m in suppl if m is not None]
        if not mols:
            return None, "No valid molecule in SDF"
        mol = mols[0]
        if mol.GetNumHeavyAtoms() == 0:
            return None, "Molecule has no heavy atoms"
        return mol, None
    except Exception as e:
        return None, f"SDF parse error: {e}"


def parse_ligand_from_mol2(path: Path):
    """Parse a ligand from MOL2 and return (mol, error_msg)."""
    if not path.exists():
        return None, "MOL2 file not found"
    try:
        mol = Chem.MolFromMol2File(str(path), removeHs=False)
        if mol is None:
            return None, "MOL2 parse failed"
        if mol.GetNumHeavyAtoms() == 0:
            return None, "Molecule has no heavy atoms"
        return mol, None
    except Exception as e:
        return None, f"MOL2 parse error: {e}"


def get_canonical_smiles(mol):
    """Get canonical SMILES."""
    try:
        return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    except Exception:
        return None


def get_ecfp4(mol, radius=2, n_bits=2048):
    """Compute an ECFP4 (Morgan) fingerprint and return indices of set bits."""
    try:
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=n_bits)
        on_bits = sorted([int(i) for i in fp.GetOnBits()])
        return on_bits
    except Exception:
        return None


def parse_ligand(row, data_dir: Path):
    """Try to parse the ligand in one record."""
    pdb_id = row["pdb_id"]
    sdf_rel = row.get("ligand_sdf_file", "")
    mol2_rel = row.get("ligand_mol2_file", "")

    sdf_path = data_dir / sdf_rel if pd.notna(sdf_rel) and sdf_rel else None
    mol2_path = data_dir / mol2_rel if pd.notna(mol2_rel) and mol2_rel else None

    # Try SDF first.
    mol, error = None, None
    if sdf_path and sdf_path.exists():
        mol, error = parse_ligand_from_sdf(sdf_path)

    # Try MOL2 if SDF fails.
    if mol is None and mol2_path and mol2_path.exists():
        mol, error = parse_ligand_from_mol2(mol2_path)
        if mol is not None:
            error = None

    # Both parsers failed.
    if mol is None:
        return {
            "pdb_id": pdb_id,
            "ligand_parse_success": False,
            "ligand_parse_error": error or "No ligand file available",
            "ligand_smiles": None,
            "ligand_ecfp4": None,
            "ligand_num_heavy_atoms": None,
        }

    smiles = get_canonical_smiles(mol)
    ecfp4 = get_ecfp4(mol)
    num_heavy = mol.GetNumHeavyAtoms()

    return {
        "pdb_id": pdb_id,
        "ligand_parse_success": True,
        "ligand_parse_error": None,
        "ligand_smiles": smiles,
        "ligand_ecfp4": json.dumps(ecfp4) if ecfp4 is not None else None,
        "ligand_num_heavy_atoms": num_heavy,
    }


def main():
    parser = argparse.ArgumentParser(description="Parse ligands for PDBbind v2020.R1")
    parser.add_argument(
        "--input",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "output" / "pdbbind_v2020_metadata_cleaned.csv",
        help="Input cleaned metadata CSV",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "output" / "pdbbind_v2020_metadata_ligand.csv",
        help="Output CSV with ligand features",
    )
    parser.add_argument(
        "--stats",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "output" / "pdbbind_v2020_metadata_ligand_stats.json",
        help="Output stats JSON",
    )
    parser.add_argument(
        "--failures",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "output" / "ligand_parse_failures.csv",
        help="Output failures CSV",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "data",
        help="Root data directory",
    )
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    total = len(df)
    print(f"[Input] Samples after cleaning: {total}")

    results = []
    for _, row in df.iterrows():
        res = parse_ligand(row, args.data_dir)
        results.append(res)

    ligand_df = pd.DataFrame(results)

    # Merge into the original table.
    merged = df.merge(
        ligand_df,
        on="pdb_id",
        how="left",
    )

    # Filter out samples that failed parsing.
    success_df = merged[merged["ligand_parse_success"] == True].copy()
    failure_df = merged[merged["ligand_parse_success"] != True].copy()

    # Output
    args.output.parent.mkdir(parents=True, exist_ok=True)
    success_df.to_csv(args.output, index=False)
    failure_df.to_csv(args.failures, index=False)

    stats = {
        "total_input": total,
        "ligand_parse_success": int(len(success_df)),
        "ligand_parse_failure": int(len(failure_df)),
        "retention_rate": round(len(success_df) / total, 4) if total > 0 else 0.0,
        "failure_reasons": dict(Counter(failure_df["ligand_parse_error"].astype(str))),
        "heavy_atoms_stats": success_df["ligand_num_heavy_atoms"].describe().to_dict() if len(success_df) > 0 else {},
    }

    with open(args.stats, "w", encoding="utf-8") as f:
        json.dump(stats, f, indent=2)

    print(f"[Output] Successfully parsed ligands: {len(success_df)}")
    print(f"[Output] Failed ligand parses: {len(failure_df)}")
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
