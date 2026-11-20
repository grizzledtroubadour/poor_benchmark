#!/usr/bin/env python3
"""
Script: 02_parse_index.py
Purpose: Parse the PDBbind v2020.R1 index (free PDBbind+ distribution) to extract protein-ligand affinity labels and metadata.

Usage:
    cd binding_affinity/
    python scripts/02_parse_index.py

Input:
    data/index/INDEX_general_PL.2020R1.lst
    data/P-L/{year_range}/{pdb_id}/{pdb_id}_protein.pdb
    data/P-L/{year_range}/{pdb_id}/{pdb_id}_ligand.sdf
    data/P-L/{year_range}/{pdb_id}/{pdb_id}_pocket.pdb

Output:
    output/pdbbind_v2020_metadata.csv
    output/pdbbind_v2020_metadata_stats.json

Example index format (no precomputed pK column):
    2tpi  2.10  1982  Kd=49uM       // 2tpi.pdf (2-mer)
    6rsa   NMR  1986  Ki=40uM       // 6rsa.pdf (UVC)
    4cts  2.90  1984  Kd<10uM       // 4cts.pdf (OAA)
"""

import os
import re
import csv
import math
import json
import argparse
from pathlib import Path
from collections import defaultdict, Counter

# Convert units to M (molar).
UNIT_FACTOR = {
    "M": 1.0,
    "mM": 1e-3,
    "uM": 1e-6,
    "um": 1e-6,
    "μM": 1e-6,
    "nM": 1e-9,
    "nm": 1e-9,
    "pM": 1e-12,
    "pm": 1e-12,
    "fM": 1e-15,
    "fm": 1e-15,
}

# Affinity string parsing: supports Kd/Ki/IC50, operators, numeric values, and units.
# Examples: Kd=49uM, Ki=0.43uM, IC50=5nM, Kd<10uM, Kd=741.3+/-53.1uM
AFFINITY_RE = re.compile(
    r"^(Kd|Ki|IC50)\s*(<=|>=|[=<>~≤≥])\s*([0-9]+\.?[0-9]*)\s*([a-zA-Zμ]+)"
)

NOTE_PATTERNS = {
    "incomplete_ligand": re.compile(r"\[Incomplete ligand\]", re.IGNORECASE),
    "covalent_complex": re.compile(r"\[Covalent complex\]", re.IGNORECASE),
    "different_protein": re.compile(r"\[Different protein in assay\]", re.IGNORECASE),
    "different_ligand": re.compile(r"\[Different ligand in assay\]", re.IGNORECASE),
    "uncommon_element": re.compile(r"\[Uncommon element:\s*([^\]]+)\]", re.IGNORECASE),
    "disulfide_bond": re.compile(r"\[Disulfide bond\]", re.IGNORECASE),
}


def parse_affinity(token: str):
    """
    Parse an affinity string.
    Return (aff_type, operator, value_M, pK), or None if parsing fails.
    """
    token = token.strip()
    # Remove the +/- uncertainty, e.g. 741.3+/-53.1uM -> 741.3uM.
    token = re.sub(r"\+/-\s*[0-9]+\.?[0-9]*", "", token)
    token = re.sub(r"\s+", "", token)

    m = AFFINITY_RE.match(token)
    if not m:
        return None

    aff_type = m.group(1)
    operator = m.group(2)
    value = float(m.group(3))
    unit = m.group(4)

    factor = UNIT_FACTOR.get(unit)
    if factor is None:
        factor = UNIT_FACTOR.get(unit.lower())
    if factor is None:
        raise ValueError(f"Unknown unit: {unit} in token '{token}'")

    value_M = value * factor
    if value_M <= 0:
        raise ValueError(f"Non-positive affinity value: {token}")

    pK = -math.log10(value_M)
    return aff_type, operator, value_M, pK


def parse_resolution(token: str):
    """Parse resolution; return None for NMR."""
    token = token.strip()
    if token.upper() == "NMR" or token == "-":
        return None
    try:
        return float(token)
    except ValueError:
        return None


def parse_ligand_name(comment: str):
    """Extract the parenthesized ligand name from comments."""
    m = re.search(r"\(([^()]+)\)", comment)
    return m.group(1).strip() if m else ""


def parse_notes(comment: str):
    """Extract flags from comments."""
    notes = {}
    for key, pattern in NOTE_PATTERNS.items():
        m = pattern.search(comment)
        if m:
            notes[key] = m.group(1).strip() if m.lastindex else True
        else:
            notes[key] = ""
    return notes


def parse_index_line(line: str):
    """
    Parse one line of the PDBbind v2020.R1 index.
    Return a dict, or None for parsing failures/comment lines.
    """
    line = line.strip()
    if not line or line.startswith("#"):
        return None

    # Separate data from comments.
    if "//" in line:
        data_part, comment = line.split("//", 1)
    else:
        data_part, comment = line, ""

    data_tokens = data_part.split()
    if len(data_tokens) < 4:
        return None

    try:
        pdb_id = data_tokens[0].lower()
        resolution = parse_resolution(data_tokens[1])
        release_year = int(data_tokens[2])
        parsed_aff = parse_affinity(data_tokens[3])
        if parsed_aff is None:
            return {"parse_error": f"Cannot parse affinity: {data_tokens[3]}", "raw": line}
        aff_type, operator, aff_value, pK = parsed_aff
    except Exception as e:
        return {"parse_error": str(e), "raw": line}

    ligand_name = parse_ligand_name(comment)
    notes = parse_notes(comment)

    rec = {
        "pdb_id": pdb_id,
        "resolution": resolution if resolution is not None else "",
        "release_year": release_year,
        "affinity_type": aff_type,
        "affinity_operator": operator,
        "affinity_value": aff_value,
        "pK": pK,
        "ligand_name": ligand_name,
        "reference": comment.strip().split()[0] if comment.strip() else "",
        "comment": comment.strip(),
        "raw": line,
    }
    rec.update(notes)
    return rec


def load_index_file(path: Path):
    """Load the index file and return parsed records."""
    records = []
    errors = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            rec = parse_index_line(line)
            if rec is None:
                continue
            if "parse_error" in rec:
                errors.append(rec)
                continue
            records.append(rec)
    return records, errors


def locate_structure_files(data_dir: Path, pdb_id: str):
    """
    Locate structure files by pdb_id under data/P-L/.
    Return relative paths as a dict, or an empty dict if no files are found.
    """
    pl_dir = data_dir / "P-L"
    if not pl_dir.exists():
        return {}

    # Search by year range.
    for year_range in pl_dir.iterdir():
        if not year_range.is_dir():
            continue
        complex_dir = year_range / pdb_id
        if complex_dir.is_dir():
            protein = complex_dir / f"{pdb_id}_protein.pdb"
            ligand_mol2 = complex_dir / f"{pdb_id}_ligand.mol2"
            ligand_sdf = complex_dir / f"{pdb_id}_ligand.sdf"
            pocket = complex_dir / f"{pdb_id}_pocket.pdb"
            return {
                "protein_file": str(protein.relative_to(data_dir)) if protein.exists() else "",
                "ligand_mol2_file": str(ligand_mol2.relative_to(data_dir)) if ligand_mol2.exists() else "",
                "ligand_sdf_file": str(ligand_sdf.relative_to(data_dir)) if ligand_sdf.exists() else "",
                "pocket_file": str(pocket.relative_to(data_dir)) if pocket.exists() else "",
                "year_range": year_range.name,
            }
    return {}


def build_metadata(data_dir: Path, output_dir: Path):
    """Parse the v2020.R1 index and generate the metadata CSV."""
    index_path = data_dir / "index" / "INDEX_general_PL.2020R1.lst"

    if not index_path.exists():
        raise FileNotFoundError(f"Index file does not exist: {index_path}")

    print(f"[Loading] {index_path}")
    records, errors = load_index_file(index_path)
    print(f"[Loading] Parsed successfully: {len(records)} records; parsing failures: {len(errors)} records")

    if errors:
        print(f"[Warning] First five records that failed parsing:")
        for e in errors[:5]:
            print(f"  - {e['raw'].strip()[:100]} | {e['parse_error']}")

    # Deduplicate by pdb_id, keeping the first occurrence.
    seen = {}
    deduped = []
    for rec in records:
        pid = rec["pdb_id"]
        if pid not in seen:
            seen[pid] = rec
            deduped.append(rec)

    print(f"[Processing] After deduplication: {len(deduped)} records")

    # Locate structure files
    print("[Processing] Locate structure files...")
    missing_structures = 0
    for rec in deduped:
        struct = locate_structure_files(data_dir, rec["pdb_id"])
        if not struct:
            missing_structures += 1
        rec.update(struct)

    # Output CSV
    output_dir.mkdir(parents=True, exist_ok=True)
    out_csv = output_dir / "pdbbind_v2020_metadata.csv"
    fieldnames = [
        "pdb_id",
        "resolution",
        "release_year",
        "affinity_type",
        "affinity_operator",
        "affinity_value",
        "pK",
        "ligand_name",
        "reference",
        "year_range",
        "protein_file",
        "ligand_mol2_file",
        "ligand_sdf_file",
        "pocket_file",
        "incomplete_ligand",
        "covalent_complex",
        "different_protein",
        "different_ligand",
        "uncommon_element",
        "disulfide_bond",
        "comment",
        "raw",
    ]
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for rec in deduped:
            writer.writerow({k: rec.get(k, "") for k in fieldnames})

    # Statistics
    stats = {
        "total_records": len(records),
        "parse_errors": len(errors),
        "unique_records": len(deduped),
        "missing_structures": missing_structures,
        "affinity_type_counts": dict(Counter(r["affinity_type"] for r in deduped)),
        "affinity_operator_counts": dict(Counter(r["affinity_operator"] for r in deduped)),
        "release_year_range": {
            "min": min(r["release_year"] for r in deduped),
            "max": max(r["release_year"] for r in deduped),
        },
        "resolution_available": sum(1 for r in deduped if r["resolution"] != ""),
        "note_counts": {
            key: sum(1 for r in deduped if r.get(key))
            for key in NOTE_PATTERNS.keys()
        },
    }

    stats_path = output_dir / "pdbbind_v2020_metadata_stats.json"
    with open(stats_path, "w", encoding="utf-8") as f:
        json.dump(stats, f, indent=2)

    print(f"\n[Output] Metadata CSV: {out_csv}")
    print(f"[Output] Statistics JSON: {stats_path}")
    print(f"[Statistics] Unique complexes: {stats['unique_records']}")
    print(f"[Statistics] Affinity type distribution: {stats['affinity_type_counts']}")
    print(f"[Statistics] Operator distribution: {stats['affinity_operator_counts']}")
    print(f"[Statistics] Year range: {stats['release_year_range']['min']} - {stats['release_year_range']['max']}")
    print(f"[Statistics] Missing structure files: {stats['missing_structures']}")


def main():
    parser = argparse.ArgumentParser(description="Parse PDBbind v2020.R1 index files")
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "data",
        help="Path to binding_affinity/data directory",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "output",
        help="Path to output directory",
    )
    args = parser.parse_args()

    build_metadata(args.data_dir, args.output_dir)


if __name__ == "__main__":
    main()
