#!/usr/bin/env python3
"""
Stage 2: Extract amino acid sequences from filtered PDB entries.
Only sequences are saved (no coordinate files) to save disk space.

长度下限 MIN_RESIDUES=31（v2 口径，见该常量处注释）：31–59aa 短链一并入池，
≤30aa 在源头排除。

Inputs:
    ss_prediction/output/stage1_filtered_entries.csv
    data/pdb/*.pdb.gz or *.cif.gz

Outputs:
    ss_prediction/output/stage2_chain_sequences.fa
    ss_prediction/output/stage2_chain_metadata.csv
    ss_prediction/output/extract_sequences.log

Usage:
    python ss_prediction/scripts/03_extract_chain_sequences.py
"""

import gzip
import csv
import logging
from pathlib import Path
from collections import OrderedDict
from concurrent.futures import ProcessPoolExecutor, as_completed

TASK_DIR = Path(__file__).resolve().parent.parent  # datasets/ss_prediction
ROOT = TASK_DIR.parent.parent                      # 仓库根（POOR）
OUTPUT_DIR = TASK_DIR / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(OUTPUT_DIR / "extract_sequences.log", mode="w"),
        logging.StreamHandler(),
    ],
)
logger = logging.getLogger(__name__)

PDB_DIR = ROOT / "data" / "pdb"
# v2 口径（2026-09-02 线性化）：MIN_RESIDUES = 31。≤30aa 短链在源头排除，
# 取代历史路线（原 MIN_RESIDUES=60 过滤 <60 短链，事后经 12–19 恢复
# 31–59 段并删除 ≤30 段——已归档 output/scripts_archive/）。31–59aa 短链
# 共 6,453 条随本步一并入池，stage2 产出由 249,018 变为 255,471 条。
MIN_RESIDUES = 31
N_WORKERS = 32

AA_MAP = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
    "MSE": "M",  # Selenomethionine
    "SEC": "C",  # Selenocysteine
    "PYL": "K",  # Pyrrolysine
}


def find_pdb_file(pdb_id: str):
    """Find local PDB file for a given PDB ID. Prefer .pdb.gz over .cif.gz."""
    pdb_path = PDB_DIR / f"{pdb_id}.pdb.gz"
    if pdb_path.exists():
        return pdb_path, "pdb"
    cif_path = PDB_DIR / f"{pdb_id}.cif.gz"
    if cif_path.exists():
        return cif_path, "cif"
    return None, None


def extract_from_pdb_gz(filepath: Path, pdb_id: str):
    """Extract all protein chain sequences from a gzipped PDB file."""
    results = []
    try:
        with gzip.open(filepath, "rt", encoding="utf-8", errors="replace") as fh:
            # chain_id -> OrderedDict of (res_seq, icode) -> res_name
            chains = {}
            current_model = 0

            for line in fh:
                if not line:
                    continue
                record = line[:6].strip()

                # NMR: only keep model 1
                if record == "MODEL":
                    current_model += 1
                    if current_model > 1:
                        break

                if record not in ("ATOM", "HETATM"):
                    continue

                # Standard PDB format
                res_name = line[17:20].strip().upper()
                if res_name not in AA_MAP:
                    continue

                chain_id = line[21:22].strip()
                if not chain_id:
                    chain_id = "A"

                res_seq = line[22:26].strip()
                icode = line[26:27].strip()
                atom_name = line[12:16].strip()

                # Only CA atoms (or N/C for robustness) to avoid double-counting
                if atom_name != "CA":
                    continue

                if chain_id not in chains:
                    chains[chain_id] = OrderedDict()

                key = (res_seq, icode)
                if key not in chains[chain_id]:
                    chains[chain_id][key] = {"res_name": res_name, "has_mse": res_name == "MSE"}

            for chain_id, residues in chains.items():
                seq = "".join(AA_MAP[r["res_name"]] for r in residues.values())
                if len(seq) >= MIN_RESIDUES:
                    results.append({
                        "pdb_id": pdb_id,
                        "chain_id": chain_id,
                        "seq": seq,
                        "num_residues": len(seq),
                        "has_mse": any(r["has_mse"] for r in residues.values()),
                        "file_type": "pdb",
                    })
    except Exception as e:
        logger.warning(f"Error parsing {filepath}: {e}")

    return results


def extract_from_cif_gz(filepath: Path, pdb_id: str):
    """Extract all protein chain sequences from a gzipped mmCIF file."""
    results = []
    try:
        with gzip.open(filepath, "rt", encoding="utf-8", errors="replace") as fh:
            in_atom_site = False
            got_loop = False
            col_idx = {}
            chains = {}
            current_model = None

            for line in fh:
                line_stripped = line.strip()
                if not line_stripped:
                    continue

                if line_stripped == "loop_":
                    got_loop = True
                    continue

                if got_loop and line_stripped.startswith("_atom_site."):
                    in_atom_site = True
                    field = line_stripped.split()[0]
                    col_idx[field] = len(col_idx)
                    continue

                if got_loop and in_atom_site and (line_stripped.startswith("#") or line_stripped.startswith("_")):
                    in_atom_site = False
                    got_loop = False
                    continue

                if not in_atom_site:
                    continue

                parts = line.split()
                if len(parts) < len(col_idx):
                    continue

                # Extract fields
                group = parts[col_idx.get("_atom_site.group_PDB", -1)] if "_atom_site.group_PDB" in col_idx else ""
                if group not in ("ATOM", "HETATM"):
                    continue

                comp_id = parts[col_idx.get("_atom_site.label_comp_id", -1)].upper() if "_atom_site.label_comp_id" in col_idx else ""
                if comp_id not in AA_MAP:
                    continue

                # Check atom name (prefer CA)
                atom_id = parts[col_idx.get("_atom_site.label_atom_id", -1)] if "_atom_site.label_atom_id" in col_idx else ""
                if atom_id != "CA":
                    continue

                chain_id = parts[col_idx.get("_atom_site.auth_asym_id", -1)] if "_atom_site.auth_asym_id" in col_idx else ""
                if not chain_id or chain_id == "?":
                    chain_id = parts[col_idx.get("_atom_site.label_asym_id", -1)] if "_atom_site.label_asym_id" in col_idx else ""

                seq_id = parts[col_idx.get("_atom_site.auth_seq_id", -1)] if "_atom_site.auth_seq_id" in col_idx else ""
                if seq_id == "?":
                    seq_id = parts[col_idx.get("_atom_site.label_seq_id", -1)] if "_atom_site.label_seq_id" in col_idx else ""

                icode = parts[col_idx.get("_atom_site.pdbx_PDB_ins_code", -1)] if "_atom_site.pdbx_PDB_ins_code" in col_idx else ""
                if icode == "?":
                    icode = ""

                # Model handling
                model_num = parts[col_idx.get("_atom_site.pdbx_PDB_model_num", -1)] if "_atom_site.pdbx_PDB_model_num" in col_idx else "1"
                if current_model is None:
                    current_model = model_num
                elif model_num != current_model:
                    break  # Only model 1

                if chain_id not in chains:
                    chains[chain_id] = OrderedDict()

                key = (seq_id, icode)
                if key not in chains[chain_id]:
                    chains[chain_id][key] = {"res_name": comp_id, "has_mse": comp_id == "MSE"}

            for chain_id, residues in chains.items():
                seq = "".join(AA_MAP[r["res_name"]] for r in residues.values())
                if len(seq) >= MIN_RESIDUES:
                    results.append({
                        "pdb_id": pdb_id,
                        "chain_id": chain_id,
                        "seq": seq,
                        "num_residues": len(seq),
                        "has_mse": any(r["has_mse"] for r in residues.values()),
                        "file_type": "cif",
                    })
    except Exception as e:
        logger.warning(f"Error parsing {filepath}: {e}")

    return results


def process_entry(row: dict):
    pdb_id = row["pdb_id"].upper()
    filepath, file_type = find_pdb_file(pdb_id)
    if filepath is None:
        return {"pdb_id": pdb_id, "error": "file_not_found"}

    if file_type == "pdb":
        return extract_from_pdb_gz(filepath, pdb_id)
    else:
        return extract_from_cif_gz(filepath, pdb_id)


def main():
    entries_path = OUTPUT_DIR / "stage1_filtered_entries.csv"
    if not entries_path.exists():
        logger.error(f"Filtered entries not found: {entries_path}")
        return

    entries = []
    with open(entries_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            entries.append(row)

    logger.info(f"Total entries to process: {len(entries)}")

    all_records = []
    errors = []

    with ProcessPoolExecutor(max_workers=N_WORKERS) as executor:
        future_to_entry = {executor.submit(process_entry, e): e for e in entries}
        for i, future in enumerate(as_completed(future_to_entry)):
            if i % 1000 == 0:
                logger.info(f"  Processed {i}/{len(entries)} entries...")
            try:
                result = future.result()
                if isinstance(result, list):
                    all_records.extend(result)
                elif isinstance(result, dict) and "error" in result:
                    errors.append(result)
            except Exception as e:
                entry = future_to_entry[future]
                errors.append({"pdb_id": entry["pdb_id"], "error": str(e)})

    logger.info(f"Done. Total chains extracted: {len(all_records)}, Errors: {len(errors)}")
    if errors:
        logger.info(f"First 5 errors: {errors[:5]}")

    # Save FASTA
    fasta_path = OUTPUT_DIR / "stage2_chain_sequences.fa"
    with open(fasta_path, "w") as f:
        for rec in all_records:
            header = f">{rec['pdb_id']}_{rec['chain_id']}"
            f.write(f"{header}\n{rec['seq']}\n")
    logger.info(f"Saved FASTA to {fasta_path} ({len(all_records)} chains)")

    # Save metadata CSV
    csv_path = OUTPUT_DIR / "stage2_chain_metadata.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["pdb_id", "chain_id", "num_residues", "has_mse", "file_type", "seq"],
        )
        writer.writeheader()
        writer.writerows(all_records)
    logger.info(f"Saved metadata to {csv_path}")

    # Quick stats
    lengths = [r["num_residues"] for r in all_records]
    logger.info(f"Chain length stats: min={min(lengths)}, max={max(lengths)}, mean={sum(lengths)/len(lengths):.1f}")

    # Per-PDB chain count stats
    from collections import Counter
    chain_counts = Counter(r["pdb_id"] for r in all_records)
    logger.info(f"Entries with chains: {len(chain_counts)}")
    logger.info(f"Max chains per entry: {max(chain_counts.values())}")
    logger.info(f"Single-chain entries: {sum(1 for c in chain_counts.values() if c == 1)}")


if __name__ == "__main__":
    main()
