#!/usr/bin/env python3
"""
Fetch sequences for proteins missing from UniProt by downloading their
PINDER holo monomer PDB structures and extracting the amino-acid sequence.

Steps:
  1. Identify proteins in corrected_seq95 without UniProt sequence
  2. Map them to holo_R_pdb / holo_L_pdb filenames via index.parquet
  3. Download missing PDB monomers from gs://pinder/2024-02/pdbs/
  4. Extract sequence from each PDB (ATOM records, canonical residue order)
  5. Append recovered sequences to data/uniprot_sequences_all.fasta
  6. Re-run representative selection and the XTAL-negative rebuild

Outputs:
  - data/pdb_structures_missing/        # downloaded PDB monomers
  - data/missing_seq_pdb_files.txt      # list of PDB files to download
  - data/missing_seq_download_log.json  # success/failure log
  - data/uniprot_sequences_all.fasta    # appended with recovered sequences
  - output/missing_seq_recovery_summary.json
"""

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.PDB import PDBParser, PPBuilder
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import warnings

warnings.filterwarnings("ignore", category=PDBConstructionWarning)

PDB_URL_PREFIX = "https://storage.googleapis.com/pinder/2024-02/pdbs"
PDB_DIR = "data/pdb_structures_missing"
FASTA_PATH = "data/uniprot_sequences_all.fasta"
INDEX_PATH = "data/index.parquet"
CORRECTED_PATH = "output/ppi_pairs_corrected_seq95.csv"
FILES_LIST_PATH = "data/missing_seq_pdb_files.txt"
LOG_PATH = "data/missing_seq_download_log.json"
SUMMARY_PATH = "output/missing_seq_recovery_summary.json"
# Optional local Swiss-Prot preload (project-root shared data dir; skipped
# automatically when absent)
SPROT_FASTA = str(
    Path(__file__).resolve().parents[3] / "data" / "uniprot" / "uniprot_sprot.fasta"
)


def load_existing_sequences(fasta_paths):
    seqs = {}
    for fasta_path in fasta_paths:
        if not os.path.exists(fasta_path):
            continue
        for record in SeqIO.parse(fasta_path, "fasta"):
            parts = record.id.split("|")
            uid = parts[1] if len(parts) >= 2 else record.id
            seqs[uid] = str(record.seq)
    return seqs


def identify_missing_proteins(corrected_path, seqs):
    df = pd.read_csv(corrected_path)
    proteins = set(df["protein_A_id"]) | set(df["protein_B_id"])
    # treat empty-sequence records (key present, value "") as missing too
    missing = sorted({p for p in proteins if p not in seqs or not seqs[p]})
    print(f"Proteins in corrected_seq95: {len(proteins)}")
    print(f"Missing UniProt sequences: {len(missing)}")
    return missing


def map_proteins_to_pdb_files(index_path, representative_ids, missing_proteins):
    idx = pd.read_parquet(
        index_path,
        columns=["id", "uniprot_R", "uniprot_L", "holo_R_pdb", "holo_L_pdb"],
    )
    idx = idx[idx["id"].isin(set(representative_ids))]

    protein_to_files = {p: [] for p in missing_proteins}
    files_set = set()

    for _, row in idx.iterrows():
        r_prot = row["uniprot_R"]
        l_prot = row["uniprot_L"]
        if r_prot in protein_to_files and row["holo_R_pdb"]:
            protein_to_files[r_prot].append(row["holo_R_pdb"])
            files_set.add(row["holo_R_pdb"])
        if l_prot in protein_to_files and row["holo_L_pdb"]:
            protein_to_files[l_prot].append(row["holo_L_pdb"])
            files_set.add(row["holo_L_pdb"])

    return protein_to_files, sorted(files_set)


def build_aria2_input(files, out_dir):
    lines = []
    for fn in files:
        url = f"{PDB_URL_PREFIX}/{fn}"
        lines.append(url)
        lines.append(f"  out={fn}")
        lines.append(f"  dir={out_dir}")
    return "\n".join(lines) + "\n"


def download_pdbs(files, pdb_dir, max_parallel=16):
    os.makedirs(pdb_dir, exist_ok=True)

    # Skip already downloaded
    to_download = [f for f in files if not os.path.exists(os.path.join(pdb_dir, f))]
    print(f"Already downloaded: {len(files) - len(to_download)}")
    print(f"To download: {len(to_download)}")

    if not to_download:
        return

    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as tf:
        tf.write(build_aria2_input(to_download, pdb_dir))
        input_file = tf.name

    try:
        subprocess.run(
            [
                "aria2c",
                "--input-file", input_file,
                "--max-concurrent-downloads", str(max_parallel),
                "--split", "4",
                "--max-connection-per-server", "4",
                "--retry-wait", "2",
                "--max-tries", "5",
                "--log-level", "warn",
                "--summary-interval", "0",
                "--console-log-level", "warn",
            ],
            check=True,
        )
    finally:
        os.unlink(input_file)


def extract_sequence_from_pdb(pdb_path):
    """Extract sequence from ATOM records. Return empty string if fails."""
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("X", pdb_path)
        ppb = PPBuilder()
        sequences = []
        for pp in ppb.build_peptides(structure, aa_only=True):
            sequences.append(str(pp.get_sequence()))
        return "".join(sequences)
    except Exception as e:
        return ""


def recover_sequences(protein_to_files, pdb_dir):
    recovered = {}
    failed_files = []
    multiple_seq_files = []

    for protein, files in protein_to_files.items():
        best_seq = ""
        best_file = ""
        for fn in files:
            pdb_path = os.path.join(pdb_dir, fn)
            if not os.path.exists(pdb_path):
                failed_files.append(fn)
                continue
            seq = extract_sequence_from_pdb(pdb_path)
            if not seq:
                failed_files.append(fn)
                continue
            if len(seq) > len(best_seq):
                best_seq = seq
                best_file = fn
        if best_seq:
            recovered[protein] = {"sequence": best_seq, "source_pdb": best_file}

    return recovered, failed_files


def append_sequences_to_fasta(recovered, fasta_path):
    os.makedirs(os.path.dirname(fasta_path) or ".", exist_ok=True)
    with open(fasta_path, "a") as f:
        for protein, info in recovered.items():
            f.write(f">{protein}\n{info['sequence']}\n")
    print(f"Appended {len(recovered)} sequences to {fasta_path}")


def re_run_pipeline():
    """Re-run representative selection and the splits rebuild with the
    updated sequence file.

    NOTE: the former 07_build_negative_pairs.py / 08_time_split_dataset.py
    steps were superseded by 09_rebuild_with_xtal_negatives.py (formerly
    23_rebuild_with_xtal_negatives.py) and archived to output/scripts_archive/.
    """
    scripts = [
        "scripts/04_build_pairs_with_corrected_ranking.py",
        "scripts/09_rebuild_with_xtal_negatives.py",
    ]
    for script in scripts:
        print(f"\n=== Running {script} ===")
        subprocess.run([sys.executable, script], check=True)


def main():
    os.makedirs("output", exist_ok=True)

    # 1. Load existing sequences
    seqs = load_existing_sequences([
        SPROT_FASTA,
        FASTA_PATH,
    ])
    print(f"Existing sequences: {len(seqs)}")

    # 2. Identify missing proteins
    missing = identify_missing_proteins(CORRECTED_PATH, seqs)

    # 3. Map to PDB files
    corrected = pd.read_csv(CORRECTED_PATH)
    representative_ids = corrected["representative_id"].tolist()
    protein_to_files, files = map_proteins_to_pdb_files(
        INDEX_PATH, representative_ids, missing
    )
    print(f"Unique PDB files needed: {len(files)}")

    with open(FILES_LIST_PATH, "w") as f:
        for fn in files:
            f.write(fn + "\n")

    # 4. Download PDBs
    download_pdbs(files, PDB_DIR)

    # 5. Extract sequences
    recovered, failed_files = recover_sequences(protein_to_files, PDB_DIR)
    print(f"Recovered sequences: {len(recovered)} / {len(missing)}")

    # 6. Append to FASTA
    append_sequences_to_fasta(recovered, FASTA_PATH)

    # 7. Re-run pipeline
    re_run_pipeline()

    # 8. Save summary
    summary = {
        "missing_proteins": len(missing),
        "pdb_files_needed": len(files),
        "pdb_files_downloaded": len([f for f in files if os.path.exists(os.path.join(PDB_DIR, f))]),
        "sequences_recovered": len(recovered),
        "recovery_rate": len(recovered) / len(missing) if missing else 0.0,
        "failed_files_count": len(set(failed_files)),
        "sample_recovered": list(recovered.keys())[:10],
    }
    with open(SUMMARY_PATH, "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
