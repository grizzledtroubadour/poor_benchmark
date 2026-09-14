#!/usr/bin/env python3
"""
Download monomer structures for all representative pairs and filter the test set
so that it only contains samples whose both chains have apo or predicted structures.

Structure selection priority per chain:
  1. apo monomer
  2. predicted monomer (AlphaFold)
  3. holo monomer (fallback)

Outputs:
  - data/pdb_structures_all/           # downloaded PDB files
  - data/structure_download_files.txt  # list of files to download
  - output/structure_selection.csv     # representative -> selected files mapping
  - output/ppi_test_structure_filtered.csv  # diagnostic: structure-filtered test set
  - output/structure_download_summary.json

NOTE: the final splits are written by 09_rebuild_with_xtal_negatives.py (which
consumes structure_selection.csv and applies the same test-set structure filter
inline). The filtered-test output here is kept as a diagnostic artifact only.
"""

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

PDB_URL_PREFIX = "https://storage.googleapis.com/pinder/2024-02/pdbs"
PDB_DIR = "data/pdb_structures_all"
INDEX_PATH = "data/index.parquet"
CORRECTED_PATH = "output/ppi_pairs_corrected_seq95.csv"
TEST_PATH = str(ROOT / "splits" / "ppi_test.csv")
OUTPUT_TEST_PATH = str(ROOT / "output" / "ppi_test_structure_filtered.csv")
FILES_LIST_PATH = "data/structure_download_files.txt"
SELECTION_PATH = "output/structure_selection.csv"
SUMMARY_PATH = "output/structure_download_summary.json"


def select_structures(corrected_path, index_path):
    corrected = pd.read_csv(corrected_path)
    idx_cols = [
        "id", "pdb_id", "uniprot_R", "uniprot_L",
        "apo_R", "apo_L", "predicted_R", "predicted_L", "holo_R", "holo_L",
        "apo_R_pdb", "apo_L_pdb", "predicted_R_pdb", "predicted_L_pdb",
        "holo_R_pdb", "holo_L_pdb",
    ]
    idx = pd.read_parquet(index_path, columns=idx_cols)
    idx = idx[idx["id"].isin(set(corrected["representative_id"]))].copy()

    df = corrected[[
        "representative_id", "pair_key", "label_str", "split",
        "protein_A_id", "protein_B_id",
    ]].merge(idx, left_on="representative_id", right_on="id", how="left")

    selections = []
    all_files = set()
    protein_has_alt = {}

    for _, row in df.iterrows():
        rep_id = row["representative_id"]
        files = []
        sources = []

        # R chain / protein_A
        if row["apo_R"] and row["apo_R_pdb"]:
            files.append(row["apo_R_pdb"])
            sources.append("apo_R")
        elif row["predicted_R"] and row["predicted_R_pdb"]:
            files.append(row["predicted_R_pdb"])
            sources.append("predicted_R")
        elif row["holo_R"] and row["holo_R_pdb"]:
            files.append(row["holo_R_pdb"])
            sources.append("holo_R")

        # L chain / protein_B
        if row["apo_L"] and row["apo_L_pdb"]:
            files.append(row["apo_L_pdb"])
            sources.append("apo_L")
        elif row["predicted_L"] and row["predicted_L_pdb"]:
            files.append(row["predicted_L_pdb"])
            sources.append("predicted_L")
        elif row["holo_L"] and row["holo_L_pdb"]:
            files.append(row["holo_L_pdb"])
            sources.append("holo_L")

        r_ok = any(s in sources for s in ["apo_R", "predicted_R"])
        l_ok = any(s in sources for s in ["apo_L", "predicted_L"])
        both_alt = r_ok and l_ok

        selections.append({
            "representative_id": rep_id,
            "pair_key": row["pair_key"],
            "label_str": row["label_str"],
            "split": row["split"],
            "protein_A_id": row["protein_A_id"],
            "protein_B_id": row["protein_B_id"],
            "selected_files": ";".join(files),
            "selected_sources": ";".join(sources),
            "R_has_apo_or_predicted": r_ok,
            "L_has_apo_or_predicted": l_ok,
            "both_chains_apo_or_predicted": both_alt,
        })

        all_files.update(files)
        protein_has_alt[row["protein_A_id"]] = protein_has_alt.get(row["protein_A_id"], False) or r_ok
        protein_has_alt[row["protein_B_id"]] = protein_has_alt.get(row["protein_B_id"], False) or l_ok

    sel_df = pd.DataFrame(selections)
    return sel_df, sorted(all_files), protein_has_alt


def download_pdbs(files, pdb_dir, max_parallel=16):
    os.makedirs(pdb_dir, exist_ok=True)
    to_download = [f for f in files if not os.path.exists(os.path.join(pdb_dir, f))]
    print(f"Already downloaded: {len(files) - len(to_download)}")
    print(f"To download: {len(to_download)}")

    if not to_download:
        return

    lines = []
    for fn in to_download:
        url = f"{PDB_URL_PREFIX}/{fn}"
        lines.append(url)
        lines.append(f"  out={fn}")
        lines.append(f"  dir={pdb_dir}")

    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as tf:
        tf.write("\n".join(lines) + "\n")
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


def filter_test_set(test_path, selection_df, protein_has_alt):
    test = pd.read_csv(test_path, low_memory=False)
    print(f"Original test set: {len(test)}")

    # Positive samples: keep only if both chains have apo/predicted
    # (current splits schema: unique_id holds the representative_id for
    # bio/xtal rows; random negatives carry neg_* ids and are handled below)
    pos = test[test["label"] == 1].copy()
    pos_merged = pos.merge(
        selection_df[["representative_id", "both_chains_apo_or_predicted"]],
        left_on="unique_id",
        right_on="representative_id",
        how="left",
    )
    pos_keep = pos_merged[pos_merged["both_chains_apo_or_predicted"] == True].copy()
    pos_keep = pos_keep.drop(columns=["representative_id", "both_chains_apo_or_predicted"])
    print(f"Test positives kept: {len(pos_keep)} / {len(pos)}")

    # Negative samples: keep only if both proteins have at least one apo/predicted structure
    neg = test[test["label"] == 0].copy()
    neg["both_alt"] = neg.apply(
        lambda r: protein_has_alt.get(r["protein_A_id"], False)
        and protein_has_alt.get(r["protein_B_id"], False),
        axis=1,
    )
    neg_keep = neg[neg["both_alt"] == True].copy()
    neg_keep = neg_keep.drop(columns=["both_alt"])
    print(f"Test negatives kept: {len(neg_keep)} / {len(neg)}")

    filtered = pd.concat([pos_keep, neg_keep], ignore_index=True)
    return filtered


def main():
    print("=== Selecting structures ===")
    sel_df, files, protein_has_alt = select_structures(CORRECTED_PATH, INDEX_PATH)
    sel_df.to_csv(SELECTION_PATH, index=False)
    print(f"Saved {SELECTION_PATH}")

    with open(FILES_LIST_PATH, "w") as f:
        for fn in files:
            f.write(fn + "\n")
    print(f"Unique structure files: {len(files)}")

    print("\n=== Downloading structures ===")
    download_pdbs(files, PDB_DIR)

    downloaded = sum(1 for f in files if os.path.exists(os.path.join(PDB_DIR, f)))
    print(f"Downloaded: {downloaded} / {len(files)}")

    print("\n=== Filtering test set ===")
    filtered_test = filter_test_set(TEST_PATH, sel_df, protein_has_alt)
    os.makedirs(os.path.dirname(OUTPUT_TEST_PATH) or ".", exist_ok=True)
    filtered_test.to_csv(OUTPUT_TEST_PATH, index=False)
    print(f"Saved {OUTPUT_TEST_PATH}: {len(filtered_test)} rows")

    # Summary
    summary = {
        "unique_structure_files": len(files),
        "downloaded_files": downloaded,
        "pairs_with_both_chains_apo_or_predicted": int(sel_df["both_chains_apo_or_predicted"].sum()),
        "proteins_with_apo_or_predicted": int(sum(protein_has_alt.values())),
        "original_test_size": int(len(pd.read_csv(TEST_PATH, low_memory=False))),
        "filtered_test_size": int(len(filtered_test)),
        "filtered_test_positives": int((filtered_test["label"] == 1).sum()),
        "filtered_test_negatives": int((filtered_test["label"] == 0).sum()),
    }
    with open(SUMMARY_PATH, "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
