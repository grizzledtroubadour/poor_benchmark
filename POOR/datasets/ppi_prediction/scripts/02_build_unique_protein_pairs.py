#!/usr/bin/env python3
"""
Build a unique protein-pair dataset from all PINDER splits.

Rules:
  - Combine train/val/test/invalid splits
  - Group by unordered protein pair {A, B}
  - If a pair has both BIO and XTAL labels, keep BIO
  - If a pair has only XTAL, keep XTAL
  - Skip pairs with missing/undefined UniProt IDs

Outputs:
  - output/ppi_unique_pairs.csv
  - output/ppi_unique_pairs_summary.json
"""

import json
import os
from pathlib import Path

import pandas as pd
from Bio import SeqIO

# Optional local Swiss-Prot preload (project-root shared data dir; skipped
# automatically when absent)
SPROT_FASTA = str(
    Path(__file__).resolve().parents[3] / "data" / "uniprot" / "uniprot_sprot.fasta"
)


def load_sequences(fasta_paths):
    seqs = {}
    for fasta_path in fasta_paths:
        if not os.path.exists(fasta_path):
            continue
        for record in SeqIO.parse(fasta_path, "fasta"):
            parts = record.id.split("|")
            uid = parts[1] if len(parts) >= 2 else record.id
            seqs[uid] = str(record.seq)
    return seqs


def build_unique_pairs(require_sequences=False):
    idx = pd.read_parquet("data/index.parquet")
    meta = pd.read_parquet("data/metadata.parquet")

    # Merge metadata
    df = idx.merge(
        meta[
            [
                "id",
                "label",
                "method",
                "resolution",
                "release_date",
                "complex_type",
                "length1",
                "length2",
                "buried_sasa",
                "intermolecular_contacts",
                "planarity",
            ]
        ],
        on="id",
        how="left",
    )

    # Filter to BIO or XTAL only
    df = df[df["label"].isin(["BIO", "XTAL"])].copy()
    print(f"BIO + XTAL rows: {len(df)}")

    # Filter invalid UniProt IDs
    df = df[
        (df["uniprot_R"].notna())
        & (df["uniprot_L"].notna())
        & (df["uniprot_R"] != "UNDEFINED")
        & (df["uniprot_L"] != "UNDEFINED")
    ].copy()
    print(f"After filtering valid UniProt IDs: {len(df)}")

    # Create unordered pair key
    df["pair_key"] = df.apply(
        lambda r: "__".join(sorted([r["uniprot_R"], r["uniprot_L"]])), axis=1
    )

    # Sort by label (BIO first) and resolution (best first)
    df["label_priority"] = df["label"].map({"BIO": 0, "XTAL": 1})
    df["resolution_f64"] = df["resolution"].astype("float64")
    df_sorted = df.sort_values(
        by=["pair_key", "label_priority", "resolution_f64"],
        ascending=[True, True, True],
    )

    # Keep first representative per pair
    representatives = df_sorted.drop_duplicates(subset=["pair_key"], keep="first")
    print(f"Unique protein pairs: {len(representatives)}")

    # Build output
    df_out = pd.DataFrame({
        "pair_key": representatives["pair_key"],
        "protein_A_id": representatives["uniprot_R"],
        "protein_B_id": representatives["uniprot_L"],
        "label": (representatives["label"] == "BIO").astype(int),
        "label_str": representatives["label"],
        "representative_id": representatives["id"],
        "pdb_id": representatives["pdb_id"],
        "split": representatives["split"],
        "complex_type": representatives["complex_type"],
        "method": representatives["method"],
        "resolution": representatives["resolution"],
        "release_date": representatives["release_date"],
        "protein_A_len": representatives["length1"],
        "protein_B_len": representatives["length2"],
        "buried_sasa": representatives["buried_sasa"],
        "intermolecular_contacts": representatives["intermolecular_contacts"],
        "planarity": representatives["planarity"],
    })

    # Add sequences if requested
    if require_sequences:
        seqs = load_sequences([
            SPROT_FASTA,
            "data/uniprot_sequences_all.fasta",
        ])
        df_out["protein_A_seq"] = df_out["protein_A_id"].map(seqs)
        df_out["protein_B_seq"] = df_out["protein_B_id"].map(seqs)
        before = len(df_out)
        df_out = df_out.dropna(subset=["protein_A_seq", "protein_B_seq"])
        df_out = df_out[
            (df_out["protein_A_seq"].str.len() > 0)
            & (df_out["protein_B_seq"].str.len() > 0)
        ]
        print(f"After requiring sequences: {len(df_out)} (dropped {before - len(df_out)})")

    os.makedirs("output", exist_ok=True)
    output_path = "output/ppi_unique_pairs.csv"
    df_out.to_csv(output_path, index=False)
    print(f"Saved {output_path}")

    # Summary
    label_counts = df_out["label_str"].value_counts()
    summary = {
        "total_pairs": int(len(df_out)),
        "bio_pairs": int(label_counts.get("BIO", 0)),
        "xtal_pairs": int(label_counts.get("XTAL", 0)),
        "n_unique_protein_A": int(df_out["protein_A_id"].nunique()),
        "n_unique_protein_B": int(df_out["protein_B_id"].nunique()),
        "split_distribution": df_out["split"].value_counts().to_dict(),
        "complex_type_distribution": df_out["complex_type"].value_counts().to_dict(),
    }
    summary_path = "output/ppi_unique_pairs_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"Saved {summary_path}")

    return df_out


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--require-sequences",
        action="store_true",
        help="Require both proteins have UniProt sequences",
    )
    args = parser.parse_args()

    build_unique_pairs(require_sequences=args.require_sequences)
