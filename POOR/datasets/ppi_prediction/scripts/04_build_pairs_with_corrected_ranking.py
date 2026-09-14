#!/usr/bin/env python3
"""
Build unique protein pairs with corrected representative selection.

Criteria (in order):
  1. Prefer samples with apo or predicted monomer structures
     - Apo with high quality > Apo with low quality > Predicted only > Holo only
  2. Prefer higher effective resolution
     - NOTE: PINDER metadata does not contain apo/predicted structure resolution.
       We use holo structure resolution as a proxy. If real apo resolution is
       required, apo/predicted PDB files must be downloaded separately.
  3. Prefer longer total sequence length
  4. Prefer newer release_date

Outputs:
  - output/ppi_pairs_corrected_uniprot.csv
  - output/ppi_pairs_corrected_seq95.csv
  - output/ppi_pairs_corrected_summary.json
"""

import json
import os
import subprocess
import tempfile
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


def run_mmseqs_cluster(input_fasta, output_tsv, min_seq_id=0.95, threads=8):
    # Re-run clustering if the output is older than the input FASTA
    if os.path.exists(output_tsv):
        out_mtime = os.path.getmtime(output_tsv)
        in_mtime = os.path.getmtime(input_fasta)
        if out_mtime >= in_mtime:
            return output_tsv
        print(f"Input FASTA newer than cluster TSV, re-clustering...")
    tmpdir = tempfile.mkdtemp(prefix="mmseqs_")
    db = os.path.join(tmpdir, "seqDB")
    clu = os.path.join(tmpdir, "cluDB")
    tmp = os.path.join(tmpdir, "tmp")
    os.makedirs(tmp, exist_ok=True)

    subprocess.run(
        ["mmseqs", "createdb", input_fasta, db],
        check=True, capture_output=True, text=True
    )
    subprocess.run(
        [
            "mmseqs", "cluster", db, clu, tmp,
            "--min-seq-id", str(min_seq_id),
            "-c", "0.8",
            "--cov-mode", "0",
            "--threads", str(threads),
        ],
        check=True, capture_output=True, text=True
    )
    subprocess.run(
        ["mmseqs", "createtsv", db, db, clu, output_tsv],
        check=True, capture_output=True, text=True
    )
    return output_tsv


def load_cluster_map(tsv_path):
    cluster_map = {}
    with open(tsv_path) as f:
        for line in f:
            rep, member = line.strip().split("\t")
            cluster_map[member] = rep
    return cluster_map


def structure_priority(row):
    """
    Higher = better.
    3: has high-quality apo
    2: has low-quality apo (but no high)
    1: has predicted only
    0: holo only
    """
    has_apo_R = bool(row["apo_R"])
    has_apo_L = bool(row["apo_L"])
    has_pred_R = bool(row["predicted_R"])
    has_pred_L = bool(row["predicted_L"])

    apo_R_high = has_apo_R and row["apo_R_quality"] == "high"
    apo_L_high = has_apo_L and row["apo_L_quality"] == "high"
    apo_R_low = has_apo_R and row["apo_R_quality"] == "low"
    apo_L_low = has_apo_L and row["apo_L_quality"] == "low"

    if apo_R_high or apo_L_high:
        return 3
    if apo_R_low or apo_L_low:
        return 2
    if has_pred_R or has_pred_L:
        return 1
    return 0


def prepare_data():
    idx = pd.read_parquet("data/index.parquet")
    meta = pd.read_parquet("data/metadata.parquet")

    df = idx.merge(
        meta[
            [
                "id",
                "label",
                "method",
                "date",
                "release_date",
                "resolution",
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

    # Filter BIO/XTAL
    df = df[df["label"].isin(["BIO", "XTAL"])].copy()

    # Filter valid UniProt IDs
    df = df[
        (df["uniprot_R"].notna()) & (df["uniprot_L"].notna())
        & (df["uniprot_R"] != "UNDEFINED") & (df["uniprot_L"] != "UNDEFINED")
    ].copy()

    # Unordered pair key
    df["pair_key"] = df.apply(
        lambda r: "__".join(sorted([r["uniprot_R"], r["uniprot_L"]])), axis=1
    )

    # Structure priority score
    df["structure_priority"] = df.apply(structure_priority, axis=1)

    # Total length
    df["total_length"] = df["length1"].fillna(0) + df["length2"].fillna(0)

    # Resolution as float64 (holo resolution proxy)
    df["resolution_f64"] = df["resolution"].astype("float64")

    # Release date parsed
    df["release_date_parsed"] = pd.to_datetime(df["release_date"], errors="coerce")
    min_date = df["release_date_parsed"].min()
    fill_date = pd.Timestamp("1900-01-01") if pd.isna(min_date) else min_date
    df["release_date_rank"] = df["release_date_parsed"].fillna(fill_date)

    # Label priority
    df["label_priority"] = df["label"].map({"BIO": 0, "XTAL": 1})

    return df


def dedup_by_key(df, key_col, sort_cols, ascending):
    df_sorted = df.sort_values(by=sort_cols, ascending=ascending)
    return df_sorted.drop_duplicates(subset=[key_col], keep="first")


def build_output(df_dedup):
    has_apo = ((df_dedup["apo_R"].fillna(False)) | (df_dedup["apo_L"].fillna(False))).astype(int)
    has_pred = ((df_dedup["predicted_R"].fillna(False)) | (df_dedup["predicted_L"].fillna(False))).astype(int)

    return pd.DataFrame({
        "pair_key": df_dedup["pair_key"],
        "protein_A_id": df_dedup["uniprot_R"],
        "protein_B_id": df_dedup["uniprot_L"],
        "label": (df_dedup["label"] == "BIO").astype(int),
        "label_str": df_dedup["label"],
        "representative_id": df_dedup["id"],
        "pdb_id": df_dedup["pdb_id"],
        "split": df_dedup["split"],
        "complex_type": df_dedup["complex_type"],
        "method": df_dedup["method"],
        "resolution": df_dedup["resolution"],
        "release_date": df_dedup["release_date"],
        "protein_A_len": df_dedup["length1"],
        "protein_B_len": df_dedup["length2"],
        "total_length": df_dedup["total_length"],
        "structure_priority": df_dedup["structure_priority"],
        "has_apo": has_apo,
        "has_predicted": has_pred,
        "buried_sasa": df_dedup["buried_sasa"],
        "intermolecular_contacts": df_dedup["intermolecular_contacts"],
        "planarity": df_dedup["planarity"],
    })


def add_sequence_clusters(df, cluster_map):
    df = df.copy()
    df["seq_cluster_key"] = df.apply(
        lambda r: "__".join(sorted([
            cluster_map.get(r["uniprot_R"], r["uniprot_R"]),
            cluster_map.get(r["uniprot_L"], r["uniprot_L"]),
        ])),
        axis=1,
    )
    return df


def main():
    os.makedirs("output", exist_ok=True)

    df = prepare_data()
    print(f"Prepared {len(df)} BIO/XTAL rows")

    # UniProt pair dedup with corrected ranking
    df_uniprot = dedup_by_key(
        df,
        key_col="pair_key",
        sort_cols=[
            "pair_key", "label_priority", "structure_priority",
            "resolution_f64", "total_length", "release_date_rank"
        ],
        ascending=[True, True, False, True, False, False],
    )
    print(f"UniProt pair dedup: {len(df_uniprot)}")

    # Load sequences and cluster
    seqs = load_sequences([SPROT_FASTA, "data/uniprot_sequences_all.fasta"])
    used_proteins = set(df_uniprot["uniprot_R"]) | set(df_uniprot["uniprot_L"])
    used_proteins = {p for p in used_proteins if p in seqs}
    fasta_path = "data/ppi_pairs_corrected_proteins.fasta"
    with open(fasta_path, "w") as f:
        for uid in sorted(used_proteins):
            f.write(f">{uid}\n{seqs[uid]}\n")
    print(f"Wrote {len(used_proteins)} sequences")

    cluster_tsv = "output/mmseqs_clusters_corrected_seq95.tsv"
    run_mmseqs_cluster(fasta_path, cluster_tsv, min_seq_id=0.95)
    cluster_map = load_cluster_map(cluster_tsv)
    print(f"Sequence clusters: {len(set(cluster_map.values()))}")

    # Sequence cluster dedup
    df_seq = add_sequence_clusters(df_uniprot, cluster_map)
    df_seq_dedup = dedup_by_key(
        df_seq,
        key_col="seq_cluster_key",
        sort_cols=[
            "seq_cluster_key", "label_priority", "structure_priority",
            "resolution_f64", "total_length", "release_date_rank"
        ],
        ascending=[True, True, False, True, False, False],
    )
    print(f"Sequence cluster dedup: {len(df_seq_dedup)}")

    # Build outputs
    out_uniprot = build_output(df_uniprot)
    out_seq = build_output(df_seq_dedup)

    out_uniprot.to_csv("output/ppi_pairs_corrected_uniprot.csv", index=False)
    out_seq.to_csv("output/ppi_pairs_corrected_seq95.csv", index=False)

    # Summary
    label_counts = out_seq["label_str"].value_counts()
    summary = {
        "uniprot_dedup": int(len(out_uniprot)),
        "seq95_dedup": int(len(out_seq)),
        "bio_pairs": int(label_counts.get("BIO", 0)),
        "xtal_pairs": int(label_counts.get("XTAL", 0)),
        "structure_priority_distribution": out_seq["structure_priority"].value_counts().to_dict(),
        "has_apo": int(out_seq["has_apo"].sum()),
        "has_predicted": int(out_seq["has_predicted"].sum()),
        "split_distribution": out_seq["split"].value_counts().to_dict(),
        "complex_type_distribution": out_seq["complex_type"].value_counts().to_dict(),
    }
    with open("output/ppi_pairs_corrected_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
