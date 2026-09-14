#!/usr/bin/env python3
"""
Deduplicate protein pairs by sequence identity using mmseqs2.

Rule:
  - Two pairs (A,B) and (A',B') are considered redundant if:
      seq_identity(A, A') > 95% AND seq_identity(B, B') > 95%
  - Chain order matters for the identity check, but the pair itself is unordered
    (i.e. {A,B} = {B,A}).
  - If a pair of sequence clusters has both BIO and XTAL representatives,
    keep BIO.

Outputs:
  - output/ppi_unique_pairs_seq95_dedup.csv
  - output/ppi_unique_pairs_seq95_dedup_summary.json
"""

import json
import os
import subprocess
import tempfile

import pandas as pd
from Bio import SeqIO


def load_sequences(fasta_path):
    seqs = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        parts = record.id.split("|")
        uid = parts[1] if len(parts) >= 2 else record.id
        seqs[uid] = str(record.seq)
    return seqs


def run_mmseqs_cluster(input_fasta, output_tsv, min_seq_id=0.95, threads=8):
    """Run mmseqs2 easy-cluster or full workflow and return cluster TSV."""
    tmpdir = tempfile.mkdtemp(prefix="mmseqs_")
    db = os.path.join(tmpdir, "seqDB")
    clu = os.path.join(tmpdir, "cluDB")
    tmp = os.path.join(tmpdir, "tmp")
    os.makedirs(tmp, exist_ok=True)

    # createdb
    subprocess.run(
        ["mmseqs", "createdb", input_fasta, db],
        check=True, capture_output=True, text=True
    )

    # cluster
    subprocess.run(
        [
            "mmseqs", "cluster", db, clu, tmp,
            "--min-seq-id", str(min_seq_id),
            "-c", "0.8",  # coverage threshold
            "--cov-mode", "0",
            "--threads", str(threads),
        ],
        check=True, capture_output=True, text=True
    )

    # createtsv
    subprocess.run(
        ["mmseqs", "createtsv", db, db, clu, output_tsv],
        check=True, capture_output=True, text=True
    )

    return output_tsv


def load_cluster_map(tsv_path):
    """Load mmseqs2 cluster TSV into dict: member -> representative."""
    cluster_map = {}
    with open(tsv_path) as f:
        for line in f:
            rep, member = line.strip().split("\t")
            cluster_map[member] = rep
    return cluster_map


def build_seq_cluster_pair_key(protein_A_id, protein_B_id, cluster_map):
    """
    Map each protein to its sequence cluster representative,
    then create an unordered pair key.
    """
    rep_a = cluster_map.get(protein_A_id)
    rep_b = cluster_map.get(protein_B_id)
    if rep_a is None or rep_b is None:
        return None
    return "__".join(sorted([rep_a, rep_b]))


def main():
    # Load unique pairs
    pairs = pd.read_csv("output/ppi_unique_pairs.csv")
    print(f"Loaded {len(pairs)} unique pairs")

    # Load sequences
    seqs = load_sequences("data/uniprot_sequences_all.fasta")
    print(f"Loaded {len(seqs)} sequences")

    # Filter pairs with sequences
    pairs["has_seq_A"] = pairs["protein_A_id"].isin(seqs)
    pairs["has_seq_B"] = pairs["protein_B_id"].isin(seqs)
    pairs_with_seq = pairs[pairs["has_seq_A"] & pairs["has_seq_B"]].copy()
    print(f"Pairs with both sequences: {len(pairs_with_seq)}")

    # Write FASTA of used proteins only
    used_proteins = set(pairs_with_seq["protein_A_id"]) | set(pairs_with_seq["protein_B_id"])
    used_fasta = "data/ppi_unique_pairs_proteins.fasta"
    with open(used_fasta, "w") as f:
        for uid in sorted(used_proteins):
            if uid in seqs:
                f.write(f">{uid}\n{seqs[uid]}\n")
    print(f"Wrote {len(used_proteins)} protein sequences to {used_fasta}")

    # Run mmseqs2 clustering at 95% identity
    cluster_tsv = "output/mmseqs_clusters_seq95.tsv"
    if not os.path.exists(cluster_tsv):
        print("Running mmseqs2 clustering at 95% identity...")
        run_mmseqs_cluster(used_fasta, cluster_tsv, min_seq_id=0.95)
    else:
        print(f"Using existing cluster file {cluster_tsv}")

    cluster_map = load_cluster_map(cluster_tsv)
    print(f"Loaded cluster map: {len(cluster_map)} members, {len(set(cluster_map.values()))} clusters")

    # Map pairs to sequence cluster representatives
    pairs_with_seq["seq_cluster_key"] = pairs_with_seq.apply(
        lambda r: build_seq_cluster_pair_key(
            r["protein_A_id"], r["protein_B_id"], cluster_map
        ),
        axis=1,
    )

    # Drop pairs that couldn't be mapped
    pairs_with_seq = pairs_with_seq.dropna(subset=["seq_cluster_key"]).copy()
    print(f"Pairs after cluster mapping: {len(pairs_with_seq)}")

    # Sort by label (BIO first) and resolution
    pairs_with_seq["label_priority"] = pairs_with_seq["label_str"].map({"BIO": 0, "XTAL": 1})
    pairs_with_seq["resolution_f64"] = pairs_with_seq["resolution"].astype("float64")
    pairs_sorted = pairs_with_seq.sort_values(
        by=["seq_cluster_key", "label_priority", "resolution_f64"],
        ascending=[True, True, True],
    )

    # Deduplicate by sequence cluster pair key
    dedup = pairs_sorted.drop_duplicates(subset=["seq_cluster_key"], keep="first")
    print(f"After 95% sequence identity dedup: {len(dedup)}")

    # Drop helper columns and save
    dedup_out = dedup.drop(
        columns=[
            "has_seq_A", "has_seq_B", "label_priority", "resolution_f64",
            "seq_cluster_key",
        ],
        errors="ignore",
    )

    output_path = "output/ppi_unique_pairs_seq95_dedup.csv"
    dedup_out.to_csv(output_path, index=False)
    print(f"Saved {output_path}")

    # Summary
    label_counts = dedup_out["label_str"].value_counts()
    summary = {
        "input_pairs": int(len(pairs)),
        "pairs_with_sequences": int(len(pairs_with_seq)),
        "sequence_clusters": int(len(set(cluster_map.values()))),
        "dedup_pairs": int(len(dedup_out)),
        "bio_pairs": int(label_counts.get("BIO", 0)),
        "xtal_pairs": int(label_counts.get("XTAL", 0)),
        "split_distribution": dedup_out["split"].value_counts().to_dict(),
        "complex_type_distribution": dedup_out["complex_type"].value_counts().to_dict(),
        "removed_pairs": int(len(pairs_with_seq) - len(dedup_out)),
    }
    summary_path = "output/ppi_unique_pairs_seq95_dedup_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"Saved {summary_path}")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
