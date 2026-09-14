#!/usr/bin/env python3
"""
Fetch UniProt sequences for all proteins in PINDER index.

Outputs:
  - data/uniprot_sequences_all.fasta
  - output/uniprot_sequence_coverage.json

Supports resumable downloads: progress is saved to a temporary FASTA every
50 batches, and reloaded on restart.
"""

import json
import os
import time
from io import StringIO
from pathlib import Path

import pandas as pd
import requests
from Bio import SeqIO
from tqdm import tqdm


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


def save_sequences(seqs, all_ids, output_path):
    """Save sequences for all PINDER IDs to a FASTA file."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as f:
        for uid in sorted(seqs.keys()):
            if uid in all_ids:
                f.write(f">{uid}\n{seqs[uid]}\n")


def fetch_uniprot_sequences(ids_to_fetch, all_ids, batch_size=100, sleep=0.3):
    """Fetch FASTA sequences from UniProt stream API in batches.

    Saves incremental progress to data/uniprot_sequences_all.tmp.fasta
    every 50 batches so the job can resume after interruption.
    """
    base_url = "https://rest.uniprot.org/uniprotkb/stream"
    tmp_fasta = "data/uniprot_sequences_all.tmp.fasta"
    fetched = {}
    checkpoint_every = 50

    for i in tqdm(range(0, len(ids_to_fetch), batch_size), desc="Fetching UniProt"):
        batch = ids_to_fetch[i : i + batch_size]
        query = "+OR+".join([f"accession:{uid}" for uid in batch])
        url = f"{base_url}?format=fasta&query={query}"
        try:
            resp = requests.get(url, timeout=120)
            if resp.status_code == 200:
                for record in SeqIO.parse(StringIO(resp.text), "fasta"):
                    parts = record.id.split("|")
                    uid = parts[1] if len(parts) >= 2 else record.id
                    fetched[uid] = str(record.seq)
            else:
                print(f"\nBatch {i//batch_size} failed: {resp.status_code} - {resp.text[:200]}")
        except Exception as e:
            print(f"\nBatch {i//batch_size} error: {e}")

        # Incremental checkpoint
        if (i // batch_size + 1) % checkpoint_every == 0:
            save_sequences(fetched, all_ids, tmp_fasta)
            print(
                f"\nCheckpoint: saved {len(fetched)} newly fetched sequences to {tmp_fasta}"
            )

        time.sleep(sleep)

    return fetched


def main():
    os.makedirs("data", exist_ok=True)
    os.makedirs("output", exist_ok=True)

    # Load index
    idx = pd.read_parquet("data/index.parquet")
    all_ids = (
        pd.concat([idx["uniprot_R"], idx["uniprot_L"]])
        .dropna()
        .unique()
        .tolist()
    )
    all_ids = [uid for uid in all_ids if uid and uid != "UNDEFINED"]
    all_ids_set = set(all_ids)
    print(f"Total valid UniProt IDs in PINDER: {len(all_ids)}")

    # Load existing sequences (final + tmp + test_val + local)
    # Optional local Swiss-Prot preload (project-root shared data dir; skipped
    # automatically when absent, sequences are then fetched from UniProt API)
    sprot_fasta = str(
        Path(__file__).resolve().parents[3] / "data" / "uniprot" / "uniprot_sprot.fasta"
    )
    existing_paths = [
        sprot_fasta,
        "data/uniprot_sequences_test_val.fasta",
        "data/uniprot_sequences_all.fasta",
        "data/uniprot_sequences_all.tmp.fasta",
    ]
    seqs = load_existing_sequences(existing_paths)
    print(f"Already available sequences: {len(seqs)}")

    ids_to_fetch = [uid for uid in all_ids if uid not in seqs]
    print(f"Need to fetch: {len(ids_to_fetch)}")

    if ids_to_fetch:
        fetched = fetch_uniprot_sequences(ids_to_fetch, all_ids_set, batch_size=100, sleep=0.3)
        seqs.update(fetched)
        print(f"Total sequences after fetching: {len(seqs)}")
        # Remove tmp file after successful completion
        if os.path.exists("data/uniprot_sequences_all.tmp.fasta"):
            os.remove("data/uniprot_sequences_all.tmp.fasta")
    else:
        print("All sequences already available.")

    # Save sequences for all PINDER IDs
    output_fasta = "data/uniprot_sequences_all.fasta"
    save_sequences(seqs, all_ids_set, output_fasta)
    print(f"Saved {output_fasta}")

    # Coverage statistics
    found = sum(1 for uid in all_ids if uid in seqs)
    coverage = {
        "total_ids": len(all_ids),
        "found": int(found),
        "missing": int(len(all_ids) - found),
        "coverage_pct": round(found / len(all_ids) * 100, 2),
        "missing_ids": [uid for uid in all_ids if uid not in seqs][:100],
    }
    with open("output/uniprot_sequence_coverage.json", "w") as f:
        json.dump(coverage, f, indent=2)
    print(
        f"Coverage: {found}/{len(all_ids)} ({coverage['coverage_pct']}%)"
    )


if __name__ == "__main__":
    main()
