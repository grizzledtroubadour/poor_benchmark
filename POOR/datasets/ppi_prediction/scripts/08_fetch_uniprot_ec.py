#!/usr/bin/env python3
"""
Fetch EC numbers for PPI UniProt IDs via UniProt REST API.

Input: ppi_prediction/data/uniprot_sequences_all.fasta
Output: ppi_prediction/output/uniprot_ec_mapping.csv

EC numbers are returned from UniProt as a semicolon-separated string.
"""

import argparse
import csv
import sys
import time
from pathlib import Path

import requests
from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_FASTA = ROOT / "data" / "uniprot_sequences_all.fasta"
DEFAULT_OUTPUT = ROOT / "output" / "uniprot_ec_mapping.csv"
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
BATCH_SIZE = 100
SLEEP_BETWEEN_BATCHES = 0.5


def load_uniprot_ids(fasta_path: Path):
    ids = []
    for rec in SeqIO.parse(str(fasta_path), "fasta"):
        ids.append(rec.id)
    return sorted(set(ids))


def fetch_ec_batch(ids: list[str]):
    query = " OR ".join(f"accession:{x}" for x in ids)
    params = {
        "query": query,
        "format": "tsv",
        "fields": "accession,ec",
        "size": len(ids),
    }
    r = requests.get(UNIPROT_SEARCH_URL, params=params, timeout=60)
    r.raise_for_status()
    lines = r.text.strip().split("\n")
    results = {}
    for line in lines[1:]:
        parts = line.split("\t")
        if len(parts) >= 2:
            uniprot_id = parts[0]
            ec_str = parts[1]
            if ec_str:
                ecs = [x.strip() for x in ec_str.split(";") if x.strip()]
                results[uniprot_id] = ecs
    return results


def main():
    parser = argparse.ArgumentParser(description="Fetch EC numbers from UniProt")
    parser.add_argument("--fasta", type=Path, default=DEFAULT_FASTA)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--batch-size", type=int, default=BATCH_SIZE)
    args = parser.parse_args()

    if not args.fasta.exists():
        print(f"FASTA not found: {args.fasta}", file=sys.stderr)
        sys.exit(1)

    uniprot_ids = load_uniprot_ids(args.fasta)
    print(f"Loaded {len(uniprot_ids)} UniProt IDs from {args.fasta}")

    args.output.parent.mkdir(parents=True, exist_ok=True)

    # Resume support
    existing = {}
    if args.output.exists():
        with open(args.output) as f:
            for r in csv.DictReader(f):
                if r["ec_numbers"]:
                    existing[r["uniprot_id"]] = [x.strip() for x in r["ec_numbers"].split(";") if x.strip()]
        print(f"Resuming: {len(existing)} IDs already fetched")

    remaining = [x for x in uniprot_ids if x not in existing]
    print(f"Remaining to fetch: {len(remaining)}")

    ec_map = dict(existing)
    total_batches = (len(remaining) + args.batch_size - 1) // args.batch_size

    with open(args.output, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["uniprot_id", "ec_numbers"])

        for i in range(total_batches):
            batch = remaining[i * args.batch_size : (i + 1) * args.batch_size]
            try:
                batch_results = fetch_ec_batch(batch)
                ec_map.update(batch_results)
                print(
                    f"Batch {i + 1}/{total_batches}: fetched {len(batch_results)} "
                    f"({sum(1 for v in batch_results.values() if v)} with EC) "
                    f"total mapped: {len(ec_map)}/{len(uniprot_ids)}"
                )
            except Exception as e:
                print(f"Batch {i + 1} failed: {e}", file=sys.stderr)
                # Continue with next batch; missing IDs can be retried later

            # Flush current results
            fh.seek(0)
            fh.truncate()
            writer.writerow(["uniprot_id", "ec_numbers"])
            for uid in uniprot_ids:
                ecs = ec_map.get(uid, [])
                writer.writerow([uid, ";".join(ecs) if ecs else ""])

            if i < total_batches - 1:
                time.sleep(SLEEP_BETWEEN_BATCHES)

    with_ec = sum(1 for ecs in ec_map.values() if ecs)
    print(f"Done. Total IDs: {len(ec_map)}, with EC: {with_ec} ({with_ec / len(ec_map) * 100:.1f}%)")
    print(f"Saved to {args.output}")


if __name__ == "__main__":
    main()
