#!/usr/bin/env python3
"""EC annotation: generate standardized EC sets per sample (ec-function-ood-annotation step 1).

Two input modes, detected from column names:
- Chain mode: `pdb`, `chain` columns -> SIFTS chain-level EC（sifts_chain_ec.tsv.gz）
- UniProt mode: `uniprot` column -> (1) SIFTS reverse mapping (pdb_chain_uniprot
  SP_PRIMARY -> chains -> union of chain ECs); (2) batch-query remaining accessions through UniProt REST
  (fields=accession,ec), writing the project-wide data/uniprot_ec_cache.json cache;
  `--offline` uses cache only, without API access.

Output: append `ec_numbers` (semicolon-separated full EC numbers) and `has_ec` to input columns.
ec_ood.py handles L3/L4 truncation; this script only produces raw EC sets.

Usage:
    python ec_annotate.py --input my_chains.csv --out my_ec_labels.csv
    python ec_annotate.py --input my_accs.csv --out my_ec_labels.csv --offline
"""
import argparse
import csv
import gzip
import json
import time
from collections import defaultdict
from pathlib import Path

import requests

POOR = Path(__file__).resolve().parents[3]
SIFTS_EC_TSV = POOR / "data" / "sifts" / "sifts_chain_ec.tsv.gz"
SIFTS_PDB_CHAIN_UNIPROT = POOR / "data" / "sifts" / "pdb_chain_uniprot.tsv.gz"
DEFAULT_CACHE = POOR / "data" / "uniprot_ec_cache.json"
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
BATCH_SIZE = 100


def load_sifts_chain_ec():
    m = defaultdict(set)
    with gzip.open(SIFTS_EC_TSV, "rt") as f:
        next(f)
        for line in f:
            p = line.strip().split("\t")
            if len(p) >= 4 and p[3] and p[3] != "-":
                m[(p[0].lower(), p[1])].add(p[3])
    return m


def load_uniprot_to_chains():
    m = defaultdict(list)
    with gzip.open(SIFTS_PDB_CHAIN_UNIPROT, "rt") as f:
        next(f)  # Comment line
        next(f)  # Header
        for line in f:
            p = line.strip().split("\t")
            if len(p) >= 3:
                m[p[2]].append((p[0].lower(), p[1]))
    return m


def fetch_ec_batch(ids):
    """Batch-query ECs through UniProt REST, as in ppi_prediction/scripts/fetch_uniprot_ec.py."""
    query = " OR ".join(f"accession:{x}" for x in ids)
    params = {"query": query, "format": "tsv", "fields": "accession,ec", "size": len(ids)}
    r = requests.get(UNIPROT_SEARCH_URL, params=params, timeout=60)
    r.raise_for_status()
    out = {}
    for line in r.text.strip().split("\n")[1:]:
        parts = line.split("\t")
        if len(parts) >= 2:
            out[parts[0]] = [x.strip() for x in parts[1].split(";") if x.strip()]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--offline", action="store_true", help="Use cache only, without UniProt API access.")
    ap.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    args = ap.parse_args()

    with open(args.input) as f:
        rows = list(csv.DictReader(f))
    fieldnames = list(rows[0].keys())
    if {"pdb", "chain"} <= set(fieldnames):
        mode = "chain"
    elif "uniprot" in fieldnames:
        mode = "uniprot"
    else:
        raise SystemExit("Unrecognized input mode: require (pdb, chain) columns or a uniprot column.")
    print(f"mode: {mode}; input rows: {len(rows)}")

    chain_ec = load_sifts_chain_ec()
    ec_of_row = {}
    if mode == "chain":
        for r in rows:
            ec_of_row[r["pdb"], r["chain"]] = chain_ec.get((r["pdb"].lower(), r["chain"]), set())
        print(f"chains with EC (SIFTS): {sum(bool(v) for v in ec_of_row.values())}/{len(rows)}")
    else:
        accs = sorted({r["uniprot"] for r in rows})
        u2c = load_uniprot_to_chains()
        for a in accs:
            s = set()
            for ck in u2c.get(a, []):
                s |= chain_ec.get(ck, set())
            ec_of_row[a] = s
        covered = sum(bool(ec_of_row[a]) for a in accs)
        print(f"accessions: {len(accs)}, with EC via SIFTS reverse mapping: {covered}")

        # Fill missing annotations through UniProt API.
        cache = json.load(open(args.cache)) if args.cache.exists() else {}
        missing = [a for a in accs if not ec_of_row[a] and a not in cache]
        print(f"cached: {sum(a in cache for a in accs)}, to query: {len(missing)}")
        if missing and not args.offline:
            for i in range(0, len(missing), BATCH_SIZE):
                batch = missing[i:i + BATCH_SIZE]
                try:
                    res = fetch_ec_batch(batch)
                except Exception as e:
                    print(f"  batch {i} query failed: {e}")
                    res = {}
                for a in batch:
                    cache[a] = res.get(a, [])
                time.sleep(0.5)
                if (i // BATCH_SIZE) % 10 == 0:
                    print(f"  {i + len(batch)}/{len(missing)}")
            args.cache.parent.mkdir(parents=True, exist_ok=True)
            json.dump(cache, open(args.cache, "w"))
        for a in accs:
            if not ec_of_row[a]:
                ec_of_row[a] = set(cache.get(a, []))
        print(f"with EC after API completion: {sum(bool(ec_of_row[a]) for a in accs)}/{len(accs)}")

    if "ec_numbers" not in fieldnames:
        fieldnames += ["ec_numbers", "has_ec"]
    for r in rows:
        key = (r["pdb"], r["chain"]) if mode == "chain" else r["uniprot"]
        ecs = ec_of_row.get(key, set())
        r["ec_numbers"] = ";".join(sorted(ecs))
        r["has_ec"] = str(bool(ecs))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    print(f"-> {args.out}")


if __name__ == "__main__":
    main()
