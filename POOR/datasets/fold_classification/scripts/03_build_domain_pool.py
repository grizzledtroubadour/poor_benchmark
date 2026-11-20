#!/usr/bin/env python3
"""Build the fold_classification domain pool.

Reconstruct a missing historical step: v2 split construction originally obtained the pool
from v1 C.A-level backups (output/splits_ca_backup/). Rebuild an equivalent pool from
primary outputs:

1. Input output/cath_s95_summary.csv: 62,901 successfully extracted S95 domains.
2. aa_seq uses authoritative domain sequences in cath-domain-seqs-S95.fa, identical to
         the v1 pool; PDB ATOM sequences omit loop residues and are not used as sequence sources.
3. Exclude five domains with length_expected >1000: 1ej6B00, 1u6gC00,
         3a6pA00, 3egwA02, 3gjxD00. As in v1, use CATH listed lengths,
         retaining cases such as 3w3uA01 with reference sequence >1000 but length_expected=999.
4. struct_file = {domain_id}.pdb; assert the file exists in pdbs/.

Output: output/domain_pool.csv, 62,896 rows; columns unique_id, aa_seq, struct_file.
Validate equivalence by row-level diff against the union of three v1 splits in output/splits_ca_backup/.

Usage:
    python scripts/03_build_domain_pool.py                      # Write output/domain_pool.csv.
    python scripts/03_build_domain_pool.py --out output/_verify_pool/domain_pool.csv
"""
import argparse
from pathlib import Path

import pandas as pd

TASK_DIR = Path(__file__).resolve().parent.parent
ROOT = TASK_DIR.parent.parent                       # Repository root (POOR)
SUMMARY_CSV = TASK_DIR / "output" / "cath_s95_summary.csv"
PDB_DIR = TASK_DIR / "pdbs"
CATH_SEQ_S95 = (
    ROOT / "data" / "CATHv44" / "sequence-data" / "cath-domain-seqs-S95.fa"
)

MAX_LEN = 1000  # Exclude domains with length_expected >1000.


def load_fasta(fasta_file: Path) -> dict:
    """Load CATH reference sequences keyed by domain_id (e.g. 12asA00)."""
    seqs = {}
    current_id = None
    current_seq = []
    with open(fasta_file) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_id is not None:
                    seqs[current_id] = "".join(current_seq)
                header = line[1:].strip()
                current_id = header.split("|")[-1].split("/")[0].strip()
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_id is not None:
            seqs[current_id] = "".join(current_seq)
    return seqs


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path,
                    default=TASK_DIR / "output" / "domain_pool.csv",
                    help="Output path (default output/domain_pool.csv)")
    args = ap.parse_args()

    summary = pd.read_csv(SUMMARY_CSV)
    fasta = load_fasta(CATH_SEQ_S95)
    print(f"Summary rows: {len(summary)}, FASTA entries: {len(fasta)}")

    dropped = summary[summary["length_expected"] > MAX_LEN]
    pool = summary[summary["length_expected"] <= MAX_LEN].copy()
    print(f"Excluded length_expected>{MAX_LEN}: {len(dropped)} "
          f"({', '.join(dropped['domain_id'])})")

    pool["unique_id"] = pool["domain_id"]
    pool["aa_seq"] = pool["unique_id"].map(fasta)
    assert pool["aa_seq"].notna().all(), "Domains missing CATH reference sequences"
    pool["struct_file"] = pool["unique_id"] + ".pdb"

    missing_pdb = [u for u in pool["unique_id"]
                   if not (PDB_DIR / f"{u}.pdb").exists()]
    assert not missing_pdb, f"Missing files in pdbs/: {missing_pdb[:5]}..."

    out = pool[["unique_id", "aa_seq", "struct_file"]].sort_values("unique_id")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    print(f"Domain pool saved: {args.out} ({len(out)} rows)")


if __name__ == "__main__":
    main()
