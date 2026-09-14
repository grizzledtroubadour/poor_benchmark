#!/usr/bin/env python3
"""Step 3 (v2): exact-duplicate removal (identical aa_seq only).

No sequence-identity clustering: tiny sequence differences may correspond to
different functions, so only completely identical sequences are merged.

Representative choice within a duplicate group: prefer a chain whose structure
file exists in output/chains_all/ (all should exist), then first occurrence.

Representative structures are COPIED to pdbs/ (complete dedup'd library for
future tasks' extra test sets).

Outputs:
  output/design_chains_dedup_exact.csv
  pdbs/{unique_id}.pdb|.cif
"""

import shutil
from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parent.parent
OUT_DIR = BASE / "output"
CHAINS_DIR = OUT_DIR / "chains_all"
PDBS_DIR = BASE / "pdbs"
PDBS_DIR.mkdir(parents=True, exist_ok=True)


def main():
    df = pd.read_csv(OUT_DIR / "design_chains_all.csv")
    print(f"input chains: {len(df)}")

    def has_struct(uid):
        return any((CHAINS_DIR / f"{uid}{e}").exists() for e in (".pdb", ".cif"))

    df["has_struct"] = df["unique_id"].map(has_struct)
    # prefer rows with structure, then keep first occurrence in file order
    df = df.sort_values("has_struct", ascending=False, kind="stable")
    dedup = df.drop_duplicates("aa_seq", keep="first").copy()
    dedup = dedup.sort_index(kind="stable")
    print(f"after exact dedup: {len(dedup)} "
          f"(merged {len(df) - len(dedup)} identical-sequence records)")
    assert dedup["has_struct"].all(), "representative without structure!"

    copied = 0
    struct_files = {}
    for uid in dedup["unique_id"]:
        for ext in (".pdb", ".cif"):
            src = CHAINS_DIR / f"{uid}{ext}"
            if src.exists():
                dst = PDBS_DIR / src.name
                if not dst.exists():
                    shutil.copy2(src, dst)
                struct_files[uid] = src.name
                copied += 1
                break
    print(f"structures copied to pdbs/: {copied}")

    dedup["struct_file"] = dedup["unique_id"].map(struct_files)
    dedup = dedup.drop(columns=["has_struct"])
    dedup.to_csv(OUT_DIR / "design_chains_dedup_exact.csv", index=False)
    print(f"design_chains_dedup_exact.csv: {len(dedup)} rows")
    print(f"  unambiguous label_candidate: "
          f"{(dedup['label_candidate'].fillna('') != '').sum()}")


if __name__ == "__main__":
    main()
