#!/usr/bin/env python3
"""Step 4: Build the fold-prediction extra test set base table.

On the dedup'd (95%) design chains:
1. Label filter: unambiguous entry-level C.A.T (label_candidate) that is
   one of the 1,163 fold-train Topology labels -> `label`.
2. train/val overlap removal (either rule triggers removal):
   a. pdb_id appears among fold train/val domain PDB ids
   b. mmseqs2 easy-search -s 7.5 vs train+val, max pident >= 95

The mmseqs search is reused by Step 5 (OOD seq axes).

Outputs:
  output/design_mmseqs_vs_trainval.m8
  output/cath_design_base.csv
  output/removed_overlap.csv
"""

import subprocess
from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parent.parent
POOR = BASE.parent.parent
OUT_DIR = BASE / "output"
FOLD = POOR / "datasets" / "fold_classification"

MMSEQS_M8 = OUT_DIR / "design_mmseqs_vs_trainval.m8"
TMP = OUT_DIR / "s04_tmp"
TMP.mkdir(exist_ok=True)

PIDENT_CUTOFF = 95.0


def write_fasta(df, path):
    with open(path, "w") as f:
        for uid, seq in zip(df["unique_id"], df["aa_seq"]):
            f.write(f">{uid}\n{seq}\n")


def main():
    dedup = pd.read_csv(OUT_DIR / "design_chains_dedup_exact.csv")
    train = pd.read_csv(FOLD / "splits" / "cath_train.csv", dtype={"label": str})
    val = pd.read_csv(FOLD / "splits" / "cath_val.csv", dtype={"label": str})
    train_labels = set(train["label"])
    tv_pdbs = (set(train["unique_id"].str[:4].str.lower())
               | set(val["unique_id"].str[:4].str.lower()))

    # ---- 1. label filter
    cand = dedup[dedup["label_candidate"].fillna("") != ""].copy()
    cand["label"] = cand["label_candidate"]
    n_uni = len(cand)
    cand = cand[cand["label"].isin(train_labels)].copy()
    print(f"label filter: {n_uni} unambiguous -> {len(cand)} with label in "
          f"train ({n_uni - len(cand)} dropped: label not in train label set)")

    # ---- 2a. PDB-id overlap
    cand["overlap_pdbid"] = cand["pdb_id"].isin(tv_pdbs)

    # ---- 2b. sequence-identity overlap (mmseqs vs train+val)
    if MMSEQS_M8.exists():
        print(f"mmseqs2: reusing {MMSEQS_M8.name}")
    else:
        q_fa = TMP / "design.fasta"
        t_fa = TMP / "trainval.fasta"
        write_fasta(cand, q_fa)
        write_fasta(pd.concat([train, val], ignore_index=True), t_fa)
        cmd = [
            "mmseqs", "easy-search", str(q_fa), str(t_fa), str(MMSEQS_M8),
            str(TMP / "mmseqs_tmp"), "-s", "7.5",
            "--format-output", "query,target,pident",
        ]
        print("mmseqs2:", " ".join(cmd))
        subprocess.run(cmd, check=True)

    hits = pd.read_csv(MMSEQS_M8, sep="\t", header=None,
                       names=["query", "target", "pident"], usecols=[0, 1, 2])
    max_pident = hits.groupby("query")["pident"].max()
    cand["max_pident_tv"] = cand["unique_id"].map(max_pident)
    cand["overlap_seq95"] = cand["max_pident_tv"].fillna(0.0) >= PIDENT_CUTOFF

    removed = cand[cand["overlap_pdbid"] | cand["overlap_seq95"]].copy()
    kept = cand[~(cand["overlap_pdbid"] | cand["overlap_seq95"])].copy()

    removed.to_csv(OUT_DIR / "removed_overlap.csv", index=False)
    print(f"overlap removal: pdb_id {cand['overlap_pdbid'].sum()}, "
          f"seq>=95% {cand['overlap_seq95'].sum()}, "
          f"total removed {len(removed)} (union)")

    kept.to_csv(OUT_DIR / "cath_design_base.csv", index=False)
    print(f"\ncath_design_base.csv: {len(kept)} rows")
    print(f"  labels covered: {kept['label'].nunique()}")
    print(f"  seq_len < 60: {(kept['seq_len'] < 60).sum()}")


if __name__ == "__main__":
    main()
