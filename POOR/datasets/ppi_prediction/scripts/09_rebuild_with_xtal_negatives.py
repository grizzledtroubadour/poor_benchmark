#!/usr/bin/env python3
"""
Rebuild PPI dataset with XTAL hard negatives + random-pairing fill.

This script replaces the old pipeline (08_time_split → 21_resample_negative_pairs)
with a single-pass rebuild that:
  1. Splits BIO + XTAL by time (70:10:20, based on BIO sample count)
  2. Hard-filters any chain >2000 aa (global, before splitting)
  3. train/val: keep only chains ≤1000 aa
  4. test: keep all (including >1000, marked OOD_ExtremeLong later),
     filter to both_chains_apo_or_predicted
  5. XTAL pairs become label=0 negatives; shortfall filled with random pairs
     from the same split's protein pool (allowing self-pairs to break the
     homodimer shortcut)
  6. Adds pair_type column (bio/xtal/random) to all splits
  7. Adds Default column to test (False for ExtremeLong)

Inputs:
  - output/ppi_pairs_corrected_seq95.csv
  - output/structure_selection.csv
  - output/ppi_unique_pairs.csv
  - data/uniprot_sequences_all.fasta

Outputs:
  - splits/ppi_{train,val,test}.csv
"""

import json
import random
import shutil
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "output"
SPLITS_DIR = ROOT / "splits"
DATA_DIR = ROOT / "data"

CORRECTED_CSV = OUTPUT_DIR / "ppi_pairs_corrected_seq95.csv"
STRUCTURE_SELECTION = OUTPUT_DIR / "structure_selection.csv"
KNOWN_PAIRS_CSV = OUTPUT_DIR / "ppi_unique_pairs.csv"
FASTA_FILE = DATA_DIR / "uniprot_sequences_all.fasta"

MAX_LEN_GLOBAL = 2000      # hard filter before splitting
MAX_LEN_TRAIN_VAL = 1000  # train/val only
TRAIN_RATIO = 0.70
VAL_RATIO = 0.10


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_sequences(fasta_path):
    seqs = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        parts = record.id.split("|")
        uid = parts[1] if len(parts) >= 2 else record.id
        seqs[uid] = str(record.seq)
    return seqs


def split_selected_files(val):
    if pd.isna(val):
        return None, None
    parts = str(val).split(";")
    if len(parts) >= 2:
        return parts[0].strip(), parts[1].strip()
    if len(parts) == 1:
        return parts[0].strip(), None
    return None, None


def load_known_pairs(pairs_csv):
    df = pd.read_csv(pairs_csv)
    bio_xtal = df[df["label_str"].isin(["BIO", "XTAL"])]
    pairs = set()
    for _, row in bio_xtal.iterrows():
        pairs.add("__".join(sorted([str(row["protein_A_id"]), str(row["protein_B_id"])])))
    return pairs


def build_sequence_file_map(structure_selection, positive_df):
    """Map protein sequence -> {A: set(files), B: set(files)} from positive pool."""
    rep_to_seqs = {}
    for _, row in positive_df.iterrows():
        rep = row.get("representative_id")
        if pd.isna(rep) or str(rep).strip() == "":
            continue
        rep_to_seqs[rep] = (row["protein_A_seq"], row["protein_B_seq"])

    seq_to_file = {}
    matched = 0
    unmatched = 0
    for _, row in structure_selection.iterrows():
        rep = row["representative_id"]
        if rep not in rep_to_seqs:
            unmatched += 1
            continue
        matched += 1
        seq_a, seq_b = rep_to_seqs[rep]
        file_a, file_b = split_selected_files(row["selected_files"])
        for seq, role, filename in [(seq_a, "A", file_a), (seq_b, "B", file_b)]:
            if filename is None:
                continue
            seq_to_file.setdefault(seq, {"A": set(), "B": set()})
            seq_to_file[seq][role].add(filename)

    print(f"  seq->file map: {len(seq_to_file)} seqs, "
          f"{matched}/{matched+unmatched} reps matched")
    return seq_to_file


def choose_file(seq_to_file, seq, prefer_role):
    if seq not in seq_to_file:
        return None
    files = seq_to_file[seq]
    if files.get(prefer_role):
        return next(iter(files[prefer_role]))
    other_role = "B" if prefer_role == "A" else "A"
    if files.get(other_role):
        return next(iter(files[other_role]))
    return None


def generate_random_negatives(protein_pool, known_pairs, n_needed, seqs, seq_to_file, seed):
    """Generate random negative pairs, allowing self-pairs."""
    random.seed(seed)
    negatives = []
    attempts = 0
    max_attempts = n_needed * 1000
    while len(negatives) < n_needed and attempts < max_attempts:
        attempts += 1
        a = random.choice(protein_pool)
        b = random.choice(protein_pool)  # allows a == b
        pair_key = "__".join(sorted([a, b]))
        if pair_key in known_pairs:
            continue
        seq_a = seqs.get(a)
        seq_b = seqs.get(b)
        if not seq_a or not seq_b:
            continue
        file_a = choose_file(seq_to_file, seq_a, "A")
        file_b = choose_file(seq_to_file, seq_b, "B")
        if file_a is None or file_b is None:
            continue
        negatives.append({
            "protein_A_id": a,
            "protein_B_id": b,
            "aa_seq1": seq_a,
            "aa_seq2": seq_b,
            "struct_file1": file_a,
            "struct_file2": file_b,
        })
    return negatives


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=== Loading inputs ===")
    seq95 = pd.read_csv(CORRECTED_CSV)
    ss = pd.read_csv(STRUCTURE_SELECTION)
    known_pairs = load_known_pairs(KNOWN_PAIRS_CSV)
    seqs = load_sequences(FASTA_FILE)
    print(f"  pairs: {len(seq95)} (BIO={sum(seq95.label_str=='BIO')}, "
          f"XTAL={sum(seq95.label_str=='XTAL')})")
    print(f"  structure_selection: {len(ss)}")
    print(f"  known pairs: {len(known_pairs)}")
    print(f"  sequences: {len(seqs)}")

    # --- Add sequences to seq95 ---
    seq95["aa_seq1"] = seq95["protein_A_id"].map(seqs)
    seq95["aa_seq2"] = seq95["protein_B_id"].map(seqs)
    seq95 = seq95.dropna(subset=["aa_seq1", "aa_seq2"]).copy()
    seq95["len1"] = seq95["aa_seq1"].str.len()
    seq95["len2"] = seq95["aa_seq2"].str.len()
    seq95["len_max"] = seq95[["len1", "len2"]].max(axis=1)

    # --- 1b-i: Global hard filter >2000 ---
    before = len(seq95)
    seq95 = seq95[seq95["len_max"] <= MAX_LEN_GLOBAL].copy()
    print(f"\n=== Global filter >{MAX_LEN_GLOBAL} ===")
    print(f"  {before} -> {len(seq95)} (dropped {before - len(seq95)})")

    # --- 1a: Time split based on BIO sample count ---
    bio = seq95[seq95["label_str"] == "BIO"].copy()
    xtal = seq95[seq95["label_str"] == "XTAL"].copy()

    bio["date_parsed"] = pd.to_datetime(bio["release_date"], errors="coerce")
    bio = bio.dropna(subset=["date_parsed"]).sort_values("date_parsed").reset_index(drop=True)
    n_bio = len(bio)
    train_cutoff = bio["date_parsed"].iloc[int(n_bio * TRAIN_RATIO)]
    val_cutoff = bio["date_parsed"].iloc[int(n_bio * (TRAIN_RATIO + VAL_RATIO))]
    print(f"\n=== Time split (BIO samples {TRAIN_RATIO}:{VAL_RATIO}:{1-TRAIN_RATIO-VAL_RATIO:.2f}) ===")
    print(f"  train cutoff: {train_cutoff}")
    print(f"  val cutoff:   {val_cutoff}")

    def assign_split(df):
        d = df.copy()
        d["date_parsed"] = pd.to_datetime(d["release_date"], errors="coerce")
        d = d.dropna(subset=["date_parsed"])
        d["split"] = "test"
        d.loc[d["date_parsed"] <= train_cutoff, "split"] = "train"
        d.loc[(d["date_parsed"] > train_cutoff) & (d["date_parsed"] <= val_cutoff), "split"] = "val"
        return d

    bio_s = assign_split(bio)
    xtal_s = assign_split(xtal)
    print(f"  BIO:  train={sum(bio_s.split=='train')}, val={sum(bio_s.split=='val')}, test={sum(bio_s.split=='test')}")
    print(f"  XTAL: train={sum(xtal_s.split=='train')}, val={sum(xtal_s.split=='val')}, test={sum(xtal_s.split=='test')}")

    # --- Merge structure selection for file lookup ---
    ss_map = ss.set_index("representative_id")
    ss_map = ss_map[["selected_files"]].copy()
    ss_map[["struct_file1", "struct_file2"]] = ss_map["selected_files"].apply(
        lambda x: pd.Series(split_selected_files(x))
    )

    # --- Build seq->file map for random negatives ---
    # Use all BIO+XTAL pairs as the positive pool for sequence lookups
    all_positive = seq95[seq95["label_str"].isin(["BIO", "XTAL"])][
        ["representative_id", "protein_A_id", "protein_B_id", "aa_seq1", "aa_seq2"]
    ].copy()
    # Rename for build_sequence_file_map compatibility
    all_positive = all_positive.rename(columns={"aa_seq1": "protein_A_seq", "aa_seq2": "protein_B_seq"})
    seq_to_file = build_sequence_file_map(ss, all_positive)

    # --- Process each split ---
    SPLITS_DIR.mkdir(parents=True, exist_ok=True)
    neg_id_counter = 0
    summary = {}

    for split_name in ["train", "val", "test"]:
        print(f"\n=== Processing {split_name} ===")
        b = bio_s[bio_s["split"] == split_name].copy()
        x = xtal_s[xtal_s["split"] == split_name].copy()

        # --- 1b-ii: train/val filter >1000 ---
        if split_name in ("train", "val"):
            b = b[b["len_max"] <= MAX_LEN_TRAIN_VAL]
            x = x[x["len_max"] <= MAX_LEN_TRAIN_VAL]

        # --- test: filter both_chains_apo_or_predicted ---
        if split_name == "test":
            b = b.merge(
                ss[["representative_id", "both_chains_apo_or_predicted"]],
                on="representative_id", how="left"
            )
            b = b[b["both_chains_apo_or_predicted"] == True].drop(
                columns=["both_chains_apo_or_predicted"]
            )
            x = x.merge(
                ss[["representative_id", "both_chains_apo_or_predicted"]],
                on="representative_id", how="left"
            )
            x = x[x["both_chains_apo_or_predicted"] == True].drop(
                columns=["both_chains_apo_or_predicted"]
            )

        n_bio = len(b)
        n_xtal = len(x)
        n_random_needed = max(n_bio - n_xtal, 0)

        print(f"  BIO={n_bio}, XTAL={n_xtal}, random_needed={n_random_needed}")

        # --- Build positive rows (BIO, label=1) ---
        b_rows = pd.DataFrame({
            "unique_id": b["representative_id"].values,
            "protein_A_id": b["protein_A_id"].values,
            "protein_B_id": b["protein_B_id"].values,
            "aa_seq1": b["aa_seq1"].values,
            "aa_seq2": b["aa_seq2"].values,
            "struct_file1": b["representative_id"].map(ss_map["struct_file1"]).values,
            "struct_file2": b["representative_id"].map(ss_map["struct_file2"]).values,
            "label": 1,
            "pair_type": "bio",
        })

        # --- Build XTAL rows (label=0) ---
        x_rows = pd.DataFrame({
            "unique_id": x["representative_id"].values,
            "protein_A_id": x["protein_A_id"].values,
            "protein_B_id": x["protein_B_id"].values,
            "aa_seq1": x["aa_seq1"].values,
            "aa_seq2": x["aa_seq2"].values,
            "struct_file1": x["representative_id"].map(ss_map["struct_file1"]).values,
            "struct_file2": x["representative_id"].map(ss_map["struct_file2"]).values,
            "label": 0,
            "pair_type": "xtal",
        })

        # --- Generate random negatives ---
        # Pool = all proteins in this split (BIO + XTAL)
        pool = sorted(
            set(b["protein_A_id"]) | set(b["protein_B_id"]) |
            set(x["protein_A_id"]) | set(x["protein_B_id"])
        )
        pool = [p for p in pool if p in seqs]  # ensure sequence available

        seed = 42 + ["train", "val", "test"].index(split_name) * 100
        rand_negs = generate_random_negatives(
            pool, known_pairs, n_random_needed, seqs, seq_to_file, seed
        )

        r_rows = pd.DataFrame(rand_negs)
        if len(r_rows) > 0:
            r_rows["unique_id"] = [f"neg_{neg_id_counter + i}" for i in range(len(r_rows))]
            neg_id_counter += len(r_rows)
            r_rows["label"] = 0
            r_rows["pair_type"] = "random"
        else:
            r_rows = pd.DataFrame(columns=[
                "unique_id", "protein_A_id", "protein_B_id",
                "aa_seq1", "aa_seq2", "struct_file1", "struct_file2",
                "label", "pair_type"
            ])

        print(f"  random generated: {len(r_rows)}")

        # --- Combine ---
        cols = ["unique_id", "protein_A_id", "protein_B_id",
                "aa_seq1", "aa_seq2", "struct_file1", "struct_file2",
                "label", "pair_type"]
        combined = pd.concat([b_rows[cols], x_rows[cols], r_rows[cols]], ignore_index=True)

        # --- test: add Default column ---
        if split_name == "test":
            combined["len_max"] = combined[["aa_seq1", "aa_seq2"]].apply(
                lambda r: max(len(r["aa_seq1"]), len(r["aa_seq2"])), axis=1
            )
            combined["Default"] = combined["len_max"] <= MAX_LEN_TRAIN_VAL
            combined = combined.drop(columns=["len_max"])
            cols = cols[:8] + ["Default", "pair_type"]
        else:
            cols = cols[:8] + ["pair_type"]

        combined = combined[cols]
        out_path = SPLITS_DIR / f"ppi_{split_name}.csv"
        combined.to_csv(out_path, index=False)
        print(f"  Saved {out_path}: {len(combined)} rows "
              f"(BIO={n_bio}, XTAL={n_xtal}, random={len(r_rows)})")

        summary[split_name] = {
            "bio": n_bio,
            "xtal": n_xtal,
            "random": len(r_rows),
            "total": len(combined),
        }

    # --- Print summary ---
    print("\n=== Summary ===")
    grand = sum(s["total"] for s in summary.values())
    for sp in ["train", "val", "test"]:
        s = summary[sp]
        print(f"  {sp}: BIO={s['bio']}, XTAL={s['xtal']}, random={s['random']}, "
              f"total={s['total']} ({s['total']/grand*100:.1f}%)")

    summary_path = OUTPUT_DIR / "xtal_rebuild_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nSaved {summary_path}")


if __name__ == "__main__":
    main()
