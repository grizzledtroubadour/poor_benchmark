#!/usr/bin/env python3
"""Stage 3.1: SSP splits (ssp_train/val/test.csv).

Rebuild splits from Stage 1/2 outputs in one pass:

1. Input output/ssp_full_data.csv (representative unique_id/aa_seq/labels);
            chain lengths come from output/stage4_cluster_representatives.csv and deposition dates from
   ``output/pdb_metadata_api.csv``。
2. Remove single chains >2000 aa (only 4O9X_A, 2,082 residues).
3. Exclude low-quality chains in QC_EXCLUDE (7UYL_K: 59% UNK). Sequence extraction skips
            UNK while structures retain it, misaligning labels; see DATA_PROCESS.md notes.
4. Put 82 chains >1000 aa directly into test, outside temporal splitting,
            with Default=False / OOD_ExtremeLong=True. Append them by unique_id at the
            end of test; later scripts/skills annotate other OOD columns (Default=False -> NaN).
5. Stably sort remaining chains by deposition_date and split by 72:8:20 quantiles:
            n_train=int(n*0.72), n_val=round(n*0.08), and the remainder in test.
            Temporal-split samples have Default=True / OOD_ExtremeLong=False.
            Cutoffs are dynamic quantiles of the current pool, not fixed dates. Published cutoffs
            are 2019-06-06 / 2021-02-09, historically computed on the >=60-aa pool before
            adding short chains at the same boundaries (see DATA_PROCESS.md history). Future
            reproductions may differ as pool composition changes.

Output: train/val columns unique_id,aa_seq,struct_file,label; test additionally includes
``Default,OOD_ExtremeLong``。``struct_file = unique_id + ".pdb"``。

Usage:
    python scripts/08_create_ssp_datasets.py                    # Write to splits/.
    python scripts/08_create_ssp_datasets.py --out-dir output/verify
"""
import argparse
from pathlib import Path

import pandas as pd

TASK_DIR = Path(__file__).resolve().parent.parent
OUTPUT_DIR = TASK_DIR / "output"

FULL_DATA = OUTPUT_DIR / "ssp_full_data.csv"
REPRESENTATIVES = OUTPUT_DIR / "stage4_cluster_representatives.csv"
METADATA = OUTPUT_DIR / "pdb_metadata_api.csv"

HARD_MAX_LEN = 2000   # Remove single chains >2000 aa.
EXTREME_LONG = 1000   # Send >1000-aa chains directly to test.

# 72:8:20 temporal split, computing cutoffs dynamically.
TRAIN_FRAC, VAL_FRAC = 0.72, 0.08

# QC exclusion: 71 of 120 structural residues in 7UYL_K are UNK (59%). AA_MAP-based
# sequence extraction skips UNK (aa_seq has only 49 residues), while structures retain it,
# breaking label/struct_label alignment.
QC_EXCLUDE = {"7UYL_K"}

BASE_COLS = ["unique_id", "aa_seq", "struct_file", "label"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", type=Path, default=TASK_DIR / "splits",
                    help="Output directory (default splits/; use output/verify for a trial)")
    args = ap.parse_args()

    full = pd.read_csv(FULL_DATA, dtype=str, keep_default_na=False)
    rep = pd.read_csv(REPRESENTATIVES, dtype=str, keep_default_na=False)
    meta = pd.read_csv(METADATA, dtype=str, keep_default_na=False)

    rep["num_residues"] = rep["num_residues"].astype(int)
    df = full.merge(rep[["chain_key", "pdb_id", "num_residues"]],
                    left_on="unique_id", right_on="chain_key", how="left")
    df = df.merge(meta[["pdb_id", "deposition_date"]], on="pdb_id", how="left")
    assert df["num_residues"].notna().all() and (df["deposition_date"] != "").all()

    # Remove single chains >2000 aa (4O9X_A, 2,082 residues).
    dropped = df[df["num_residues"] > HARD_MAX_LEN]
    df = df[df["num_residues"] <= HARD_MAX_LEN].copy()
    print(f"dropped len>{HARD_MAX_LEN}: {len(dropped)} "
          f"({', '.join(dropped['unique_id'])})")

    # Exclude low-quality chains.
    qc = df[df["unique_id"].isin(QC_EXCLUDE)]
    df = df[~df["unique_id"].isin(QC_EXCLUDE)].copy()
    print(f"dropped QC_EXCLUDE: {len(qc)} ({', '.join(qc['unique_id'])})")

    extreme = df[df["num_residues"] > EXTREME_LONG].sort_values("unique_id")
    rest = df[df["num_residues"] <= EXTREME_LONG].sort_values(
        "deposition_date", kind="stable")
    print(f"extreme-long (>{EXTREME_LONG}): {len(extreme)}, time-split pool: {len(rest)}")

    # 72:8:20 temporal split by counts after stable deposition_date sorting.
    n = len(rest)
    n_train = int(n * TRAIN_FRAC)
    n_val = round(n * VAL_FRAC)
    train = rest.iloc[:n_train]
    val = rest.iloc[n_train:n_train + n_val]
    test_time = rest.iloc[n_train + n_val:]

    def base(frame):
        out = frame[["unique_id", "aa_seq"]].copy()
        out["struct_file"] = out["unique_id"] + ".pdb"
        out["label"] = frame["labels"]
        return out[BASE_COLS]

    test = pd.concat([base(test_time), base(extreme)])
    test["Default"] = ["True"] * len(test_time) + ["False"] * len(extreme)
    test["OOD_ExtremeLong"] = ["False"] * len(test_time) + ["True"] * len(extreme)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    for name, frame in [("ssp_train", base(train)), ("ssp_val", base(val)),
                        ("ssp_test", test)]:
        path = args.out_dir / f"{name}.csv"
        frame.to_csv(path, index=False, lineterminator="\n")
        print(f"  {path}  rows={len(frame)}")
    print(f"cutoffs (computed): train <= {train['deposition_date'].max()}, "
          f"val <= {val['deposition_date'].max()}")
    print(f"date range: train {train['deposition_date'].min()} ~ {train['deposition_date'].max()}, "
          f"val {val['deposition_date'].min()} ~ {val['deposition_date'].max()}, "
          f"test {test_time['deposition_date'].min()} ~ {test_time['deposition_date'].max()}")


if __name__ == "__main__":
    main()
