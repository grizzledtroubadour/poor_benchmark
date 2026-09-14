#!/usr/bin/env python3
"""Move the 3 train rows with aa_seq < 60 into test as OOD_ExtremeShort.

These rows (P05959 55aa, B1KRG2 59aa, P0C7A9 46aa) slipped into train because
the extreme-length pre-split in 02_split_dataset.py judged by the metadata
`seq_length` (UniProt domain length) while `aa_seq` is the structure-extracted
sequence — the two calibers disagree on these 3 samples.

New test rows follow the pH extreme-length convention: Default=False,
OOD_ExtremeShort=True, OOD_ExtremeLong=False, all other OOD columns NaN,
InD=False. Structure files already exist in pdbs/ and are untouched.

Usage: python3 07_move_short_train_rows_to_test.py [--apply]
Without --apply: dry-run report only.
"""
import os, sys
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
TD = os.path.join(ROOT, "datasets", "enzyme_optimal_ph")
TRAIN = os.path.join(TD, "splits", "optimal_ph_prediction_train.csv")
TEST = os.path.join(TD, "splits", "optimal_ph_prediction_test.csv")
MIN_LEN = 60

def main():
    apply = "--apply" in sys.argv
    train = pd.read_csv(TRAIN, dtype=str, keep_default_na=False)
    test = pd.read_csv(TEST, dtype=str, keep_default_na=False)

    short_mask = train["aa_seq"].str.len() < MIN_LEN
    movers = train[short_mask]
    print(f"train rows with aa_seq < {MIN_LEN}: {len(movers)}")
    for _, r in movers.iterrows():
        print(f"  {r['unique_id']}  len={len(r['aa_seq'])}  label={r['label']}")
    if len(movers) == 0:
        print("nothing to do")
        return
    if not apply:
        print("dry-run; pass --apply to write")
        return

    train_new = train[~short_mask]
    core = ["unique_id", "aa_seq", "struct_file", "label"]
    rows = []
    for _, r in movers.iterrows():
        row = {c: "" for c in test.columns}
        for c in core:
            row[c] = r[c]
        row["Default"] = "False"
        row["OOD_ExtremeShort"] = "True"
        row["OOD_ExtremeLong"] = "False"
        row["InD"] = "False"
        rows.append(row)
    test_new = pd.concat([test, pd.DataFrame(rows, columns=test.columns)], ignore_index=True)

    train_new.to_csv(TRAIN, index=False, lineterminator="\n")
    test_new.to_csv(TEST, index=False, lineterminator="\n")
    print(f"train {len(train)} -> {len(train_new)}; test {len(test)} -> {len(test_new)}")

if __name__ == "__main__":
    main()
