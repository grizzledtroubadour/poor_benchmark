#!/usr/bin/env python3
"""
Create length-based OOD split for kcat prediction with UniProt-level non-overlap.

Pipeline:
1. Separate all extreme-length samples (< NORMAL_MIN or > NORMAL_MAX).
2. Random 7:1:2 split of normal-length UniProts into train/val/test (no UniProt overlap).
3. Compute reactant ECFP4 max Tanimoto similarity of normal val/test rows to train rows.
   Move entire UniProts with max similarity < SIMILARITY_THRESHOLD to train.
4. Compute reactant ECFP4 max Tanimoto similarity of extreme-length rows to the expanded train.
   Discard entire UniProts with max similarity < SIMILARITY_THRESHOLD.
   Put the remaining extreme-length rows into test as OOD samples.

Input: enzyme_kinetics_prediction/data/kcat_with_reactant_ecfp4.csv
Output: enzyme_kinetics_prediction/splits/kcat_{train,val,test}.csv
        enzyme_kinetics_prediction/splits/kcat_split_summary.json
"""
import json
import pandas as pd
import numpy as np
from pathlib import Path
from rdkit import RDLogger

RDLogger.DisableLog('rdApp.*')

TASK_DIR = Path(__file__).resolve().parent.parent  # datasets/enzyme_kinetics_prediction
DATA_DIR = TASK_DIR / "data"
SPLITS_DIR = TASK_DIR / "splits"
OUT_DIR = TASK_DIR / "output"
SPLITS_DIR.mkdir(parents=True, exist_ok=True)
OUT_DIR.mkdir(parents=True, exist_ok=True)

FP_NBITS = 2048
SIMILARITY_THRESHOLD = 0.6
RANDOM_SEED = 42
TRAIN_RATIO = 0.7
VAL_RATIO = 0.1
NORMAL_MIN = 60
NORMAL_MAX = 1000


def drop_unnamed_index_cols(df):
    """Drop pandas artifact columns like Unnamed: 0, Unnamed: 0.1."""
    unnamed = [c for c in df.columns if str(c).startswith('Unnamed')]
    if unnamed:
        df = df.drop(columns=unnamed)
    return df


def compute_max_similarity(query_fp, train_fp):
    """Compute max Tanimoto similarity for each query against the train set."""
    q = query_fp.astype(np.float32)
    t = train_fp.astype(np.float32)
    intersection = q @ t.T
    q_sum = q.sum(axis=1, keepdims=True)
    t_sum = t.sum(axis=1, keepdims=True)
    union = q_sum + t_sum.T - intersection
    sim = np.divide(intersection, union, out=np.zeros_like(intersection), where=union != 0)
    return sim.max(axis=1)


def compute_uniprot_max_similarity(df, train_fp, train_df):
    """Compute per-UniProt max similarity of df rows to train_fp."""
    if len(df) == 0 or len(train_df) == 0:
        return {}
    fp_cols = [f'ecfp4_bit_{i}' for i in range(FP_NBITS)]
    query_fp = df[fp_cols].values.astype(np.uint8)
    sims = compute_max_similarity(query_fp, train_fp)
    df = df.copy()
    df['__sim__'] = sims
    return df.groupby('uniprot')['__sim__'].max().to_dict()


def main():
    print("=" * 80)
    print("Loading full dataset with ECFP4 fingerprints")
    print("=" * 80)

    df = pd.read_csv(DATA_DIR / 'kcat_with_reactant_ecfp4.csv', dtype={'ecfp4': str})
    df = drop_unnamed_index_cols(df)
    print(f"Total rows: {len(df):,}")
    print(f"Unique UniProts: {df['uniprot'].nunique():,}")

    # Compute sequence lengths
    seq_len = df['sequence'].str.len()
    is_extreme = (seq_len < NORMAL_MIN) | (seq_len > NORMAL_MAX)

    print(f"\nLength distribution:")
    print(f"  Normal ({NORMAL_MIN}-{NORMAL_MAX}): {(~is_extreme).sum():,} rows, {df.loc[~is_extreme, 'uniprot'].nunique():,} UniProts")
    print(f"  Extreme (<{NORMAL_MIN} or >{NORMAL_MAX}): {is_extreme.sum():,} rows, {df.loc[is_extreme, 'uniprot'].nunique():,} UniProts")
    print(f"    - Short (<{NORMAL_MIN}): {(seq_len < NORMAL_MIN).sum():,} rows")
    print(f"    - Long (>{NORMAL_MAX}): {(seq_len > NORMAL_MAX).sum():,} rows")

    # Step 1: Separate extreme and normal samples
    print("\n" + "=" * 80)
    print("Step 1: Separate extreme-length samples")
    print("=" * 80)
    normal_df = df[~is_extreme].copy()
    extreme_df = df[is_extreme].copy()

    # Step 2: Random 7:1:2 split of normal-length UniProts
    print("\n" + "=" * 80)
    print("Step 2: Random 7:1:2 split of normal-length UniProts")
    print("=" * 80)
    normal_uniprots = normal_df['uniprot'].unique()
    np.random.seed(RANDOM_SEED)
    shuffled_uniprots = np.random.permutation(normal_uniprots)

    n_total = len(shuffled_uniprots)
    n_train = int(n_total * TRAIN_RATIO)
    n_val = int(n_total * VAL_RATIO)

    train_uniprots = set(shuffled_uniprots[:n_train])
    val_uniprots = set(shuffled_uniprots[n_train:n_train + n_val])
    test_uniprots = set(shuffled_uniprots[n_train + n_val:])

    print(f"  Initial normal UniProts: train={len(train_uniprots):,}, val={len(val_uniprots):,}, test={len(test_uniprots):,}")

    train_df = normal_df[normal_df['uniprot'].isin(train_uniprots)].copy()
    val_df = normal_df[normal_df['uniprot'].isin(val_uniprots)].copy()
    test_normal_df = normal_df[normal_df['uniprot'].isin(test_uniprots)].copy()

    print(f"  Initial normal rows: train={len(train_df):,}, val={len(val_df):,}, test={len(test_normal_df):,}")

    # Step 3: Compute similarity of normal val/test rows to train, move low-sim UniProts to train
    print("\n" + "=" * 80)
    print(f"Step 3: Filter normal val/test by reactant similarity to train (threshold={SIMILARITY_THRESHOLD})")
    print("=" * 80)

    fp_cols = [f'ecfp4_bit_{i}' for i in range(FP_NBITS)]
    train_fp = train_df[fp_cols].values.astype(np.uint8)

    val_uniprot_sim = compute_uniprot_max_similarity(val_df, train_fp, train_df)
    test_normal_uniprot_sim = compute_uniprot_max_similarity(test_normal_df, train_fp, train_df)

    val_move_uniprots = {up for up, sim in val_uniprot_sim.items() if sim < SIMILARITY_THRESHOLD}
    test_move_uniprots = {up for up, sim in test_normal_uniprot_sim.items() if sim < SIMILARITY_THRESHOLD}

    print(f"  Val UniProts to move to train (sim < {SIMILARITY_THRESHOLD}): {len(val_move_uniprots):,}")
    print(f"  Test normal UniProts to move to train (sim < {SIMILARITY_THRESHOLD}): {len(test_move_uniprots):,}")

    # Move entire UniProts to train
    train_df = pd.concat([
        train_df,
        val_df[val_df['uniprot'].isin(val_move_uniprots)],
        test_normal_df[test_normal_df['uniprot'].isin(test_move_uniprots)]
    ], ignore_index=True)
    val_df = val_df[~val_df['uniprot'].isin(val_move_uniprots)].reset_index(drop=True)
    test_normal_df = test_normal_df[~test_normal_df['uniprot'].isin(test_move_uniprots)].reset_index(drop=True)

    print(f"\nAfter moving low-sim UniProts:")
    print(f"  Train rows: {len(train_df):,}, UniProts: {train_df['uniprot'].nunique():,}")
    print(f"  Val rows:   {len(val_df):,}, UniProts: {val_df['uniprot'].nunique():,}")
    print(f"  Test normal rows: {len(test_normal_df):,}, UniProts: {test_normal_df['uniprot'].nunique():,}")

    # Step 4: Compute similarity of extreme-length rows to expanded train, remove low-sim UniProts
    print("\n" + "=" * 80)
    print(f"Step 4: Filter extreme-length samples by similarity to expanded train")
    print("=" * 80)

    train_fp = train_df[fp_cols].values.astype(np.uint8)
    extreme_uniprot_sim = compute_uniprot_max_similarity(extreme_df, train_fp, train_df)

    extreme_keep_uniprots = {up for up, sim in extreme_uniprot_sim.items() if sim >= SIMILARITY_THRESHOLD}
    extreme_remove_uniprots = {up for up, sim in extreme_uniprot_sim.items() if sim < SIMILARITY_THRESHOLD}

    extreme_test_df = extreme_df[extreme_df['uniprot'].isin(extreme_keep_uniprots)].copy()
    extreme_removed_df = extreme_df[extreme_df['uniprot'].isin(extreme_remove_uniprots)].copy()

    print(f"  Extreme UniProts kept for test (sim >= {SIMILARITY_THRESHOLD}): {len(extreme_keep_uniprots):,}")
    print(f"  Extreme UniProts removed (sim < {SIMILARITY_THRESHOLD}): {len(extreme_remove_uniprots):,}")

    # Final test = normal test + extreme test
    test_df = pd.concat([test_normal_df, extreme_test_df], ignore_index=True)

    # Save splits
    print("\n" + "=" * 80)
    print("Saving final splits")
    print("=" * 80)
    train_df.to_csv(SPLITS_DIR / 'kcat_train.csv', index=False)
    val_df.to_csv(SPLITS_DIR / 'kcat_val.csv', index=False)
    test_df.to_csv(SPLITS_DIR / 'kcat_test.csv', index=False)

    print(f"  Train: {SPLITS_DIR / 'kcat_train.csv'} ({len(train_df):,} rows, {train_df['uniprot'].nunique():,} UniProts)")
    print(f"  Val:   {SPLITS_DIR / 'kcat_val.csv'} ({len(val_df):,} rows, {val_df['uniprot'].nunique():,} UniProts)")
    print(f"  Test:  {SPLITS_DIR / 'kcat_test.csv'} ({len(test_df):,} rows, {test_df['uniprot'].nunique():,} UniProts)")

    # Non-overlap check
    train_ups = set(train_df['uniprot'].unique())
    val_ups = set(val_df['uniprot'].unique())
    test_ups = set(test_df['uniprot'].unique())
    print(f"\nUniProt overlap check:")
    print(f"  Train ∩ Val:   {len(train_ups & val_ups)}")
    print(f"  Train ∩ Test:  {len(train_ups & test_ups)}")
    print(f"  Val ∩ Test:    {len(val_ups & test_ups)}")

    # Summary
    train_lengths = train_df['sequence'].str.len()
    val_lengths = val_df['sequence'].str.len()
    test_lengths = test_df['sequence'].str.len()

    summary = {
        'total_rows': int(len(df)),
        'normal_length_rows': int((~is_extreme).sum()),
        'extreme_length_rows': int(is_extreme.sum()),
        'short_rows': int((seq_len < NORMAL_MIN).sum()),
        'long_rows': int((seq_len > NORMAL_MAX).sum()),
        'train_rows': int(len(train_df)),
        'val_rows': int(len(val_df)),
        'test_rows': int(len(test_df)),
        'train_uniprots': int(len(train_ups)),
        'val_uniprots': int(len(val_ups)),
        'test_uniprots': int(len(test_ups)),
        'train_normal_rows': int(((train_lengths >= NORMAL_MIN) & (train_lengths <= NORMAL_MAX)).sum()),
        'train_extreme_rows': int(((train_lengths < NORMAL_MIN) | (train_lengths > NORMAL_MAX)).sum()),
        'val_normal_rows': int(((val_lengths >= NORMAL_MIN) & (val_lengths <= NORMAL_MAX)).sum()),
        'val_extreme_rows': int(((val_lengths < NORMAL_MIN) | (val_lengths > NORMAL_MAX)).sum()),
        'test_normal_rows': int(((test_lengths >= NORMAL_MIN) & (test_lengths <= NORMAL_MAX)).sum()),
        'test_extreme_rows': int(((test_lengths < NORMAL_MIN) | (test_lengths > NORMAL_MAX)).sum()),
        'test_short_rows': int((test_lengths < NORMAL_MIN).sum()),
        'test_long_rows': int((test_lengths > NORMAL_MAX).sum()),
        'normal_val_uniprots_moved_to_train': int(len(val_move_uniprots)),
        'normal_test_uniprots_moved_to_train': int(len(test_move_uniprots)),
        'extreme_uniprots_removed': int(len(extreme_remove_uniprots)),
        'extreme_uniprots_kept_for_test': int(len(extreme_keep_uniprots)),
        'similarity_threshold': SIMILARITY_THRESHOLD,
        'train_ratio': TRAIN_RATIO,
        'val_ratio': VAL_RATIO,
        'random_seed': RANDOM_SEED,
    }

    with open(OUT_DIR / 'kcat_split_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"\nSaved summary to {OUT_DIR / 'kcat_split_summary.json'}")

    print("\n" + "=" * 80)
    print("Done!")
    print("=" * 80)


if __name__ == '__main__':
    main()
