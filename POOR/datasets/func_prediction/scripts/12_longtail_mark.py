#!/usr/bin/env python3
"""Flag predicted-label long tails: True if any test label has training frequency <=10.

Definition:
1. Count labels in splits/{task}_train.csv (semicolon-separated label column).
2. Select labels with frequency <=10 as the long-tail set.
3. A test label intersecting that set (any semantics) sets OOD_LongTail=True.
4. Default=False extreme-length rows are excluded and left empty (NaN: not computed),
   following project convention; finalize_ood_columns.py --mode nan is a fallback check.

I/O: read/write splits/{task}_train.csv and splits/{task}_test.csv directly.
Historical deviation: the original pipeline ran before extreme-length chains were merged,
reading output/datasets_split_raw/{task}_train.csv and {task}_test_ood.csv (not retained).
The current version runs on splits/ and flags only Default=True rows; values match historical results (DATA_PROCESS.md 4.8).
"""

import pandas as pd
from collections import Counter
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent  # datasets/func_prediction/
SPLITS_DIR = BASE / 'splits'
FIXED_THRESHOLD = 10
TASKS = ['ec', 'go_mf', 'go_cc', 'go_bp']


def split_labels(labels_str):
    return {l.strip() for l in str(labels_str).split(';') if l.strip()}


def get_longtail_labels(train_path, threshold=FIXED_THRESHOLD):
    """Get labels with training frequency <=threshold."""
    train_df = pd.read_csv(train_path)

    label_counts = Counter()
    for labels_str in train_df['label']:
        label_counts.update(split_labels(labels_str))

    longtail = {l for l, c in label_counts.items() if c <= threshold}
    return longtail, len(label_counts)


def mark_longtail(task):
    print(f"\n{'='*60}")
    print(f"{task.upper()} OOD_LongTail Marking (threshold ≤ {FIXED_THRESHOLD})")
    print(f"{'='*60}")

    train_path = SPLITS_DIR / f'{task}_train.csv'
    test_path = SPLITS_DIR / f'{task}_test.csv'

    longtail_labels, total_labels = get_longtail_labels(train_path)
    print(f"  Total train labels: {total_labels}")
    print(f"  Long-tail labels (≤ {FIXED_THRESHOLD}): {len(longtail_labels)} "
          f"({len(longtail_labels)/total_labels*100:.1f}%)")

    # Load test and mark only Default=True rows.
    test_df = pd.read_csv(test_path, dtype=str, keep_default_na=False)
    mask_default = test_df['Default'] == 'True'

    new_col = pd.Series('', index=test_df.index, dtype=object)
    new_col[mask_default] = test_df.loc[mask_default, 'label'].apply(
        lambda s: 'True' if split_labels(s) & longtail_labels else 'False')
    test_df['OOD_LongTail'] = new_col
    test_df.to_csv(test_path, index=False, lineterminator='\n')

    nf = (new_col == 'True').sum()
    n_default = mask_default.sum()
    print(f"\n  Test samples with long-tail labels: {nf}/{n_default} ({nf/n_default*100:.2f}%)")
    print(f"  Saved to: {test_path}")

    return nf, n_default


if __name__ == '__main__':
    for task in TASKS:
        mark_longtail(task)
