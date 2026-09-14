#!/usr/bin/env python3
"""
Step 2: Integrate extreme-length chains into test CSVs.

For each task (EC, GO-MF, GO-CC, GO-BP):
1. Load extreme-length chains from Step 1 output
2. Clean labels (EC L4->L3 merge, GO frequency filter)
3. Filter: labels must appear in train set
4. Filter GO test: reviewed or has experimental GO evidence
5. Add to test CSV with OOD_ExtremeShort / OOD_ExtremeLong columns
6. Set Default=False for extreme-length samples
7. Other OOD columns set to False

Output: Updated splits/{task}_test.csv
"""

import json
import shutil
from collections import Counter
from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parent.parent  # func_prediction/
POOR_ROOT = BASE.parent.parent  # POOR/
SPLITS_DIR = BASE / 'splits'
OUTPUT_DIR = BASE / 'output'
BACKUP_DIR = OUTPUT_DIR / 'splits_pre_extreme_backup'

# EC cleaning helpers (from clean_datasets.py)
def parse_ec_level(ec_number):
    parts = ec_number.split('.')
    return sum(1 for p in parts if p != '-')

def get_ec_prefix(ec_number, level=3):
    parts = ec_number.split('.')
    non_dash = [p for p in parts if p != '-']
    return '.'.join(non_dash[:level])

def normalize_ec(ec_number):
    parts = ec_number.split('.')
    while len(parts) < 4:
        parts.append('-')
    return '.'.join(parts)

def clean_ec_labels(ec_str, train_ec_set):
    """Clean EC labels: keep L3+, merge L4->L3, filter to train set."""
    if not ec_str:
        return ''
    ecs = [normalize_ec(e.strip()) for e in ec_str.split(';') if e.strip()]
    
    # Step 1: Filter to L3+ only
    filtered = [ec for ec in ecs if parse_ec_level(ec) >= 3]
    
    # Step 2: Merge L4 -> L3
    merged = []
    seen = set()
    for ec in filtered:
        level = parse_ec_level(ec)
        if level == 4:
            prefix = get_ec_prefix(ec, 3)
            new_ec = f"{prefix}.-"
        else:
            new_ec = ec
        if new_ec not in seen:
            seen.add(new_ec)
            merged.append(new_ec)
    
    # Step 3: Filter to train set
    result = [ec for ec in merged if ec in train_ec_set]
    return ';'.join(result) if result else ''

def clean_go_labels(go_str, valid_go_set):
    """Clean GO labels: filter to valid set (frequency 5-5000 in train)."""
    if not go_str:
        return ''
    gos = [g.strip() for g in go_str.split(';') if g.strip()]
    result = [g for g in gos if g in valid_go_set]
    return ';'.join(result) if result else ''

def load_train_labels(splits_dir, task):
    """Load all labels from train CSV."""
    df = pd.read_csv(splits_dir / f'{task}_train.csv')
    labels = set()
    for l in df['label']:
        for item in str(l).split(';'):
            item = item.strip()
            if item:
                labels.add(item)
    return labels

def load_train_go_freq(splits_dir, task):
    """Count GO term frequency in train set."""
    df = pd.read_csv(splits_dir / f'{task}_train.csv')
    counts = Counter()
    for labels in df['label']:
        for item in str(labels).split(';'):
            item = item.strip()
            if item:
                counts[item] += 1
    return counts

# OOD columns that exist in current test CSVs
EC_TEST_OOD_COLS = [
    'seq_Redundancy_90', 'seq_Redundancy_80', 'seq_Redundancy_70',
    'seq_Redundancy_60', 'seq_Redundancy_50', 'seq_Redundancy_40',
    'seq_Redundancy_30',
    'TM-score_0.9', 'TM-score_0.8', 'TM-score_0.7',
    'TM-score_0.6', 'TM-score_0.5', 'TM-score_0.4', 'TM-score_0.3',
    'OOD_IDR', 'OOD_Orphan', 'OOD_NewEC_L4', 'OOD_LongTail', 'OOD_Combinatorial',
]
GO_TEST_OOD_COLS = [
    'seq_Redundancy_90', 'seq_Redundancy_80', 'seq_Redundancy_70',
    'seq_Redundancy_60', 'seq_Redundancy_50', 'seq_Redundancy_40',
    'seq_Redundancy_30',
    'TM-score_0.9', 'TM-score_0.8', 'TM-score_0.7',
    'TM-score_0.6', 'TM-score_0.5', 'TM-score_0.4', 'TM-score_0.3',
    'OOD_IDR', 'OOD_Orphan', 'OOD_LongTail', 'OOD_Combinatorial',
]

EXP_EVIDENCE = {'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'TAS', 'IC', 'HTP', 'HDA'}


def main():
    print("=" * 60)
    print("Step 2: Integrate extreme-length chains into test CSVs")
    print("=" * 60)

    # Load extreme-length chains
    ext_df = pd.read_csv(OUTPUT_DIR / 'extreme_length_chains.csv')
    print(f"Extreme-length chains: {len(ext_df)}")
    print(f"  ExtremeShort: {len(ext_df[ext_df['extreme_type']=='ExtremeShort'])}")
    print(f"  ExtremeLong:  {len(ext_df[ext_df['extreme_type']=='ExtremeLong'])}")
    print()

    # Backup current splits
    BACKUP_DIR.mkdir(exist_ok=True)
    for f in SPLITS_DIR.glob('*_test.csv'):
        shutil.copy2(f, BACKUP_DIR / f.name)
    print(f"Backed up test CSVs to {BACKUP_DIR}")
    print()

    # Process EC task
    print("=" * 40)
    print("EC task")
    print("=" * 40)
    integrate_task('ec', ext_df, is_ec=True)

    # Process GO tasks
    for task, label_col in [('go_mf', 'go_mf_labels'), ('go_cc', 'go_cc_labels'), ('go_bp', 'go_bp_labels')]:
        print()
        print("=" * 40)
        print(f"{task.upper()} task")
        print("=" * 40)
        integrate_task(task, ext_df, is_ec=False, go_label_col=label_col)


def integrate_task(task, ext_df, is_ec=False, go_label_col=None):
    """Integrate extreme-length chains for a specific task."""
    # Load existing test CSV
    test_path = SPLITS_DIR / f'{task}_test.csv'
    test_df = pd.read_csv(test_path)
    print(f"Current test: {len(test_df)} samples")

    # Determine OOD columns for this task
    ood_cols = EC_TEST_OOD_COLS if is_ec else GO_TEST_OOD_COLS

    # Get train labels
    train_labels = load_train_labels(SPLITS_DIR, task)
    print(f"Train labels: {len(train_labels)}")

    # For GO: get valid GO terms (frequency 5-5000 in train)
    valid_go_set = None
    if not is_ec:
        go_freq = load_train_go_freq(SPLITS_DIR, task)
        valid_go_set = {g for g, c in go_freq.items() if 5 <= c <= 5000}
        print(f"Valid GO terms (freq 5-5000): {len(valid_go_set)}")

    # Filter extreme chains for this task
    if is_ec:
        # EC task: chains with ec_labels
        task_ext = ext_df[ext_df['ec_labels'].notna() & (ext_df['ec_labels'] != '')].copy()
        print(f"Extreme chains with EC labels: {len(task_ext)}")

        # Clean EC labels
        task_ext['cleaned_labels'] = task_ext['ec_labels'].apply(
            lambda x: clean_ec_labels(x, train_labels)
        )
    else:
        # GO task: chains with corresponding GO labels
        task_ext = ext_df[ext_df[go_label_col].notna() & (ext_df[go_label_col] != '')].copy()
        print(f"Extreme chains with {go_label_col}: {len(task_ext)}")

        # Clean GO labels
        task_ext['cleaned_labels'] = task_ext[go_label_col].apply(
            lambda x: clean_go_labels(x, valid_go_set)
        )

    # Filter: must have at least one valid label
    task_ext = task_ext[task_ext['cleaned_labels'] != ''].copy()
    print(f"After label cleaning (must be in train): {len(task_ext)}")

    if len(task_ext) == 0:
        print("No extreme chains to add for this task, skipping")
        # Still add OOD_ExtremeShort / OOD_ExtremeLong columns (all False)
        test_df['OOD_ExtremeShort'] = False
        test_df['OOD_ExtremeLong'] = False
        test_df.to_csv(test_path, index=False)
        print(f"Updated {test_path} (no new samples, added columns)")
        return

    # GO test filtering: reviewed or has experimental GO evidence
    if not is_ec:
        before = len(task_ext)
        task_ext = task_ext[
            (task_ext['is_reviewed'] == True) | (task_ext['has_exp_go'] == True)
        ].copy()
        print(f"After GO test filter (reviewed or exp GO): {len(task_ext)} (removed {before - len(task_ext)})")

    # Deduplicate by chain_key (in case same chain appears for multiple tasks)
    task_ext = task_ext.drop_duplicates(subset='chain_key', keep='first')
    print(f"After dedup: {len(task_ext)}")

    # Also check chain_key not already in test
    existing_keys = set(test_df['unique_id'])
    task_ext = task_ext[~task_ext['chain_key'].isin(existing_keys)].copy()
    print(f"After removing existing test keys: {len(task_ext)}")

    # Build new rows
    new_rows = []
    for _, row in task_ext.iterrows():
        is_short = row['extreme_type'] == 'ExtremeShort'
        new_row = {
            'unique_id': row['chain_key'],
            'aa_seq': row['aa_seq'],
            'struct_file': f"{row['chain_key']}.pdb" if not row.get('is_cif', False) else f"{row['chain_key']}.cif",
            'label': row['cleaned_labels'],
            # All existing OOD columns set to False
            **{col: False for col in ood_cols},
            'Default': False,
            'OOD_ExtremeShort': is_short,
            'OOD_ExtremeLong': not is_short,
        }
        new_rows.append(new_row)

    # Add OOD_ExtremeShort / OOD_ExtremeLong to existing test rows (all False)
    test_df['OOD_ExtremeShort'] = False
    test_df['OOD_ExtremeLong'] = False

    # Append new rows
    new_df = pd.DataFrame(new_rows)
    # Ensure column order matches
    all_cols = list(test_df.columns)
    for col in all_cols:
        if col not in new_df.columns:
            new_df[col] = False
    new_df = new_df[all_cols]
    combined = pd.concat([test_df, new_df], ignore_index=True)

    # Save
    combined.to_csv(test_path, index=False)
    print(f"Updated {test_path}: {len(test_df)} -> {len(combined)} samples")
    print(f"  Added ExtremeShort: {sum(1 for r in new_rows if r['OOD_ExtremeShort'])}")
    print(f"  Added ExtremeLong:  {sum(1 for r in new_rows if r['OOD_ExtremeLong'])}")
    print(f"  Default=True: {combined['Default'].sum()}")
    print(f"  Default=False: {(~combined['Default']).sum()}")


if __name__ == '__main__':
    main()
