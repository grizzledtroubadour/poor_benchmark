#!/usr/bin/env python3
"""
Dataset cleaning script for EC and GO annotations.

EC cleaning rules:
1. Keep only level 3 and 4 EC numbers (remove level 1 and 2)
2. If a level-3 prefix (e.g., 1.1.1) exists in the dataset (as 1.1.1.-), 
   merge all level-4 numbers under that prefix (e.g., 1.1.1.1 -> 1.1.1.-)
3. Remove EC numbers not present in the training set

GO cleaning rules:
1. Keep only GO terms with 5-5000 training samples (default min_count=5, max_count=5000)
2. Filter by namespace (MF/CC/BP)
"""

import json
import pandas as pd
from collections import Counter, defaultdict
from pathlib import Path


def parse_ec_level(ec_number: str) -> int:
    """Parse EC number and return its level (1-4)."""
    parts = ec_number.split('.')
    return sum(1 for p in parts if p != '-')


def get_ec_prefix(ec_number: str, level: int = 3) -> str:
    """Get the first N non-dash parts of an EC number."""
    parts = ec_number.split('.')
    non_dash = [p for p in parts if p != '-']
    return '.'.join(non_dash[:level])


def normalize_ec(ec_number: str) -> str:
    """Normalize EC to standard 4-part format."""
    parts = ec_number.split('.')
    while len(parts) < 4:
        parts.append('-')
    return '.'.join(parts)


def load_ec_annotations(path: str) -> dict:
    """Load EC annotations from SIFTS JSON."""
    with open(path) as f:
        data = json.load(f)
    
    ec_map = {}
    for chain_key, info in data.items():
        ec_numbers = []
        for ann in info.get('annotations', []):
            ec = ann.get('ec_number', '')
            if ec:
                ec_numbers.append(normalize_ec(ec))
        if ec_numbers:
            ec_map[chain_key] = ec_numbers
    return ec_map


def load_go_annotations(path: str) -> dict:
    """Load GO annotations from SIFTS JSON."""
    with open(path) as f:
        data = json.load(f)
    
    go_map = {}
    for chain_key, info in data.items():
        go_terms = []
        for ann in info.get('annotations', []):
            go_id = ann.get('go_id', '')
            if go_id:
                go_terms.append(go_id)
        if go_terms:
            go_map[chain_key] = go_terms
    return go_map


def load_go_namespace_map(path: str) -> dict:
    """Load GO ID to namespace mapping."""
    with open(path) as f:
        data = json.load(f)
    return data


def clean_ec(ec_map: dict, split_dfs: dict) -> dict:
    """
    Clean EC annotations.
    Rules:
    1. Keep only level 3 and 4 EC numbers.
    2. Merge ALL level-4 to level-3 (global merge).
    3. Remove EC numbers not present in training set.
    Returns cleaned EC map and statistics.
    """
    train_keys = set(split_dfs['train']['chain_key'])
    val_keys = set(split_dfs['val']['chain_key'])
    test_keys = set(split_dfs['test']['chain_key'])
    all_keys = train_keys | val_keys | test_keys
    
    # Step 1: Filter to level 3 and 4 only
    filtered_ec = {}
    level_counts = Counter()
    
    for chain_key in all_keys:
        if chain_key not in ec_map:
            continue
        
        ecs = ec_map[chain_key]
        filtered = []
        for ec in ecs:
            level = parse_ec_level(ec)
            level_counts[level] += 1
            if level >= 3:
                filtered.append(ec)
        
        if filtered:
            filtered_ec[chain_key] = filtered
    
    print(f"  After level filtering (L3+): {len(filtered_ec)} chains")
    print(f"  Level distribution: {dict(sorted(level_counts.items()))}")
    
    # Step 2: Merge ALL level-4 to level-3
    merged_ec = {}
    merge_count = 0
    for chain_key, ecs in filtered_ec.items():
        new_ecs = []
        seen = set()
        for ec in ecs:
            level = parse_ec_level(ec)
            prefix = get_ec_prefix(ec, 3)
            
            if level == 4:
                # Always merge L4 to L3
                new_ec = f"{prefix}.-"
                merge_count += 1
            else:
                new_ec = ec
            
            if new_ec not in seen:
                seen.add(new_ec)
                new_ecs.append(new_ec)
        
        if new_ecs:
            merged_ec[chain_key] = new_ecs
    
    print(f"  Merged {merge_count} L4->L3 annotations")
    print(f"  After merging: {len(merged_ec)} chains")
    
    # Step 3: Remove EC numbers not in training set
    train_ecs = set()
    for chain_key in train_keys:
        if chain_key in merged_ec:
            train_ecs.update(merged_ec[chain_key])
    
    print(f"  Unique EC numbers in train: {len(train_ecs)}")
    
    cleaned_ec = {}
    removed_non_train = 0
    for chain_key in all_keys:
        if chain_key not in merged_ec:
            continue
        
        new_ecs = [ec for ec in merged_ec[chain_key] if ec in train_ecs]
        if new_ecs:
            cleaned_ec[chain_key] = new_ecs
        else:
            removed_non_train += 1
    
    print(f"  Removed {removed_non_train} chains with no train ECs")
    print(f"  Final: {len(cleaned_ec)} chains")
    
    return cleaned_ec


def clean_go(go_map: dict, ns_map: dict, split_dfs: dict, namespace: str, min_count: int = 5, max_count: int = 5000) -> dict:
    """
    Clean GO annotations for a specific namespace.
    Returns cleaned GO map and statistics.
    """
    train_keys = set(split_dfs['train']['chain_key'])
    val_keys = set(split_dfs['val']['chain_key'])
    test_keys = set(split_dfs['test']['chain_key'])
    all_keys = train_keys | val_keys | test_keys
    
    # Filter GO terms by namespace
    ns_go_map = {}
    for chain_key in all_keys:
        if chain_key not in go_map:
            continue
        
        terms = [go_id for go_id in go_map[chain_key] 
                 if go_id in ns_map and ns_map[go_id] == namespace]
        if terms:
            ns_go_map[chain_key] = terms
    
    print(f"  After namespace filtering ({namespace}): {len(ns_go_map)} chains")
    
    # Count GO terms in training set
    train_counts = Counter()
    for chain_key in train_keys:
        if chain_key in ns_go_map:
            for go_id in ns_go_map[chain_key]:
                train_counts[go_id] += 1
    
    print(f"  Unique GO terms in train: {len(train_counts)}")
    
    # Filter by count range
    valid_terms = {go_id for go_id, cnt in train_counts.items() 
                   if min_count <= cnt <= max_count}
    
    print(f"  GO terms in range [{min_count}, {max_count}]: {len(valid_terms)}")
    print(f"  Removed {len(train_counts) - len(valid_terms)} terms outside range")
    
    # Remove invalid terms from all splits
    cleaned_go = {}
    removed_chains = 0
    for chain_key in all_keys:
        if chain_key not in ns_go_map:
            continue
        
        new_terms = [go_id for go_id in ns_go_map[chain_key] if go_id in valid_terms]
        if new_terms:
            cleaned_go[chain_key] = new_terms
        else:
            removed_chains += 1
    
    print(f"  Removed {removed_chains} chains with no valid GO terms")
    print(f"  Final: {len(cleaned_go)} chains")
    
    return cleaned_go


def generate_output(cleaned_labels: dict, split_dfs: dict, output_dir: Path, task_name: str):
    """Generate cleaned output CSVs."""
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create label DataFrame
    rows = []
    for chain_key, labels in cleaned_labels.items():
        rows.append({
            'chain_key': chain_key,
            'labels': ';'.join(labels),
            'num_labels': len(labels)
        })
    label_df = pd.DataFrame(rows)
    
    # Save full label list
    label_df.to_csv(output_dir / f'{task_name}_labels.csv', index=False)
    
    # Save per-split chain keys
    for split_name, df in split_dfs.items():
        split_keys = set(df['chain_key'])
        valid_keys = split_keys & set(cleaned_labels.keys())
        
        split_label_df = label_df[label_df['chain_key'].isin(valid_keys)].copy()
        split_label_df.to_csv(output_dir / f'{task_name}_{split_name}.csv', index=False)
        
        # Also save the metadata subset
        meta_df = df[df['chain_key'].isin(valid_keys)].copy()
        meta_df.to_csv(output_dir / f'{task_name}_{split_name}_meta.csv', index=False)
        
        print(f"    {split_name}: {len(valid_keys)} samples")
    
    return label_df


def compute_statistics(label_df: pd.DataFrame, split_dfs: dict, cleaned_labels: dict, task_name: str) -> dict:
    """Compute and return statistics."""
    stats = {
        'task': task_name,
        'total_samples': len(label_df),
        'total_unique_labels': len(set(label for labels in cleaned_labels.values() for label in labels)),
        'splits': {}
    }
    
    for split_name, df in split_dfs.items():
        split_keys = set(df['chain_key']) & set(cleaned_labels.keys())
        split_labels = Counter()
        for key in split_keys:
            for label in cleaned_labels[key]:
                split_labels[label] += 1
        
        stats['splits'][split_name] = {
            'num_samples': len(split_keys),
            'num_labels': len(split_labels),
            'label_counts': dict(split_labels.most_common(20))
        }
    
    return stats


def main():
    base_dir = Path(__file__).resolve().parent.parent  # datasets/func_prediction/
    data_dir = base_dir.parent.parent / 'data' / 'sifts'  # POOR/data/sifts/
    # 历史执行中清洗结果覆盖写回划分目录，因此输入/输出同为 output/datasets_split_raw/
    split_dir = base_dir / 'output' / 'datasets_split_raw'
    output_dir = base_dir / 'output' / 'datasets_split_raw'
    output_dir.mkdir(exist_ok=True)
    
    print("=" * 60)
    print("Dataset Cleaning")
    print("=" * 60)
    
    # Load annotations
    print("\nLoading EC annotations...")
    ec_map = load_ec_annotations(data_dir / 'sifts_ec_annotations.json')
    print(f"  Loaded {len(ec_map)} chains with EC annotations")
    
    print("\nLoading GO annotations...")
    go_map = load_go_annotations(data_dir / 'sifts_go_annotations.json')
    print(f"  Loaded {len(go_map)} chains with GO annotations")
    
    print("\nLoading GO namespace mapping...")
    ns_map = load_go_namespace_map(base_dir / 'output' / 'go_id_to_namespace.json')
    print(f"  Loaded {len(ns_map)} GO term mappings")
    
    all_stats = {}
    
    # ===== EC Cleaning =====
    print("\n" + "=" * 60)
    print("EC Cleaning")
    print("=" * 60)
    
    ec_splits = {
        'train': pd.read_csv(split_dir / 'ec_train.csv'),
        'val': pd.read_csv(split_dir / 'ec_val.csv'),
        'test': pd.read_csv(split_dir / 'ec_test.csv')
    }
    
    print(f"\nOriginal splits: train={len(ec_splits['train'])}, val={len(ec_splits['val'])}, test={len(ec_splits['test'])}")
    
    cleaned_ec = clean_ec(ec_map, ec_splits)
    ec_label_df = generate_output(cleaned_ec, ec_splits, output_dir, 'ec')
    ec_stats = compute_statistics(ec_label_df, ec_splits, cleaned_ec, 'EC')
    all_stats['EC'] = ec_stats
    
    # ===== GO Cleaning =====
    go_tasks = {
        'go_mf': 'mf',
        'go_cc': 'cc', 
        'go_bp': 'bp'
    }
    
    for task_name, ns in go_tasks.items():
        print("\n" + "=" * 60)
        print(f"GO-{ns.upper()} Cleaning")
        print("=" * 60)
        
        go_splits = {
            'train': pd.read_csv(split_dir / f'{task_name}_train.csv'),
            'val': pd.read_csv(split_dir / f'{task_name}_val.csv'),
            'test': pd.read_csv(split_dir / f'{task_name}_test.csv')
        }
        
        print(f"\nOriginal splits: train={len(go_splits['train'])}, val={len(go_splits['val'])}, test={len(go_splits['test'])}")
        
        cleaned_go = clean_go(go_map, ns_map, go_splits, ns)
        go_label_df = generate_output(cleaned_go, go_splits, output_dir, task_name)
        go_stats = compute_statistics(go_label_df, go_splits, cleaned_go, f'GO-{ns.upper()}')
        all_stats[f'GO-{ns.upper()}'] = go_stats
    
    # Save statistics
    with open(output_dir / 'cleaning_stats.json', 'w') as f:
        json.dump(all_stats, f, indent=2)
    
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    
    for task, stats in all_stats.items():
        print(f"\n{task}:")
        print(f"  Total samples: {stats['total_samples']}")
        print(f"  Total unique labels: {stats['total_unique_labels']}")
        for split, split_stats in stats['splits'].items():
            print(f"  {split}: {split_stats['num_samples']} samples, {split_stats['num_labels']} labels")
    
    print(f"\nCleaned datasets saved to: {output_dir}")


if __name__ == '__main__':
    main()
