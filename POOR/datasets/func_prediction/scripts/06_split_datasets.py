#!/usr/bin/env python3
"""
按任务（EC / GO-MF / GO-CC / GO-BP）和时间划分数据集。

划分规则：
- 时间切分：按 deposition_date 排序后按分位数切分
- 划分比例：EC = 72:8:20；GO-MF / GO-CC / GO-BP = 63:7:30（训练:验证:测试）
- EC 测试集：不过滤
- GO 测试集：只保留 reviewed 条目 或 携带任意实验证据 GO 标注的条目
  （实验证据：EXP, IDA, IPI, IMP, IGI, IEP, TAS, IC, HTP, HDA）
"""

import csv
import json
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import ijson
import pandas as pd

TASK_DIR = Path(__file__).resolve().parent.parent  # datasets/func_prediction/
POOR_ROOT = TASK_DIR.parent.parent  # POOR/
OUT_DIR = TASK_DIR / 'output'
EXP_EVIDENCE = {'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'TAS', 'IC', 'HTP', 'HDA'}

GO_NAMESPACE_MAP = {
    'molecular_function': 'mf',
    'cellular_component': 'cc',
    'biological_process': 'bp',
}


def parse_go_obo(obopath):
    """解析 go-basic.obo，返回 go_id -> namespace 映射。"""
    go2ns = {}
    current_id = None
    current_ns = None
    with open(obopath) as f:
        for line in f:
            line = line.strip()
            if line == '[Term]':
                if current_id and current_ns:
                    go2ns[current_id] = current_ns
                current_id = None
                current_ns = None
            elif line.startswith('id: GO:'):
                current_id = 'GO:' + line.split('GO:')[1].strip()
            elif line.startswith('namespace:'):
                current_ns = line.split('namespace:')[1].strip()
        # Last term
        if current_id and current_ns:
            go2ns[current_id] = current_ns
    return go2ns


def load_go_annotations_stream(go_json_path, target_keys):
    """用 ijson 流式解析 GO annotations，只保留 target_keys。"""
    chain_go = defaultdict(lambda: {'mf': False, 'cc': False, 'bp': False,
                                     'exp_mf': False, 'exp_cc': False, 'exp_bp': False,
                                     'has_any_exp': False})
    
    with open(go_json_path, 'rb') as f:
        # ijson 解析顶层对象
        for chain_key, record in ijson.kvitems(f, ''):
            if chain_key not in target_keys:
                continue
            for anno in record.get('annotations', []):
                go_id = anno.get('go_id', '')
                evidence = anno.get('evidence', '')
                is_exp = evidence in EXP_EVIDENCE
                
                ns = go2ns.get(go_id, '')
                if is_exp:
                    chain_go[chain_key]['has_any_exp'] = True
                if ns == 'molecular_function':
                    chain_go[chain_key]['mf'] = True
                    if is_exp:
                        chain_go[chain_key]['exp_mf'] = True
                elif ns == 'cellular_component':
                    chain_go[chain_key]['cc'] = True
                    if is_exp:
                        chain_go[chain_key]['exp_cc'] = True
                elif ns == 'biological_process':
                    chain_go[chain_key]['bp'] = True
                    if is_exp:
                        chain_go[chain_key]['exp_bp'] = True
    
    return chain_go


def time_split(df, train_ratio=0.63, val_ratio=0.07, test_ratio=0.30):
    """按 deposition_date 排序后切分。"""
    assert abs(train_ratio + val_ratio + test_ratio - 1.0) < 1e-6
    df = df.sort_values('deposition_date').reset_index(drop=True)
    n = len(df)
    train_end = int(n * train_ratio)
    val_end = int(n * (train_ratio + val_ratio))
    
    train = df.iloc[:train_end].copy()
    val = df.iloc[train_end:val_end].copy()
    test = df.iloc[val_end:].copy()
    
    return train, val, test


def main():
    print("Loading final representatives...")
    df = pd.read_csv(OUT_DIR / 'final_nonredundant_representatives.csv')
    
    print("Loading deposition dates...")
    with open(OUT_DIR / 'deposition_dates.json') as f:
        dates = json.load(f)
    df['deposition_date'] = df['chain_key'].map(dates)
    df = df[df['deposition_date'].notna()].copy()
    print(f"  Chains with deposition date: {len(df)} / 47852")
    
    print("Parsing GO OBO...")
    global go2ns
    # go-basic.obo 需从 Gene Ontology 下载（https://purl.obolibrary.org/obo/go/go-basic.obo）
    go2ns = parse_go_obo(OUT_DIR / 'go-basic.obo')
    print(f"  GO terms mapped: {len(go2ns)}")
    
    print("Loading GO annotations (streaming)...")
    target_keys = set(df['chain_key'])
    chain_go = load_go_annotations_stream(
        POOR_ROOT / 'data/sifts/sifts_go_annotations.json',
        target_keys
    )
    print(f"  Chains with GO annotations: {len(chain_go)}")
    
    # Add GO flags to df
    for ns in ['mf', 'cc', 'bp']:
        df[f'has_go_{ns}'] = df['chain_key'].map(lambda k: chain_go.get(k, {}).get(ns, False))
        df[f'has_exp_go_{ns}'] = df['chain_key'].map(lambda k: chain_go.get(k, {}).get(f'exp_{ns}', False))
    
    # Datasets
    datasets = {}
    
    # EC dataset
    ec_df = df[df['ec_count'] > 0].copy()
    datasets['ec'] = ec_df
    print(f"\nEC dataset: {len(ec_df)} chains")
    
    # GO sub-datasets
    for ns, label in [('mf', 'GO-MF'), ('cc', 'GO-CC'), ('bp', 'GO-BP')]:
        go_df = df[df[f'has_go_{ns}']].copy()
        datasets[f'go_{ns}'] = go_df
        print(f"{label} dataset: {len(go_df)} chains")
    
    # Split and filter
    out_dir = OUT_DIR / 'datasets_split_raw'
    out_dir.mkdir(exist_ok=True)
    
    summary = []
    
    for name, d in datasets.items():
        if len(d) == 0:
            print(f"\n{name}: empty, skipping")
            continue
        
        if name == 'ec':
            train, val, test = time_split(d, train_ratio=0.72, val_ratio=0.08, test_ratio=0.20)
        else:
            train, val, test = time_split(d)
        
        # Test set filtering
        if name == 'ec':
            # EC: no filtering
            test_filtered = test.copy()
            filter_reason = "no filter"
        else:
            # GO: reviewed OR any experimental GO evidence (including HDA)
            test_filtered = test[
                (test['is_reviewed'] == True) | 
                (test['chain_key'].map(lambda k: chain_go.get(k, {}).get('has_any_exp', False)))
            ].copy()
            filter_reason = "reviewed or exp_go"
        
        # Save
        train.to_csv(out_dir / f'{name}_train.csv', index=False)
        val.to_csv(out_dir / f'{name}_val.csv', index=False)
        test_filtered.to_csv(out_dir / f'{name}_test.csv', index=False)
        
        print(f"\n{name}:")
        print(f"  Train: {len(train)}  ({train['deposition_date'].min()} ~ {train['deposition_date'].max()})")
        print(f"  Val:   {len(val)}  ({val['deposition_date'].min()} ~ {val['deposition_date'].max()})")
        print(f"  Test:  {len(test)} -> {len(test_filtered)} after filtering ({filter_reason})")
        if len(test_filtered) > 0:
            print(f"         ({test_filtered['deposition_date'].min()} ~ {test_filtered['deposition_date'].max()})")
        
        summary.append({
            'dataset': name,
            'total': len(d),
            'train': len(train),
            'val': len(val),
            'test_raw': len(test),
            'test_filtered': len(test_filtered),
            'filter': filter_reason,
            'train_date_max': train['deposition_date'].max() if len(train) > 0 else None,
            'val_date_min': val['deposition_date'].min() if len(val) > 0 else None,
            'val_date_max': val['deposition_date'].max() if len(val) > 0 else None,
            'test_date_min': test['deposition_date'].min() if len(test) > 0 else None,
        })
    
    # Save summary
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(out_dir / 'split_summary_new.csv', index=False)
    print(f"\nSummary saved to {out_dir / 'split_summary_new.csv'}")


if __name__ == '__main__':
    main()
