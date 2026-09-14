#!/usr/bin/env python3
"""UniProt 级代表链筛选（DATA_PROCESS.md 2.2 的重实现，补缺口脚本）。

按 UniProt ID 分组，每组仅保留一条最优 PDB chain。候选链按以下优先级逐层比较：
  1. 提取模式        A 档 > B 档（C 档不回退进入候选池）
  2. 来源质量        Swiss-Prot > TrEMBL（组内恒同，实际无区分力）
  3. EC 优先         有 EC > 无 EC
  4. 实验 GO 丰富度  降序
  5. GO 丰富度       降序
  6. 分辨率          升序（缺失排最后）
  7. R-free          升序（缺失排最后）
  8. 映射长度        降序
  9. PDB ID + chain ID 字母序升序（确定性 tiebreak）

候选池即时构建自以下来源（均为现行产物）：
  - data/sifts/annotated_chains.csv        初筛链元数据（mapped_length ∈ [60,1000] 且有 EC/GO）
  - data/sifts/sifts_uniprot_mapping.json  SIFTS 残基段 → A/B 档判定（存在 pdb_beg/pdb_end 均有效的段 → A 档）
  - data/sifts/sifts_go_annotations.json   GO 证据代码 → 实验 GO 数（ijson 流式，10 种实验证据代码）
  - output/uniprot_metadata.json           reviewed 状态（Swiss-Prot > TrEMBL）

输出：output/representative_chains.csv（供 03_extract_representative_structures.py 使用；
历史执行为 55,124 条代表链）。

已知口径偏差（详见 DATA_PROCESS.md 注意事项）：
  - 历史执行在候选入池前对 B 档链做边界推断 identity ≥ 0.90 过滤（419,370/439,427 通过），
    其产物 data/intermediate/b_boundary_inference.jsonl 未随仓库迁移。本脚本默认 B 档全量入池；
    若重做了边界推断，可用 --b-keep-list 传入通过的 chain_key 清单（每行一个或含 chain_key 列的 CSV）
    以复现历史过滤。
"""

import argparse
import csv
import json
from pathlib import Path

import ijson
import pandas as pd

TASK_DIR = Path(__file__).resolve().parent.parent  # datasets/func_prediction/
POOR_ROOT = TASK_DIR.parent.parent                 # POOR/
SIFTS_DIR = POOR_ROOT / 'data' / 'sifts'
OUT_DIR = TASK_DIR / 'output'

EXP_EVIDENCE = {'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'TAS', 'IC', 'HTP', 'HDA'}


def load_filtered_chains():
    """初筛：mapped_length ∈ [60, 1000] 且有 EC 或 GO 标注（对应文档 1.5）。"""
    df = pd.read_csv(SIFTS_DIR / 'annotated_chains.csv', low_memory=False)
    df = df[(df['has_ec'] | df['has_go']) &
            (df['mapped_length'] >= 60) & (df['mapped_length'] <= 1000)].copy()
    # 剔除 chain_id / pdb_id 缺失的脏行（无法定位结构，不参与筛选）
    df = df[df['pdb_id'].notna() & df['chain_id'].notna()].copy()
    df['chain_key'] = df['pdb_id'] + '_' + df['chain_id']
    return df


def classify_extract_mode(chain_keys):
    """A 档：SIFTS 段中存在 pdb_beg/pdb_end 均有效者；否则 B 档。"""
    with open(SIFTS_DIR / 'sifts_uniprot_mapping.json') as f:
        mapping = json.load(f)
    modes = {}
    for ck in chain_keys:
        info = mapping.get(ck)
        mode = 'B'
        if info:
            for seg in info.get('segments', []):
                if seg.get('pdb_beg') and seg.get('pdb_end'):
                    mode = 'A'
                    break
        modes[ck] = mode
    return modes


def count_go_evidence(chain_keys):
    """流式解析 SIFTS GO 注释，统计每条链的实验证据 GO 数与 GO 总数。"""
    go_count, exp_go_count = {}, {}
    with open(SIFTS_DIR / 'sifts_go_annotations.json', 'rb') as f:
        for chain_key, record in ijson.kvitems(f, ''):
            if chain_key not in chain_keys:
                continue
            seen, seen_exp = set(), set()
            for ann in record.get('annotations', []):
                go_id = ann.get('go_id', '')
                if not go_id:
                    continue
                seen.add(go_id)
                if ann.get('evidence', '') in EXP_EVIDENCE:
                    seen_exp.add(go_id)
            go_count[chain_key] = len(seen)
            exp_go_count[chain_key] = len(seen_exp)
    return go_count, exp_go_count


def load_b_keep_list(path):
    """读取 B 档边界推断通过清单（每行一个 chain_key，或含 chain_key 列的 CSV）。"""
    with open(path) as f:
        first = f.readline().strip()
    if first == 'chain_key' or 'chain_key' in first.split(','):
        return set(pd.read_csv(path)['chain_key'])
    with open(path) as f:
        return {line.strip() for line in f if line.strip()}


def main():
    ap = argparse.ArgumentParser(description='UniProt-level representative chain selection')
    ap.add_argument('--out', type=Path, default=OUT_DIR / 'representative_chains.csv',
                    help='输出 CSV（默认 output/representative_chains.csv）')
    ap.add_argument('--b-keep-list', type=Path, default=None,
                    help='可选：B 档边界推断 identity≥0.90 通过的 chain_key 清单（复现历史过滤）')
    args = ap.parse_args()

    print('[1/4] Loading filtered candidate chains...')
    df = load_filtered_chains()
    print(f'  candidates: {len(df):,} chains, {df["uniprot"].nunique():,} UniProts')

    print('[2/4] Classifying extract mode (A/B) from SIFTS segments...')
    chain_keys = set(df['chain_key'])
    df['extract_mode'] = df['chain_key'].map(classify_extract_mode(chain_keys))
    if args.b_keep_list:
        keep = load_b_keep_list(args.b_keep_list)
        before = len(df)
        df = df[(df['extract_mode'] == 'A') | (df['chain_key'].isin(keep))].copy()
        print(f'  B-keep filter: {before:,} -> {len(df):,} chains')
    print(f'  A: {(df["extract_mode"] == "A").sum():,}, B: {(df["extract_mode"] == "B").sum():,}')

    print('[3/4] Counting GO evidence (streaming) + reviewed status...')
    go_count, exp_go_count = count_go_evidence(chain_keys)
    df['go_count'] = df['chain_key'].map(go_count).fillna(0).astype(int)
    df['exp_go_count'] = df['chain_key'].map(exp_go_count).fillna(0).astype(int)
    with open(OUT_DIR / 'uniprot_metadata.json') as f:
        uniprot_meta = json.load(f)
    df['is_reviewed'] = df['uniprot'].map(
        lambda u: bool(uniprot_meta.get(u, {}).get('reviewed', False)))

    print('[4/4] Selecting one representative chain per UniProt...')
    df['_mode_pri'] = df['extract_mode'].map({'A': 0, 'B': 1})
    df['_has_ec'] = df['has_ec'].astype(int)
    df['_reviewed'] = df['is_reviewed'].astype(int)
    df = df.sort_values(
        ['uniprot', '_mode_pri', '_reviewed', '_has_ec', 'exp_go_count', 'go_count',
         'resolution', 'r_free', 'mapped_length', 'pdb_id', 'chain_id'],
        ascending=[True, True, False, False, False, False,
                   True, True, False, True, True],
        na_position='last',
        kind='mergesort',  # stable sort for deterministic tiebreak
    )
    rep = df.drop_duplicates(subset='uniprot', keep='first').copy()

    out_cols = ['uniprot', 'pdb_id', 'chain_id', 'chain_key', 'extract_mode',
                'is_reviewed', 'ec_count', 'exp_go_count', 'go_count',
                'resolution', 'r_free', 'mapped_length', 'method']
    rep = rep[out_cols].sort_values('uniprot').reset_index(drop=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    rep.to_csv(args.out, index=False)

    print(f'\nRepresentative chains: {len(rep):,} '
          f'(A: {(rep["extract_mode"] == "A").sum():,}, B: {(rep["extract_mode"] == "B").sum():,})')
    print(f'  Swiss-Prot: {rep["is_reviewed"].sum():,}, with EC: {(rep["ec_count"] > 0).sum():,}, '
          f'with exp GO: {(rep["exp_go_count"] > 0).sum():,}')
    print(f'Saved to {args.out}')


if __name__ == '__main__':
    main()
