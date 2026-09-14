#!/usr/bin/env python3
"""
Convert current kcat splits into POOD-Benchmark compatible format.

Adds standard columns (unique_id, file_type, aa_seq, labels), computes
protein-level OOD labels on the test set, and ensures splits/ only contains
the three final files.

Implemented OOD columns on test set:
- Default
- OOD_ExtremeLength, OOD_ExtremeLong, OOD_ExtremeShort
- seq_Redundancy_90 ~ seq_Redundancy_30   [superseded]
- OOD_Orphan                              [superseded]
- TM-score_0.9 ~ TM-score_0.3             [superseded]
- OOD_IDR                                 [superseded]

[SUPERSEDED 说明] 本脚本为首轮（历史）一次性流水线。其中三段重计算口径已被
skill 统一实现取代，重算时必须改用对应 skill 脚本，不得再用本脚本的内联实现：
  - 序列同源（mmseqs2 seq_Redundancy_* + OOD_Orphan）：
    .skills/homology-ood-annotation/scripts/seq_homology_ood.py
  - 结构同源（Foldseek TM-score_*）：
    .skills/homology-ood-annotation/scripts/struct_homology_ood.py
  - IDR（metapredict OOD_IDR / idr_ratio）：
    .skills/ood-annotation-toolkit/scripts/idr_ood.py
现行仍有效的功能：加 POOD 标准列、长度 OOD（Default/OOD_Extreme*）。
列名 schema（pdb_file/labels/file_type → struct_file/label）其后经全项目统一迁移，
以当前 splits 表头为准（见 DATA_PROCESS.md 注意事项）。
"""
import hashlib
import json
import shutil
import subprocess
import tempfile
from itertools import groupby
from operator import itemgetter
from pathlib import Path

import metapredict
import numpy as np
import pandas as pd

# Paths
TASK_DIR = Path(__file__).resolve().parent.parent  # datasets/enzyme_kinetics_prediction
SPLITS_DIR = TASK_DIR / "splits"
OUTPUT_DIR = TASK_DIR / "output"
PDBS_DIR = TASK_DIR / "pdbs"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Parameters
NORMAL_MIN = 60
NORMAL_MAX = 1000
SEQ_REDUNDANCY_THRESHOLDS = [90, 80, 70, 60, 50, 40, 30]
TM_SCORE_THRESHOLDS = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3]
IDR_DISORDER_THRESHOLD = 0.5
IDR_FRACTION_THRESHOLD = 0.30
ORPHAN_EVALUE_THRESHOLD = 1e-3


def load_splits():
    train = pd.read_csv(SPLITS_DIR / 'kcat_train.csv', dtype=str, low_memory=False)
    val = pd.read_csv(SPLITS_DIR / 'kcat_val.csv', dtype=str, low_memory=False)
    test = pd.read_csv(SPLITS_DIR / 'kcat_test.csv', dtype=str, low_memory=False)
    return train, val, test


def add_standard_columns(df):
    """Add POOD standard columns: unique_id, file_type, aa_seq, labels.

    unique_id = AF-<uniprot>-F1-model_v6（蛋白级标识，与 struct_file 同名）。
    kcat 为（酶, 反应物）对级数据，同一 UniProt 多底物多行；为保证 split 内
    唯一（README 约定 "unique within a task"），对有重复的 unique_id 追加
    反应物内容哈希后缀 -<md5(reactant_smiles)[:8]>；已验证
    (uniprot, reactant_smiles) 组合无重复，故后缀保证唯一。
    """
    df = df.copy()
    uid = df['uniprot'].apply(lambda u: f"AF-{u}-F1-model_v6")
    dup_mask = uid.duplicated(keep=False)
    uid.loc[dup_mask] = [
        f"{u}-{hashlib.md5(smi.encode('utf-8')).hexdigest()[:8]}"
        for u, smi in zip(uid.loc[dup_mask], df.loc[dup_mask, 'reactant_smiles'])
    ]
    assert uid.is_unique, "unique_id 加后缀后仍不唯一"
    df['unique_id'] = uid
    df['file_type'] = 'pdb'
    df['aa_seq'] = df['sequence']
    df['labels'] = df['log10_value'].astype(str)
    return df


def dedup_by_sequence(df):
    return df.drop_duplicates(subset=['sequence'], keep='first').reset_index(drop=True).copy()


def dedup_by_uniprot(df):
    return df.drop_duplicates(subset=['uniprot'], keep='first').reset_index(drop=True).copy()


def write_fasta(df, path, prefix=''):
    with open(path, 'w') as f:
        for idx, row in df.iterrows():
            seq_id = f"{prefix}{row['uniprot']}_{idx}"
            f.write(f">{seq_id}\n{row['sequence']}\n")


def parse_fasta_id(seq_id):
    parts = seq_id.rsplit('_', 1)
    return parts[0], int(parts[1])


def run_mmseqs_search(query_fasta, target_fasta, output_path, tmp_dir):
    cmd = [
        'mmseqs', 'easy-search',
        str(query_fasta), str(target_fasta), str(output_path), str(tmp_dir),
        '--format-mode', '4',
        '--format-output', 'query,target,fident,evalue',
        '--threads', '8',
    ]
    subprocess.run(cmd, check=True, capture_output=True, text=True)


def compute_sequence_redundancy(trainval_unique, test_unique, work_dir):
    """
    Search test sequences against train+val. For each test sequence, compute
    max identity and best hit e-value. Return dicts keyed by test index.
    """
    print("\n" + "=" * 80)
    print("Computing sequence redundancy (MMseqs2 search test vs train+val)")
    print("=" * 80)

    train_fasta = work_dir / 'seq_trainval.fasta'
    test_fasta = work_dir / 'seq_test.fasta'
    write_fasta(trainval_unique, train_fasta, prefix='tv_')
    write_fasta(test_unique, test_fasta, prefix='te_')

    output_path = work_dir / 'seq_search.tsv'
    tmp_dir = work_dir / 'tmp_seq_search'
    run_mmseqs_search(test_fasta, train_fasta, output_path, tmp_dir)

    cols = ['query', 'target', 'fident', 'evalue']
    if output_path.stat().st_size == 0:
        print("  No hits found.")
        return {}, {}

    hits = pd.read_csv(output_path, sep='\t', skiprows=1, names=cols)
    hits['fident'] = hits['fident'].astype(float)
    hits['evalue'] = hits['evalue'].astype(float)
    hits['test_idx'] = hits['query'].apply(lambda x: int(x.rsplit('_', 1)[1]))

    max_identity = hits.groupby('test_idx')['fident'].max()
    best_evalue = hits.groupby('test_idx')['evalue'].min()

    print(f"  Test sequences with hit: {len(max_identity)}")
    print(f"  Max identity range: {max_identity.min():.4f} ~ {max_identity.max():.4f}")

    return max_identity.to_dict(), best_evalue.to_dict()


def compute_orphan_labels(max_identity, best_evalue, test_unique):
    """OOD_Orphan: no significant hit (no hit OR best e-value > 1e-3)."""
    labels = []
    for idx in range(len(test_unique)):
        if idx not in max_identity:
            labels.append(True)
        elif best_evalue.get(idx, 1.0) > ORPHAN_EVALUE_THRESHOLD:
            labels.append(True)
        else:
            labels.append(False)
    return pd.Series(labels, index=test_unique.index)


def compute_seq_redundancy_labels(max_identity, test_unique):
    """seq_Redundancy_XX: True if max identity < XX%."""
    label_dict = {}
    for thresh in SEQ_REDUNDANCY_THRESHOLDS:
        col = f'seq_Redundancy_{thresh}'
        labels = []
        for idx in range(len(test_unique)):
            identity = max_identity.get(idx, 0.0)
            labels.append(identity < thresh / 100.0)
        label_dict[col] = pd.Series(labels, index=test_unique.index)
    return label_dict


def run_foldseek_search(test_pdbs, trainval_pdbs, work_dir):
    """
    Create Foldseek DB from train+val PDBs and search test PDBs.
    Foldseek createdb works on directories, so we symlink PDBs into temp dirs.
    Returns path to result TSV with query, target, alntmscore.
    """
    print("\n" + "=" * 80)
    print("Computing structural similarity (Foldseek test vs train+val)")
    print("=" * 80)

    # Create temp dirs with symlinks
    trainval_dir = work_dir / 'trainval_pdbs'
    test_dir = work_dir / 'test_pdbs'
    trainval_dir.mkdir()
    test_dir.mkdir()

    for p in trainval_pdbs:
        (trainval_dir / p.name).symlink_to(p.resolve())
    for p in test_pdbs:
        (test_dir / p.name).symlink_to(p.resolve())

    # Create DB for train+val
    trainval_db = work_dir / 'trainval_db'
    subprocess.run(['foldseek', 'createdb', str(trainval_dir), str(trainval_db), '--threads', '8'],
                   check=True, capture_output=True, text=True)

    # Create DB for test
    test_db = work_dir / 'test_db'
    subprocess.run(['foldseek', 'createdb', str(test_dir), str(test_db), '--threads', '8'],
                   check=True, capture_output=True, text=True)

    # Search
    result_db = work_dir / 'foldseek_result'
    tmp_dir = work_dir / 'tmp_foldseek'
    subprocess.run(['foldseek', 'search', str(test_db), str(trainval_db), str(result_db), str(tmp_dir),
                    '-a', '--threads', '8'],
                   check=True, capture_output=True, text=True)

    # Convert to TSV
    result_tsv = work_dir / 'foldseek_result.tsv'
    subprocess.run(['foldseek', 'convertalis', str(test_db), str(trainval_db), str(result_db), str(result_tsv),
                    '--format-output', 'query,target,alntmscore', '--threads', '8'],
                   check=True, capture_output=True, text=True)

    return result_tsv


def compute_tm_score_labels(test_unique, trainval_unique, work_dir):
    """TM-score_X.X: True if max TM-score < threshold."""
    # Map uniprot to PDB path
    test_pdbs = []
    test_uniprot_order = []
    for _, row in test_unique.iterrows():
        up = row['uniprot']
        pdb_path = PDBS_DIR / f"AF-{up}-F1-model_v6.pdb"
        if pdb_path.exists():
            test_pdbs.append(pdb_path)
            test_uniprot_order.append(up)

    trainval_pdbs = []
    trainval_uniprots = []
    for _, row in trainval_unique.iterrows():
        up = row['uniprot']
        pdb_path = PDBS_DIR / f"AF-{up}-F1-model_v6.pdb"
        if pdb_path.exists():
            trainval_pdbs.append(pdb_path)
            trainval_uniprots.append(up)

    print(f"  Test PDBs available: {len(test_pdbs)}/{len(test_unique)}")
    print(f"  Train+Val PDBs available: {len(trainval_pdbs)}/{len(trainval_unique)}")

    if len(test_pdbs) == 0 or len(trainval_pdbs) == 0:
        print("  Warning: insufficient PDB files for Foldseek. Returning all NaN.")
        return {f'TM-score_{t:.1f}': pd.Series(np.nan, index=test_unique.index) for t in TM_SCORE_THRESHOLDS}

    result_tsv = run_foldseek_search(test_pdbs, trainval_pdbs, work_dir)

    # Parse results
    cols = ['query', 'target', 'alntmscore']
    if result_tsv.stat().st_size == 0:
        print("  No Foldseek hits found.")
        max_tm = {up: 0.0 for up in test_uniprot_order}
    else:
        hits = pd.read_csv(result_tsv, sep='\t', names=cols)
        hits['alntmscore'] = hits['alntmscore'].astype(float)
        # query is PDB filename stem
        hits['query_uniprot'] = hits['query'].str.replace('AF-', '', regex=False).str.replace('-F1-model_v6', '', regex=False)
        max_tm_series = hits.groupby('query_uniprot')['alntmscore'].max()
        max_tm = max_tm_series.to_dict()

    label_dict = {}
    for thresh in TM_SCORE_THRESHOLDS:
        col = f'TM-score_{thresh:.1f}'
        labels = []
        for _, row in test_unique.iterrows():
            up = row['uniprot']
            tm = max_tm.get(up, 0.0)
            labels.append(tm < thresh)
        label_dict[col] = pd.Series(labels, index=test_unique.index)

    return label_dict


def compute_idr_regions(disorder_scores):
    """Find contiguous regions with score > threshold; return total length."""
    positions = np.where(disorder_scores > IDR_DISORDER_THRESHOLD)[0]
    if len(positions) == 0:
        return 0

    # Group consecutive positions
    total_len = 0
    for _, group in groupby(enumerate(positions), lambda x: x[0] - x[1]):
        group_list = list(group)
        total_len += len(group_list)

    return total_len


def compute_idr_labels(test_unique):
    """OOD_IDR: contiguous disorder regions > 30% of sequence length."""
    print("\n" + "=" * 80)
    print("Computing IDR labels (metapredict)")
    print("=" * 80)

    labels = []
    ood_count = 0
    for idx, row in test_unique.iterrows():
        seq = row['sequence']
        disorder = metapredict.predict_disorder(seq)
        idr_len = compute_idr_regions(disorder)
        idr_frac = idr_len / len(seq)
        is_idr = idr_frac > IDR_FRACTION_THRESHOLD
        labels.append(is_idr)
        if is_idr:
            ood_count += 1

        if (idx + 1) % 500 == 0 or (idx + 1) == len(test_unique):
            print(f"  Processed {idx + 1}/{len(test_unique)}, IDR OOD so far: {ood_count}")

    print(f"  IDR OOD total: {ood_count}")
    return pd.Series(labels, index=test_unique.index)


def map_unique_labels_to_rows(unique_labels, unique_df, full_df):
    """Map labels from unique-sequence dataframe to full dataframe rows by sequence."""
    seq_to_indices = {}
    for idx, row in full_df.iterrows():
        seq_to_indices.setdefault(row['sequence'], []).append(idx)

    full_labels = pd.Series(False, index=full_df.index)
    for unique_idx, is_ood in unique_labels.items():
        if is_ood:
            seq = unique_df.iloc[unique_idx]['sequence']
            for full_idx in seq_to_indices.get(seq, []):
                full_labels.loc[full_idx] = True
    return full_labels


def map_unique_dict_to_rows(label_dict, unique_df, full_df):
    """Map a dict of label Series from unique df to full df."""
    seq_to_indices = {}
    for idx, row in full_df.iterrows():
        seq_to_indices.setdefault(row['sequence'], []).append(idx)

    full_dict = {}
    for col, unique_labels in label_dict.items():
        full_labels = pd.Series(False, index=full_df.index)
        for unique_idx, is_ood in unique_labels.items():
            if is_ood:
                seq = unique_df.iloc[unique_idx]['sequence']
                for full_idx in seq_to_indices.get(seq, []):
                    full_labels.loc[full_idx] = True
        full_dict[col] = full_labels
    return full_dict


def main():
    print("=" * 80)
    print("Creating POOD-Benchmark compatible splits")
    print("=" * 80)

    train, val, test = load_splits()
    print(f"Loaded: train={len(train):,}, val={len(val):,}, test={len(test):,}")

    # Add standard columns
    train = add_standard_columns(train)
    val = add_standard_columns(val)
    test = add_standard_columns(test)

    # Compute length-based OOD labels for test
    test['seq_len'] = test['aa_seq'].str.len().astype(int)
    test['OOD_ExtremeLength'] = (test['seq_len'] < NORMAL_MIN) | (test['seq_len'] > NORMAL_MAX)
    test['OOD_ExtremeLong'] = test['seq_len'] > NORMAL_MAX
    test['OOD_ExtremeShort'] = test['seq_len'] < NORMAL_MIN
    test['Default'] = ~test['OOD_ExtremeLength']

    # Deduplicate for sequence-based computations
    trainval = pd.concat([train, val], ignore_index=True)
    trainval_unique = dedup_by_sequence(trainval)
    test_unique = dedup_by_sequence(test)

    print(f"Unique sequences: train+val={len(trainval_unique):,}, test={len(test_unique):,}")

    with tempfile.TemporaryDirectory(prefix='pood_ood_', dir='.') as tmp:
        work_dir = Path(tmp)

        # Sequence redundancy and orphan
        # [superseded] 重算改用 .skills/homology-ood-annotation/scripts/seq_homology_ood.py
        max_identity, best_evalue = compute_sequence_redundancy(trainval_unique, test_unique, work_dir)
        orphan_unique = compute_orphan_labels(max_identity, best_evalue, test_unique)
        seq_red_labels_unique = compute_seq_redundancy_labels(max_identity, test_unique)

        # Structural similarity
        # [superseded] 重算改用 .skills/homology-ood-annotation/scripts/struct_homology_ood.py
        trainval_unique_up = dedup_by_uniprot(trainval)
        test_unique_up = dedup_by_uniprot(test)
        tm_labels_unique = compute_tm_score_labels(test_unique_up, trainval_unique_up, work_dir)

        # IDR
        # [superseded] 重算改用 .skills/ood-annotation-toolkit/scripts/idr_ood.py
        idr_unique = compute_idr_labels(test_unique)

    # Map unique labels back to full test rows
    test['OOD_Orphan'] = map_unique_labels_to_rows(orphan_unique, test_unique, test)

    seq_red_full = map_unique_dict_to_rows(seq_red_labels_unique, test_unique, test)
    for col, labels in seq_red_full.items():
        test[col] = labels

    # TM-score labels are at UniProt level; map by uniprot
    up_to_indices = {}
    for idx, row in test.iterrows():
        up_to_indices.setdefault(row['uniprot'], []).append(idx)

    for col, unique_labels in tm_labels_unique.items():
        full_labels = pd.Series(False, index=test.index)
        for unique_idx, is_ood in unique_labels.items():
            if is_ood:
                up = test_unique_up.iloc[unique_idx]['uniprot']
                for full_idx in up_to_indices.get(up, []):
                    full_labels.loc[full_idx] = True
        test[col] = full_labels

    test['OOD_IDR'] = map_unique_labels_to_rows(idr_unique, test_unique, test)

    # Default=False（极端长度）行：除 OOD_ExtremeShort/Long 外的 OOD 列直接置 NaN
    # （"未计算"语义，全项目统一惯例，流水线位置写入；
    #   finalize_ood_columns.py --mode nan 仅作兜底校验）
    nd_mask = ~test['Default']
    if nd_mask.any():
        for col in ([f'seq_Redundancy_{t}' for t in SEQ_REDUNDANCY_THRESHOLDS] +
                    ['OOD_Orphan'] +
                    [f'TM-score_{t:.1f}' for t in TM_SCORE_THRESHOLDS] +
                    ['OOD_IDR']):
            test[col] = test[col].astype(object)
            test.loc[nd_mask, col] = pd.NA
        print(f"  Default=False rows -> NaN: {int(nd_mask.sum())} rows")

    # Reorder columns: standard first, then original columns, then OOD columns at the end
    standard_cols = ['unique_id', 'file_type', 'aa_seq', 'labels']
    ood_cols = (['Default', 'OOD_ExtremeLength', 'OOD_ExtremeLong', 'OOD_ExtremeShort'] +
                [f'seq_Redundancy_{t}' for t in SEQ_REDUNDANCY_THRESHOLDS] +
                ['OOD_Orphan'] +
                [f'TM-score_{t:.1f}' for t in TM_SCORE_THRESHOLDS] +
                ['OOD_IDR'])

    other_cols = [c for c in test.columns if c not in standard_cols + ood_cols]
    test = test[standard_cols + other_cols + ood_cols]

    # Train/val only have standard + original columns, with standard first
    train_other = [c for c in train.columns if c not in standard_cols + ood_cols]
    val_other = [c for c in val.columns if c not in standard_cols + ood_cols]
    train = train[standard_cols + train_other]
    val = val[standard_cols + val_other]

    # Save final splits
    print("\n" + "=" * 80)
    print("Saving final POOD-compatible splits")
    print("=" * 80)
    train.to_csv(SPLITS_DIR / 'kcat_train.csv', index=False)
    val.to_csv(SPLITS_DIR / 'kcat_val.csv', index=False)
    test.to_csv(SPLITS_DIR / 'kcat_test.csv', index=False)
    print(f"  Train: {len(train):,} rows, {len(train.columns)} cols")
    print(f"  Val:   {len(val):,} rows, {len(val.columns)} cols")
    print(f"  Test:  {len(test):,} rows, {len(test.columns)} cols")

    # Save OOD subset files to output/ (not splits/)
    for col in ood_cols:
        subset = test[test[col] == True].copy()
        if len(subset) > 0:
            subset.to_csv(OUTPUT_DIR / f'kcat_test_ood_{col}.csv', index=False)
            print(f"  Saved {col} subset: {len(subset):,} rows -> output/kcat_test_ood_{col}.csv")

    # Summary
    summary = {
        'train_rows': int(len(train)),
        'val_rows': int(len(val)),
        'test_rows': int(len(test)),
        'ood_counts': {col: int(test[col].sum()) for col in ood_cols},
        'thresholds': {
            'seq_redundancy': SEQ_REDUNDANCY_THRESHOLDS,
            'tm_score': TM_SCORE_THRESHOLDS,
            'orphan_evalue': ORPHAN_EVALUE_THRESHOLD,
            'idr_disorder_score': IDR_DISORDER_THRESHOLD,
            'idr_fraction': IDR_FRACTION_THRESHOLD,
        }
    }
    summary_path = OUTPUT_DIR / 'kcat_pood_summary.json'
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"\nSaved summary to: {summary_path}")
    print(json.dumps(summary, indent=2))

    print("\n" + "=" * 80)
    print("Done!")
    print("=" * 80)


if __name__ == '__main__':
    main()
