#!/usr/bin/env python3
"""标记 OOD_LongTail（预测标签长尾）：测试样本任一标签在训练集中频次 ≤ 10 → True。

口径：
1. 统计 splits/{task}_train.csv 中每个标签（label 列，';' 分隔）的出现次数
2. 取频次 ≤ 10 的标签作为长尾标签集合
3. 测试样本标签与长尾集合有交集（any 语义）→ OOD_LongTail=True
4. Default=False（极端长度）行不参与计算，直接留空（NaN，"未计算"语义，
   全项目统一惯例）；finalize_ood_columns.py --mode nan 仅作兜底校验

I/O：直接读写 splits/{task}_train.csv / splits/{task}_test.csv。
历史偏差说明：旧流程在极端长度链并入之前运行，读写的是中间文件
output/datasets_split_raw/{task}_train.csv 与 {task}_test_ood.csv（均未保留）；
现行版本改为在 splits/ 上运行且只标记 Default=True 行，数值与历史一致（见 DATA_PROCESS.md 4.8）。
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
    """获取训练集中频次 ≤ threshold 的标签。"""
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

    # Load test and mark (仅 Default=True 行)
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
