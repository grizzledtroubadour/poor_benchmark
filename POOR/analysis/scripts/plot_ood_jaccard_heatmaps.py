#!/usr/bin/env python3
"""Generate per-task OOD-axis Jaccard heatmaps."""
import os
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT_DIR = os.path.join(ROOT, "analysis/output/ood_jaccard_heatmaps")
os.makedirs(OUT_DIR, exist_ok=True)


def get_ood_cols(df):
    return [c for c in df.columns if c.startswith('OOD_') or c.startswith('seq_Redundancy_') or c.startswith('TM-score_')]


def select_representative_axes(cols):
    """Keep only one threshold per continuous OOD family.
    - seq_Redundancy -> keep Seq<70
    - TM-score -> keep TM<0.7
    - keep all categorical axes as-is
    """
    keep = []
    for c in cols:
        if 'seq_Redundancy_' in c:
            if c.endswith('_70'):
                keep.append(c)
        elif 'TM-score_' in c:
            if c.endswith('_0.7'):
                keep.append(c)
        else:
            keep.append(c)
    return keep


def clean_label(col):
    """Shorten labels for display."""
    col = col.replace('seq_Redundancy_70', 'Seq<70')
    col = col.replace('TM-score_0.7', 'TM<0.7')
    col = col.replace('OOD_', '')
    col = col.replace('ExtremeShort', 'ExtShort')
    col = col.replace('ExtremeLong', 'ExtLong')
    col = col.replace('SuperfamilyHoldout', 'SFholdout')
    col = col.replace('FoldHoldout', 'FoldHold')
    col = col.replace('Combinatorial', 'Comb')
    col = col.replace('LongTail', 'LongT')
    return col


def order_cols(cols):
    """Order representative OOD columns by category."""
    priority = {}
    for c in cols:
        if 'ExtremeShort' in c:
            priority[c] = (0, 0)
        elif 'ExtremeLong' in c:
            priority[c] = (0, 1)
        elif 'seq_Redundancy_70' in c:
            priority[c] = (1, 0)
        elif 'TM-score_0.7' in c:
            priority[c] = (2, 0)
        elif 'Orphan' in c:
            priority[c] = (3, 0)
        elif 'IDR' in c:
            priority[c] = (3, 1)
        elif 'NewEC' in c:
            priority[c] = (4, 0)
        elif 'LongTail' in c:
            priority[c] = (4, 1)
        elif 'FoldHoldout' in c or 'SuperfamilyHoldout' in c:
            priority[c] = (4, 2)
        elif 'Combinatorial' in c:
            priority[c] = (4, 3)
        else:
            priority[c] = (5, 0)
    return sorted(cols, key=lambda x: priority[x])


def compute_jaccard(df, col_a, col_b):
    a = df[col_a].fillna(False).astype(bool)
    b = df[col_b].fillna(False).astype(bool)
    inter = (a & b).sum()
    union = (a | b).sum()
    if union == 0:
        return np.nan
    return inter / union


task_files = [
    ('kcat', 'datasets/enzyme_kinetics_prediction/splits/kcat_test.csv'),
    ('optimal_ph', 'datasets/enzyme_optimal_ph/splits/optimal_ph_prediction_test.csv'),
    ('cath', 'datasets/fold_classification/splits/cath_test.csv'),
    ('ec', 'datasets/func_prediction/splits/ec_test.csv'),
    ('go_bp', 'datasets/func_prediction/splits/go_bp_test.csv'),
    ('go_cc', 'datasets/func_prediction/splits/go_cc_test.csv'),
    ('go_mf', 'datasets/func_prediction/splits/go_mf_test.csv'),
    ('lba', 'datasets/ligand_binding_affinity/splits/ligand_binding_affinity_test.csv'),
    ('lbs', 'datasets/ligand_binding_site/splits/ligand_binding_site_test.csv'),
    ('ppi', 'datasets/ppi_prediction/splits/ppi_test.csv'),
    ('ppis', 'datasets/ppis_prediction/splits/ppis_test.csv'),
    ('ssp', 'datasets/ss_prediction/splits/ssp_test.csv'),
]

summary = []
for task_name, path in task_files:
    df = pd.read_csv(os.path.join(ROOT, path), low_memory=False)
    cols = get_ood_cols(df)
    cols = select_representative_axes(cols)
    cols = order_cols(cols)
    labels = [clean_label(c) for c in cols]
    n = len(cols)
    mat = np.zeros((n, n))
    for i, ca in enumerate(cols):
        for j, cb in enumerate(cols):
            if i == j:
                mat[i, j] = 1.0
            else:
                mat[i, j] = compute_jaccard(df, ca, cb)

    mean_j = np.nanmean(mat[np.triu_indices_from(mat, k=1)])
    max_j = np.nanmax(mat[np.triu_indices_from(mat, k=1)])
    summary.append({
        'task': task_name,
        'n_ood_axes': n,
        'mean_jaccard': round(mean_j, 4),
        'max_jaccard': round(max_j, 4),
    })

    # Plot: mask diagonal so 1.0 values don't dominate the color scale
    mask = np.eye(n, dtype=bool)
    # also mask lower triangle to show only upper triangle (optional, cleaner)
    # mask = mask | np.tril(np.ones((n,n), dtype=bool), k=-1)

    fig, ax = plt.subplots(figsize=(max(6, n * 0.55), max(5, n * 0.5)))
    # Fixed color scale cap at 0.3; values above are annotated with asterisk
    vmax = 0.3

    # Cap visualization values at vmax so colorbar stops at 0.3
    plot_mat = mat.copy()
    plot_mat = np.where(plot_mat > vmax, vmax, plot_mat)

    # Build annotation array: show true value if <= vmax, else value + '*'
    annot_arr = np.empty_like(mat, dtype=object)
    for i in range(n):
        for j in range(n):
            if i == j:
                annot_arr[i, j] = ''
            elif mat[i, j] > vmax:
                annot_arr[i, j] = f'{mat[i,j]:.2f}*'
            else:
                annot_arr[i, j] = f'{mat[i,j]:.2f}'

    sns.heatmap(
        plot_mat,
        mask=mask,
        xticklabels=labels,
        yticklabels=labels,
        cmap='YlOrRd',
        vmin=0,
        vmax=vmax,
        square=True,
        linewidths=0.3,
        annot=annot_arr,
        fmt='',
        annot_kws={'size': 6},
        cbar_kws={'shrink': 0.8, 'label': 'Jaccard', 'ticks': [0, 0.1, 0.2, 0.3]},
        ax=ax,
    )
    ax.set_title(f'{task_name.upper()} OOD-axis Overlap (Jaccard)\nmean={mean_j:.3f}, max={max_j:.3f}', fontsize=12)
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor', fontsize=7)
    plt.setp(ax.get_yticklabels(), rotation=0, fontsize=7)
    plt.tight_layout()
    out_path = os.path.join(OUT_DIR, f'{task_name}_ood_jaccard.png')
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved {out_path}')

with open(os.path.join(OUT_DIR, 'summary.json'), 'w') as f:
    json.dump(summary, f, indent=2)

print('\nSummary:')
for s in summary:
    print(s)
