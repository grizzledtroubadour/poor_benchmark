#!/usr/bin/env python3
"""Generate two main figures for OOD-axis orthogonality:
1. Cross-task averaged Jaccard heatmap (representative axes, no extreme length axes).
2. Line plot of Jaccard per task for selected axis pairs.
"""
import os
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
# Keep text as editable TrueType in vector outputs (Illustrator-friendly)
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
import seaborn as sns

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT_DIR = os.path.join(ROOT, "analysis/output/ood_jaccard_heatmaps")
os.makedirs(OUT_DIR, exist_ok=True)


def get_ood_cols(df):
    return [c for c in df.columns if c.startswith('OOD_') or c.startswith('seq_Redundancy_') or c.startswith('TM-score_')]


# Axis family grouping rules for cross-task aggregation.
# For each task, pick one representative column per canonical family:
# - seq_Redundancy_40
# - TM-score_0.4
# - Orphan
# - IDR
# - NewEC*: pick the variant with largest True count (incl. CATH's OOD_NewEC EC7 prior)
# - LongTail*: group all label-frequency longtails (excluding EC-level ones if separate)
# - FoldHoldout / SuperfamilyHoldout / Combinatorial kept as-is if present

def canonical_family(col):
    if 'seq_Redundancy_40' in col:
        return 'Seq<40'
    if 'TM-score_0.4' in col:
        return 'TM<0.4'
    if col == 'OOD_Orphan':
        return 'Orphan'
    if col == 'OOD_IDR':
        return 'IDR'
    if col.startswith('OOD_NewEC'):
        # OOD_NewEC (CATH EC7 pre-split prior) and OOD_NewEC_L4/L3
        return 'New EC'
    if col == 'OOD_LongTail':
        return 'Long-tail (label-frequency)'
    if col.startswith('OOD_LongTail_EC'):
        return 'Long-tail (EC-level)'
    if col.startswith(('OOD_LongTail_pH', 'OOD_LongTail_KcatBin', 'OOD_LongTail_AffBin')):
        # regression value-bin long-tails: merge into label-frequency long-tail
        return 'Long-tail (label-frequency)'
    if col == 'OOD_FoldHoldout':
        return 'Fold holdout'
    if col == 'OOD_SuperfamilyHoldout':
        return 'Superfamily holdout'
    if col == 'OOD_Combinatorial':
        return 'Combinatorial'
    return None


def select_representative_axes(df):
    """Return a dict family -> selected column name, choosing the variant with the most True counts."""
    cols = get_ood_cols(df)
    families = {}
    for c in cols:
        fam = canonical_family(c)
        if fam is None:
            continue
        if fam not in families:
            families[fam] = []
        families[fam].append(c)
    selected = {}
    for fam, fam_cols in families.items():
        # pick the column with the largest number of True samples
        best = max(fam_cols, key=lambda c: df[c].fillna(False).sum())
        selected[fam] = best
    return selected


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

# Fixed canonical axis order for cross-task plots (figure 2a uses the stricter 40/0.4 rungs)
canonical_axes = [
    'Seq<40',
    'TM<0.4',
    'Orphan',
    'IDR',
    'New EC',
    'Long-tail (label-frequency)',
    'Long-tail (EC-level)',
    'Fold holdout',
    'Superfamily holdout',
    'Combinatorial',
]

# Collect per-task Jaccard matrices in canonical space
task_mats = {}
task_axes = {}
for task_name, path in task_files:
    df = pd.read_csv(os.path.join(ROOT, path), low_memory=False)
    selected = select_representative_axes(df)
    # invert to family list
    present_axes = [fam for fam in canonical_axes if fam in selected]
    n = len(present_axes)
    mat = np.full((n, n), np.nan)
    for i, fam_i in enumerate(present_axes):
        for j, fam_j in enumerate(present_axes):
            if i == j:
                mat[i, j] = 1.0
            else:
                mat[i, j] = compute_jaccard(df, selected[fam_i], selected[fam_j])
    task_mats[task_name] = mat
    task_axes[task_name] = present_axes

# Build 3D array: task x canonical_axis x canonical_axis
m = len(canonical_axes)
stacked = np.full((len(task_files), m, m), np.nan)
for t, (task_name, _) in enumerate(task_files):
    present = task_axes[task_name]
    for i, fam_i in enumerate(present):
        for j, fam_j in enumerate(present):
            idx_i = canonical_axes.index(fam_i)
            idx_j = canonical_axes.index(fam_j)
            stacked[t, idx_i, idx_j] = task_mats[task_name][i, j]

mean_mat = np.nanmean(stacked, axis=0)

# Figure 1: Cross-task averaged heatmap
fig, ax = plt.subplots(figsize=(7, 6))
mask = np.eye(m, dtype=bool)
plot_mat = mean_mat.copy()
vmax = 0.3
plot_mat = np.where(plot_mat > vmax, vmax, plot_mat)

annot_arr = np.empty_like(mean_mat, dtype=object)
for i in range(m):
    for j in range(m):
        if i == j:
            annot_arr[i, j] = ''
        elif np.isnan(mean_mat[i, j]):
            annot_arr[i, j] = ''
        elif mean_mat[i, j] > vmax:
            annot_arr[i, j] = f'{mean_mat[i,j]:.2f}*'
        else:
            annot_arr[i, j] = f'{mean_mat[i,j]:.2f}'

sns.heatmap(
    plot_mat,
    mask=mask,
    xticklabels=canonical_axes,
    yticklabels=canonical_axes,
    cmap='YlOrRd',
    vmin=0,
    vmax=vmax,
    square=True,
    linewidths=0.3,
    annot=annot_arr,
    fmt='',
    annot_kws={'size': 8},
    cbar_kws={'shrink': 0.8, 'label': 'Jaccard', 'ticks': [0, 0.1, 0.2, 0.3]},
    ax=ax,
)
ax.set_title('Cross-task Average OOD-axis Overlap (Jaccard)', fontsize=13)
plt.setp(ax.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor', fontsize=9)
plt.setp(ax.get_yticklabels(), rotation=0, fontsize=9)
plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, 'figure2a_cross_task_avg_jaccard.png'), dpi=300, bbox_inches='tight')
fig.savefig(os.path.join(OUT_DIR, 'figure2a_cross_task_avg_jaccard.pdf'), bbox_inches='tight')
plt.close(fig)
print('Saved figure2a_cross_task_avg_jaccard.png/pdf')

# Figure 2b: Lollipop-style dot plot per axis pair (y-axis = axis pairs)
# Four axes: Seq<40, TM<0.4, IDR, New EC -> all six pairwise combinations.
selected_axes = ['Seq<40', 'TM<0.4', 'IDR', 'New EC']
selected_pairs = [
    ('Seq<40', 'TM<0.4'),
    ('Seq<40', 'IDR'),
    ('Seq<40', 'New EC'),
    ('TM<0.4', 'IDR'),
    ('TM<0.4', 'New EC'),
    ('IDR', 'New EC'),
]

pair_records = []
for task_name, path in task_files:
    df = pd.read_csv(os.path.join(ROOT, path), low_memory=False)
    selected = select_representative_axes(df)
    if not all(a in selected for a in selected_axes):
        continue
    for a, b in selected_pairs:
        j = compute_jaccard(df, selected[a], selected[b])
        pair_records.append({
            'task': task_name,
            'pair': f'{a} -- {b}',
            'jaccard': j,
        })

pair_df = pd.DataFrame(pair_records)

# Task ordering by mean Jaccard across selected pairs (low -> high)
task_mean_jaccard = pair_df.groupby('task')['jaccard'].mean().sort_values(ascending=True)
task_order = list(task_mean_jaccard.index)

# Fixed color per task using a high-distinguishability palette
# Chosen to maximize pairwise perceptual distance; includes distinct hues and lightnesses.
task_palette = [
    '#0072B2',  # blue
    '#D55E00',  # vermillion
    '#009E73',  # bluish green
    '#CC79A7',  # reddish purple
    '#F0E442',  # yellow
    '#56B4E9',  # sky blue
    '#E69F00',  # orange
    '#000000',  # black
    '#999999',  # gray
    '#117733',  # green
    '#882255',  # dark pink
    '#44AA99',  # teal
]
task_color = {t: task_palette[i % len(task_palette)] for i, t in enumerate(task_order)}

pair_order = [f'{a} -- {b}' for a, b in selected_pairs]

# Near-square figure, serif font for LaTeX consistency
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['DejaVu Serif'],
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 8,
})
fig, ax = plt.subplots(figsize=(7, 6))

# Draw faint vertical grid and reference line
ax.axvline(x=0.3, color='gray', linestyle='--', linewidth=0.8, alpha=0.6)
ax.text(0.31, -0.35, 'J = 0.3', ha='left', va='center', fontsize=8, color='gray')

# Plot dots per axis pair
for y_idx, pair_str in enumerate(pair_order):
    sub = pair_df[pair_df['pair'] == pair_str].copy()
    sub['task'] = pd.Categorical(sub['task'], categories=task_order, ordered=True)
    sub = sub.sort_values('task')
    # slight vertical jitter for readability
    y_pos = y_idx + np.linspace(-0.18, 0.18, len(sub))
    for (_, row), y in zip(sub.iterrows(), y_pos):
        ax.plot(
            row['jaccard'], y,
            marker='o',
            markersize=7,
            color=task_color[row['task']],
            markeredgecolor='white',
            markeredgewidth=0.6,
        )
    # light horizontal guide line per axis pair
    ax.axhline(y=y_idx, color='lightgray', linestyle='-', linewidth=0.4, alpha=0.5)

ax.set_yticks(range(len(pair_order)))
ax.set_yticklabels(pair_order)
ax.set_xlabel('Jaccard', fontsize=11)
ax.set_ylabel('OOD-axis pair', fontsize=11)
ax.set_title('OOD-axis Pair Overlap Across Tasks', fontsize=12, pad=10)
ax.set_xlim(-0.02, max(pair_df['jaccard'].max() * 1.05, 0.35))
ax.tick_params(axis='y', labelsize=9)
ax.tick_params(axis='x', labelsize=9)
ax.grid(axis='x', linestyle=':', alpha=0.4)
sns.despine(ax=ax, top=True, right=True, left=False, bottom=False)

# Legend: task -> color
handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=task_color[t],
                      markeredgecolor='white', markersize=8, label=t) for t in task_order]
ax.legend(handles=handles, title='Task', fontsize=8, title_fontsize=9,
          loc='upper right', frameon=False, ncol=2)

plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, 'figure2b_pair_overlap_per_task.png'), dpi=300, bbox_inches='tight')
fig.savefig(os.path.join(OUT_DIR, 'figure2b_pair_overlap_per_task.pdf'), bbox_inches='tight')
plt.close(fig)
print('Saved figure2b_pair_overlap_per_task.png/pdf')

# Save numerical data
with open(os.path.join(OUT_DIR, 'cross_task_jaccard_data.json'), 'w') as f:
    json.dump({
        'labels': canonical_axes,
        'mean_matrix': mean_mat.tolist(),
        'pair_records': pair_records,
    }, f, indent=2)
print('Saved cross_task_jaccard_data.json')
