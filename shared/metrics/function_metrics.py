"""
Evaluation metrics for the POOR benchmark and protein function prediction.

Covers all POOR benchmark task types:
  - Multi-label classification (EC, GO): Fmax (CAFA protein-centric), micro/macro AUPR
  - Binary classification (PPI): AUROC, AUPR, F1
  - Single-label classification (fold/CATH): Accuracy, macro-F1
  - Regression (kcat, pH, LBA): RMSE, MAE, Pearson r, Spearman r
  - Per-residue prediction (SSP, PPIS, LBS): residue-level AUPR, F1
"""
import numpy as np
import pandas as pd
import torch
from typing import Tuple, Optional, List
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import (
    average_precision_score,
    roc_auc_score,
    f1_score,
    accuracy_score,
)


# ============================================================
# Multi-label classification metrics (CAFA-style)
# ============================================================

def f1_score_max(pred: torch.Tensor, target: torch.Tensor) -> torch.Tensor:
    """
    faster CAFA-style protein-centric Fmax (PyTorch tensor implementation).

    Enumerates all possible thresholds and picks the one with max F1.

    Args:
        pred: predictions of shape (B, N) — tensor
        target: binary targets of shape (B, N) — tensor

    Returns:
        fmax score (scalar tensor)
    """
    order = pred.argsort(descending=True, dim=1)
    target = target.gather(1, order)
    precision = target.cumsum(1) / torch.ones_like(target).cumsum(1)
    recall = target.cumsum(1) / (target.sum(1, keepdim=True) + 1e-10)
    is_start = torch.zeros_like(target).bool()
    is_start[:, 0] = 1
    is_start = torch.scatter(is_start, 1, order, is_start)

    all_order = pred.flatten().argsort(descending=True)
    order = order + torch.arange(order.shape[0], device=order.device).unsqueeze(1) * order.shape[1]
    order = order.flatten()
    inv_order = torch.zeros_like(order)
    inv_order[order] = torch.arange(order.shape[0], device=order.device)
    is_start = is_start.flatten()[all_order]
    all_order = inv_order[all_order]
    precision = precision.flatten()
    recall = recall.flatten()
    all_precision = precision[all_order] - \
                    torch.where(is_start, torch.zeros_like(precision), precision[all_order - 1])
    all_precision = all_precision.cumsum(0) / is_start.cumsum(0)
    all_recall = recall[all_order] - \
                 torch.where(is_start, torch.zeros_like(recall), recall[all_order - 1])
    all_recall = all_recall.cumsum(0) / pred.shape[0]
    all_f1 = 2 * all_precision * all_recall / (all_precision + all_recall + 1e-10)
    return all_f1.max()


def fmax_score(
    probs: np.ndarray,
    labels: np.ndarray,
    thresholds: Optional[np.ndarray] = None,
) -> Tuple[float, float]:
    """
    CAFA-style protein-centric Fmax (numpy implementation).

    For each threshold *t*, computes micro-averaged precision and recall
    across all proteins, then takes the maximum F1 over all thresholds.

    Args:
        probs:   [N, C] predicted probabilities (numpy or torch-compatible)
        labels:  [N, C] multi-hot ground truth
        thresholds: optional 1-D array of thresholds to sweep.
                    If *None*, a set of unique probability values
                    (capped at 500) is used.

    Returns:
        (fmax, best_threshold)  — both floats
    """
    probs = np.asarray(probs, dtype=np.float64)
    labels = np.asarray(labels, dtype=np.float64)

    # Handle 1-D arrays (single binary task)
    if probs.ndim == 1:
        probs = probs.reshape(-1, 1)
    if labels.ndim == 1:
        labels = labels.reshape(-1, 1)

    if thresholds is None:
        thresholds = np.unique(probs)
        if len(thresholds) > 500:
            thresholds = np.linspace(0.01, 0.99, 500)
    elif np.isscalar(thresholds):
        thresholds = np.array([thresholds])

    n_proteins = probs.shape[0]
    thresholds = np.sort(np.asarray(thresholds, dtype=np.float64))

    # Count tp/fp/fn per (protein, threshold) via binary search + prefix
    # sums: each (protein, class) entry falls into one threshold bin, so the
    # cost is O(N*C log T) instead of a Python loop over T full-matrix
    # comparisons. Counts are integers, so results are bit-identical to the
    # explicit loop. preds = probs >= t means an entry with prob p counts
    # toward every threshold th[k] <= p, i.e. bins k < searchsorted(p).
    labels_bool = labels > 0.5
    n_thresholds = len(thresholds)
    kmax = np.searchsorted(thresholds, probs, side="right")  # [N, C] in 0..T
    rows = np.repeat(np.arange(n_proteins), probs.shape[1])
    flat_pos = labels_bool.ravel()
    flat_k = kmax.ravel()

    tp_diff = np.zeros((n_thresholds + 1, n_proteins), dtype=np.int64)
    fp_diff = np.zeros((n_thresholds + 1, n_proteins), dtype=np.int64)
    for diff, mask in ((tp_diff, flat_pos), (fp_diff, ~flat_pos)):
        r = rows[mask]
        k = flat_k[mask]
        np.add.at(diff, (np.zeros_like(k), r), 1)   # counts from bin 0 ...
        np.add.at(diff, (k, r), -1)                  # ... up to bin k-1
    tp = np.cumsum(tp_diff, axis=0)[:n_thresholds].astype(np.float64)  # [T, N]
    fp = np.cumsum(fp_diff, axis=0)[:n_thresholds].astype(np.float64)  # [T, N]
    fn = labels_bool.sum(axis=1).astype(np.float64)[None, :] - tp      # [T, N]

    # CAFA: per-protein precision/recall, then average over proteins.
    # Rows are contiguous [N] vectors, so the row sums reduce in the same
    # order as the original per-threshold loop (bit-identical values).
    denom_p = tp + fp
    precision_i = np.where(denom_p > 0, tp / (denom_p + 1e-10), 0.0)
    denom_r = tp + fn
    recall_i = np.where(denom_r > 0, tp / (denom_r + 1e-10), 0.0)

    precision = precision_i.sum(axis=1) / n_proteins  # [T]
    recall = recall_i.sum(axis=1) / n_proteins        # [T]
    f1 = 2 * precision * recall / (precision + recall + 1e-10)

    # Strict > keeps the earliest threshold on ties, as the loop did.
    best_fmax = 0.0
    best_t = 0.5
    for j in range(n_thresholds):
        if f1[j] > best_fmax:
            best_fmax = float(f1[j])
            best_t = float(thresholds[j])

    return best_fmax, best_t


def macro_fmax_score(
    probs: np.ndarray,
    labels: np.ndarray,
    thresholds: Optional[np.ndarray] = None,
) -> Tuple[float, float]:
    """
    Label-centric macro Fmax.

    For each candidate threshold, compute precision/recall per label,
    then average F1 over labels with at least one positive sample.

    Args:
        probs: [N, C] predicted probabilities
        labels: [N, C] multi-hot ground truth
        thresholds: optional 1-D array of thresholds to sweep

    Returns:
        (macro_fmax, best_threshold)
    """
    probs = np.asarray(probs, dtype=np.float64)
    labels = np.asarray(labels, dtype=np.float64)

    if thresholds is None:
        thresholds = np.unique(probs)
        if len(thresholds) > 500:
            thresholds = np.linspace(0.01, 0.99, 500)

    best_fmax = 0.0
    best_t = 0.5

    active_labels = labels.sum(axis=0) > 0
    if active_labels.sum() == 0:
        return 0.0, 0.5

    probs_active = probs[:, active_labels]
    labels_active = labels[:, active_labels]

    for t in thresholds:
        preds = (probs_active >= t).astype(np.float64)
        tp = (preds * labels_active).sum(axis=0)
        fp = (preds * (1 - labels_active)).sum(axis=0)
        fn = ((1 - preds) * labels_active).sum(axis=0)

        precision = np.where((tp + fp) > 0, tp / (tp + fp + 1e-10), 0.0)
        recall = np.where((tp + fn) > 0, tp / (tp + fn + 1e-10), 0.0)

        f1_per_label = np.where(
            (precision + recall) > 0,
            2 * precision * recall / (precision + recall + 1e-10),
            0.0,
        )
        macro_f1 = f1_per_label.mean()

        if macro_f1 > best_fmax:
            best_fmax = macro_f1
            best_t = t

    return best_fmax, best_t


def micro_aupr(probs: np.ndarray, labels: np.ndarray) -> float:
    """Micro-AUPR: flatten all (protein, label) pairs, compute PR-AUC."""
    probs = np.asarray(probs, dtype=np.float64).ravel()
    labels = np.asarray(labels, dtype=np.float64).ravel()
    if labels.sum() == 0:
        return 0.0
    # Reshape to 2D so sklearn treats this as binary (single-label) classification
    probs = probs.reshape(-1, 1)
    labels = labels.reshape(-1, 1)
    return float(average_precision_score(labels, probs))


def macro_aupr(probs: np.ndarray, labels: np.ndarray) -> float:
    """Macro-AUPR: per-label AUPR, averaged (over labels with >=1 positive)."""
    probs = np.asarray(probs, dtype=np.float64)
    labels = np.asarray(labels, dtype=np.float64)
    if probs.ndim == 1:
        probs = probs.reshape(-1, 1)
        labels = labels.reshape(-1, 1)
    n_labels = labels.shape[1]

    aups = []
    for c in range(n_labels):
        if labels[:, c].sum() == 0:
            continue
        aupr = average_precision_score(labels[:, c], probs[:, c])
        aups.append(aupr)

    if len(aups) == 0:
        return 0.0
    return float(np.mean(aups))


def select_threshold(val_probs: np.ndarray, val_labels: np.ndarray) -> Tuple[float, float]:
    """Select the best Fmax threshold on the validation set."""
    fmax, best_t = fmax_score(val_probs, val_labels)
    return best_t, fmax


# ============================================================
# Regression metrics
# ============================================================

def rmse(preds: np.ndarray, targets: np.ndarray) -> float:
    """Root Mean Squared Error."""
    preds = np.asarray(preds, dtype=np.float64).ravel()
    targets = np.asarray(targets, dtype=np.float64).ravel()
    return float(np.sqrt(np.mean((preds - targets) ** 2)))


def mae(preds: np.ndarray, targets: np.ndarray) -> float:
    """Mean Absolute Error."""
    preds = np.asarray(preds, dtype=np.float64).ravel()
    targets = np.asarray(targets, dtype=np.float64).ravel()
    return float(np.mean(np.abs(preds - targets)))


def pearson_r(preds: np.ndarray, targets: np.ndarray) -> float:
    """Pearson correlation coefficient."""
    preds = np.asarray(preds, dtype=np.float64).ravel()
    targets = np.asarray(targets, dtype=np.float64).ravel()
    if len(preds) < 2:
        return 0.0
    r, _ = pearsonr(preds, targets)
    return float(r) if not np.isnan(r) else 0.0


def spearman_r(preds: np.ndarray, targets: np.ndarray) -> float:
    """Spearman rank correlation coefficient."""
    preds = np.asarray(preds, dtype=np.float64).ravel()
    targets = np.asarray(targets, dtype=np.float64).ravel()
    if len(preds) < 2:
        return 0.0
    r, _ = spearmanr(preds, targets)
    return float(r) if not np.isnan(r) else 0.0


# ============================================================
# Binary classification metrics (e.g. PPI)
# ============================================================

def auroc(probs: np.ndarray, labels: np.ndarray) -> float:
    """Area Under the ROC Curve.

    Returns NaN when the labels contain a single class (AUROC undefined),
    so downstream aggregation can drop the value instead of treating a
    meaningless 0.0 as a valid score.
    """
    probs = np.asarray(probs, dtype=np.float64).ravel()
    labels = np.asarray(labels, dtype=np.float64).ravel().astype(int)
    if len(np.unique(labels)) < 2:
        return float("nan")
    return float(roc_auc_score(labels, probs))


def binary_aupr(probs: np.ndarray, labels: np.ndarray) -> float:
    """Area Under the Precision-Recall Curve for binary classification."""
    probs = np.asarray(probs, dtype=np.float64).ravel()
    labels = np.asarray(labels, dtype=np.float64).ravel()
    if labels.sum() == 0:
        return 0.0
    return float(average_precision_score(labels, probs))


def binary_f1(
    probs: np.ndarray,
    labels: np.ndarray,
    threshold: float = 0.5,
) -> float:
    """Binary F1 at a fixed threshold."""
    preds = (np.asarray(probs, dtype=np.float64).ravel() >= threshold).astype(int)
    labels = np.asarray(labels, dtype=np.float64).ravel().astype(int)
    return float(f1_score(labels, preds, zero_division=0))


# ============================================================
# Single-label classification metrics (e.g. fold/CATH)
# ============================================================

def accuracy(probs_or_logits: np.ndarray, labels: np.ndarray) -> float:
    """
    Top-1 accuracy for single-label classification.

    Args:
        probs_or_logits: [N, C] logits or probabilities
        labels: [N] integer class indices
    """
    arr = np.asarray(probs_or_logits, dtype=np.float64)
    if arr.ndim == 1:
        preds = (arr >= 0.5).astype(int)
    else:
        preds = arr.argmax(axis=1)
    labels = np.asarray(labels, dtype=np.int64).ravel()
    return float(accuracy_score(labels, preds))


def macro_f1(
    probs_or_logits: np.ndarray,
    labels: np.ndarray,
    num_classes: Optional[int] = None,
) -> float:
    """
    Macro-averaged F1 for single-label classification.

    Args:
        probs_or_logits: [N, C] logits or probabilities
        labels: [N] integer class indices
        num_classes: optional, for explicit class count
    """
    arr = np.asarray(probs_or_logits, dtype=np.float64)
    if arr.ndim == 1:
        preds = (arr >= 0.5).astype(int)
    else:
        preds = arr.argmax(axis=1)
    labels = np.asarray(labels, dtype=np.int64).ravel()
    return float(f1_score(labels, preds, average="macro", zero_division=0))


# ============================================================
# Per-residue prediction metrics (e.g. SSP, PPIS, LBS)
# ============================================================

def residue_aupr(
    probs: np.ndarray,
    labels: np.ndarray,
    mask: Optional[np.ndarray] = None,
) -> float:
    """
    Per-residue AUPR.

    Args:
        probs:   [N] or [N, C] predicted probabilities for residue-level binary tasks
        labels:  [N] or [N, C] ground-truth labels (0/1, -1 for masked)
        mask:    optional [N] boolean mask; if None, positions with label == -1 are masked out
    """
    probs = np.asarray(probs, dtype=np.float64).ravel()
    labels = np.asarray(labels, dtype=np.float64).ravel()

    if mask is None:
        mask = labels >= 0
    else:
        mask = np.asarray(mask, dtype=bool).ravel()

    probs = probs[mask]
    labels = labels[mask]
    if labels.sum() == 0:
        return 0.0
    return float(average_precision_score(labels, probs))


def residue_f1(
    probs: np.ndarray,
    labels: np.ndarray,
    mask: Optional[np.ndarray] = None,
    threshold: float = 0.5,
) -> float:
    """
    Per-residue F1 at a fixed threshold.

    Args:
        probs:   [N] predicted probabilities
        labels:  [N] ground-truth labels (0/1, -1 for masked)
        mask:    optional boolean mask
        threshold: probability threshold for positive prediction
    """
    probs = np.asarray(probs, dtype=np.float64).ravel()
    labels = np.asarray(labels, dtype=np.float64).ravel()

    if mask is None:
        mask = labels >= 0
    else:
        mask = np.asarray(mask, dtype=bool).ravel()

    preds = (probs[mask] >= threshold).astype(int)
    gt = labels[mask].astype(int)
    return float(f1_score(gt, preds, zero_division=0))


# ============================================================
# Label parsing utility
# ============================================================

def parse_labels(label_str: str, label_vocab: List[str]) -> Tuple[np.ndarray, List[str]]:
    """
    Parse a semicolon-separated label string into a multi-hot vector.

    Args:
        label_str: e.g. "3.2.1.-;2.7.7.-" or "GO:0006412;GO:0005524"
        label_vocab: list of labels (defines the index mapping)

    Returns:
        (multi_hot_vector, unseen_labels)
    """
    label_to_idx = {l: i for i, l in enumerate(label_vocab)}
    multi_hot = np.zeros(len(label_vocab), dtype=np.float32)
    unseen = []

    if pd.isna(label_str) or label_str == "":
        return multi_hot, unseen

    for label in label_str.split(";"):
        label = label.strip()
        if label in label_to_idx:
            multi_hot[label_to_idx[label]] = 1.0
        else:
            unseen.append(label)

    return multi_hot, unseen
