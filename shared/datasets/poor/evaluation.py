"""
Unified OOD evaluation for the POOR (POOD) benchmark.

Supports all task types through automatic dispatch:
  - regression      (kcat, optimal_pH, ligand_binding_affinity)
  - binary          (ppi_prediction)
  - multilabel      (EC, GO-MF/BP/CC)
  - single_label    (fold_classification)
  - residual        (ss_prediction, ppis_prediction, ligand_binding_site)
"""
import json
import os
import re
import numpy as np
import pandas as pd
from typing import Optional, Union, Tuple, List, Dict

from shared.metrics.function_metrics import (
    fmax_score,
    micro_aupr,
    macro_aupr,
    rmse,
    mae,
    pearson_r,
    spearman_r,
    auroc,
    binary_aupr,
    binary_f1,
    accuracy,
    macro_f1,
    residue_aupr,
    residue_f1,
)


# ============================================================
# Task-type detection
# ============================================================

def detect_task_type(
    test_df: pd.DataFrame,
    label_col: str = "label",
    num_classes: Optional[int] = None,
) -> str:
    """
    Infer the task type from the test DataFrame's *label* column.

    Returns one of: ``"regression"``, ``"binary"``, ``"multilabel"``,
    ``"single_label"``, ``"residual"``.
    """
    sample_labels = test_df[label_col].dropna().head(200).tolist()
    if len(sample_labels) == 0:
        return "regression"

    # Residual: labels look like "[0 1 0 ...]" (bracketed arrays)
    first = str(sample_labels[0]).strip()
    if first.startswith("[") and first.endswith("]"):
        return "residual"

    # Multi-label: semicolon-separated tokens (EC, GO)
    if any(";" in str(l) for l in sample_labels):
        return "multilabel"

    # Try float conversion
    float_count = 0
    int_count = 0
    str_count = 0
    unique_vals = set()
    for l in sample_labels:
        try:
            v = float(l)
            float_count += 1
            unique_vals.add(v)
            if v == int(v):
                int_count += 1
        except (ValueError, TypeError):
            str_count += 1

    if float_count == len(sample_labels):
        # Binary: labels are exactly 0/1 values
        if unique_vals <= {0.0, 1.0}:
            return "binary"
        # Single-label: integer class indices with >2 classes (e.g. fold)
        if int_count == len(sample_labels) and len(unique_vals) > 2:
            return "single_label"
        # Any other numeric labels (continuous values or integer targets
        # with few distinct values, e.g. kcat 0.6989, pH 7.0) → regression
        return "regression"

    # Non-numeric class labels (e.g. CATH fold codes "1.10.10")
    return "single_label"


# ============================================================
# Label parsing helpers
# ============================================================

def _parse_residual_labels(label_series: pd.Series) -> Tuple[np.ndarray, np.ndarray]:
    """
    Parse bracketed per-residue label strings into a flat (probs, labels) format.

    Returns:
        flat_labels: 1-D int array (all residues concatenated, -1 = masked)
        boundaries:  list of (start, end) offsets per protein
    """
    flat_labels = []
    boundaries = []
    for val in label_series:
        s = str(val).strip().strip("[]")
        parts = s.split()
        arr = np.array([int(float(p)) if p.strip() not in ("", "nan") else -1 for p in parts])
        start = len(flat_labels)
        flat_labels.extend(arr.tolist())
        end = len(flat_labels)
        boundaries.append((start, end))
    return np.array(flat_labels, dtype=np.int64), boundaries


# ============================================================
# Subset evaluation
# ============================================================

def evaluate_subset(
    predictions: np.ndarray,
    labels: np.ndarray,
    task_type: str,
    threshold: Optional[float] = None,
) -> Dict:
    """
    Evaluate a single subset with metrics appropriate for the task type.

    Args:
        predictions: for multilabel/binary/residual — probabilities [N] or [N, C]
                     for regression — predicted values [N]
                     for single_label — logits/probs [N, C]
        labels: ground truth
        task_type: one of "regression", "binary", "multilabel", "single_label", "residual"
        threshold: optional fixed threshold for Fmax/F1

    Returns:
        dict of metric_name → value, always includes "n_samples"
    """
    n = len(predictions) if not hasattr(predictions, 'shape') else predictions.shape[0]
    result = {"n_samples": n}

    if task_type == "regression":
        preds = np.asarray(predictions, dtype=np.float64).ravel()
        tgts = np.asarray(labels, dtype=np.float64).ravel()
        result["rmse"] = rmse(preds, tgts)
        result["mae"] = mae(preds, tgts)
        result["pearson_r"] = pearson_r(preds, tgts)
        result["spearman_r"] = spearman_r(preds, tgts)

    elif task_type == "binary":
        probs = np.asarray(predictions, dtype=np.float64).ravel()
        labs = np.asarray(labels, dtype=np.float64).ravel()
        result["auroc"] = auroc(probs, labs)
        result["aupr"] = binary_aupr(probs, labs)
        result["f1"] = binary_f1(probs, labs, threshold=0.5)

    elif task_type == "multilabel":
        probs = np.asarray(predictions, dtype=np.float64)
        labs = np.asarray(labels, dtype=np.float64)
        if probs.ndim == 1:
            probs = probs.reshape(-1, 1)
        if labs.ndim == 1:
            labs = labs.reshape(-1, 1)
        if threshold is not None:
            # Fixed threshold (e.g. selected on the validation set): no
            # free sweep on the evaluation subset, avoiding test leakage.
            fmax, best_t = fmax_score(probs, labs, thresholds=np.array([threshold]))
        else:
            fmax, best_t = fmax_score(probs, labs)
        result["fmax"] = fmax
        result["best_threshold"] = best_t
        result["micro_aupr"] = micro_aupr(probs, labs)
        result["macro_aupr"] = macro_aupr(probs, labs)

    elif task_type == "single_label":
        probs_or_logits = np.asarray(predictions, dtype=np.float64)
        labs = np.asarray(labels, dtype=np.int64).ravel()
        result["accuracy"] = accuracy(probs_or_logits, labs)
        result["macro_f1"] = macro_f1(probs_or_logits, labs)

    elif task_type == "residual":
        labs = np.asarray(labels).ravel().astype(np.int64)
        mask = labs >= 0
        unique_labs = np.unique(labs[mask]) if mask.any() else np.array([])
        probs = np.asarray(predictions, dtype=np.float64)
        if unique_labs.size > 2:
            # Per-residue multi-class (e.g. SSP 8-state): argmax accuracy +
            # macro-F1 over valid residues, matching the KPLM
            # residual_classification metric (per-residue argmax accuracy).
            if probs.ndim < 2:
                raise ValueError(
                    "Multi-class residual labels require per-class predictions "
                    f"[N, L, C] or [n_residues, C], got shape {probs.shape}"
                )
            n_classes = probs.shape[-1]
            probs = probs.reshape(-1, n_classes)
            valid_labs = labs[mask]
            result["accuracy"] = accuracy(probs[mask], valid_labs)
            result["macro_f1"] = macro_f1(probs[mask], valid_labs)
        else:
            probs = probs.ravel()
            result["aupr"] = residue_aupr(probs, labs, mask=mask)
            result["f1"] = residue_f1(probs, labs, mask=mask, threshold=0.5)

    else:
        raise ValueError(f"Unknown task_type: {task_type}")

    return result


# ============================================================
# OOD-column evaluation
# ============================================================

def _get_ood_columns(test_df: pd.DataFrame) -> List[str]:
    """Identify all OOD annotation columns in the test DataFrame."""
    ood_cols = []
    for col in test_df.columns:
        if (
            col == "Default"
            or col == "InD"
            or col.startswith("OOD_")
            or col.startswith("seq_Redundancy")
            or col.startswith("TM-score")
        ):
            ood_cols.append(col)
    return ood_cols


def evaluate_by_ood_columns(
    predictions: np.ndarray,
    labels: np.ndarray,
    test_df: pd.DataFrame,
    task_type: Optional[str] = None,
    threshold: Optional[float] = None,
) -> pd.DataFrame:
    """
    Evaluate on *Default=True* (in-distribution) and each OOD column subset.

    Args:
        predictions: model predictions (probs, logits, or regression values)
        labels: ground truth labels (same length as test_df rows)
        test_df: test DataFrame with OOD annotation columns
        task_type: if None, auto-detected from test_df
        threshold: optional fixed threshold for Fmax/F1

    Returns:
        DataFrame with columns: subset, n_samples, <metrics...>
    """
    if task_type is None:
        task_type = detect_task_type(test_df)

    predictions = np.asarray(predictions)
    labels = np.asarray(labels)

    # Align test_df rows with predictions
    assert len(test_df) == len(predictions), (
        f"Mismatch: test_df has {len(test_df)} rows but predictions has {len(predictions)} rows"
    )

    ood_cols = _get_ood_columns(test_df)
    results = []

    for col in ood_cols:
        if col not in test_df.columns:
            continue

        # Handle boolean columns; NaN/missing → False
        mask = test_df[col].fillna(False).astype(bool).values
        if mask.sum() == 0:
            continue

        subset_preds = predictions[mask]
        subset_labels = labels[mask]
        metrics = evaluate_subset(subset_preds, subset_labels, task_type, threshold=threshold)

        row = {"subset": col}
        row.update(metrics)
        results.append(row)

    return pd.DataFrame(results)


def evaluate_model_predictions(
    predictions: np.ndarray,
    labels: np.ndarray,
    test_df: pd.DataFrame,
    task_type: Optional[str] = None,
    threshold: Optional[float] = None,
    output_dir: Optional[str] = None,
    save_predictions: bool = True,
) -> Tuple[pd.DataFrame, str]:
    """
    One-stop evaluation entry point.

    Auto-detects task type, evaluates ID (Default=True) and all OOD subsets,
    prints a summary table, and optionally saves results to disk.

    Args:
        predictions: model predictions
        labels: ground truth
        test_df: test DataFrame with OOD columns
        task_type: if None, auto-detected
        threshold: optional fixed threshold for Fmax/F1
        output_dir: if provided, save results CSV + JSON here
        save_predictions: if True, save raw predictions to CSV

    Returns:
        (results_df, task_type)
    """
    if task_type is None:
        task_type = detect_task_type(test_df)

    print(f"Detected task type: {task_type}")

    results_df = evaluate_by_ood_columns(
        predictions, labels, test_df, task_type=task_type, threshold=threshold
    )

    # Print summary
    print("\n" + "=" * 80)
    print(f"OOD Evaluation Summary (task_type={task_type})")
    print("=" * 80)
    if len(results_df) > 0:
        cols_to_show = ["subset", "n_samples"] + [
            c for c in results_df.columns
            if c not in ("subset", "n_samples")
        ]
        print(results_df[cols_to_show].to_string(index=False, float_format="%.4f"))
    print("=" * 80)

    # Save to disk
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

        results_csv = os.path.join(output_dir, "test_results.csv")
        results_df.to_csv(results_csv, index=False)
        print(f"Results saved to {results_csv}")

        results_json = os.path.join(output_dir, "test_results.json")
        results_json_data = {
            "task_type": task_type,
            "n_test_samples": len(predictions),
            "subsets": results_df.to_dict(orient="records"),
        }
        with open(results_json, "w") as f:
            json.dump(results_json_data, f, indent=2, default=str)
        print(f"JSON saved to {results_json}")

        if save_predictions:
            pred_csv = os.path.join(output_dir, "test_predictions.csv")
            pred_df = pd.DataFrame({"prediction": predictions.ravel()})
            if len(labels) == len(predictions):
                pred_df["label"] = np.asarray(labels).ravel()
            # Include OOD subset columns for traceability
            for col in _get_ood_columns(test_df):
                pred_df[col] = test_df[col].values if col in test_df.columns else None
            pred_df.to_csv(pred_csv, index=False)
            print(f"Predictions saved to {pred_csv}")

    return results_df, task_type
