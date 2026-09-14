"""Shared evaluation metrics for the POOR benchmark and protein function prediction."""
from shared.metrics.function_metrics import (
    # Multi-label classification (CAFA-style)
    f1_score_max,
    fmax_score,
    macro_fmax_score,
    micro_aupr,
    macro_aupr,
    select_threshold,
    # Regression
    rmse,
    mae,
    pearson_r,
    spearman_r,
    # Binary classification
    auroc,
    binary_aupr,
    binary_f1,
    # Single-label classification
    accuracy,
    macro_f1,
    # Per-residue prediction
    residue_aupr,
    residue_f1,
    # Label parsing
    parse_labels,
)

__all__ = [
    # Multi-label
    "f1_score_max",
    "fmax_score",
    "macro_fmax_score",
    "micro_aupr",
    "macro_aupr",
    "select_threshold",
    # Regression
    "rmse",
    "mae",
    "pearson_r",
    "spearman_r",
    # Binary
    "auroc",
    "binary_aupr",
    "binary_f1",
    # Single-label
    "accuracy",
    "macro_f1",
    # Per-residue
    "residue_aupr",
    "residue_f1",
    # Parsing
    "parse_labels",
]
