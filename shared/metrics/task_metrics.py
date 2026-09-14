"""Lightweight task-level metric adapters.

This module provides the small function surface consumed by
``shared.metrics.prodpr`` (binary_auroc / f1_score_max_value / mae),
implemented on top of :mod:`shared.metrics.function_metrics`.
"""
from __future__ import annotations

import numpy as np

from shared.metrics.function_metrics import auroc, fmax_score, mae  # noqa: F401

# Fixed threshold grid for the CAFA-style Fmax sweep. Using a shared grid keeps
# ProDPR strata computations deterministic and fast (identical thresholds for
# model and null scores, subsets and references).
_FMAX_THRESHOLDS = np.linspace(0.01, 0.99, 500)


def f1_score_max_value(scores: np.ndarray, targets: np.ndarray) -> float:
    """Protein-level CAFA Fmax over a fixed threshold grid.

    Returns NaN for empty inputs so callers can treat it as "undefined".
    """
    scores = np.asarray(scores, dtype=np.float64)
    targets = np.asarray(targets, dtype=np.float64)
    if scores.size == 0 or targets.size == 0 or scores.shape[0] == 0:
        return float("nan")
    value, _ = fmax_score(scores, targets, thresholds=_FMAX_THRESHOLDS)
    return float(value)


def binary_auroc(scores: np.ndarray, targets: np.ndarray) -> float:
    """Micro-averaged binary AUROC; NaN when undefined."""
    scores = np.asarray(scores, dtype=np.float64)
    targets = np.asarray(targets, dtype=np.float64)
    if scores.size == 0 or targets.size == 0 or scores.shape[0] == 0:
        return float("nan")
    try:
        return float(auroc(scores, targets))
    except Exception:
        return float("nan")
