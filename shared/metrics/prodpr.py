#!/usr/bin/env python3
"""Compute ProDPR from per-sample predictions and OOD annotations.

This script implements the ProDPR metric described in the paper.
ProDPR (Protein Distributional Performance Retention) is a finite-sample
estimator of an OOD-space retention functional. It evaluates the expected
fraction of effective predictive utility retained across a high-dimensional
protein OOD descriptor space.

Each protein is mapped to a task-specific OOD descriptor that combines:
  - Continuous novelty coordinates: sequence novelty, structural novelty,
    length novelty, and label-cardinality novelty.
  - Discrete biological scenario annotations: e.g., IDR, orphan, long-tail,
    new-function, combinatorial-function.

The OOD descriptor space is partitioned into biologically interpretable
strata (joint sequence--structure novelty regions, biological scenario
patterns, or interactions between continuous novelty and discrete scenarios).
Only strata with sufficient support are retained.

For each valid stratum, the retained utility is computed as:

  r_t(f, ell) = (M(f, E_{t,ell}) - M(f_null, E_{t,ell})) /
                (M(f, B_{t,ell}) - M(f_null, B_{t,ell}) + eps)

where E_{t,ell} is the evaluated stratum and B_{t,ell} is a matched
counterfactual reference subset that balances confounding factors
(sequence/structure novelty bins, length bins, and label-cardinality bins).

The empirical ProDPR is the weighted average over all valid strata:

  ProDPR_t(f) = sum_{ell in Omega_t} mu_t(ell) * r_t(f, ell)

with macro weights mu_t(ell) = 1 / |Omega_t| by default. The global ProDPR
across a benchmark is the task-level macro average.

The implementation is intentionally sample-level. Aggregated summary CSVs are
not enough to compute strict ProDPR because DPR needs a task-specific null
baseline and matched controls for each stratum.

Typical GO/EC multi-label use:

  python KPLM/tasks/scripts/calculate_prodpr.py \
    --metadata-csv data/datasets_split/go/go_cc_test.csv \
    --train-csv data/datasets_split/go/go_cc_train.csv \
    --pred-path outputs/go_cc_logits.npy \
    --output-prefix data/result_summary/saprot_go_cc_prodpr

If --true-path is not provided, y_true is constructed from metadata-csv labels
using the label vocabulary inferred from --train-csv.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
import zlib
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from shared.metrics.task_metrics import binary_auroc, f1_score_max_value, mae

# Default substring that auto-discovered scenario axes must contain: the
# new-function family (OOD_NewEC_L3, OOD_NewEC_L4, ...). Pass an empty
# filter to keep every OOD-prefixed column instead.
DEFAULT_SCENARIO_FILTER = "NewEC"


@dataclass
class DprResult:
    dpr: Optional[float]
    model_subset: Optional[float]
    null_subset: Optional[float]
    model_reference: Optional[float]
    null_reference: Optional[float]
    subset_size: int
    reference_size: int


@dataclass
class ProDprSummary:
    prodpr: Optional[float]
    worst_stratum: Optional[float]
    n_min: int
    weight_mode: str
    metric: str
    null_baseline: str
    scenario_cols: List[str]
    strata_count: int
    reference_size: int
    output_prefix: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute ProDPR from per-sample predictions, labels, and OOD metadata."
    )
    parser.add_argument("--metadata-csv", required=True, help="Test metadata CSV aligned with predictions.")
    parser.add_argument("--pred-path", required=True, help="Prediction score/logit matrix: .npy, .npz, .csv, or .tsv.")
    parser.add_argument("--true-path", default=None, help="Optional true-label matrix. If omitted, use metadata labels.")
    parser.add_argument("--train-csv", default=None, help="Training CSV used for label vocabulary and label-prior null baseline. Required when --null-baseline=label_prior.")
    parser.add_argument("--array-key", default=None, help="Array key for .npz inputs.")
    parser.add_argument("--output-prefix", required=True, help="Prefix for output JSON/CSV files.")

    parser.add_argument(
        "--metric",
        default="f1_max",
        choices=["f1_max", "accuracy", "macro_f1", "residue_accuracy", "residue_auroc", "auroc", "mae", "spearman"],
        help="Task utility metric.",
    )
    parser.add_argument("--null-baseline", default="label_prior", choices=["label_prior", "zero"], help="Null predictor.")
    parser.add_argument("--epsilon", type=float, default=1e-8, help="Numerical stability constant in DPR.")
    parser.add_argument("--n-min", type=int, default=20, help="Minimum support for valid strata.")
    parser.add_argument("--seed", type=int, default=0, help="Random seed for matched-control sampling.")
    parser.add_argument("--weight-mode", default="macro", choices=["macro", "empirical"], help="Group weighting mode (macro used by default).")

    parser.add_argument("--label-col", default="labels", help="Metadata label column for building y_true.")
    parser.add_argument(
        "--label-vocab-path",
        default=None,
        help=(
            "Optional label vocabulary file (one 'index<TAB>label' line per "
            "entry, the same file the training dataset uses, e.g. "
            "ec_vocab.txt). When given, y_true columns follow this file's "
            "index order instead of an alphabetically sorted vocabulary "
            "re-inferred from the train csv."
        ),
    )
    parser.add_argument("--seq-prefix", default="seq_Redundancy_", help="Prefix for sequence novelty cutoff flags.")
    parser.add_argument("--str-prefix", default="TM-score_", help="Prefix for structure novelty cutoff flags.")
    parser.add_argument("--scenario-prefix", default="OOD_", help="Prefix for task-specific OOD scenario flags.")
    parser.add_argument("--seq-cols", default=None, help="Comma-separated override for sequence novelty columns.")
    parser.add_argument("--str-cols", default=None, help="Comma-separated override for structure novelty columns.")
    parser.add_argument(
        "--scenario-cols",
        default=None,
        help=(
            "Comma-separated override for scenario columns. If omitted, "
            "columns are auto-discovered via --scenario-prefix and "
            "--scenario-filter."
        ),
    )
    parser.add_argument(
        "--scenario-filter",
        default=DEFAULT_SCENARIO_FILTER,
        help=(
            "Substring that auto-discovered scenario columns must contain. "
            "Default 'NewEC' selects new-function axes (OOD_NewEC_L3, "
            "OOD_NewEC_L4, ...); pass an empty string to keep every "
            "prefixed column."
        ),
    )
    parser.add_argument(
        "--flag-direction",
        default="le",
        choices=["le", "ge"],
        help="Meaning of cutoff flags. le: True means similarity <= cutoff; ge: True means similarity >= cutoff.",
    )
    parser.add_argument(
        "--reference-pool",
        default="auto",
        choices=["auto", "default", "ind"],
        help=(
            "Reference pool definition. 'auto' (default): no-scenario on the "
            "selected axes, preferring the highest similarity bins, "
            "restricted to Default=True rows; supports --control matched "
            "sampling. 'default': all rows flagged by the --default-col "
            "column, always used whole. 'ind': all InD=True rows, always "
            "used whole."
        ),
    )
    parser.add_argument(
        "--control",
        default="matched",
        choices=["matched", "pool"],
        help=(
            "Control construction for each stratum. 'matched' (default): "
            "sample pool controls matching the stratum's background "
            "dimensions (unselected axes); only honoured with "
            "--reference-pool auto. The 'default'/'ind' pools "
            "are canonical fixed reference sets and always use the whole "
            "pool. 'pool': always use the full reference pool."
        ),
    )
    parser.add_argument(
        "--control-draws",
        type=int,
        default=50,
        help=(
            "Number of matched-control draws averaged per stratum when "
            "--control matched is active with background dimensions. "
            "Averaging the reference metric over K draws reduces Monte "
            "Carlo noise by ~1/sqrt(K). Ignored for --control pool."
        ),
    )
    parser.add_argument(
        "--default-col",
        default="Default",
        help=(
            "Boolean metadata column marking rows with fully computed OOD "
            "annotations. When the column exists, the reference pool "
            "(denominator) is restricted to rows where it is True; strata "
            "(numerator) are unaffected. NaN counts as False. Pass an empty "
            "string to disable the restriction."
        ),
    )
    return parser.parse_args()


def warn(msg: str) -> None:
    print(f"[prodpr] WARNING: {msg}", file=sys.stderr)


def split_label_values(value: object) -> List[str]:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return []
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return []
    return [item.strip() for item in re.split(r"[;,]", text) if item.strip()]


def infer_label_mapping(*csv_paths: Optional[str], label_col: str = "labels") -> Dict[str, int]:
    labels: List[str] = []
    for csv_path in csv_paths:
        if not csv_path:
            continue
        df = pd.read_csv(csv_path, usecols=lambda col: col == label_col)
        if label_col not in df.columns:
            continue
        for value in df[label_col]:
            labels.extend(split_label_values(value))
    unique = sorted(set(labels))
    return {label: idx for idx, label in enumerate(unique)}


def load_vocab_file(path: str) -> Dict[str, int]:
    """Load a dataset label vocabulary file (``index<TAB>label`` per line,
    e.g. ec_vocab.txt), returning {label: index} in file index order."""
    entries: Dict[int, str] = {}
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 2:
                raise ValueError(f"Malformed line in label vocab {path!r}: {line!r}")
            idx, label = int(parts[0]), parts[1]
            entries[idx] = label
    return {label: idx for idx, label in sorted(entries.items())}


def resolve_label_mapping(args: argparse.Namespace) -> Dict[str, int]:
    """Label vocabulary for y_true construction.

    Prefers the dataset vocab file (--label-vocab-path, the same file the
    training pipeline uses) over re-inferring from the train csv. The
    re-inferred vocabulary sorts labels alphabetically, which silently
    permutes y_true columns whenever the dataset vocab uses a different
    order — with identical column counts no shape check can catch that.
    """
    vocab_path = getattr(args, "label_vocab_path", None)
    if vocab_path:
        mapping = load_vocab_file(vocab_path)
        if args.train_csv:
            inferred = infer_label_mapping(args.train_csv, label_col=args.label_col)
            if inferred and set(inferred) != set(mapping):
                only_file = set(mapping) - set(inferred)
                only_train = set(inferred) - set(mapping)
                warn(
                    f"Label vocab file and train csv label sets differ: "
                    f"{len(only_file)} labels only in the vocab file, "
                    f"{len(only_train)} only in train. Train-only labels are "
                    f"dropped from y_true; vocab-file order is used."
                )
        return mapping
    return infer_label_mapping(args.train_csv, label_col=args.label_col)


def labels_to_matrix(values: Sequence[object], label_to_idx: Dict[str, int]) -> np.ndarray:
    y = np.zeros((len(values), len(label_to_idx)), dtype=np.float32)
    for row_idx, value in enumerate(values):
        for label in split_label_values(value):
            col_idx = label_to_idx.get(label)
            if col_idx is not None:
                y[row_idx, col_idx] = 1.0
    return y


def labels_to_vector(values: Sequence[object], label_to_idx: Dict[str, int]) -> np.ndarray:
    y = np.full((len(values),), -1, dtype=np.int64)
    for row_idx, value in enumerate(values):
        label = str(value).strip()
        if label in label_to_idx:
            y[row_idx] = label_to_idx[label]
    if (y < 0).any():
        missing = int((y < 0).sum())
        warn(f"{missing} labels were not found in the training label vocabulary; they are set to class 0.")
        y[y < 0] = 0
    return y


def load_numeric_matrix(path: str, array_key: Optional[str] = None) -> np.ndarray:
    source = Path(path)
    suffix = source.suffix.lower()
    if suffix == ".npy":
        arr = np.load(source)
    elif suffix == ".npz":
        data = np.load(source)
        key = array_key or next(iter(data.files))
        arr = data[key]
    elif suffix in {".csv", ".tsv"}:
        sep = "\t" if suffix == ".tsv" else ","
        df = pd.read_csv(source, sep=sep)
        numeric_df = df.select_dtypes(include=[np.number])
        if numeric_df.empty:
            raise ValueError(f"No numeric columns found in {path}")
        arr = numeric_df.to_numpy()
    else:
        raise ValueError(f"Unsupported matrix file extension: {source.suffix}")
    arr = np.asarray(arr, dtype=np.float32)
    if arr.ndim == 1:
        arr = arr.reshape(-1, 1)
    if arr.ndim != 2:
        raise ValueError(f"Expected a 2D matrix in {path}, got shape {arr.shape}")
    return arr


def build_y_true(args: argparse.Namespace, metadata: pd.DataFrame, y_pred: np.ndarray) -> Tuple[np.ndarray, Dict[str, int]]:
    if args.true_path:
        y_true = load_numeric_matrix(args.true_path, args.array_key)
        label_to_idx = resolve_label_mapping(args)
        if args.metric in {"mae", "spearman"}:
            label_to_idx = {}
        elif args.metric in {"accuracy", "macro_f1", "auroc"}:
            y_true = y_true.reshape(-1).astype(np.int64)
            if args.metric == "auroc" and y_pred.ndim == 2 and y_pred.shape[1] not in (1, 2):
                raise ValueError(
                    f"auroc expects binary scores ([N] / [N,1] / [N,2]), got {y_pred.shape}"
                )
            if args.metric != "auroc" and label_to_idx and len(label_to_idx) != y_pred.shape[1]:
                warn(
                    "Train label vocabulary size does not match prediction columns; "
                    "the label-prior null baseline cannot be built from the train "
                    "split and build_null_scores will raise."
                )
                label_to_idx = {}
        elif label_to_idx and len(label_to_idx) != y_true.shape[1]:
            warn(
                "Train label vocabulary size does not match --true-path columns; "
                "the label-prior null baseline cannot be built from the train "
                "split and build_null_scores will raise."
            )
            label_to_idx = {}
    else:
        if args.label_col not in metadata.columns:
            raise ValueError(f"--true-path not provided and metadata lacks label column {args.label_col!r}")
        label_to_idx = resolve_label_mapping(args)
        if not label_to_idx:
            label_to_idx = infer_label_mapping(args.metadata_csv, label_col=args.label_col)
        if not label_to_idx:
            raise ValueError("Could not infer any labels. Provide --true-path or a valid --train-csv/--label-col.")
        if args.metric in {"accuracy", "macro_f1", "auroc"}:
            y_true = labels_to_vector(metadata[args.label_col].tolist(), label_to_idx)
        elif args.metric in {"mae", "spearman"}:
            y_true = pd.to_numeric(
                metadata[args.label_col], errors="raise"
            ).to_numpy(dtype=np.float32).reshape(-1, 1)
            label_to_idx = {}
        else:
            y_true = labels_to_matrix(metadata[args.label_col].tolist(), label_to_idx)

    if y_true.shape[0] != y_pred.shape[0]:
        raise ValueError(f"Row mismatch: y_true has {y_true.shape[0]} rows, y_pred has {y_pred.shape[0]}")
    if args.metric == "f1_max" and y_true.shape[1] != y_pred.shape[1]:
        raise ValueError(f"Class mismatch: y_true has {y_true.shape[1]} columns, y_pred has {y_pred.shape[1]}")
    if args.metric in {"accuracy", "macro_f1"} and y_pred.ndim != 2:
        raise ValueError(f"Expected 2D class-score matrix for {args.metric}, got {y_pred.shape}")
    if args.metric == "auroc" and y_pred.ndim == 2 and y_pred.shape[1] not in (1, 2):
        raise ValueError(f"auroc expects binary scores ([N] / [N,1] / [N,2]), got {y_pred.shape}")
    if args.metric in {"accuracy", "macro_f1", "auroc"}:
        return y_true.reshape(-1).astype(np.int64), label_to_idx
    return y_true.astype(np.float32), label_to_idx


def require_train_csv(args: argparse.Namespace) -> None:
    """label_prior null 的先验必须来自训练集（参考集）。

    拒绝静默回退到测试集标签频率：那会让 null 分子分母都混入测试标签，
    属于评估泄漏。要么传训练集，要么显式改用 zero null。
    """
    if args.null_baseline == "label_prior" and not args.train_csv and getattr(args, "metric", None) != "spearman":
        raise ValueError(
            "null_baseline='label_prior' requires the training split (--train-csv) "
            "to estimate the label prior; pass --train-csv or switch to "
            "--null-baseline zero"
        )


def build_null_scores(
    args: argparse.Namespace,
    metadata: pd.DataFrame,
    y_true: np.ndarray,
    label_to_idx: Dict[str, int],
) -> np.ndarray:
    require_train_csv(args)
    if args.metric == "spearman":
        # Null rho = 0 is data-free; compute_dpr substitutes 0.0 directly and
        # never reads these scores. A zero placeholder keeps the shape valid.
        return np.zeros((y_true.shape[0], 1), dtype=np.float32)

    if args.metric in {"accuracy", "macro_f1"}:
        num_classes = len(label_to_idx) if label_to_idx else int(np.max(y_true)) + 1
        if args.null_baseline == "zero":
            prior = np.zeros(num_classes, dtype=np.float32)
        elif args.train_csv and label_to_idx:
            train_df = pd.read_csv(args.train_csv)
            train_y = labels_to_vector(train_df[args.label_col].tolist(), label_to_idx)
            counts = np.bincount(train_y.astype(np.int64), minlength=num_classes).astype(np.float32)
            prior = counts / max(float(counts.sum()), 1.0)
        else:
            raise ValueError(
                f"Cannot build label-prior null from train csv {args.train_csv!r}: "
                "no usable label vocabulary (train/eval vocabulary mismatch)"
            )
        return np.tile(prior.reshape(1, -1), (y_true.shape[0], 1))

    if args.metric == "mae":
        if args.null_baseline == "zero":
            prior = np.zeros((1,), dtype=np.float32)
        elif args.train_csv:
            train_df = pd.read_csv(args.train_csv)
            train_values = pd.to_numeric(
                train_df[args.label_col], errors="coerce"
            ).dropna()
            if train_values.empty:
                raise ValueError(
                    f"Cannot build label-prior null from train csv {args.train_csv!r}: "
                    f"no numeric values in column {args.label_col!r}"
                )
            prior = np.asarray([float(train_values.median())], dtype=np.float32)
        else:
            raise ValueError(
                "null_baseline='label_prior' requires the training split (--train-csv)"
            )
        return np.tile(prior.reshape(1, -1), (y_true.shape[0], 1))

    if args.metric == "auroc":
        # Binary tasks: the null is a constant positive-class probability
        # (train positive rate), giving AUROC = 0.5 by construction.
        if args.null_baseline == "zero":
            prior = np.zeros((1,), dtype=np.float32)
        elif args.train_csv and label_to_idx:
            train_df = pd.read_csv(args.train_csv)
            train_y = labels_to_vector(train_df[args.label_col].tolist(), label_to_idx)
            prior = np.asarray([float(np.mean(train_y == 1))], dtype=np.float32)
        else:
            raise ValueError(
                f"Cannot build label-prior null from train csv {args.train_csv!r}: "
                "no usable label vocabulary (train/eval vocabulary mismatch)"
            )
        return np.tile(prior.reshape(1, -1), (y_true.shape[0], 1))

    if args.null_baseline == "zero":
        prior = np.zeros(y_true.shape[1], dtype=np.float32)
    elif args.train_csv and label_to_idx:
        train_df = pd.read_csv(args.train_csv)
        train_y = labels_to_matrix(train_df[args.label_col].tolist(), label_to_idx)
        prior = train_y.mean(axis=0).astype(np.float32)
    else:
        raise ValueError(
            f"Cannot build label-prior null from train csv {args.train_csv!r}: "
            "no usable label vocabulary (train/eval vocabulary mismatch)"
        )
    return np.tile(prior.reshape(1, -1), (y_true.shape[0], 1))


def accuracy_score(scores: np.ndarray, targets: np.ndarray) -> Optional[float]:
    if scores.shape[0] == 0:
        return None
    preds = np.argmax(scores, axis=-1).reshape(-1)
    truth = np.asarray(targets).reshape(-1).astype(np.int64)
    return float((preds == truth).mean()) if truth.size else None


def macro_f1_score(scores: np.ndarray, targets: np.ndarray) -> Optional[float]:
    if scores.shape[0] == 0:
        return None
    preds = np.argmax(scores, axis=-1).reshape(-1).astype(np.int64)
    truth = np.asarray(targets).reshape(-1).astype(np.int64)
    # Fixed class set from the score matrix width: averaging over the same
    # classes in every stratum keeps macro-F1 comparable across strata.
    num_classes = scores.shape[-1] if scores.ndim > 1 else int(np.max(truth)) + 1
    f1s = []
    for cls in range(num_classes):
        tp = np.logical_and(preds == cls, truth == cls).sum()
        fp = np.logical_and(preds == cls, truth != cls).sum()
        fn = np.logical_and(preds != cls, truth == cls).sum()
        denom = 2 * tp + fp + fn
        f1s.append(0.0 if denom == 0 else (2 * tp) / denom)
    return float(np.mean(f1s)) if f1s else None


def auroc_score(scores: np.ndarray, targets: np.ndarray) -> Optional[float]:
    if scores.shape[0] == 0:
        return None
    scores = np.asarray(scores)
    # Binary AUROC over one score per sample: use the positive-class column
    # when given a 2-column score matrix.
    if scores.ndim == 2 and scores.shape[1] > 1:
        if scores.shape[1] != 2:
            raise ValueError(f"auroc expects binary scores, got shape {scores.shape}")
        scores = scores[:, 1]
    value = binary_auroc(scores.reshape(-1), targets)
    return None if np.isnan(value) else value


def residue_indices_for_proteins(lengths: np.ndarray, protein_indices: Sequence[int]) -> np.ndarray:
    starts = np.concatenate([[0], np.cumsum(lengths[:-1])])
    residue_indices = []
    for protein_idx in np.asarray(protein_indices, dtype=np.int64):
        start = int(starts[protein_idx])
        stop = start + int(lengths[protein_idx])
        residue_indices.extend(range(start, stop))
    return np.asarray(residue_indices, dtype=np.int64)


def residue_accuracy_score(scores: object, targets: object, indices: Sequence[int]) -> Optional[float]:
    if not isinstance(scores, dict):
        raise ValueError("residue_accuracy requires an .npz pred-path with logits, labels, and lengths")
    lengths = np.asarray(scores["lengths"], dtype=np.int64)
    residue_idx = residue_indices_for_proteins(lengths, indices)
    if residue_idx.size == 0:
        return None
    labels = np.asarray(scores["labels"], dtype=np.int64)[residue_idx]
    if "prior" in scores:
        preds = np.full(labels.shape, int(np.argmax(scores["prior"])), dtype=np.int64)
    else:
        logits = np.asarray(scores["logits"], dtype=np.float32)[residue_idx]
        preds = np.argmax(logits, axis=-1).astype(np.int64)
    return float((preds == labels).mean()) if labels.size else None


def residue_auroc_score(scores: object, targets: object, indices: Sequence[int]) -> Optional[float]:
    if not isinstance(scores, dict):
        raise ValueError("residue_auroc requires an .npz pred-path with logits, labels, and lengths")
    lengths = np.asarray(scores["lengths"], dtype=np.int64)
    residue_idx = residue_indices_for_proteins(lengths, indices)
    if residue_idx.size == 0:
        return None
    labels = np.asarray(scores["labels"], dtype=np.int64)[residue_idx]
    if "prior" in scores:
        prior = np.asarray(scores["prior"], dtype=np.float32).reshape(1, -1)
        logits = np.repeat(prior, labels.size, axis=0)
    else:
        logits = np.asarray(scores["logits"], dtype=np.float32)[residue_idx]
    # AUROC is a binary rank statistic over one score per residue: use the
    # positive-class column. Passing the full [M, C] matrix would ravel to
    # M*C scores against M labels and fail inside roc_auc_score.
    if logits.ndim == 2 and logits.shape[1] > 1:
        if logits.shape[1] != 2:
            raise ValueError(
                f"residue_auroc expects binary residue logits, got {logits.shape[1]} classes"
            )
        scores_1d = logits[:, 1]
    else:
        scores_1d = logits.reshape(-1)
    value = binary_auroc(scores_1d, labels)
    return None if np.isnan(value) else value


def spearman_value(scores: np.ndarray, targets: np.ndarray) -> Optional[float]:
    x = np.asarray(scores, dtype=float).reshape(-1)
    y = np.asarray(targets, dtype=float).reshape(-1)
    if x.shape != y.shape or not (np.isfinite(x).all() and np.isfinite(y).all()):
        raise ValueError("Spearman requires aligned finite predictions and labels")
    if len(x) < 2 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return None
    return float(np.corrcoef(pd.Series(x).rank(method="average"), pd.Series(y).rank(method="average"))[0, 1])


def metric_value(metric: str, scores: object, targets: object, indices: Sequence[int]) -> Optional[float]:
    idx = np.asarray(indices, dtype=np.int64)
    if idx.size == 0:
        return None
    if metric == "residue_accuracy":
        return residue_accuracy_score(scores, targets, idx)
    if metric == "residue_auroc":
        return residue_auroc_score(scores, targets, idx)
    score_array = np.asarray(scores)
    target_array = np.asarray(targets)
    if metric == "f1_max":
        value = f1_score_max_value(score_array[idx], target_array[idx])
        return None if np.isnan(value) else value
    if metric == "accuracy":
        return accuracy_score(score_array[idx], target_array[idx])
    if metric == "macro_f1":
        return macro_f1_score(score_array[idx], target_array[idx])
    if metric == "auroc":
        return auroc_score(score_array[idx], target_array[idx])
    if metric == "mae":
        value = mae(score_array[idx], target_array[idx])
        return None if np.isnan(value) else value
    if metric == "spearman":
        return spearman_value(score_array[idx], target_array[idx])
    raise ValueError(f"Unsupported metric: {metric}")


def _dpr_value(
    metric: str,
    model_subset: Optional[float],
    null_subset: Optional[float],
    model_reference: Optional[float],
    null_reference: Optional[float],
    epsilon: float,
) -> Optional[float]:
    if None in (model_subset, null_subset, model_reference, null_reference):
        return None
    direction = -1.0 if metric == "mae" else 1.0
    ref_gap = direction * (
        float(model_reference) - float(null_reference)
    )
    # Guard against near-zero / negative denominators: when the model
    # does not beat the null baseline on the ID reference set, DPR is
    # undefined (a ~epsilon denominator would explode to ±1e8 and poison
    # the macro average). Scale the threshold relative to the metric's
    # magnitude and treat negative gaps (model worse than null) the same.
    ref_scale = max(
        abs(float(model_reference)), abs(float(null_reference)), 1e-8
    )
    if ref_gap <= 0.01 * ref_scale:
        return None
    numerator = direction * (
        float(model_subset) - float(null_subset)
    )
    return numerator / (ref_gap + epsilon)


def _mean_metric(values: List[Optional[float]]) -> Optional[float]:
    valid = [float(v) for v in values if v is not None and not math.isnan(float(v))]
    if not valid:
        return None
    return sum(valid) / len(valid)


def compute_dpr_multi(
    metric: str,
    y_pred: np.ndarray,
    y_true: np.ndarray,
    y_null: np.ndarray,
    subset_indices: Sequence[int],
    reference_draws: Sequence[Sequence[int]],
    epsilon: float,
) -> DprResult:
    """DPR with the reference metric averaged over multiple control draws.

    Each draw is one matched-control sample; averaging the reference metric
    over draws approximates its expectation and removes single-draw Monte
    Carlo noise. A single deterministic draw (e.g. the full pool) reduces
    to the plain DPR.
    """
    subset_indices = np.asarray(subset_indices, dtype=np.int64)
    draws = [np.asarray(d, dtype=np.int64) for d in reference_draws if len(d) > 0]
    model_subset = metric_value(metric, y_pred, y_true, subset_indices)
    if metric == "spearman":
        # A non-informative predictor has zero rank correlation by definition.
        null_subset: Optional[float] = 0.0
        null_reference: Optional[float] = 0.0
    else:
        null_subset = metric_value(metric, y_null, y_true, subset_indices)
        null_reference = _mean_metric(
            [metric_value(metric, y_null, y_true, d) for d in draws]
        )
    model_reference = _mean_metric(
        [metric_value(metric, y_pred, y_true, d) for d in draws]
    )
    dpr = _dpr_value(
        metric, model_subset, null_subset, model_reference, null_reference, epsilon
    )
    return DprResult(
        dpr=dpr,
        model_subset=model_subset,
        null_subset=null_subset,
        model_reference=model_reference,
        null_reference=null_reference,
        subset_size=int(subset_indices.size),
        reference_size=int(draws[0].size) if draws else 0,
    )


def parse_column_override(value: Optional[str]) -> Optional[List[str]]:
    if value is None:
        return None
    return [item.strip() for item in value.split(",") if item.strip()]


def infer_columns(df: pd.DataFrame, prefix: str, override: Optional[str]) -> List[str]:
    columns = parse_column_override(override)
    if columns is not None:
        missing = [col for col in columns if col not in df.columns]
        if missing:
            raise ValueError(f"Missing requested columns: {missing}")
        return columns
    return [col for col in df.columns if col.startswith(prefix)]


def infer_scenario_columns(
    df: pd.DataFrame,
    override: Optional[str],
    prefix: str = "OOD_",
    name_filter: str = DEFAULT_SCENARIO_FILTER,
) -> List[str]:
    columns = parse_column_override(override)
    if columns is not None:
        missing = [col for col in columns if col not in df.columns]
        if missing:
            raise ValueError(f"Missing requested scenario columns: {missing}")
        return columns

    selected = [
        col for col in df.columns
        if col.startswith(prefix) and (not name_filter or name_filter in col)
    ]
    if not selected:
        warn(
            f"No scenario columns with prefix {prefix!r}"
            + (f" containing {name_filter!r}" if name_filter else "")
            + " found; strata degrade to novelty bins only."
        )
    return selected


def as_bool_array(series: pd.Series) -> np.ndarray:
    if pd.api.types.is_bool_dtype(series):
        return series.to_numpy(dtype=bool)
    if pd.api.types.is_numeric_dtype(series):
        # Float/int 0-1 flags (NaN counts as False). A float column appears
        # whenever a boolean column with NaNs is written to csv and re-read.
        return series.fillna(0).to_numpy(dtype=float) > 0.5
    text = series.astype(str).str.strip().str.lower()
    return text.isin({"1", "1.0", "true", "t", "yes", "y"}).to_numpy(dtype=bool)


def cutoff_value(column: str) -> float:
    matches = re.findall(r"[-+]?\d*\.?\d+", column)
    if not matches:
        return float("nan")
    return float(matches[-1])


def cutoff_label(value: float) -> str:
    return (f"{value:g}").replace(".", "p")


def assign_cutoff_bins(
    df: pd.DataFrame,
    columns: List[str],
    axis_name: str,
    flag_direction: str,
) -> Tuple[pd.Series, pd.Series]:
    if not columns:
        return pd.Series([f"{axis_name}_all"] * len(df), index=df.index), pd.Series([0.0] * len(df), index=df.index)

    col_values = [(col, cutoff_value(col)) for col in columns]
    col_values = [(col, val) for col, val in col_values if not math.isnan(val)]
    if not col_values:
        return pd.Series([f"{axis_name}_all"] * len(df), index=df.index), pd.Series([0.0] * len(df), index=df.index)

    values = np.array([val for _, val in col_values], dtype=np.float64)
    bool_matrix = np.vstack([as_bool_array(df[col]) for col, _ in col_values]).T

    labels: List[str] = []
    ranks: List[float] = []
    max_value = float(np.max(values))
    min_value = float(np.min(values))
    for row in bool_matrix:
        true_values = values[row]
        if flag_direction == "le":
            if true_values.size == 0:
                labels.append(f"{axis_name}_gt_{cutoff_label(max_value)}")
                ranks.append(max_value + 1.0)
            else:
                chosen = float(np.min(true_values))
                labels.append(f"{axis_name}_le_{cutoff_label(chosen)}")
                ranks.append(chosen)
        else:
            if true_values.size == 0:
                labels.append(f"{axis_name}_lt_{cutoff_label(min_value)}")
                ranks.append(min_value - 1.0)
            else:
                chosen = float(np.max(true_values))
                labels.append(f"{axis_name}_ge_{cutoff_label(chosen)}")
                ranks.append(chosen)
    return pd.Series(labels, index=df.index), pd.Series(ranks, index=df.index)


def scenario_matrix(df: pd.DataFrame, scenario_cols: List[str]) -> np.ndarray:
    if not scenario_cols:
        return np.zeros((len(df), 0), dtype=bool)
    return np.vstack([as_bool_array(df[col]) for col in scenario_cols]).T


def mask_to_indices(mask: np.ndarray) -> np.ndarray:
    return np.flatnonzero(np.asarray(mask, dtype=bool))


def choose_reference_indices(
    args: argparse.Namespace,
    df: pd.DataFrame,
    no_scenario_mask: np.ndarray,
    default_mask: Optional[np.ndarray] = None,
) -> np.ndarray:
    base = np.ones(len(df), dtype=bool) if default_mask is None else np.asarray(default_mask, dtype=bool)
    if not base.any():
        warn("--default-col excluded every row; ignoring it for the reference pool.")
        base = np.ones(len(df), dtype=bool)

    pool_mode = getattr(args, "reference_pool", "auto")
    if pool_mode in {"default", "ind"}:
        if pool_mode == "default":
            col = getattr(args, "default_col", "Default")
            if not col:
                raise ValueError(
                    "--reference-pool 'default' conflicts with an empty "
                    "--default-col: the pool column is the one --default-col "
                    "would have named. Set --default-col or pick another pool."
                )
        else:
            col = "InD"
        if col not in df.columns:
            raise ValueError(f"--reference-pool {pool_mode!r} requires a {col!r} column in metadata")
        indices = mask_to_indices(as_bool_array(df[col]) & base)
        if indices.size == 0:
            raise ValueError(f"--reference-pool {pool_mode!r} produced an empty reference subset")
        return indices

    candidate = np.asarray(no_scenario_mask, dtype=bool) & base
    if not candidate.any():
        warn("No no-scenario proteins found; using all default-flagged proteins as reference candidate.")
        candidate = base.copy()

    seq_rank = df["_prodpr_seq_rank"].to_numpy(dtype=float) if "_prodpr_seq_rank" in df.columns else None
    str_rank = df["_prodpr_str_rank"].to_numpy(dtype=float) if "_prodpr_str_rank" in df.columns else None

    attempts: List[np.ndarray] = []
    if seq_rank is not None and str_rank is not None:
        attempts.append(candidate & (seq_rank == np.nanmax(seq_rank[candidate])) & (str_rank == np.nanmax(str_rank[candidate])))
    if seq_rank is not None:
        attempts.append(candidate & (seq_rank == np.nanmax(seq_rank[candidate])))
    if str_rank is not None:
        attempts.append(candidate & (str_rank == np.nanmax(str_rank[candidate])))
    attempts.append(candidate)
    attempts.append(base)

    for mask in attempts:
        indices = mask_to_indices(mask)
        if indices.size >= args.n_min:
            return indices
    indices = mask_to_indices(attempts[-1])
    warn(f"Reference subset has fewer than n_min={args.n_min} samples; using {indices.size} samples.")
    return indices


def pattern_name(pattern: Sequence[bool], scenario_cols: Sequence[str]) -> str:
    active = [col for col, flag in zip(scenario_cols, pattern) if flag]
    return "+".join(active) if active else "ID"


def _stratum_rng(seed: int, key: str) -> np.random.Generator:
    """Derive an independent RNG for one stratum.

    Seeding per stratum (instead of sharing one global RNG across strata)
    decouples the control draws: adding/removing/reordering strata no longer
    reshuffles every other stratum's controls. crc32 gives a stable hash
    across processes (unlike Python's salted hash()).
    """
    digest = zlib.crc32(key.encode("utf-8"))
    return np.random.default_rng(np.random.SeedSequence([int(seed), int(digest)]))


def _profile_key(key: object) -> Tuple[str, ...]:
    if isinstance(key, tuple):
        return tuple(str(v) for v in key)
    return (str(key),)


def build_pool_profile_index(
    df: pd.DataFrame,
    ref_pool_indices: np.ndarray,
    match_cols: List[str],
) -> Dict[Tuple[str, ...], np.ndarray]:
    """Map each background profile in the reference pool to its row indices.

    Built once per run so that per-stratum, per-draw control sampling only
    needs dict lookups instead of re-scanning the full frame.
    """
    index: Dict[Tuple[str, ...], np.ndarray] = {}
    if ref_pool_indices.size == 0:
        return index
    for key, sub in df.loc[ref_pool_indices].groupby(match_cols, sort=False, dropna=False):
        index[_profile_key(key)] = sub.index.to_numpy(dtype=np.int64)
    return index


def sample_bg_matched_control(
    df: pd.DataFrame,
    group_indices: np.ndarray,
    ref_pool_indices: np.ndarray,
    rng: np.random.Generator,
    match_cols: List[str],
    pool_index: Dict[Tuple[str, ...], np.ndarray],
) -> np.ndarray:
    """Sample controls whose background dimensions match the stratum members.

    Background dimensions are the descriptor axes NOT selected for evaluation:
    the joint pattern of unselected scenario columns plus the bins of any
    disabled novelty axes. For each background profile present in the stratum,
    sample an equal number of controls from pool members sharing that exact
    profile, so the only difference between stratum members and their
    controls is the selected OOD axes. Profiles with no pool candidates fall
    back to the whole pool. Pools no larger than the request are taken in
    full (the expectation of without-replacement sampling).
    """
    if ref_pool_indices.size == 0:
        return np.array([], dtype=np.int64)

    selected: List[int] = []
    for key, sub_group in df.loc[group_indices].groupby(match_cols, sort=False, dropna=False):
        candidates = pool_index.get(_profile_key(key))
        if candidates is None or candidates.size == 0:
            candidates = ref_pool_indices
        needed = len(sub_group)
        if candidates.size <= needed:
            # Small pool: taking every candidate equals the expectation of
            # without-replacement sampling and avoids with-replacement bias.
            selected.extend(int(x) for x in candidates)
        else:
            sampled = rng.choice(candidates, size=needed, replace=False)
            selected.extend(int(x) for x in sampled)
    return np.asarray(selected, dtype=np.int64)


def compute_prodpr_unified(
    args: argparse.Namespace,
    df: pd.DataFrame,
    scenario_cols: List[str],
    scenario_flags: np.ndarray,
    y_pred: np.ndarray,
    y_true: np.ndarray,
    y_null: np.ndarray,
    baseline_pool: np.ndarray,
) -> Tuple[Optional[float], Optional[float], pd.DataFrame]:
    """Compute ProDPR over a single unified cross of all OOD axes.

    Every sample is mapped to one cell of ``seq_bin x str_bin x scenario
    pattern``; valid cells (>= n_min) become cross strata. An axis (a
    scenario column or a novelty bin) whose samples all fall outside valid
    cells falls back to a single-axis stratum built from its still-unused
    samples; a fallback that cannot reach n_min samples is dropped. Fallback
    candidates are retried smallest-first, re-counting after every claim so
    small axes secure their samples before larger ones. Each retained sample
    belongs to at most one stratum, so the strata form a partition of the
    covered test set. The zero-shift cell (all selected axes at their
    reference state) is excluded from the strata: it is the reference
    population, not an OOD region.

    Controls: when the descriptor uses only a subset of the available axes,
    each stratum's control is sampled from the reference pool so that the
    background dimensions — the unselected scenario columns' joint pattern
    plus any disabled novelty axes' bins — match the stratum members exactly;
    the ratio then isolates the selected OOD factors. When every axis is
    selected (no background dimensions remain), the control is the full
    reference pool itself.
    """
    rows: List[Dict[str, object]] = []
    total = len(df)

    pattern_labels: List[str] = []
    label_bits: Dict[str, str] = {}
    for row_flags in scenario_flags:
        active = any(row_flags.tolist())
        label = pattern_name(row_flags.tolist(), scenario_cols) if active else "none"
        pattern_labels.append(label)
        label_bits.setdefault(
            label, "".join("1" if x else "0" for x in row_flags.tolist())
        )
    label_bits.setdefault(
        "none", "0" * scenario_flags.shape[1] if scenario_flags.shape[1] else ""
    )
    df["_prodpr_pattern"] = pattern_labels

    cross_axes = ["_prodpr_seq_bin"]
    has_str_bin = "_prodpr_str_bin" in df.columns
    if has_str_bin:
        cross_axes.append("_prodpr_str_bin")
    cross_axes.append("_prodpr_pattern")

    used = np.zeros(total, dtype=bool)
    bg_cols = [
        c for c in ("_prodpr_bg_seq_bin", "_prodpr_bg_str_bin", "_prodpr_bg_pattern")
        if c in df.columns
    ]
    control_mode = getattr(args, "control", "matched")
    n_draws = max(1, int(getattr(args, "control_draws", 50)))
    # Matched-control sampling is only meaningful for the 'auto' pool: the
    # 'default'/'ind' pools are canonical fixed reference sets whose
    # whole-pool metric is the intended denominator ('ind' rows are
    # all-False on every OOD axis, so background matching would have nothing
    # to match and silently degrade to whole-pool fallback anyway).
    pool_mode = getattr(args, "reference_pool", "auto")
    matching_allowed = control_mode == "matched" and pool_mode == "auto"
    if bg_cols and control_mode == "matched" and not matching_allowed:
        print(
            f"[prodpr] --control matched ignored: reference pool "
            f"{pool_mode!r} always uses the whole pool as the denominator "
            f"(matched sampling is only available with --reference-pool auto)."
        )
    pool_profile_index = (
        build_pool_profile_index(df, baseline_pool, bg_cols)
        if bg_cols and matching_allowed
        else None
    )

    def make_controls(subset_indices: np.ndarray, key: str) -> List[np.ndarray]:
        """Control draws for one stratum.

        Matched mode with background dimensions: K independent draws from a
        stratum-local RNG (seeded by the stratum key, so strata never
        reshuffle each other); the reference metric is averaged over draws.
        Otherwise: the full reference pool, deterministically.
        """
        if pool_profile_index is not None:
            rng = _stratum_rng(getattr(args, "seed", 0), key)
            return [
                sample_bg_matched_control(
                    df, subset_indices, baseline_pool, rng, bg_cols, pool_profile_index
                )
                for _ in range(n_draws)
            ]
        return [baseline_pool]

    # Zero-shift cell: every *selected* axis sits at its reference state
    # (highest similarity bin on selected novelty axes, "none" scenario
    # pattern). Its retention is ~1 by construction — it is the reference
    # population, not an OOD stratum — so it is excluded from the strata but
    # still marked used (it must not reappear via axis fallback).
    if "_prodpr_seq_rank" in df.columns:
        seq_rank = df["_prodpr_seq_rank"].to_numpy(dtype=float)
        seq_top = seq_rank == np.nanmax(seq_rank)
    else:
        seq_top = np.ones(total, dtype=bool)
    if "_prodpr_str_rank" in df.columns:
        str_rank = df["_prodpr_str_rank"].to_numpy(dtype=float)
        str_top = str_rank == np.nanmax(str_rank)
    else:
        str_top = np.ones(total, dtype=bool)
    zero_shift = seq_top & str_top & (df["_prodpr_pattern"].astype(str).to_numpy() == "none")

    def add_row(
        stratum_type: str,
        seq_bin: Optional[str],
        str_bin: Optional[str],
        scenario_pattern: Optional[str],
        pattern_bits: Optional[str],
        indices: np.ndarray,
        result: DprResult,
    ) -> None:
        row: Dict[str, object] = {
            "stratum_type": stratum_type,
            "seq_bin": seq_bin,
            "str_bin": str_bin,
            "scenario_pattern": scenario_pattern,
            "pattern_bits": pattern_bits,
            **asdict(result),
        }
        if args.weight_mode == "empirical":
            row["weight"] = indices.size / total
        else:
            row["weight"] = 1.0
        rows.append(row)

    # 1. Cross strata: cells of seq_bin x str_bin x scenario pattern. With
    # background (unselected) axes present, controls are sampled to match the
    # members' background pattern; otherwise the full reference pool is used.
    grouped = df.groupby(cross_axes, sort=False, dropna=False)
    for key, group in grouped:
        subset_indices = group.index.to_numpy(dtype=np.int64)
        if subset_indices.size < args.n_min:
            continue
        if zero_shift[subset_indices].all():
            used[subset_indices] = True  # reference corner: not an OOD stratum
            continue
        key_tuple = key if isinstance(key, tuple) else (key,)
        if has_str_bin:
            seq_bin, str_bin, pattern = key_tuple[0], key_tuple[1], key_tuple[2]
        else:
            seq_bin, str_bin, pattern = key_tuple[0], None, key_tuple[1]
        control_draws = make_controls(
            subset_indices, f"cross|{seq_bin}|{str_bin}|{pattern}"
        )
        if not control_draws or all(d.size == 0 for d in control_draws):
            continue
        result = compute_dpr_multi(
            args.metric, y_pred, y_true, y_null,
            subset_indices, control_draws, args.epsilon
        )
        add_row(
            "cross", seq_bin, str_bin, pattern, label_bits.get(pattern),
            subset_indices, result,
        )
        used[subset_indices] = True

    # 2. Fallback strata: axes completely outside the valid cross cells are
    # retried as single-axis strata over their still-unused samples. The
    # axis with the fewest available samples is retried first; available
    # counts are re-computed after every claim, so small axes secure their
    # samples before larger ones can take them. Ties break by candidate
    # order (novelty bins before scenario columns).
    claimed = np.zeros(total, dtype=bool)

    candidates: List[Tuple[Optional[str], Optional[str], Optional[str], np.ndarray]] = []
    for bin_col, rank_col in (
        ("_prodpr_seq_bin", "_prodpr_seq_rank"),
        ("_prodpr_str_bin", "_prodpr_str_rank"),
    ):
        if bin_col not in df.columns or rank_col not in df.columns:
            continue  # disabled axes (no rank) are background dimensions, not fallback candidates
        order = df.groupby(bin_col, sort=False)[rank_col].first().sort_values()
        for bin_label in order.index:
            candidates.append(
                (bin_col, bin_label, None, (df[bin_col] == bin_label).to_numpy(dtype=bool))
            )
    for j, col in enumerate(scenario_cols):
        candidates.append((None, None, col, scenario_flags[:, j]))

    # Only axes with no sample inside a valid cross cell may fall back.
    pending: List[Tuple[Optional[str], Optional[str], Optional[str], np.ndarray]] = [
        (bin_col, bin_label, scenario_col, np.asarray(axis_mask, dtype=bool))
        for bin_col, bin_label, scenario_col, axis_mask in candidates
        if not used[np.asarray(axis_mask, dtype=bool)].any()
    ]

    while pending:
        avail_counts = [
            int(np.flatnonzero(axis_mask & ~claimed).size)
            for _, _, _, axis_mask in pending
        ]
        pick = min(range(len(pending)), key=lambda i: (avail_counts[i], i))
        bin_col, bin_label, scenario_col, axis_mask = pending.pop(pick)
        subset_indices = np.flatnonzero(axis_mask & ~claimed)
        if subset_indices.size < args.n_min:
            continue  # single axis below threshold: drop the axis, keep going
        seq_bin = bin_label if bin_col == "_prodpr_seq_bin" else None
        str_bin = bin_label if bin_col == "_prodpr_str_bin" else None
        axis_label = bin_label if bin_label is not None else scenario_col
        control_draws = make_controls(subset_indices, f"fallback|{axis_label}")
        if not control_draws or all(d.size == 0 for d in control_draws):
            continue
        result = compute_dpr_multi(
            args.metric, y_pred, y_true, y_null,
            subset_indices, control_draws, args.epsilon
        )
        add_row(
            "fallback", seq_bin, str_bin, scenario_col, None,
            subset_indices, result,
        )
        claimed[subset_indices] = True

    if not rows:
        return None, None, pd.DataFrame()

    # Macro weighting: equal importance to all valid strata (default per paper)
    if args.weight_mode == "macro":
        for row in rows:
            row["weight"] = 1.0 / len(rows)

    table = pd.DataFrame(rows)
    valid_dprs = [float(x) for x in table["dpr"].tolist() if x is not None and not math.isnan(float(x))]
    prodpr = sum(valid_dprs) / len(valid_dprs) if valid_dprs else None
    worst_stratum = min(valid_dprs) if valid_dprs else None

    return prodpr, worst_stratum, table


def prepare_metadata(args: argparse.Namespace, metadata: pd.DataFrame) -> Tuple[pd.DataFrame, List[str], np.ndarray]:
    df = metadata.copy().reset_index(drop=True)
    seq_cols = infer_columns(df, args.seq_prefix, args.seq_cols)
    str_cols = infer_columns(df, args.str_prefix, args.str_cols)
    scenario_cols = infer_scenario_columns(
        df,
        args.scenario_cols,
        args.scenario_prefix,
        getattr(args, "scenario_filter", DEFAULT_SCENARIO_FILTER),
    )

    seq_bin, seq_rank = assign_cutoff_bins(df, seq_cols, "seq", args.flag_direction)
    df["_prodpr_seq_bin"] = seq_bin
    if seq_cols:
        df["_prodpr_seq_rank"] = seq_rank
    else:
        # Axis disabled (e.g. --seq-cols ""): demote it to a background
        # matching dimension, binned from all prefix-matching columns present
        # in the metadata.
        bg_seq_cols = [c for c in df.columns if c.startswith(args.seq_prefix)]
        if bg_seq_cols:
            df["_prodpr_bg_seq_bin"] = assign_cutoff_bins(df, bg_seq_cols, "seq", args.flag_direction)[0]
    if str_cols:
        str_bin, str_rank = assign_cutoff_bins(df, str_cols, "str", args.flag_direction)
        df["_prodpr_str_bin"] = str_bin
        df["_prodpr_str_rank"] = str_rank
    else:
        bg_str_cols = [c for c in df.columns if c.startswith(args.str_prefix)]
        if bg_str_cols:
            df["_prodpr_bg_str_bin"] = assign_cutoff_bins(df, bg_str_cols, "str", args.flag_direction)[0]

    flags = scenario_matrix(df, scenario_cols)

    # Background pattern: joint values of scenario columns that are present in
    # the metadata but NOT selected as descriptor axes. Controls are matched
    # on the background dimensions (unselected scenario pattern + disabled
    # novelty bins) so the retention ratio isolates the selected factors.
    unselected = [
        col for col in df.columns
        if col.startswith(args.scenario_prefix) and col not in scenario_cols
    ]
    if unselected:
        bg_flags = scenario_matrix(df, unselected)
        df["_prodpr_bg_pattern"] = [
            pattern_name(row.tolist(), unselected) if any(row.tolist()) else "none"
            for row in bg_flags
        ]
    return df, scenario_cols, flags


def write_outputs(
    args: argparse.Namespace,
    summary: ProDprSummary,
    strata_table: pd.DataFrame,
) -> None:
    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    summary_path = prefix.with_name(prefix.name + "_summary.json")
    summary_csv_path = prefix.with_name(prefix.name + "_summary.csv")
    strata_path = prefix.with_name(prefix.name + "_strata.csv")

    with summary_path.open("w", encoding="utf-8") as f:
        json.dump(asdict(summary), f, indent=2, ensure_ascii=False)
    pd.DataFrame([asdict(summary)]).to_csv(summary_csv_path, index=False)
    strata_table.to_csv(strata_path, index=False)

    print(f"Saved summary JSON: {summary_path}")
    print(f"Saved summary CSV:  {summary_csv_path}")
    print(f"Saved strata details: {strata_path}")


def load_residue_prediction_bundle(path: str) -> Dict[str, np.ndarray]:
    data = np.load(path)
    required = {"logits", "labels", "lengths"}
    missing = required - set(data.files)
    if missing:
        raise ValueError(f"Residue prediction npz is missing keys: {sorted(missing)}")
    return {
        "logits": np.asarray(data["logits"], dtype=np.float32),
        "labels": np.asarray(data["labels"], dtype=np.int64),
        "lengths": np.asarray(data["lengths"], dtype=np.int64),
    }


def build_residue_null_scores(args: argparse.Namespace, bundle: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
    num_classes = int(bundle["logits"].shape[1])
    if args.null_baseline == "zero":
        prior = np.zeros(num_classes, dtype=np.float32)
    else:
        require_train_csv(args)
        try:
            train_df = pd.read_csv(args.train_csv)
            parsed = []
            for value in train_df[args.label_col]:
                parsed.extend(int(item) for item in str(value).strip("[]").replace("\n", " ").split())
        except Exception as exc:
            raise ValueError(
                f"Cannot build label-prior null from train csv {args.train_csv!r}: {exc}"
            ) from exc
        if not parsed:
            raise ValueError(
                f"Cannot build label-prior null from train csv {args.train_csv!r}: "
                f"no residue labels found in column {args.label_col!r}"
            )
        labels_for_prior = np.asarray(parsed, dtype=np.int64)
        labels_for_prior = labels_for_prior[(labels_for_prior >= 0) & (labels_for_prior < num_classes)]
        counts = np.bincount(labels_for_prior, minlength=num_classes).astype(np.float32)
        prior = counts / max(float(counts.sum()), 1.0)
    return {"prior": prior, "labels": bundle["labels"], "lengths": bundle["lengths"]}


def default_row_mask(args: argparse.Namespace, df: pd.DataFrame) -> Optional[np.ndarray]:
    """Boolean mask of rows with fully computed OOD annotations.

    Returns None when the filter is disabled or the column is absent. NaN is
    treated as False — rows with uncomputed annotations never enter the
    reference pool (denominator); strata (numerator) are unaffected.
    """
    col = getattr(args, "default_col", "Default")
    if not col or col not in df.columns:
        return None
    return as_bool_array(df[col])


def align_metadata_rows(args: argparse.Namespace, metadata: pd.DataFrame, expected_rows: int) -> pd.DataFrame:
    if len(metadata) == expected_rows:
        return metadata

    run_dir = Path(args.pred_path).resolve().parent
    aligned_path = run_dir / "test_metadata.csv"
    if aligned_path.exists():
        aligned = pd.read_csv(aligned_path).reset_index(drop=True)
        if len(aligned) == expected_rows:
            warn(
                f"Metadata has {len(metadata)} rows but predictions have {expected_rows}; "
                f"using aligned metadata {aligned_path}."
            )
            return aligned

    indices_path = run_dir / "test_row_indices.npy"
    if indices_path.exists():
        indices = np.asarray(np.load(indices_path), dtype=np.int64).reshape(-1)
        if len(indices) == expected_rows and (
            len(indices) == 0 or (indices.min() >= 0 and indices.max() < len(metadata))
        ):
            warn(
                f"Metadata has {len(metadata)} rows but predictions have {expected_rows}; "
                f"applying row indices from {indices_path}."
            )
            return metadata.iloc[indices].reset_index(drop=True)

    raise ValueError(
        f"Row mismatch: metadata has {len(metadata)} rows, predictions have {expected_rows}. "
        f"Expected {aligned_path} or {indices_path} for alignment."
    )


def main() -> None:
    args = parse_args()
    require_train_csv(args)

    metadata = pd.read_csv(args.metadata_csv).reset_index(drop=True)
    if args.label_col not in metadata.columns:
        alternate = next(
            (column for column in ("label", "labels") if column in metadata.columns),
            None,
        )
        if alternate is not None:
            warn(
                f"Metadata has no {args.label_col!r} column; using {alternate!r} "
                "for labels and train-prior estimation."
            )
            args.label_col = alternate
    if args.metric in {"residue_accuracy", "residue_auroc"}:
        y_pred = load_residue_prediction_bundle(args.pred_path)
        metadata = align_metadata_rows(args, metadata, len(y_pred["lengths"]))
        y_true = np.ones((len(metadata), 1), dtype=np.float32)
        y_null = build_residue_null_scores(args, y_pred)
    else:
        y_pred = load_numeric_matrix(args.pred_path, args.array_key)
        metadata = align_metadata_rows(args, metadata, y_pred.shape[0])
        y_true, label_to_idx = build_y_true(args, metadata, y_pred)
        y_null = build_null_scores(args, metadata, y_true, label_to_idx)
    df, scenario_cols, scenario_flags = prepare_metadata(args, metadata)

    no_scenario_mask = ~scenario_flags.any(axis=1) if scenario_flags.shape[1] else np.ones(len(df), dtype=bool)
    default_mask = default_row_mask(args, df)
    if default_mask is not None:
        print(f"[prodpr] --default-col: reference pool restricted to {int(default_mask.sum())}/{len(default_mask)} rows")
    baseline_pool = choose_reference_indices(args, df, no_scenario_mask, default_mask=default_mask)

    prodpr, worst_stratum, strata_table = compute_prodpr_unified(
        args, df, scenario_cols, scenario_flags, y_pred, y_true, y_null, baseline_pool
    )

    summary = ProDprSummary(
        prodpr=prodpr,
        worst_stratum=worst_stratum,
        n_min=args.n_min,
        weight_mode=args.weight_mode,
        metric=args.metric,
        null_baseline=args.null_baseline,
        scenario_cols=scenario_cols,
        strata_count=int(len(strata_table)),
        reference_size=int(baseline_pool.size),
        output_prefix=args.output_prefix,
    )
    write_outputs(args, summary, strata_table)


if __name__ == "__main__":
    main()