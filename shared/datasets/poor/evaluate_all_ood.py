#!/usr/bin/env python3
"""
Batch OOD evaluation across all POOR benchmark tasks.

Scans the POOR benchmark data directory, loads any saved predictions,
and produces a comprehensive evaluation summary across all tasks and OOD subsets.

Usage:
    python -m shared.datasets.poor.evaluate_all_ood \
        --data_root data/poor_benchmark \
        --predictions_root results/predictions \
        --output_dir results/ood_summary

    # If predictions are stored alongside model outputs:
    python -m shared.datasets.poor.evaluate_all_ood \
        --data_root data/poor_benchmark \
        --predictions_root PLM-SAE/outputs/poor_function \
        --output_dir results/ood_summary
"""
import os
import sys
import argparse
import json
import glob
import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple

from shared.datasets.poor.evaluation import (
    _parse_residual_labels,
    detect_task_type,
    evaluate_by_ood_columns,
    evaluate_model_predictions,
)


# ============================================================
# Task directory mapping
# ============================================================

# Maps task directory name to (task_prefix, label_col, expected_task_type)
TASK_MAP = {
    "enzyme_kinetics_prediction": ("kcat", "label", "regression"),
    "enzyme_optimal_ph": ("optimal_ph_prediction", "label", "regression"),
    "ligand_binding_affinity": ("ligand_binding_affinity", "label", "regression"),
    "fold_classification": ("cath", "label", "single_label"),
    "func_prediction": (None, "label", "multilabel"),  # has sub-tasks (ec, go_mf, go_cc, go_bp)
    "ppi_prediction": ("ppi", "label", "binary"),
    "ss_prediction": ("ssp", "label", "residual"),
    "ppis_prediction": ("ppis", "label", "residual"),
    "ligand_binding_site": ("ligand_binding_site", "label", "residual"),
}


def find_prediction_files(predictions_root: str, task_prefix: str) -> List[str]:
    """
    Search for prediction files matching a task prefix.

    Looks for .npy, .npz, and .csv files containing 'predict' in the name
    within the predictions_root. Only files whose name contains the task
    prefix are considered, so predictions of other tasks are never matched.
    Returns candidates in sorted (deterministic) order.
    """
    candidates = []

    # Direct file match patterns
    patterns = [
        os.path.join(predictions_root, f"*{task_prefix}*predict*.npy"),
        os.path.join(predictions_root, f"*{task_prefix}*predict*.npz"),
        os.path.join(predictions_root, f"*{task_prefix}*predict*.csv"),
        os.path.join(predictions_root, "**", f"*{task_prefix}*predict*.npy"),
        os.path.join(predictions_root, "**", f"*{task_prefix}*predict*.npz"),
        os.path.join(predictions_root, "**", f"*{task_prefix}*predict*.csv"),
        os.path.join(predictions_root, "**", f"*{task_prefix}*test_predict*"),
    ]

    for pattern in patterns:
        candidates.extend(glob.glob(pattern, recursive=True))

    # Deduplicate, deterministic order
    return sorted(set(candidates))


def _find_vocab_file(task_dir: str, task_prefix: str) -> Optional[str]:
    """Locate the label vocab file for a task (``{prefix}_vocab.txt`` preferred)."""
    preferred = os.path.join(task_dir, f"{task_prefix}_vocab.txt")
    if os.path.exists(preferred):
        return preferred
    matches = sorted(glob.glob(os.path.join(task_dir, "*_vocab.txt")))
    return matches[0] if matches else None


def _load_vocab_file(vocab_path: str) -> Dict[str, int]:
    """Load a label vocab file (``index\\tlabel`` per line) → {label: index}."""
    vocab = {}
    with open(vocab_path, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) == 2:
                vocab[parts[1]] = int(parts[0])
    return vocab


def _infer_vocab_from_train(
    splits_dir: str, task_prefix: str, label_col: str, split_tokens: bool
) -> Optional[Dict[str, int]]:
    """Fallback vocab: sorted unique labels (or ';' tokens) from the train split."""
    if not os.path.isdir(splits_dir):
        return None
    train_files = sorted(
        f for f in os.listdir(splits_dir) if "train" in f and f.endswith(".csv")
    )
    if not train_files:
        return None
    train_df = pd.read_csv(os.path.join(splits_dir, train_files[0]), usecols=[label_col])
    values = set()
    for val in train_df[label_col].dropna():
        if split_tokens:
            values.update(tok.strip() for tok in str(val).split(";") if tok.strip())
        else:
            values.add(str(val))
    return {label: idx for idx, label in enumerate(sorted(values))}


def parse_test_labels(
    test_df: pd.DataFrame,
    task_type_hint: str,
    task_dir: str,
    splits_dir: str,
    task_prefix: str,
    label_col: str = "label",
) -> np.ndarray:
    """
    Parse the label column of a POOR test CSV according to the task type:

      - regression/binary: float values
      - multilabel: ';'-separated tokens → multi-hot via the task vocab
      - single_label: int class indices, or string classes mapped via vocab
      - residual: bracketed per-residue strings → flat int array (-1 = masked)
    """
    labels_raw = test_df[label_col]

    if task_type_hint in ("regression", "binary"):
        return labels_raw.astype(float).values

    if task_type_hint == "multilabel":
        vocab_file = _find_vocab_file(task_dir, task_prefix)
        if vocab_file is not None:
            vocab = _load_vocab_file(vocab_file)
        else:
            vocab = _infer_vocab_from_train(splits_dir, task_prefix, label_col, split_tokens=True)
            if vocab is None:
                raise ValueError(f"No label vocab found for task '{task_prefix}'")
            print(f"  [WARN] No vocab file for '{task_prefix}'; inferred {len(vocab)} "
                  "labels from the train split (column order may differ from the model's)")
        n_classes = max(vocab.values()) + 1
        multi_hot = np.zeros((len(test_df), n_classes), dtype=np.float64)
        for i, val in enumerate(labels_raw):
            if pd.isna(val):
                continue
            for tok in str(val).split(";"):
                tok = tok.strip()
                if tok in vocab:
                    multi_hot[i, vocab[tok]] = 1.0
        return multi_hot

    if task_type_hint == "single_label":
        # Integer class indices need no vocab
        try:
            return labels_raw.astype(int).values
        except (ValueError, TypeError):
            pass
        vocab_file = _find_vocab_file(task_dir, task_prefix)
        vocab = None
        if vocab_file is not None:
            candidate = _load_vocab_file(vocab_file)
            if all(str(v) in candidate for v in labels_raw.dropna()):
                vocab = candidate
        if vocab is None:
            vocab = _infer_vocab_from_train(splits_dir, task_prefix, label_col, split_tokens=False)
            if vocab is None:
                raise ValueError(f"No class vocab found for task '{task_prefix}'")
            print(f"  [WARN] Vocab file does not cover test labels for '{task_prefix}'; "
                  f"inferred {len(vocab)} classes from the train split "
                  "(column order may differ from the model's)")
        mapped = [vocab.get(str(v), -1) for v in labels_raw]
        if any(m < 0 for m in mapped):
            raise ValueError(f"Test labels of task '{task_prefix}' not covered by the class vocab")
        return np.array(mapped, dtype=np.int64)

    if task_type_hint == "residual":
        flat_labels, _ = _parse_residual_labels(labels_raw)
        return flat_labels

    raise ValueError(f"Unknown task type hint: {task_type_hint}")


def evaluate_task(
    task_dir: str,
    task_prefix: str,
    task_type_hint: str,
    predictions_root: str,
    output_dir: str,
    label_col: str = "label",
) -> Optional[Dict]:
    """
    Evaluate a single task if predictions are available.

    Returns:
        dict with task_name, task_type, n_subsets, results, or None if no predictions found.
    """
    splits_dir = os.path.join(task_dir, "splits")

    # Find test CSV
    test_csv = None
    for name in [f"{task_prefix}_test.csv", f"{task_prefix}_test.csv"]:
        path = os.path.join(splits_dir, name)
        if os.path.exists(path):
            test_csv = path
            break

    if test_csv is None:
        # Try listing files
        if os.path.isdir(splits_dir):
            test_files = [f for f in os.listdir(splits_dir) if "test" in f and f.endswith(".csv")]
            if test_files:
                test_csv = os.path.join(splits_dir, test_files[0])

    if test_csv is None:
        print(f"  [SKIP] No test CSV found in {splits_dir}")
        return None

    # Find predictions
    pred_files = find_prediction_files(predictions_root, task_prefix)
    if not pred_files:
        print(f"  [SKIP] No prediction files found for task '{task_prefix}'")
        return None

    if len(pred_files) > 1:
        # Deterministic choice: most recently modified file, ties broken by name
        pred_file = max(pred_files, key=lambda p: (os.path.getmtime(p), p))
        print(f"  Multiple prediction files found for '{task_prefix}', "
              f"using latest: {pred_file}")
        for p in pred_files:
            if p != pred_file:
                print(f"    (skipped: {p})")
    else:
        pred_file = pred_files[0]
    print(f"  Test CSV: {test_csv}")
    print(f"  Predictions: {pred_file}")

    # Load predictions
    if pred_file.endswith(".npy"):
        predictions = np.load(pred_file, allow_pickle=True)
    elif pred_file.endswith(".npz"):
        data = np.load(pred_file, allow_pickle=True)
        predictions = data[list(data.keys())[0]]
    elif pred_file.endswith(".csv"):
        df = pd.read_csv(pred_file)
        pred_cols = [c for c in df.columns if c.startswith("prediction") or c == "prob" or c == "score"]
        predictions = df[pred_cols[0]].values if pred_cols else df.iloc[:, 0].values
    else:
        print(f"  [SKIP] Unknown prediction file format: {pred_file}")
        return None

    # Load test CSV
    test_df = pd.read_csv(test_csv)

    # Extract labels
    if label_col not in test_df.columns:
        print(f"  [SKIP] No '{label_col}' column in {test_csv}")
        return None
    labels = parse_test_labels(
        test_df, task_type_hint, task_dir, splits_dir, task_prefix, label_col
    )

    # Align: predictions must cover every test row (no silent truncation)
    assert len(predictions) == len(test_df), (
        f"Task '{task_prefix}': predictions have {len(predictions)} rows but "
        f"{test_csv} has {len(test_df)} rows"
    )

    # Run evaluation
    task_output = os.path.join(output_dir, task_prefix)
    os.makedirs(task_output, exist_ok=True)

    results_df, detected_type = evaluate_model_predictions(
        predictions=predictions,
        labels=labels,
        test_df=test_df,
        task_type=task_type_hint,
        output_dir=task_output,
        save_predictions=False,
    )

    return {
        "task": task_prefix,
        "task_type": detected_type,
        "n_test_samples": len(test_df),
        "n_subsets": len(results_df),
        "results": results_df.to_dict(orient="records"),
    }


def evaluate_all_ood(
    data_root: str = "data/poor_benchmark",
    predictions_root: str = "results",
    output_dir: str = "results/ood_summary",
) -> pd.DataFrame:
    """
    Run OOD evaluation across all POOR benchmark tasks.

    Args:
        data_root: root directory of POOR benchmark data
        predictions_root: root directory containing per-task prediction files
        output_dir: where to save the summary

    Returns:
        DataFrame with all results concatenated
    """
    os.makedirs(output_dir, exist_ok=True)

    all_results = []
    failed_tasks = []

    for task_dir_name, (task_prefix, label_col, task_type_hint) in TASK_MAP.items():
        task_dir = os.path.join(data_root, task_dir_name)
        if not os.path.isdir(task_dir):
            continue

        print(f"\n{'='*60}")
        print(f"Task: {task_dir_name}")

        # Handle func_prediction with multiple sub-tasks
        sub_tasks = ["ec", "go_mf", "go_cc", "go_bp"] if task_dir_name == "func_prediction" else [task_prefix]
        for sub_task in sub_tasks:
            hint = "multilabel" if task_dir_name == "func_prediction" else task_type_hint
            sub_output = os.path.join(output_dir, sub_task) if task_dir_name == "func_prediction" else output_dir
            try:
                if task_dir_name == "func_prediction":
                    print(f"\n  Sub-task: {sub_task}")
                result = evaluate_task(
                    task_dir=task_dir,
                    task_prefix=sub_task,
                    task_type_hint=hint,
                    predictions_root=predictions_root,
                    output_dir=sub_output,
                    label_col=label_col,
                )
                if result:
                    all_results.append(result)
            except Exception as e:
                failed_tasks.append((sub_task, str(e)))
                print(f"  [ERROR] Task '{sub_task}' failed: {e}")

    # Build summary
    if not all_results:
        print("\nNo predictions found for any task. Nothing to evaluate.")
        if failed_tasks:
            print("Failed tasks:")
            for name, err in failed_tasks:
                print(f"  - {name}: {err}")
        return pd.DataFrame()

    # Flatten results into a summary table
    summary_rows = []
    for task_result in all_results:
        for subset_result in task_result["results"]:
            row = {
                "task": task_result["task"],
                "task_type": task_result["task_type"],
                "subset": subset_result.get("subset"),
                "n_samples": subset_result.get("n_samples"),
            }
            # Add all metrics
            for k, v in subset_result.items():
                if k not in ("subset", "n_samples"):
                    row[k] = v
            summary_rows.append(row)

    summary_df = pd.DataFrame(summary_rows)
    summary_path = os.path.join(output_dir, "ood_summary.csv")
    summary_df.to_csv(summary_path, index=False)
    print(f"\n{'='*60}")
    print(f"Summary saved to {summary_path}")
    print(f"Total: {len(summary_df)} subset evaluations across {len(all_results)} tasks")
    if failed_tasks:
        print(f"Failed tasks ({len(failed_tasks)}):")
        for name, err in failed_tasks:
            print(f"  - {name}: {err}")

    # Also save as JSON
    summary_json = os.path.join(output_dir, "ood_summary.json")
    with open(summary_json, "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"JSON saved to {summary_json}")

    return summary_df


def main():
    parser = argparse.ArgumentParser(description="Batch OOD Evaluation for POOR Benchmark")
    parser.add_argument("--data_root", type=str, default="data/poor_benchmark")
    parser.add_argument("--predictions_root", type=str, required=True,
                        help="Root directory containing per-task prediction files")
    parser.add_argument("--output_dir", type=str, default="results/ood_summary")
    args = parser.parse_args()

    evaluate_all_ood(
        data_root=args.data_root,
        predictions_root=args.predictions_root,
        output_dir=args.output_dir,
    )


if __name__ == "__main__":
    main()
