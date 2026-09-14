#!/usr/bin/env python3
"""
Generic evaluation CLI for the POOR benchmark.

Given a set of model predictions and the corresponding test CSV, this script
auto-detects the task type and runs OOD-aware evaluation on all subsets.

Usage:
    python -m shared.datasets.poor.run_evaluation \
        --predictions results/predictions.npy \
        --test_csv data/poor_benchmark/func_prediction/splits/ec_test.csv \
        --output_dir results/eval

    # Or load labels from the test CSV directly:
    python -m shared.datasets.poor.run_evaluation \
        --predictions results/predictions.npy \
        --test_csv data/poor_benchmark/ligand_binding_affinity/splits/ligand_binding_affinity_test.csv \
        --output_dir results/eval

    # For multi-label tasks, also provide labels:
    python -m shared.datasets.poor.run_evaluation \
        --predictions results/predictions.npy \
        --labels results/labels.npy \
        --test_csv data/poor_benchmark/func_prediction/splits/ec_test.csv \
        --output_dir results/eval
"""
import os
import sys
import argparse
import json
import numpy as np
import pandas as pd

from shared.datasets.poor.evaluation import (
    detect_task_type,
    evaluate_model_predictions,
)


def parse_args():
    parser = argparse.ArgumentParser(description="POOR Benchmark Unified Evaluation")
    parser.add_argument(
        "--predictions", type=str, required=True,
        help="Path to .npy file containing model predictions (probs/logits/values)"
    )
    parser.add_argument(
        "--labels", type=str, default=None,
        help="Path to .npy file containing ground-truth labels (if not in test_csv)"
    )
    parser.add_argument(
        "--test_csv", type=str, required=True,
        help="Path to the POOR benchmark test CSV with OOD annotation columns"
    )
    parser.add_argument(
        "--task_type", type=str, default=None,
        choices=[None, "regression", "binary", "multilabel", "single_label", "residual"],
        help="Force task type (auto-detected if not specified)"
    )
    parser.add_argument(
        "--threshold", type=float, default=None,
        help="Fixed threshold for Fmax/F1 (auto-selected if not specified)"
    )
    parser.add_argument(
        "--output_dir", type=str, default="./results/eval",
        help="Directory to save evaluation results"
    )
    return parser.parse_args()


def load_predictions(path: str) -> np.ndarray:
    """Load predictions from .npy, .npz, or .csv."""
    if path.endswith(".npy"):
        return np.load(path, allow_pickle=True)
    elif path.endswith(".npz"):
        data = np.load(path, allow_pickle=True)
        # Try common keys
        for key in ["predictions", "preds", "probs", "logits", "arr_0"]:
            if key in data:
                return data[key]
        return data[list(data.keys())[0]]
    elif path.endswith(".csv"):
        df = pd.read_csv(path)
        # Try common column names
        for col in ["prediction", "predictions", "probs", "pred", "score"]:
            if col in df.columns:
                return df[col].values
        return df.iloc[:, 0].values
    else:
        raise ValueError(f"Unknown file format: {path}")


def load_labels(path: str) -> np.ndarray:
    """Load labels from .npy or .csv."""
    if path.endswith(".npy"):
        return np.load(path, allow_pickle=True)
    elif path.endswith(".csv"):
        df = pd.read_csv(path)
        for col in ["label", "labels", "y", "target"]:
            if col in df.columns:
                return df[col].values
        return df.iloc[:, 0].values
    else:
        raise ValueError(f"Unknown file format: {path}")


def main():
    args = parse_args()

    # Load predictions
    predictions = load_predictions(args.predictions)
    print(f"Loaded predictions: shape={predictions.shape}, dtype={predictions.dtype}")

    # Load labels
    if args.labels:
        labels = load_labels(args.labels)
        print(f"Loaded labels: shape={labels.shape}")
    else:
        # For regression/binary tasks, labels come from the test CSV
        test_df = pd.read_csv(args.test_csv)
        if "label" in test_df.columns:
            # For regression/binary: parse as float
            try:
                labels = test_df["label"].astype(float).values
            except (ValueError, TypeError):
                labels = test_df["label"].values
        else:
            print("ERROR: No 'label' column in test CSV and no --labels provided.")
            sys.exit(1)
        print(f"Extracted labels from test CSV: {len(labels)} samples")

    # Load test CSV
    test_df = pd.read_csv(args.test_csv)
    print(f"Test CSV: {len(test_df)} rows, columns: {list(test_df.columns[:10])}...")

    # Run evaluation
    results_df, task_type = evaluate_model_predictions(
        predictions=predictions,
        labels=labels,
        test_df=test_df,
        task_type=args.task_type,
        threshold=args.threshold,
        output_dir=args.output_dir,
        save_predictions=True,
    )

    print(f"\n{'='*60}")
    print(f"Evaluation complete. Task type: {task_type}")
    print(f"Results saved to: {args.output_dir}")
    print(f"  - test_results.csv")
    print(f"  - test_results.json")
    print(f"  - test_predictions.csv")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
