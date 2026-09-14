"""
POOR benchmark dataset and evaluation utilities.

Uses lazy imports so that ``evaluation`` can be imported without the
heavy ``protein_dataset`` dependency chain (which requires KPLM's src).

CLI entry points:
  - ``python -m shared.datasets.poor.run_evaluation`` — single-task eval
  - ``python -m shared.datasets.poor.evaluate_all_ood`` — batch eval across all tasks
"""

__all__ = [
    "PoorProteinDataset",
    "detect_task_type",
    "evaluate_subset",
    "evaluate_by_ood_columns",
    "evaluate_model_predictions",
    "run_evaluation_main",
    "evaluate_all_ood",
]


def __getattr__(name):
    if name == "PoorProteinDataset":
        from .protein_dataset import PoorProteinDataset
        return PoorProteinDataset

    if name in ("detect_task_type", "evaluate_subset", "evaluate_by_ood_columns", "evaluate_model_predictions"):
        from . import evaluation as _eval
        return getattr(_eval, name)

    if name == "run_evaluation_main":
        from .run_evaluation import main
        return main

    if name == "evaluate_all_ood":
        from .evaluate_all_ood import evaluate_all_ood
        return evaluate_all_ood

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
