#!/usr/bin/env python3
"""Compute per-split experimental vs predicted (AlphaFold) structure proportions."""
import os
import json
import pandas as pd
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def is_predicted(filename: str) -> bool:
    """Return True if filename indicates an AlphaFold/predicted structure."""
    if not isinstance(filename, str):
        return False
    name = os.path.basename(filename).lower()
    return name.startswith("af-") or name.startswith("af__") or "alphafold" in name


def load_csv(path: str) -> pd.DataFrame:
    return pd.read_csv(path, low_memory=False)


# Accumulator: task -> split -> source -> count
results = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
task_totals = defaultdict(lambda: defaultdict(int))

# 1. enzyme_kinetics_prediction (kcat) — 100 % AlphaFold v6
for split in ["train", "val", "test"]:
    df = load_csv(f"{ROOT}/datasets/enzyme_kinetics_prediction/splits/kcat_{split}.csv")
    pred = df["struct_file"].apply(is_predicted).sum()
    total = len(df)
    results["kcat"][split]["Predicted"] = int(pred)
    results["kcat"][split]["Experimental"] = int(total - pred)
    task_totals["kcat"][split] = total

# 2. enzyme_optimal_ph — mixed AF / PDB, use metadata
ph_meta = load_csv(f"{ROOT}/datasets/enzyme_optimal_ph/data/metadata.csv")
ph_meta = ph_meta[["uniprot_id", "structure_source"]].drop_duplicates()
for split in ["train", "val", "test"]:
    df = load_csv(f"{ROOT}/datasets/enzyme_optimal_ph/splits/optimal_ph_prediction_{split}.csv")
    merged = df.merge(ph_meta, left_on="unique_id", right_on="uniprot_id", how="left")
    # structure_source: AF = predicted, PDB = experimental
    pred = (merged["structure_source"] == "AF").sum()
    exp = (merged["structure_source"] == "PDB").sum()
    total = len(df)
    results["optimal_ph"][split]["Predicted"] = int(pred)
    results["optimal_ph"][split]["Experimental"] = int(exp)
    results["optimal_ph"][split]["Unknown"] = int(total - pred - exp)
    task_totals["optimal_ph"][split] = total

# 3. fold_classification — CATH S95 domains from RCSB PDB => experimental
for split in ["train", "val", "test"]:
    df = load_csv(f"{ROOT}/datasets/fold_classification/splits/cath_{split}.csv")
    total = len(df)
    results["fold_classification"][split]["Experimental"] = total
    task_totals["fold_classification"][split] = total

# 4. func_prediction (EC / GO) — SIFTS PDB => experimental
for split in ["train", "val", "test"]:
    for sub in ["ec", "go_bp", "go_cc", "go_mf"]:
        df = load_csv(f"{ROOT}/datasets/func_prediction/splits/{sub}_{split}.csv")
        total = len(df)
        key = f"func_prediction_{sub}"
        results[key][split]["Experimental"] = total
        task_totals[key][split] = total

# 5. ligand_binding_affinity — PDBbind => experimental
for split in ["train", "val", "test"]:
    df = load_csv(f"{ROOT}/datasets/ligand_binding_affinity/splits/ligand_binding_affinity_{split}.csv")
    total = len(df)
    results["ligand_binding_affinity"][split]["Experimental"] = total
    task_totals["ligand_binding_affinity"][split] = total

# 6. ligand_binding_site — RCSB PDB local mirror => experimental
for split in ["train", "val", "test"]:
    df = load_csv(f"{ROOT}/datasets/ligand_binding_site/splits/ligand_binding_site_{split}.csv")
    total = len(df)
    results["ligand_binding_site"][split]["Experimental"] = total
    task_totals["ligand_binding_site"][split] = total

# 7. ppi_prediction — PINDER; pair-level, predicted if any chain is AF
for split in ["train", "val", "test"]:
    df = load_csv(f"{ROOT}/datasets/ppi_prediction/splits/ppi_{split}.csv")
    pred1 = df["struct_file1"].apply(is_predicted)
    pred2 = df["struct_file2"].apply(is_predicted)
    any_pred = (pred1 | pred2).sum()
    both_pred = (pred1 & pred2).sum()
    both_exp = (~pred1 & ~pred2).sum()
    total = len(df)
    results["ppi_prediction"][split]["Any predicted"] = int(any_pred)
    results["ppi_prediction"][split]["Both experimental"] = int(both_exp)
    results["ppi_prediction"][split]["Both predicted"] = int(both_pred)
    task_totals["ppi_prediction"][split] = total

# 8. ppis_prediction — DIPS-PLUS derived from PDB => experimental
for split in ["train", "val", "test"]:
    df = load_csv(f"{ROOT}/datasets/ppis_prediction/splits/ppis_{split}.csv")
    total = len(df)
    results["ppis_prediction"][split]["Experimental"] = total
    task_totals["ppis_prediction"][split] = total

# 9. ss_prediction — RCSB PDB (X-ray/NMR/EM filtered) => experimental
for split in ["train", "val", "test"]:
    df = load_csv(f"{ROOT}/datasets/ss_prediction/splits/ssp_{split}.csv")
    total = len(df)
    results["ss_prediction"][split]["Experimental"] = total
    task_totals["ss_prediction"][split] = total

# Convert to JSON-serializable structure
out = {
    "per_task_per_split": {
        task: {
            split: dict(counts)
            for split, counts in splits.items()
        }
        for task, splits in results.items()
    },
    "totals": {
        task: dict(splits)
        for task, splits in task_totals.items()
    },
}

os.makedirs(f"{ROOT}/analysis/output", exist_ok=True)
with open(f"{ROOT}/analysis/output/structure_source_stats.json", "w") as f:
    json.dump(out, f, indent=2)

# Print human-readable summary
print(f"{'Task':<25} {'Split':<6} {'Total':>8} {'Experimental':>12} {'Predicted':>12} {'Other':>10}")
print("-" * 75)
for task, splits in results.items():
    for split, counts in splits.items():
        total = task_totals[task][split]
        exp = counts.get("Experimental", 0) + counts.get("Both experimental", 0)
        pred = counts.get("Predicted", 0) + counts.get("Any predicted", 0)
        other = total - exp - pred
        print(f"{task:<25} {split:<6} {total:>8} {exp:>12} ({100*exp/total:5.1f}%) {pred:>12} ({100*pred/total:5.1f}%) {other:>10}")
