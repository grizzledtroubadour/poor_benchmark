#!/usr/bin/env python3
"""
EnzyBase12k metadata analyzer
对 Zenodo 18405148 的 metadata.csv 进行全面分析
"""
import os
import sys
import json
import re
from collections import Counter
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
OUTPUT_DIR = ROOT / "output"
FIG_DIR = OUTPUT_DIR / "figures"
LOG_DIR = ROOT / "logs"
for d in [OUTPUT_DIR, FIG_DIR, LOG_DIR]:
    d.mkdir(parents=True, exist_ok=True)

META_PATH = DATA_DIR / "metadata.csv"

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
print(f"[INFO] Loading {META_PATH}")
df = pd.read_csv(META_PATH)
print(f"[INFO] Shape: {df.shape}")

# ---------------------------------------------------------------------------
# Basic stats
# ---------------------------------------------------------------------------
basic = {
    "n_rows": int(len(df)),
    "n_columns": int(len(df.columns)),
    "columns": list(df.columns),
    "column_dtypes": {c: str(df[c].dtype) for c in df.columns},
    "missing_counts": {c: int(df[c].isna().sum()) for c in df.columns},
    "missing_ratios": {c: round(float(df[c].isna().sum()) / len(df), 4) for c in df.columns},
}

# ---------------------------------------------------------------------------
# pH optimum analysis
# ---------------------------------------------------------------------------
ph = pd.to_numeric(df["ph_optimum"], errors="coerce")
ph_stats = {
    "count": int(ph.count()),
    "missing": int(ph.isna().sum()),
    "mean": round(float(ph.mean()), 3),
    "median": round(float(ph.median()), 3),
    "std": round(float(ph.std()), 3),
    "min": round(float(ph.min()), 3),
    "max": round(float(ph.max()), 3),
    "q1": round(float(ph.quantile(0.25)), 3),
    "q3": round(float(ph.quantile(0.75)), 3),
    "p01": round(float(ph.quantile(0.01)), 3),
    "p99": round(float(ph.quantile(0.99)), 3),
}

# ---------------------------------------------------------------------------
# Categorical columns
# ---------------------------------------------------------------------------
cat_summary = {}
for col in ["structure_source", "structure_mode", "cut_to_uniprot_domain"]:
    cat_summary[col] = df[col].value_counts(dropna=False).to_dict()

# ---------------------------------------------------------------------------
# PDB resolution
# ---------------------------------------------------------------------------
res = pd.to_numeric(df["resolution"], errors="coerce")
has_pdb = df["structure_source"] == "PDB"
res_stats = {
    "n_with_resolution": int(res.notna().sum()),
    "n_pdb_entries": int(has_pdb.sum()),
    "mean": round(float(res.mean()), 3),
    "median": round(float(res.median()), 3),
    "std": round(float(res.std()), 3),
    "min": round(float(res.min()), 3),
    "max": round(float(res.max()), 3),
    "q1": round(float(res.quantile(0.25)), 3),
    "q3": round(float(res.quantile(0.75)), 3),
}

# ---------------------------------------------------------------------------
# PDB IDs
# ---------------------------------------------------------------------------
pdb_ids = df["pdb_id_final"].dropna().astype(str)
pdb_id_counts = pdb_ids.value_counts().head(20).to_dict()

# ---------------------------------------------------------------------------
# EC analysis
# ---------------------------------------------------------------------------
# EC column is e.g. "1.1.1.1", may be missing
ec_series = df["ec_id"].dropna().astype(str)
# Extract top-level EC class
ec_top = ec_series.str.extract(r"^(\d+)")[0].value_counts().to_dict()
# Full EC counts
ec_counts = ec_series.value_counts().head(20).to_dict()

# ---------------------------------------------------------------------------
# Organism analysis
# ---------------------------------------------------------------------------
org_counts = df["organism"].value_counts().head(20).to_dict()

# ---------------------------------------------------------------------------
# Sequence length
# ---------------------------------------------------------------------------
seq_len = pd.to_numeric(df["seq_length"], errors="coerce")
seq_stats = {
    "count": int(seq_len.count()),
    "mean": round(float(seq_len.mean()), 1),
    "median": round(float(seq_len.median()), 1),
    "std": round(float(seq_len.std()), 1),
    "min": int(seq_len.min()),
    "max": int(seq_len.max()),
    "q1": round(float(seq_len.quantile(0.25)), 1),
    "q3": round(float(seq_len.quantile(0.75)), 1),
}

# ---------------------------------------------------------------------------
# Duplicates / uniqueness
# ---------------------------------------------------------------------------
uniprot_unique = df["uniprot_id"].nunique()
uniprot_duplicates = df["uniprot_id"].duplicated().sum()
duplicate_uniprots = df["uniprot_id"][df["uniprot_id"].duplicated(keep=False)].value_counts().head(10).to_dict()

# ---------------------------------------------------------------------------
# Cross tabulations
# ---------------------------------------------------------------------------
ph_by_source = df.groupby("structure_source")["ph_optimum"].agg(["count", "mean", "median", "std"]).round(3).to_dict()
ph_by_mode = df.groupby("structure_mode")["ph_optimum"].agg(["count", "mean", "median", "std"]).round(3).to_dict()

# ---------------------------------------------------------------------------
# Write JSON report
# ---------------------------------------------------------------------------
report = {
    "basic": basic,
    "ph_optimum": ph_stats,
    "resolution": res_stats,
    "seq_length": seq_stats,
    "categorical": cat_summary,
    "pdb_id_top": pdb_id_counts,
    "ec_top_level": ec_top,
    "ec_top_full": ec_counts,
    "organism_top": org_counts,
    "uniqueness": {
        "n_unique_uniprot": int(uniprot_unique),
        "n_duplicate_uniprot_rows": int(uniprot_duplicates),
        "duplicate_uniprots_top": duplicate_uniprots,
    },
    "ph_by_structure_source": ph_by_source,
    "ph_by_structure_mode": ph_by_mode,
}

report_path = OUTPUT_DIR / "metadata_analysis.json"
with open(report_path, "w") as f:
    json.dump(report, f, indent=2, default=str)
print(f"[INFO] Report saved to {report_path}")

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------
def savefig(name):
    plt.tight_layout()
    path = FIG_DIR / name
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"[INFO] Figure saved to {path}")

# pH distribution
plt.figure(figsize=(8, 5))
plt.hist(ph.dropna(), bins=60, color="steelblue", edgecolor="white")
plt.axvline(ph.median(), color="red", linestyle="--", label=f"median={ph.median():.2f}")
plt.axvline(ph.mean(), color="orange", linestyle="--", label=f"mean={ph.mean():.2f}")
plt.xlabel("pH optimum")
plt.ylabel("Count")
plt.title("Distribution of pH Optimum (n={})".format(ph.count()))
plt.legend()
savefig("ph_optimum_distribution.png")

# pH by structure source (boxplot)
plt.figure(figsize=(6, 5))
df.boxplot(column="ph_optimum", by="structure_source", ax=plt.gca())
plt.suptitle("")
plt.title("pH Optimum by Structure Source")
plt.ylabel("pH optimum")
savefig("ph_by_structure_source.png")

# pH by structure mode (boxplot)
plt.figure(figsize=(8, 5))
df.boxplot(column="ph_optimum", by="structure_mode", ax=plt.gca())
plt.suptitle("")
plt.title("pH Optimum by Structure Mode")
plt.ylabel("pH optimum")
plt.xticks(rotation=15, ha="right")
savefig("ph_by_structure_mode.png")

# Structure source pie
plt.figure(figsize=(6, 6))
source_counts = df["structure_source"].value_counts()
plt.pie(source_counts, labels=source_counts.index, autopct="%1.1f%%", startangle=90)
plt.title("Structure Source Distribution")
savefig("structure_source_pie.png")

# Structure mode bar
plt.figure(figsize=(8, 5))
mode_counts = df["structure_mode"].value_counts()
mode_counts.plot(kind="bar", color="teal")
plt.title("Structure Mode Distribution")
plt.ylabel("Count")
plt.xticks(rotation=15, ha="right")
savefig("structure_mode_bar.png")

# Resolution distribution (PDB only)
plt.figure(figsize=(8, 5))
res_pdb = res[has_pdb].dropna()
plt.hist(res_pdb, bins=50, color="green", edgecolor="white")
plt.axvline(res_pdb.median(), color="red", linestyle="--", label=f"median={res_pdb.median():.2f}")
plt.xlabel("Resolution (Å)")
plt.ylabel("Count")
plt.title(f"Resolution Distribution for PDB Structures (n={len(res_pdb)})")
plt.legend()
savefig("resolution_distribution.png")

# Sequence length distribution
plt.figure(figsize=(8, 5))
plt.hist(seq_len.dropna(), bins=80, color="purple", edgecolor="white")
plt.axvline(seq_len.median(), color="red", linestyle="--", label=f"median={seq_len.median():.0f}")
plt.xlabel("Sequence Length")
plt.ylabel("Count")
plt.title(f"Sequence Length Distribution (n={seq_len.count()})")
plt.legend()
savefig("seq_length_distribution.png")

# EC top-level class bar
plt.figure(figsize=(8, 5))
ec_top_series = pd.Series(ec_top).sort_index()
ec_top_series.plot(kind="bar", color="coral")
plt.xlabel("EC Top-level Class")
plt.ylabel("Count")
plt.title("Distribution of EC Top-level Classes")
plt.xticks(rotation=0)
savefig("ec_top_level_distribution.png")

# Organism top 15
plt.figure(figsize=(10, 6))
org_top15 = df["organism"].value_counts().head(15)
org_top15.plot(kind="barh", color="olive")
plt.xlabel("Count")
plt.title("Top 15 Organisms")
plt.gca().invert_yaxis()
savefig("organism_top15.png")

print("[INFO] Metadata analysis complete.")
