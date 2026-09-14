"""Recompute OOD coverage matrix from current split CSVs (post schema migration).

Emits:
  analysis/output/ood_coverage_matrix.csv
  analysis/output/ood_coverage_matrix.md
  analysis/output/ood_overlap_matrix.csv  (pairwise Jaccard overlap of selected OOD axes)
"""
import os, glob, json
import pandas as pd
import numpy as np
from collections import OrderedDict

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUTDIR = os.path.join(ROOT, "analysis/output")
os.makedirs(OUTDIR, exist_ok=True)

TASKS = OrderedDict([
    ("kcat", "datasets/enzyme_kinetics_prediction/splits/kcat_test.csv"),
    ("pH", "datasets/enzyme_optimal_ph/splits/optimal_ph_prediction_test.csv"),
    ("CATH", "datasets/fold_classification/splits/cath_test.csv"),
    ("EC", "datasets/func_prediction/splits/ec_test.csv"),
    ("GO-BP", "datasets/func_prediction/splits/go_bp_test.csv"),
    ("GO-CC", "datasets/func_prediction/splits/go_cc_test.csv"),
    ("GO-MF", "datasets/func_prediction/splits/go_mf_test.csv"),
    ("LBA", "datasets/ligand_binding_affinity/splits/ligand_binding_affinity_test.csv"),
    ("LBS", "datasets/ligand_binding_site/splits/ligand_binding_site_test.csv"),
    ("PPI", "datasets/ppi_prediction/splits/ppi_test.csv"),
    ("PPIS", "datasets/ppis_prediction/splits/ppis_test.csv"),
    ("SSP", "datasets/ss_prediction/splits/ssp_test.csv"),
])

# Each row: (display name, list of candidate column names to sum/union, dtype-of-absence marker)
SCENARIOS = OrderedDict([
    ("Test set size", None),
    ("Default (ID) subset", ["Default"]),
    ("InD subset (no OOD flag)", ["InD"]),
    ("Extreme short (<60 aa)", ["OOD_ExtremeShort"]),
    ("Extreme long (>1000 aa)", ["OOD_ExtremeLong"]),
    ("Seq identity < 90%", ["seq_Redundancy_90"]),
    ("Seq identity < 80%", ["seq_Redundancy_80"]),
    ("Seq identity < 70%", ["seq_Redundancy_70"]),
    ("Seq identity < 60%", ["seq_Redundancy_60"]),
    ("Seq identity < 50%", ["seq_Redundancy_50"]),
    ("Seq identity < 40%", ["seq_Redundancy_40"]),
    ("Seq identity < 30%", ["seq_Redundancy_30"]),
    ("Orphan (no significant hit)", ["OOD_Orphan"]),
    ("TM-score < 0.9", ["TM-score_0.9"]),
    ("TM-score < 0.8", ["TM-score_0.8"]),
    ("TM-score < 0.7", ["TM-score_0.7"]),
    ("TM-score < 0.6", ["TM-score_0.6"]),
    ("TM-score < 0.5", ["TM-score_0.5"]),
    ("TM-score < 0.4", ["TM-score_0.4"]),
    ("TM-score < 0.3", ["TM-score_0.3"]),
    ("IDR-containing", ["OOD_IDR"]),
    ("New EC L4", ["OOD_NewEC_L4", "OOD_NewEC"]),
    ("New EC L3", ["OOD_NewEC_L3"]),
    ("Long-tail (label-frequency)", ["OOD_LongTail"]),
    ("Long-tail EC L3 ≤5", ["OOD_LongTail_EC_L3_le5"]),
    ("Long-tail EC L3 ≤10", ["OOD_LongTail_EC_L3_le10"]),
    ("Long-tail EC L4 ≤5", ["OOD_LongTail_EC_L4_le5"]),
    ("Long-tail EC L4 ≤10", ["OOD_LongTail_EC_L4_le10"]),
    ("Long-tail pH bin ≤50", ["OOD_LongTail_pHbin_le50"]),
    ("Long-tail pH bin ≤100", ["OOD_LongTail_pHbin_le100"]),
    ("Long-tail Aff bin ≤50", ["OOD_LongTail_AffBin_le50"]),
    ("Long-tail Aff bin ≤100", ["OOD_LongTail_AffBin_le100"]),
    ("Long-tail kcat bin ≤50", ["OOD_LongTail_KcatBin_le50"]),
    ("Long-tail kcat bin ≤100", ["OOD_LongTail_KcatBin_le100"]),
    ("Combinatorial (new label pair)", ["OOD_Combinatorial"]),
    ("Fold holdout", ["OOD_FoldHoldout"]),
    ("Superfamily holdout", ["OOD_SuperfamilyHoldout"]),
    ("Has EC annotation", ["has_ec_annotation"]),
])


def load():
    dfs = {}
    for task, path in TASKS.items():
        df = pd.read_csv(os.path.join(ROOT, path), low_memory=False)
        # ensure booleans
        for c in df.columns:
            if df[c].dtype == object:
                try:
                    df[c] = df[c].map(lambda x: {"True": True, "False": False}.get(x, x))
                except Exception:
                    pass
        dfs[task] = df
    return dfs


def count(df, cols):
    if cols is None:
        return len(df)
    present = [c for c in cols if c in df.columns]
    if not present:
        return None
    m = df[present[0]].fillna(False).astype(bool)
    for c in present[1:]:
        m = m | df[c].fillna(False).astype(bool)
    return int(m.sum())


def build_matrix(dfs):
    rows = []
    for scen, cols in SCENARIOS.items():
        row = {"OOD Scenario": scen}
        for task in TASKS:
            df = dfs[task]
            n = len(df) if scen == "Test set size" else count(df, cols)
            if n is None:
                row[task] = "-"
            else:
                pct = n / len(df) * 100
                if scen == "Test set size":
                    row[task] = f"{n}"
                else:
                    row[task] = f"{n} ({pct:.1f}%)"
        rows.append(row)
    return pd.DataFrame(rows)


def jaccard(a, b):
    if len(a) == 0 and len(b) == 0:
        return 0.0
    return float((a & b).sum()) / float((a | b).sum())


def pairwise_overlap(dfs):
    # For tasks with Default + OOD axes, compute Jaccard between selected axes.
    # We focus on the core continuous/discrete OOD columns.
    overlap_rows = []
    axes = ["OOD_ExtremeShort", "OOD_ExtremeLong", "OOD_Orphan",
            "seq_Redundancy_30", "TM-score_0.5", "OOD_IDR",
            "OOD_NewEC_L4", "OOD_NewEC", "OOD_LongTail"]
    for task, df in dfs.items():
        avail = [a for a in axes if a in df.columns]
        # skip if too few
        if len(avail) < 2:
            continue
        for i in range(len(avail)):
            for j in range(i+1, len(avail)):
                ai = df[avail[i]].fillna(False).astype(bool)
                aj = df[avail[j]].fillna(False).astype(bool)
                overlap_rows.append({
                    "task": task,
                    "axis_a": avail[i], "axis_b": avail[j],
                    "jaccard": round(jaccard(ai, aj), 3),
                    "a_only": int((ai & ~aj).sum()),
                    "b_only": int((~ai & aj).sum()),
                    "a_and_b": int((ai & aj).sum()),
                })
    return pd.DataFrame(overlap_rows)


def main():
    dfs = load()
    mat = build_matrix(dfs)
    mat.to_csv(os.path.join(OUTDIR, "ood_coverage_matrix.csv"), index=False)
    # markdown
    md = mat.to_markdown(index=False)
    with open(os.path.join(OUTDIR, "ood_coverage_matrix.md"), "w") as f:
        f.write(md + "\n")
    overlap = pairwise_overlap(dfs)
    overlap.to_csv(os.path.join(OUTDIR, "ood_overlap_matrix.csv"), index=False)
    print(md)
    print(f"\n\nOverlap matrix rows: {len(overlap)}; file: analysis/output/ood_overlap_matrix.csv")
    # summary high-overlap pairs
    high = overlap[overlap["jaccard"] > 0.5].sort_values("jaccard", ascending=False)
    print("\nHigh-overlap pairs (Jaccard > 0.5):")
    print(high.to_string(index=False))


if __name__ == "__main__":
    main()
