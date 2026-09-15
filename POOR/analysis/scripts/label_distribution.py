"""Label long-tail distribution stats for benchmark_analysis.md §4 (csv-only).

Classification tasks (EC / GO-BP / GO-CC / GO-MF / fold / SSP):
  - primary stats: individual-label frequency in TRAIN split
    (multi-label EC/GO labels are ';'-separated; per-residue SSP: unique tokens).
    n_labels here matches the class counts in paper Table 1.
  - secondary stats ("combo_*"): unique label *combinations* (whole label string
    as one unit), kept for reference.
  - head20/tail50 coverage: top-20% / bottom-50% of labels' share of label occurrences
  - singleton ratio: labels occurring exactly once (train)
  - Zipf coefficient: slope of log(freq) ~ log(rank) fit

Class imbalance: PPI (sample-level pos%), PPIS (residue-level pos%).
Regression tasks: train/test mean+-std and range of the float label.

Output: analysis/output/label_distribution.json
"""
import os, json, math
from collections import OrderedDict, Counter
import pandas as pd
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUTDIR = os.path.join(ROOT, "analysis/output")

CLS_TASKS = OrderedDict([
    ("EC", "datasets/func_prediction/splits/ec_train.csv"),
    ("GO-BP", "datasets/func_prediction/splits/go_bp_train.csv"),
    ("GO-CC", "datasets/func_prediction/splits/go_cc_train.csv"),
    ("GO-MF", "datasets/func_prediction/splits/go_mf_train.csv"),
    ("Fold", "datasets/fold_classification/splits/cath_train.csv"),
    ("SSP", "datasets/ss_prediction/splits/ssp_train.csv"),
])

REG_TASKS = OrderedDict([
    ("LBA", ("datasets/ligand_binding_affinity/splits/ligand_binding_affinity_{}.csv")),
    ("kcat", ("datasets/enzyme_kinetics_prediction/splits/kcat_{}.csv")),
    ("pH", ("datasets/enzyme_optimal_ph/splits/optimal_ph_prediction_{}.csv")),
])


def tokenize(v, task):
    if task == "SSP":
        return [t.strip("[]") for t in v.split()]
    return v.split(";")  # individual labels; single-label tasks yield one token


def _stats(cnt):
    freqs = np.array(sorted(cnt.values(), reverse=True), dtype=float)
    total = freqs.sum()
    n = len(freqs)
    k_head = max(1, int(math.ceil(n * 0.2)))
    head20 = freqs[:k_head].sum() / total
    k_tail = int(n * 0.5)
    tail50 = freqs[n - k_tail:].sum() / total if k_tail > 0 else 0.0
    singletons = float((freqs == 1).mean())
    ranks = np.arange(1, n + 1)
    slope = np.polyfit(np.log(ranks), np.log(freqs), 1)[0]
    return OrderedDict([
        ("n_labels", n),
        ("head20_coverage", round(float(head20), 4)),
        ("tail50_coverage", round(float(tail50), 4)),
        ("singleton_ratio", round(singletons, 4)),
        ("zipf_slope", round(float(slope), 3)),
    ])


def long_tail(path, task):
    df = pd.read_csv(os.path.join(ROOT, path), dtype=str, keep_default_na=False)
    cnt = Counter()
    for v in df["label"]:
        for t in tokenize(v, task):
            if t and t != "-1":
                cnt[t] += 1
    out = _stats(cnt)
    if task != "SSP":
        # secondary: whole label combination as one unit (previous methodology)
        combo_cnt = Counter(v for v in df["label"] if v)
        for k, v in _stats(combo_cnt).items():
            out["combo_" + k] = v
    return out


def main():
    out = OrderedDict()
    for task, path in CLS_TASKS.items():
        out[task] = long_tail(path, task)

    # class imbalance
    ppi_tr = pd.read_csv(os.path.join(ROOT, "datasets/ppi_prediction/splits/ppi_train.csv"),
                         dtype=str, keep_default_na=False)
    ppi_te = pd.read_csv(os.path.join(ROOT, "datasets/ppi_prediction/splits/ppi_test.csv"),
                         dtype=str, keep_default_na=False)
    out["PPI"] = OrderedDict([
        ("train_pos_pct", round(float((ppi_tr["label"] == "1").mean()) * 100, 1)),
        ("test_pos_pct", round(float((ppi_te["label"] == "1").mean()) * 100, 1)),
    ])
    def residue_pos(path):
        df = pd.read_csv(os.path.join(ROOT, path), dtype=str, keep_default_na=False)
        pos = tot = 0
        for v in df["label"]:
            toks = [t.strip("[]") for t in v.split()]
            pos += sum(1 for t in toks if t == "1")
            tot += sum(1 for t in toks if t in ("0", "1"))
        return round(pos / tot * 100, 1)
    out["PPIS"] = OrderedDict([
        ("train_pos_pct", residue_pos("datasets/ppis_prediction/splits/ppis_train.csv")),
        ("test_pos_pct", residue_pos("datasets/ppis_prediction/splits/ppis_test.csv")),
    ])

    # regression label stats
    for task, pat in REG_TASKS.items():
        tr = pd.read_csv(os.path.join(ROOT, pat.format("train")))["label"].astype(float)
        te = pd.read_csv(os.path.join(ROOT, pat.format("test")))["label"].astype(float)
        out[task] = OrderedDict([
            ("train_mean", round(float(tr.mean()), 2)),
            ("train_std", round(float(tr.std()), 2)),
            ("test_mean", round(float(te.mean()), 2)),
            ("test_std", round(float(te.std()), 2)),
            ("range", [round(float(min(tr.min(), te.min())), 2),
                       round(float(max(tr.max(), te.max())), 2)]),
        ])

    fp = os.path.join(OUTDIR, "label_distribution.json")
    with open(fp, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    print(json.dumps(out, indent=2, ensure_ascii=False))
    print("->", fp)


if __name__ == "__main__":
    main()
