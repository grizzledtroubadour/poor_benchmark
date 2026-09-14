"""Dataset scale statistics for benchmark_analysis.md §1 (csv-only, no structure parsing).

For each of the 12 task views:
  - train/val/test row counts and total
  - number of classes (unique individual labels; ';'-separated multi-labels split)
  - sequence length min/max/median over all splits (PPI: per-chain over aa_seq1/aa_seq2)
  - whether train/val contain extreme-length samples (<60 or >1000 aa; PPI: either chain)

Output: analysis/output/dataset_scale_stats.json
"""
import os, json
from collections import OrderedDict
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUTDIR = os.path.join(ROOT, "analysis/output")
os.makedirs(OUTDIR, exist_ok=True)

TASKS = OrderedDict([
    ("fold_classification", "datasets/fold_classification/splits/cath_{split}.csv"),
    ("func_prediction (EC)", "datasets/func_prediction/splits/ec_{split}.csv"),
    ("func_prediction (GO-BP)", "datasets/func_prediction/splits/go_bp_{split}.csv"),
    ("func_prediction (GO-CC)", "datasets/func_prediction/splits/go_cc_{split}.csv"),
    ("func_prediction (GO-MF)", "datasets/func_prediction/splits/go_mf_{split}.csv"),
    ("ppi_prediction", "datasets/ppi_prediction/splits/ppi_{split}.csv"),
    ("ppis_prediction", "datasets/ppis_prediction/splits/ppis_{split}.csv"),
    ("ss_prediction", "datasets/ss_prediction/splits/ssp_{split}.csv"),
    ("ligand_binding_site", "datasets/ligand_binding_site/splits/ligand_binding_site_{split}.csv"),
    ("ligand_binding_affinity", "datasets/ligand_binding_affinity/splits/ligand_binding_affinity_{split}.csv"),
    ("enzyme_kinetics_prediction", "datasets/enzyme_kinetics_prediction/splits/kcat_{split}.csv"),
    ("enzyme_optimal_ph", "datasets/enzyme_optimal_ph/splits/optimal_ph_prediction_{split}.csv"),
])

PER_RESIDUE = {"ppis_prediction", "ss_prediction", "ligand_binding_site"}
REGRESSION = {"ligand_binding_affinity", "enzyme_kinetics_prediction", "enzyme_optimal_ph"}


def seq_cols(df):
    if "aa_seq1" in df.columns:
        return ["aa_seq1", "aa_seq2"]
    return ["aa_seq"]


def main():
    stats = OrderedDict()
    tot_train = tot_val = tot_test = 0
    for task, pat in TASKS.items():
        dfs = {}
        for split in ("train", "val", "test"):
            dfs[split] = pd.read_csv(os.path.join(ROOT, pat.format(split=split)),
                                     dtype=str, keep_default_na=False)
        n = {s: len(dfs[s]) for s in dfs}
        # sequence lengths over all splits, per chain
        # ('|' separates chains of a complex, e.g. LBA multi-chain entries)
        lens = []
        for s, df in dfs.items():
            for c in seq_cols(df):
                for v in df[c]:
                    for part in v.split("|"):
                        lens.append(len(part))
        all_lens = pd.Series(lens, dtype=int)
        # extreme-length in train/val (sample-level: any chain extreme)
        extreme_tv = 0
        extreme_tv_rows = 0
        for s in ("train", "val"):
            df = dfs[s]
            mask = None
            for c in seq_cols(df):
                m = (df[c].str.len() < 60) | (df[c].str.len() > 1000)
                mask = m if mask is None else (mask | m)
            extreme_tv += int(mask.sum())
            extreme_tv_rows += len(df)
        # classes
        if task in REGRESSION:
            classes = None
        elif task in PER_RESIDUE:
            classes = None
        elif task == "ppi_prediction":
            classes = 2
        else:
            labels = set()
            for df in dfs.values():
                for v in df["label"]:
                    for t in v.split(";"):
                        if t:
                            labels.add(t)
            classes = len(labels)
        stats[task] = OrderedDict([
            ("train", n["train"]), ("val", n["val"]), ("test", n["test"]),
            ("total", n["train"] + n["val"] + n["test"]),
            ("classes", classes),
            ("seq_len_min", int(all_lens.min())),
            ("seq_len_max", int(all_lens.max())),
            ("seq_len_median", float(all_lens.median())),
            ("extreme_in_train_val", extreme_tv),
            ("train_val_rows", extreme_tv_rows),
        ])
        tot_train += n["train"]; tot_val += n["val"]; tot_test += n["test"]
    stats["_grand_total"] = OrderedDict([
        ("train", tot_train), ("val", tot_val), ("test", tot_test),
        ("total", tot_train + tot_val + tot_test),
    ])
    out = os.path.join(OUTDIR, "dataset_scale_stats.json")
    with open(out, "w") as f:
        json.dump(stats, f, indent=2, ensure_ascii=False)
    print(json.dumps(stats, indent=2, ensure_ascii=False))
    print("->", out)


if __name__ == "__main__":
    main()
