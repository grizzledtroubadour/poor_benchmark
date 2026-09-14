#!/usr/bin/env python3
"""Step 7: OOD annotation for the EC design extra test set -> splits/ec_design.csv.

Column set/order mirrors func_prediction/splits/ec_test.csv exactly:
  unique_id, aa_seq, struct_file, label,
  seq_Redundancy_90..30, TM-score_0.9..0.3,
  OOD_IDR, OOD_Orphan, OOD_NewEC_L4, OOD_LongTail, OOD_Combinatorial,
  Default, OOD_ExtremeShort, OOD_ExtremeLong

Conventions (func_prediction DATA_PROCESS sections 4.3-4.10):
  - mmseqs easy-search -s 7 vs EC train+val (reuses Step 6 .m8)
  - foldseek easy-search (alntmscore) vs EC train+val structures
  - OOD_IDR: metapredict v3, idr_ratio > 0.3 (non-standard aa filtered)
  - OOD_NewEC_L4: any 4-level EC (4 numeric parts) of the chain absent from
    train+val L4 set (SIFTS raw annotations); chains without L4 -> False
  - OOD_LongTail: any label with EC-train frequency <= 10
  - OOD_Combinatorial: any label pair never co-occurring in train+val
  - ExtremeShort (<60) / ExtremeLong (>1000): other OOD columns = False
    (func 4.10 convention for non-normal-split samples)
  - Default = True for ALL rows (extra test set, no pre-split picked-out samples)
"""

import json
import os
import subprocess
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent.parent
POOR = BASE.parent.parent
OUT_DIR = BASE / "output"
SPLITS_DIR = BASE / "splits"
SPLITS_DIR.mkdir(exist_ok=True)
PDBS_DIR = BASE / "pdbs"
TMP = OUT_DIR / "s07_tmp"
TMP.mkdir(exist_ok=True)

FUNC = POOR / "datasets" / "func_prediction"
FUNC_PDBS = FUNC / "pdbs"
SIFTS_EC = POOR / "data" / "sifts" / "sifts_ec_annotations.json"

MMSEQS_M8 = OUT_DIR / "design_ec_mmseqs_vs_trainval.m8"
FOLDSEEK_M8 = OUT_DIR / "design_ec_foldseek_vs_trainval.m8"

SEQ_THRESHOLDS = [90, 80, 70, 60, 50, 40, 30]
TM_THRESHOLDS = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3]
SHORT_LEN, LONG_LEN = 60, 1000
LONGTAIL_FREQ = 10
IDR_THRESHOLD = 0.3

# ------------------------------------------------------------------- load
base = pd.read_csv(OUT_DIR / "ec_design_base.csv")
train = pd.read_csv(FUNC / "splits" / "ec_train.csv")
val = pd.read_csv(FUNC / "splits" / "ec_val.csv")
print(f"ec design base: {len(base)} rows; train={len(train)}, val={len(val)}")

test = base.copy()

# ------------------------------------------------- ExtremeShort/Long flag
test["OOD_ExtremeShort"] = test["seq_len"] < SHORT_LEN
test["OOD_ExtremeLong"] = test["seq_len"] > LONG_LEN
test["is_normal_len"] = ~(test["OOD_ExtremeShort"] | test["OOD_ExtremeLong"])
print(f"ExtremeShort: {test['OOD_ExtremeShort'].sum()}, "
      f"ExtremeLong: {test['OOD_ExtremeLong'].sum()}")

# ---------------------------------------------------------- mmseqs2 axes
hits = pd.read_csv(MMSEQS_M8, sep="\t", header=None,
                   names=["query", "target", "pident"], usecols=[0, 1, 2])
max_pident = hits.groupby("query")["pident"].max()
print(f"mmseqs: queries with hit "
      f"{test['unique_id'].isin(set(max_pident.index)).sum()}/{len(test)}")

for thr in SEQ_THRESHOLDS:
    col = f"seq_Redundancy_{thr}"
    # .to_numpy() to avoid index-misalignment (reindexed by unique_id)
    test[col] = (max_pident.reindex(test["unique_id"]).fillna(0.0)
                 .to_numpy() < thr)
test["OOD_Orphan"] = ~test["unique_id"].isin(set(max_pident.index))

# ---------------------------------------------------------- foldseek axes
if FOLDSEEK_M8.exists():
    print(f"foldseek: reusing {FOLDSEEK_M8.name}")
else:
    q_dir = TMP / "fs_query"
    t_dir = TMP / "fs_target"
    for d in (q_dir, t_dir):
        d.mkdir(exist_ok=True)
    for _, r in test.iterrows():
        src = (PDBS_DIR / r["struct_file"]).resolve()
        dst = q_dir / r["struct_file"]
        if not dst.exists():
            os.symlink(src, dst)
    for uid, sf in pd.concat([train[["unique_id", "struct_file"]],
                              val[["unique_id", "struct_file"]]]).itertuples(
            index=False, name=None):
        src = (FUNC_PDBS / sf).resolve()
        dst = t_dir / sf
        if not dst.exists():
            os.symlink(src, dst)
    cmd = ["foldseek", "easy-search", str(q_dir), str(t_dir),
           str(FOLDSEEK_M8), str(TMP / "foldseek_tmp"),
           "--format-output", "query,target,alntmscore"]
    print("foldseek:", " ".join(cmd))
    subprocess.run(cmd, check=True)

fs = pd.read_csv(FOLDSEEK_M8, sep="\t", header=None,
                 names=["query", "target", "alntmscore"], usecols=[0, 1, 2])
fs["query"] = (fs["query"].astype(str)
               .str.replace(r"\.(pdb|cif)$", "", regex=True)
               .str.split("/").str[-1])
max_tm = fs.groupby("query")["alntmscore"].max()
print(f"foldseek: queries with hit {len(max_tm)}/{len(test)}")

for thr in TM_THRESHOLDS:
    col = f"TM-score_{thr}"
    test[col] = (max_tm.reindex(test["unique_id"]).fillna(0.0)
                 .to_numpy() < thr)

# ------------------------------------------------------------------- IDR
import metapredict as meta_pred

clean = test["aa_seq"].str.replace(r"[^ACDEFGHIKLMNPQRSTVWY]", "", regex=True)
scores = meta_pred.predict_disorder_batch(list(clean), show_progress_bar=False)
ratios = []
for sc, uid in zip(scores, test["unique_id"]):
    arr = np.asarray(sc[1], dtype=float)
    ratios.append(float((arr > 0.5).mean()) if len(arr) else np.nan)
test["idr_ratio"] = ratios
test["OOD_IDR"] = test["idr_ratio"] > IDR_THRESHOLD

# ----------------------------------------------------------- OOD_NewEC_L4
sifts = json.load(open(SIFTS_EC))


def full_l4(ec):
    parts = ec.split(".")
    return len(parts) == 4 and parts[3] not in ("-", "")


def l4_set_of(uids):
    out = set()
    for u in uids:
        ann = sifts.get(u)
        if ann:
            for a in ann.get("annotations", []):
                ec = a.get("ec_number") or ""
                if full_l4(ec):
                    out.add(ec)
    return out


trainval_l4 = l4_set_of(list(train["unique_id"]) + list(val["unique_id"]))
print(f"train+val L4 EC set: {len(trainval_l4)}")


def newec_l4(ec_l4_set):
    l4s = {e for e in str(ec_l4_set).split(";") if full_l4(e)}
    if not l4s:
        return False  # no L4 -> does not participate
    return any(e not in trainval_l4 for e in l4s)


test["OOD_NewEC_L4"] = test["ec_l4_set"].map(newec_l4)

# ------------------------------------------------------------ OOD_LongTail
freq = {}
for l in train["label"]:
    for x in str(l).split(";"):
        freq[x] = freq.get(x, 0) + 1
rare = {x for x, n in freq.items() if n <= LONGTAIL_FREQ}
print(f"LongTail rare labels (train freq <= {LONGTAIL_FREQ}): {len(rare)}")
test["OOD_LongTail"] = test["label"].map(
    lambda l: any(x in rare for x in str(l).split(";")))

# ------------------------------------------------------ OOD_Combinatorial
tv_pairs = set()
for l in pd.concat([train["label"], val["label"]]):
    xs = sorted(set(str(l).split(";")))
    for a, b in combinations(xs, 2):
        tv_pairs.add((a, b))


def combinatorial(label):
    xs = sorted(set(str(label).split(";")))
    if len(xs) < 2:
        return False
    return any((a, b) not in tv_pairs for a, b in combinations(xs, 2))


test["OOD_Combinatorial"] = test["label"].map(combinatorial)

# --------------------------------------- ES/EL rows: other OOD = False
dep_cols = ([f"seq_Redundancy_{t}" for t in SEQ_THRESHOLDS]
            + [f"TM-score_{t}" for t in TM_THRESHOLDS]
            + ["OOD_IDR", "OOD_Orphan", "OOD_NewEC_L4", "OOD_LongTail",
               "OOD_Combinatorial"])
es_el = ~test["is_normal_len"]
for c in dep_cols:
    test.loc[es_el, c] = False

# --------------------------------------------------------------- Default
# extra test set：无"划分前预挑出"样本，按 Default 本意整表 True（2026-08-27 改判）
test["Default"] = True

# ------------------------------------------------------------------ write
final_cols = (["unique_id", "aa_seq", "struct_file", "label"]
              + [f"seq_Redundancy_{t}" for t in SEQ_THRESHOLDS]
              + [f"TM-score_{t}" for t in TM_THRESHOLDS]
              + ["OOD_IDR", "OOD_Orphan", "OOD_NewEC_L4", "OOD_LongTail",
                 "OOD_Combinatorial", "Default",
                 "OOD_ExtremeShort", "OOD_ExtremeLong"])
out = test[final_cols].copy()
bool_cols = [c for c in final_cols if c.startswith(("seq_Redundancy", "TM-score", "OOD_", "Default"))]
for c in bool_cols:
    out[c] = out[c].map({True: "True", False: "False"})
out.to_csv(SPLITS_DIR / "ec_design.csv", index=False)
print(f"\nWrote {SPLITS_DIR / 'ec_design.csv'} ({len(out)} rows)")

print("\n=== OOD summary (ec_design) ===")
for c in bool_cols:
    n_true = (test[c] == True).sum()  # noqa: E712
    print(f"  {c}: True={n_true}")
print(f"  idr_ratio: median={test['idr_ratio'].median():.3f}, "
      f"max={test['idr_ratio'].max():.3f}")
print("Done.")
