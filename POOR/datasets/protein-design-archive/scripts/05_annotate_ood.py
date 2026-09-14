#!/usr/bin/env python3
"""Step 5: OOD annotation for the design extra test set -> splits/cath_design.csv.

Column conventions mirror fold_classification/scripts/annotate_test_ood.py:
  - seq_Redundancy_* / TM-score_* / OOD_Orphan / OOD_IDR: dependent axes
    computed vs fold train+val.
  - OOD_SuperfamilyHoldout: chain's candidate superfamilies (PDA cath_full
    codes whose C.A.T == label) — holdout iff NONE appears in train/val.
  - OOD_NewEC: SIFTS EC7 lookup by pdb_id (prior semantics; designs are
    expected to be all False). If any True -> treated as prior row
    (dependent axes False).
  - OOD_ExtremeShort: seq_len < 60; ES rows keep idr_ratio, all other OOD
    columns are empty (NaN), matching the fold convention.
  - OOD_LongTail: label train frequency <= 10 (NaN for ES rows).
  - Default = True for ALL rows (extra test set, but no pre-split picked-out samples).
  - InD = all OOD flags False.

Heavy searches are checkpointed:
  output/design_mmseqs_vs_trainval.m8   (from Step 4)
  output/design_foldseek_vs_trainval.m8
"""

import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent.parent
POOR = BASE.parent.parent
OUT_DIR = BASE / "output"
SPLITS_DIR = BASE / "splits"
SPLITS_DIR.mkdir(exist_ok=True)
PDBS_DIR = BASE / "pdbs"
TMP = OUT_DIR / "s05_tmp"
TMP.mkdir(exist_ok=True)

FOLD = POOR / "datasets" / "fold_classification"
FOLD_PDBS = FOLD / "pdbs"
SIFTS_EC = POOR / "data" / "sifts" / "sifts_ec_annotations.json"

MMSEQS_M8 = OUT_DIR / "design_mmseqs_vs_trainval.m8"
FOLDSEEK_M8 = OUT_DIR / "design_foldseek_vs_trainval.m8"

SEQ_THRESHOLDS = [90, 80, 70, 60, 50, 40, 30]
TM_THRESHOLDS = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3]
SHORT_LEN = 60
LONGTAIL_FREQ = 10
IDR_THRESHOLD = 0.3

# ------------------------------------------------------------------- load
base = pd.read_csv(OUT_DIR / "cath_design_base.csv")
train = pd.read_csv(FOLD / "splits" / "cath_train.csv", dtype={"label": str})
val = pd.read_csv(FOLD / "splits" / "cath_val.csv", dtype={"label": str})
meta = pd.read_csv(FOLD / "output" / "cath_s95_summary.csv",
                   usecols=["domain_id", "C", "A", "T", "H"])
meta["superfamily"] = (meta["C"].astype(int).astype(str) + "."
                       + meta["A"].astype(int).astype(str) + "."
                       + meta["T"].astype(int).astype(str) + "."
                       + meta["H"].astype(int).astype(str))
sf_map = dict(zip(meta["domain_id"], meta["superfamily"]))
trainval_sfs = set(train["unique_id"].map(sf_map)) | set(val["unique_id"].map(sf_map))

test = base.copy()
print(f"design test base: {len(test)} rows")

# ------------------------------------------------------- SIFTS EC7 prior
import json
sifts = json.load(open(SIFTS_EC))


def has_ec7(pdb_id, chain_id):
    for key in (f"{pdb_id.upper()}_{chain_id}", f"{pdb_id.upper()}_A",
                f"{pdb_id.upper()}_0"):
        ann = sifts.get(key)
        if ann:
            for a in ann.get("annotations", []):
                ec = a.get("ec_number") or ""
                if ec.startswith("7."):
                    return True
    return False


test["OOD_NewEC"] = [has_ec7(p, c) for p, c in
                     zip(test["pdb_id"], test["chain_id"])]
n_prior = int((test["OOD_NewEC"] == True).sum())  # noqa: E712
print(f"OOD_NewEC (EC7) prior rows: {n_prior}")
test["is_prior"] = test["OOD_NewEC"] == True  # noqa: E712

# ------------------------------------------------------ ExtremeShort flag
test["seq_len"] = test["aa_seq"].str.len()
test["OOD_ExtremeShort"] = (~test["is_prior"]) & (test["seq_len"] < SHORT_LEN)
test["is_dependent"] = (~test["is_prior"]) & (~test["OOD_ExtremeShort"])
print(f"ExtremeShort: {test['OOD_ExtremeShort'].sum()}, "
      f"dependent rows: {test['is_dependent'].sum()}")

# ---------------------------------------------------------- mmseqs2 axes
hits = pd.read_csv(MMSEQS_M8, sep="\t", header=None,
                   names=["query", "target", "pident"], usecols=[0, 1, 2])
max_pident = hits.groupby("query")["pident"].max()
print(f"mmseqs: queries with hit {len(max_pident)}/{len(test)}")

for thr in SEQ_THRESHOLDS:
    col = f"seq_Redundancy_{thr}"
    v = max_pident.reindex(test["unique_id"]).fillna(0.0) < thr
    test[col] = np.where(test["is_dependent"], v.to_numpy(),
                         np.where(test["is_prior"], False, np.nan))
test["OOD_Orphan"] = np.where(
    test["is_dependent"],
    ~test["unique_id"].isin(set(max_pident.index)),
    np.where(test["is_prior"], False, np.nan))

# ---------------------------------------------------------- foldseek axes
if FOLDSEEK_M8.exists():
    print(f"foldseek: reusing {FOLDSEEK_M8.name}")
else:
    q_dir = TMP / "fs_query"
    t_dir = TMP / "fs_target"
    for d in (q_dir, t_dir):
        d.mkdir(exist_ok=True)
    for _, r in test.iterrows():
        src = PDBS_DIR / r["struct_file"]
        dst = q_dir / r["struct_file"]
        if not dst.exists():
            os.symlink(src.resolve(), dst)
    for uid in pd.concat([train["unique_id"], val["unique_id"]]):
        src = FOLD_PDBS / f"{uid}.pdb"
        dst = t_dir / f"{uid}.pdb"
        if not dst.exists():
            os.symlink(src.resolve(), dst)
    cmd = [
        "foldseek", "easy-search", str(q_dir), str(t_dir), str(FOLDSEEK_M8),
        str(TMP / "foldseek_tmp"),
        "--format-output", "query,target,alntmscore",
    ]
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
    v = max_tm.reindex(test["unique_id"]).fillna(0.0) < thr
    test[col] = np.where(test["is_dependent"], v.to_numpy(),
                         np.where(test["is_prior"], False, np.nan))

# ------------------------------------------------------------------- IDR
import metapredict as meta_pred

clean = test["aa_seq"].str.replace(r"[^ACDEFGHIKLMNPQRSTVWY]", "", regex=True)
ratios = {}
uids = list(test["unique_id"])
seqs = list(clean)
valid_idx = [i for i, s in enumerate(seqs) if len(s) > 0]
scores = meta_pred.predict_disorder_batch([seqs[i] for i in valid_idx],
                                          show_progress_bar=False)
for i, sc in zip(valid_idx, scores):
    arr = np.asarray(sc[1], dtype=float)
    ratios[uids[i]] = float((arr > 0.5).mean()) if len(arr) else np.nan
for i, u in enumerate(uids):
    if u not in ratios:
        ratios[u] = np.nan
test["idr_ratio"] = test["unique_id"].map(ratios)
test["OOD_IDR"] = np.where(
    test["is_dependent"], test["idr_ratio"] > IDR_THRESHOLD,
    np.where(test["is_prior"], False, np.nan))

# ------------------------------------------------- SuperfamilyHoldout
def candidate_sfs(row):
    sfs = set()
    for code in str(row["cath_codes"]).split(";"):
        parts = code.split(".")
        if len(parts) >= 4 and ".".join(parts[:3]) == row["label"]:
            sfs.add(".".join(parts[:4]))
    return sfs


test["candidate_superfamilies"] = test.apply(candidate_sfs, axis=1)
sf_holdout = test["candidate_superfamilies"].map(
    lambda s: len(s & trainval_sfs) == 0)
test["OOD_SuperfamilyHoldout"] = np.where(
    test["is_dependent"], sf_holdout.to_numpy(),
    np.where(test["is_prior"], False, np.nan))
print(f"SF-holdout: {(test['OOD_SuperfamilyHoldout'] == True).sum()}")  # noqa: E712

# --------------------------------------------------------------- LongTail
freq = train["label"].value_counts()
rare_topos = set(freq[freq <= LONGTAIL_FREQ].index)
lt = test["label"].isin(rare_topos)
test["OOD_LongTail"] = np.where(~test["OOD_ExtremeShort"], lt.to_numpy(), np.nan)

# ES rows: all other OOD columns empty (fold convention)
es_mask = test["OOD_ExtremeShort"] == True  # noqa: E712
test["OOD_NewEC"] = test["OOD_NewEC"].astype(object)
test.loc[es_mask, "OOD_NewEC"] = np.nan

# ---------------------------------------------------------- Default/Pure
# extra test set：无"划分前预挑出"样本，按 Default 本意整表 True
# （2026-08-27 由 False 改判，与 fold "ES 行也是 Default=True" 惯例一致）
test["Default"] = True
seq_cols = [f"seq_Redundancy_{t}" for t in SEQ_THRESHOLDS]
tm_cols = [f"TM-score_{t}" for t in TM_THRESHOLDS]
ood_flag_cols = (seq_cols + tm_cols
                 + ["OOD_Orphan", "OOD_SuperfamilyHoldout", "OOD_NewEC",
                    "OOD_ExtremeShort", "OOD_LongTail", "OOD_IDR"]
                 + [c for c in test.columns if c.startswith("OOD_LongTail_EC_")])
flags = test[ood_flag_cols].apply(lambda s: pd.to_numeric(s, errors="coerce").fillna(0))
test["InD"] = flags.sum(axis=1) == 0
print(f"InD: {test['InD'].sum()}")

# ------------------------------------------------------------------ write
final_cols = (["unique_id", "aa_seq", "struct_file", "label", "Default",
               "InD"] + seq_cols + tm_cols
              + ["OOD_Orphan", "OOD_SuperfamilyHoldout", "OOD_NewEC",
                 "OOD_ExtremeShort", "OOD_LongTail", "idr_ratio", "OOD_IDR"])
out = test[final_cols].copy()
bool_cols = (["Default", "InD"] + seq_cols + tm_cols
             + ["OOD_Orphan", "OOD_SuperfamilyHoldout", "OOD_NewEC",
                "OOD_ExtremeShort", "OOD_LongTail", "OOD_IDR"])
for c in bool_cols:
    s = pd.to_numeric(out[c], errors="coerce")
    mapped = pd.Series(pd.NA, index=s.index, dtype=object)
    mapped[s == 1] = "True"
    mapped[s == 0] = "False"
    out[c] = mapped
out.to_csv(SPLITS_DIR / "cath_design.csv", index=False)
print(f"\nWrote {SPLITS_DIR / 'cath_design.csv'} ({len(out)} rows)")

print("\n=== OOD summary (cath_design) ===")
for c in bool_cols + ["idr_ratio"]:
    if c == "idr_ratio":
        print(f"  {c}: NaN={out[c].isna().sum()}")
        continue
    s = test[c]
    print(f"  {c}: True={(s == True).sum()}, False={(s == False).sum()}, "  # noqa: E712
          f"NaN={s.isna().sum()}")
print("Done.")
