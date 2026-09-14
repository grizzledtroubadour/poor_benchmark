#!/usr/bin/env python3
"""Rebuild fold_classification splits with C.A.T (Topology) labels.

v2 (2026-08): prior OOD redesigned.
  - label = "C.A.T" string (e.g. "3.40.50"), ~1,163 classes.
  - Stratified 7:1:2 split at the Topology level; every Topology present in
    both train and test.
  - Singleton Topologies (1 sample) are dropped.
  - NewEC prior: domains carrying an **EC7** (translocase) L1 class per SIFTS
    are held out of train/val entirely (`OOD_NewEC`); EC7 chosen to minimize
    distribution shift on the test set (~3% of test).
  - Pure-EC7 zero-pool Topologies (all domains are EC7) are dropped — they
    could never appear in train.
  - SuperfamilyHoldout is NOT a prior anymore: it is annotated post-split by
    05_annotate_test_ood.py (test samples whose superfamily is absent from
    train/val). All non-EC7 domains join the split pool.

Outputs:
  splits/cath_train.csv, splits/cath_val.csv, splits/cath_test.csv
  (test carries the OOD_NewEC flag; 05_annotate_test_ood.py adds the rest)
  output/dropped_domains.csv  -- removed singleton + pure-EC7 domains
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent.parent
SUMMARY_CSV = BASE / "output" / "cath_s95_summary.csv"
SPLITS_DIR = BASE / "splits"
POOL_CSV = BASE / "output" / "domain_pool.csv"
BACKUP_DIR = BASE / "output" / "splits_ca_backup"
OUT_DIR = BASE / "output"
SIFTS_EC_JSON = BASE.parent.parent / "data" / "sifts" / "sifts_ec_annotations.json"

SEED = 42
TEST_FRAC = 0.2
VAL_FRAC = 0.1
HOLDOUT_EC_L1 = {"7"}  # EC7 translocases

# ---------------------------------------------------------------- load
print("Loading data ...")
summary = pd.read_csv(SUMMARY_CSV)
core_cols = ["unique_id", "aa_seq", "struct_file"]
if POOL_CSV.exists():
    # 默认：03_build_domain_pool.py 从 summary+CATH 参考序列重建的样本池
    print(f"  pool: {POOL_CSV}")
    pool_all = pd.read_csv(POOL_CSV)[core_cols]
else:
    # 回退：v1 C.A 粒度划分备份（与 03 产物行级等价，保留作历史对照）
    print(f"  pool: fallback to {BACKUP_DIR} (v1 C.A splits)")
    old_train = pd.read_csv(BACKUP_DIR / "cath_train.csv", dtype={"label": str})
    old_val = pd.read_csv(BACKUP_DIR / "cath_val.csv", dtype={"label": str})
    old_test = pd.read_csv(BACKUP_DIR / "cath_test.csv", dtype={"label": str})
    pool_all = pd.concat(
        [old_train[core_cols], old_val[core_cols], old_test[core_cols]],
        ignore_index=True,
    )
assert pool_all["unique_id"].is_unique, "duplicate unique_id across pool"

meta = summary[["domain_id", "pdb_id", "chain", "C", "A", "T", "H"]].copy()
df = pool_all.merge(meta, left_on="unique_id", right_on="domain_id", how="left")
assert df["C"].notna().all(), "unique_id missing from cath_s95_summary.csv"

df["label"] = (
    df["C"].astype(int).astype(str)
    + "."
    + df["A"].astype(int).astype(str)
    + "."
    + df["T"].astype(int).astype(str)
)
df["topology"] = df["label"]
df["superfamily"] = df["label"] + "." + df["H"].astype(int).astype(str)
print(f"  pooled samples: {len(df)}")

# ------------------------------------------------- drop singleton topologies
topo_counts = df["topology"].value_counts()
singleton_topos = set(topo_counts[topo_counts == 1].index)
dropped_singletons = df[df["topology"].isin(singleton_topos)]
df = df[~df["topology"].isin(singleton_topos)].copy()
print(f"  dropped singleton topologies: {len(singleton_topos)} "
      f"({len(dropped_singletons)} domains)")

# -------------------------------------------- EC7 annotation (NewEC prior)
print(f"Loading SIFTS EC annotations: {SIFTS_EC_JSON}")
with open(SIFTS_EC_JSON) as fh:
    ec_ann = json.load(fh)


def get_ec_l1_set(pdb_id, chain):
    chains = [chain] if chain != "0" else ["0", "A", " "]
    for ch in chains:
        rec = ec_ann.get(f"{pdb_id.upper()}_{ch}")
        if rec:
            return {
                a["ec_number"].split(".")[0]
                for a in rec["annotations"]
                if a.get("ec_number") and a["ec_number"] != "-"
            }
    return set()


df["ec_l1"] = [
    get_ec_l1_set(p, c) for p, c in zip(df["pdb_id"], df["chain"])
]
n_ec = (df["ec_l1"].apply(len) > 0).sum()
print(f"  EC annotated: {n_ec} / {len(df)} ({n_ec / len(df) * 100:.1f}%)")

df["is_newec"] = df["ec_l1"].apply(lambda s: bool(s & HOLDOUT_EC_L1))
priors = df[df["is_newec"]].copy()
pool = df[~df["is_newec"]].copy()
print(f"  NewEC (EC7) prior samples: {len(priors)} "
      f"across {priors['topology'].nunique()} topologies")
print(f"  split pool: {len(pool)}")

# --------------------------- drop pure-EC7 zero-pool topologies
zero_pool_topos = set(priors["topology"]) - set(pool["topology"])
dropped_newec = priors[priors["topology"].isin(zero_pool_topos)].copy()
priors = priors[~priors["topology"].isin(zero_pool_topos)]
print(f"  dropped pure-EC7 zero-pool topologies: {len(zero_pool_topos)} "
      f"({len(dropped_newec)} domains)")
print(f"  priors kept in test: {len(priors)}")

# ------------------------------------------------------ stratified split
rng = np.random.default_rng(SEED)
train_idx, val_idx, test_idx = [], [], []
n_forced_train = 0
for topo, grp in pool.groupby("topology"):
    idx = grp.index.to_numpy().copy()
    rng.shuffle(idx)
    n = len(idx)
    if n == 1:
        # sole pool sample must be in train; test coverage comes from priors
        train_idx.extend(idx)
        n_forced_train += 1
        continue
    n_test = max(1, int(round(TEST_FRAC * n)))
    n_val = max(0, int(round(VAL_FRAC * n)))
    if n - n_test - n_val < 1:
        n_val = max(0, n - n_test - 1)
    n_train = n - n_test - n_val
    train_idx.extend(idx[:n_train])
    val_idx.extend(idx[n_train:n_train + n_val])
    test_idx.extend(idx[n_train + n_val:])

train_df = pool.loc[train_idx]
val_df = pool.loc[val_idx]
test_df = pd.concat([pool.loc[test_idx], priors], ignore_index=True)
print(f"  n=1 topologies forced to train: {n_forced_train}")

# ------------------------------------------------------------------- save
OUT_DIR.mkdir(parents=True, exist_ok=True)
dropped = pd.concat([dropped_singletons, dropped_newec], ignore_index=True)
dropped[["unique_id", "topology", "superfamily", "is_newec"]].to_csv(
    OUT_DIR / "dropped_domains.csv", index=False
)

out_cols = ["unique_id", "aa_seq", "struct_file", "label"]
train_df[out_cols].to_csv(SPLITS_DIR / "cath_train.csv", index=False)
val_df[out_cols].to_csv(SPLITS_DIR / "cath_val.csv", index=False)

test_out = test_df[["unique_id", "aa_seq", "struct_file", "label",
                    "is_newec"]].copy()
test_out = test_out.rename(columns={"is_newec": "OOD_NewEC"})
test_out.to_csv(SPLITS_DIR / "cath_test.csv", index=False)

# ------------------------------------------------------------ validation
print("\n=== Validation ===")
tr_ids, va_ids, te_ids = (
    set(train_df["unique_id"]), set(val_df["unique_id"]), set(test_df["unique_id"])
)
assert not (tr_ids & va_ids) and not (tr_ids & te_ids) and not (va_ids & te_ids), \
    "unique_id leakage across splits"
print(f"  sizes: train={len(train_df)}, val={len(val_df)}, test={len(test_df)} "
      f"(total={len(train_df) + len(val_df) + len(test_df)})")
print(f"  unique_id overlap: none")

all_topos = set(df["topology"]) - singleton_topos - zero_pool_topos
tr_topos, te_topos = set(train_df["topology"]), set(test_df["topology"])
assert all_topos <= tr_topos, f"{len(all_topos - tr_topos)} topologies missing in train"
assert all_topos <= te_topos, f"{len(all_topos - te_topos)} topologies missing in test"
print(f"  topologies: {len(all_topos)}, all present in train and test")

tr_freq = train_df["topology"].value_counts()
print(f"  train topology freq: min={tr_freq.min()}, median={tr_freq.median()}, "
      f"max={tr_freq.max()}")

import re
pat = re.compile(r"^\d+\.\d+\.\d+$")
assert train_df["label"].str.match(pat).all() and test_df["label"].str.match(pat).all()
print(f"  label format OK (C.A.T strings)")

# EC7 must not leak into train/val
assert not train_df["is_newec"].any() and not val_df["is_newec"].any()
print(f"  no EC7 (NewEC) samples in train/val")
print("\nDone.")
