#!/usr/bin/env python3
"""Annotate cath_test.csv with OOD columns (Topology-label splits, v2).

v2 conventions:
  - Only ONE prior: `OOD_NewEC` (EC7-holdout, set aside before splitting).
      Prior rows FULLY participate in all dependent OOD axes since
                  2026-08-27 (historical backfill: output/scripts_archive/fill_prior_ood.py;
                  previously placeholder False under a purely dependent definition).
  - `OOD_SuperfamilyHoldout` is a POST-SPLIT annotation: test samples whose
      superfamily (C.A.T.H) does not appear in train/val. These are normal
      split samples and fully participate in all dependent OOD axes.
      - ExtremeShort rows (len(aa_seq) <60, including prior-held-out rows): since 2026-08-27,
                  annotate length faithfully for all rows independently of Default, as in PPI;
                  previously required non-prior and <60 aa:
      Default = True (they went through the normal split) and FULLY
                  participate in all dependent OOD axes since 2026-08-27; historical backfill:
                  output/scripts_archive/fill_es_ood.py. Previously all were left NaN.
      - InD = all OOD flags False (formerly Pure_ID). The final union also includes
                  OOD_LongTail_EC_* and is recomputed by
      .skills/ood-annotation-toolkit/scripts/finalize_ood_columns.py --mode ind
                  after all OOD columns are available.
  - All other rows: fully computed.

Stages write checkpoints under output/ and are skipped on re-run:
  output/ood_mmseqs_test_vs_trainval.m8
  output/ood_foldseek_test_vs_trainval.m8
  output/idr_ratio_full.csv   (incremental per-sequence cache)
"""

import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent.parent
SPLITS_DIR = BASE / "splits"
SUMMARY_CSV = BASE / "output" / "cath_s95_summary.csv"
BACKUP_TEST = BASE / "output" / "splits_ca_backup" / "cath_test.csv"
PDB_DIR = BASE / "pdbs"
OUT_DIR = BASE / "output"
TMP_DIR = OUT_DIR / "ood_tmp"

MMSEQS_M8 = OUT_DIR / "ood_mmseqs_test_vs_trainval.m8"
FOLDSEEK_M8 = OUT_DIR / "ood_foldseek_test_vs_trainval.m8"
IDR_CSV = OUT_DIR / "idr_ratio_full.csv"

SEQ_THRESHOLDS = [90, 80, 70, 60, 50, 40, 30]
TM_THRESHOLDS = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3]
SHORT_LEN = 60
LONGTAIL_FREQ = 10
IDR_THRESHOLD = 0.3

OUT_DIR.mkdir(parents=True, exist_ok=True)
TMP_DIR.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------------- load
print("Loading splits ...")
train = pd.read_csv(SPLITS_DIR / "cath_train.csv", dtype={"label": str})
val = pd.read_csv(SPLITS_DIR / "cath_val.csv", dtype={"label": str})
test = pd.read_csv(SPLITS_DIR / "cath_test.csv", dtype={"label": str})
print(f"  train={len(train)}, val={len(val)}, test={len(test)}")

# superfamily (C.A.T.H) for every domain, from the CATH summary
meta = pd.read_csv(SUMMARY_CSV, usecols=["domain_id", "C", "A", "T", "H"])
meta["superfamily"] = (
    meta["C"].astype(int).astype(str) + "."
    + meta["A"].astype(int).astype(str) + "."
    + meta["T"].astype(int).astype(str) + "."
    + meta["H"].astype(int).astype(str)
)
sf_map = dict(zip(meta["domain_id"], meta["superfamily"]))
for df in (train, val, test):
    df["superfamily"] = df["unique_id"].map(sf_map)
assert test["superfamily"].notna().all()

test["seq_len"] = test["aa_seq"].str.len()
test["is_prior"] = test["OOD_NewEC"] == True
print(f"  prior rows (OOD_NewEC): {test['is_prior'].sum()}")

# ------------------------------------------------- ExtremeShort / Default
test["OOD_ExtremeShort"] = test["seq_len"] < SHORT_LEN  # Include prior-held-out rows (since 2026-08-27).
# Default marks samples from ordinary splitting. Short sequences were not pre-reserved;
# they undergo the same stratified 7:1:2 split, so extreme-short samples remain Default=True,
# as in PPI. Default=False is reserved for explicitly pre-held-out scenarios.
test["Default"] = True
n_es = test["OOD_ExtremeShort"].sum()
print(f"  ExtremeShort: {n_es} ({n_es / len(test) * 100:.1f}%)")

# rows that participate in dependent OOD axes
# Since 2026-08-27, evaluate dependent OOD for all rows, including ES and NewEC-prior rows;
# mmseqs/foldseek queries already cover all test samples.
test["is_dependent"] = True
print(f"  dependent rows: {test['is_dependent'].sum()}")


def write_fasta(df, path):
    with open(path, "w") as fh:
        for uid, seq in zip(df["unique_id"], df["aa_seq"]):
            fh.write(f">{uid}\n{seq}\n")


# ---------------------------------------------------- mmseqs2 (seq axes)
if MMSEQS_M8.exists():
    print(f"mmseqs2: reusing {MMSEQS_M8.name}")
else:
    q_fa = TMP_DIR / "test.fasta"
    t_fa = TMP_DIR / "trainval.fasta"
    write_fasta(test, q_fa)
    write_fasta(pd.concat([train, val], ignore_index=True), t_fa)
    cmd = [
        "mmseqs", "easy-search", str(q_fa), str(t_fa), str(MMSEQS_M8),
        str(TMP_DIR / "mmseqs_tmp"),
        "-s", "7.5",
        "--format-output", "query,target,pident",
    ]
    print("mmseqs2:", " ".join(cmd))
    subprocess.run(cmd, check=True)

print("Parsing mmseqs2 results ...")
hits = pd.read_csv(
    MMSEQS_M8, sep="\t", header=None, names=["query", "target", "pident"],
    usecols=[0, 1, 2],
)
max_pident = hits.groupby("query")["pident"].max()
print(f"  queries with hit: {len(max_pident)} / {len(test)}")

for thr in SEQ_THRESHOLDS:
    col = f"seq_Redundancy_{thr}"
    val_ser = max_pident.reindex(test["unique_id"]).fillna(0.0) < thr
    test[col] = np.where(test["is_dependent"], val_ser.to_numpy(), np.where(
        test["is_prior"], False, np.nan))
test["OOD_Orphan"] = np.where(
    test["is_dependent"],
    ~test["unique_id"].isin(set(max_pident.index)),
    np.where(test["is_prior"], False, np.nan),
)

# ----------------------------------------------------- foldseek (TM axes)
if FOLDSEEK_M8.exists():
    print(f"foldseek: reusing {FOLDSEEK_M8.name}")
else:
    # foldseek needs per-structure PDB files; build query/target dirs of
    # symlinks to avoid copying ~10 GB.
    q_dir = TMP_DIR / "fs_query"
    t_dir = TMP_DIR / "fs_target"
    for d in (q_dir, t_dir):
        d.mkdir(exist_ok=True)
    for uid in test["unique_id"]:
        src = PDB_DIR / f"{uid}.pdb"
        dst = q_dir / f"{uid}.pdb"
        if not dst.exists():
            os.symlink(src, dst)
    for uid in pd.concat([train["unique_id"], val["unique_id"]]):
        src = PDB_DIR / f"{uid}.pdb"
        dst = t_dir / f"{uid}.pdb"
        if not dst.exists():
            os.symlink(src, dst)
    cmd = [
        "foldseek", "easy-search", str(q_dir), str(t_dir), str(FOLDSEEK_M8),
        str(TMP_DIR / "foldseek_tmp"),
        "--format-output", "query,target,alntmscore",
    ]
    print("foldseek:", " ".join(cmd))
    subprocess.run(cmd, check=True)

print("Parsing foldseek results ...")
fs = pd.read_csv(
    FOLDSEEK_M8, sep="\t", header=None, names=["query", "target", "alntmscore"],
    usecols=[0, 1, 2],
)
# query/target may carry paths or .pdb suffixes depending on foldseek version
fs["query"] = fs["query"].astype(str).str.replace(r"\.pdb$", "", regex=True)
fs["query"] = fs["query"].str.split("/").str[-1]
max_tm = fs.groupby("query")["alntmscore"].max()
print(f"  queries with hit: {len(max_tm)} / {len(test)}")

for thr in TM_THRESHOLDS:
    col = f"TM-score_{thr}"
    val_ser = max_tm.reindex(test["unique_id"]).fillna(0.0) < thr
    test[col] = np.where(test["is_dependent"], val_ser.to_numpy(), np.where(
        test["is_prior"], False, np.nan))

# ------------------------------------------------------------------- IDR
# per-sequence cache: old test idr_ratio + incremental idr_ratio_full.csv
idr_cache = {}
old_test = pd.read_csv(BACKUP_TEST, dtype={"label": str})
idr_cache.update(dict(zip(old_test["unique_id"], old_test["idr_ratio"])))
if IDR_CSV.exists():
    idr_df = pd.read_csv(IDR_CSV)
    idr_cache.update(dict(zip(idr_df["unique_id"], idr_df["idr_ratio"])))

missing = [
    u for u in test["unique_id"]
    if u not in idr_cache or pd.isna(idr_cache.get(u))
]
print(f"IDR: {len(test) - len(missing)} cached, {len(missing)} to compute")
if missing:
    import metapredict as meta
    seqs = test.set_index("unique_id").loc[missing, "aa_seq"]
    # old pipeline filtered non-standard residues before prediction
    clean = seqs.str.replace(r"[^ACDEFGHIKLMNPQRSTVWY]", "", regex=True)
    valid = clean[clean.str.len() > 0]
    ratios = {u: np.nan for u in missing if u not in valid.index}
    uids = list(valid.index)
    seq_list = list(valid.values)
    print(f"  running metapredict on {len(seq_list)} sequences ...")
    CHUNK = 512
    for i in range(0, len(seq_list), CHUNK):
        chunk_uids = uids[i:i + CHUNK]
        chunk_seqs = seq_list[i:i + CHUNK]
        scores = meta.predict_disorder_batch(chunk_seqs, show_progress_bar=False)
        # each element is [sequence, per-residue scores]
        for uid, sc in zip(chunk_uids, scores):
            arr = np.asarray(sc[1], dtype=float)
            ratios[uid] = float((arr > 0.5).mean()) if len(arr) else np.nan
        print(f"    {min(i + CHUNK, len(seq_list))}/{len(seq_list)}")
    idr_cache.update(ratios)
    # persist the incremental cache
    idr_out = pd.DataFrame(
        {"unique_id": list(idr_cache.keys()),
         "idr_ratio": list(idr_cache.values())}
    )
    idr_out.to_csv(IDR_CSV, index=False)

test["idr_ratio"] = test["unique_id"].map(idr_cache)
test["OOD_IDR"] = np.where(
    test["is_dependent"],
    test["idr_ratio"] > IDR_THRESHOLD,
    np.where(test["is_prior"], False, np.nan),
)

# ----------------------------------------- SuperfamilyHoldout (post-split)
trainval_sfs = set(train["superfamily"]) | set(val["superfamily"])
sf_holdout = ~test["superfamily"].isin(trainval_sfs)
print(f"SF-holdout (post-split): {sf_holdout.sum()} samples "
      f"({sf_holdout.sum() / len(test) * 100:.1f}%), "
      f"{test.loc[sf_holdout, 'superfamily'].nunique()} superfamilies, "
      f"{test.loc[sf_holdout, 'label'].nunique()} topologies")
test["OOD_SuperfamilyHoldout"] = np.where(
    test["is_dependent"], sf_holdout.to_numpy(),
    np.where(test["is_prior"], False, np.nan),
)

# --------------------------------------------------------------- LongTail
freq = train["label"].value_counts()
rare_topos = set(freq[freq <= LONGTAIL_FREQ].index)
print(f"LongTail: {len(rare_topos)} rare topologies (train freq <= {LONGTAIL_FREQ})")
lt = test["label"].isin(rare_topos)
# Since 2026-08-27, compute for all rows including ES.
test["OOD_LongTail"] = lt.to_numpy()


# ----------------------------------------------------------------- InD
# InD: all OOD flags False; samples belong to no OOD scenario.
# Formerly Pure_ID. This union excludes later-added OOD_LongTail_EC_* columns;
# the final definition is recomputed by .skills/ood-annotation-toolkit/scripts/finalize_ood_columns.py
# --mode ind after all OOD columns are ready.
seq_cols = [f"seq_Redundancy_{t}" for t in SEQ_THRESHOLDS]
tm_cols = [f"TM-score_{t}" for t in TM_THRESHOLDS]
ood_flag_cols = (seq_cols + tm_cols
                 + ["OOD_Orphan", "OOD_SuperfamilyHoldout", "OOD_NewEC",
                    "OOD_ExtremeShort", "OOD_LongTail", "OOD_IDR"]
                 + [c for c in test.columns if c.startswith("OOD_LongTail_EC_")])
flags = test[ood_flag_cols].apply(
    lambda s: pd.to_numeric(s, errors="coerce").fillna(0)
)
test["InD"] = (flags.sum(axis=1) == 0)
print(f"InD (no OOD at all): {test['InD'].sum()} "
      f"({test['InD'].sum() / len(test) * 100:.1f}%)")

# ------------------------------------------------------------------ write
final_cols = (
    ["unique_id", "aa_seq", "struct_file", "label", "Default", "InD"]
    + seq_cols + tm_cols
    + ["OOD_Orphan", "OOD_SuperfamilyHoldout", "OOD_NewEC",
       "OOD_ExtremeShort", "OOD_LongTail", "idr_ratio", "OOD_IDR"]
)

out = test[final_cols].copy()
bool_cols = (["Default", "InD"] + seq_cols + tm_cols
             + ["OOD_Orphan", "OOD_SuperfamilyHoldout", "OOD_NewEC",
                "OOD_ExtremeShort", "OOD_LongTail", "OOD_IDR"])
for c in bool_cols:
    s = pd.to_numeric(out[c], errors="coerce")  # -> 1.0 / 0.0 / NaN
    mapped = pd.Series(pd.NA, index=s.index, dtype=object)
    mapped[s == 1] = "True"
    mapped[s == 0] = "False"
    out[c] = mapped
out.to_csv(SPLITS_DIR / "cath_test.csv", index=False)
print(f"\nWrote {SPLITS_DIR / 'cath_test.csv'} ({len(out)} rows, "
      f"{len(final_cols)} cols)")

# ------------------------------------------------------------- summary
print("\n=== OOD summary (new test) ===")
for c in bool_cols + ["idr_ratio"]:
    if c == "idr_ratio":
        print(f"  {c}: NaN={out[c].isna().sum()}")
        continue
    s = test[c] if c in test else out[c]
    print(f"  {c}: True={(s == True).sum()}, False={(s == False).sum()}, "
          f"NaN={s.isna().sum()}")
print("Done.")
