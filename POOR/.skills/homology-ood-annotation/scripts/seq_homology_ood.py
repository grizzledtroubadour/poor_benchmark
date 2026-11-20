#!/usr/bin/env python3
"""Sequence-homology OOD computation (homology-ood-annotation step 1).

Canonical definition, based on fold_classification/scripts/annotate_test_ood.py:
- Search: query=test sequences, target=train+val sequences.
- Tool: mmseqs2 `easy-search -s 7.5` (default e-value <=1e-3).
  The repaired lbs implementation uses -s 7; see the historical table in SKILL.md.
- Per test sample, take maximum reference pident; no hit means 0.0.
  seq_Redundancy_<T> = max pident <T; OOD_Orphan means no hit.
  This ensures OOD_Orphan ⊆ seq_Redundancy_30 ⊆ ... ⊆ seq_Redundancy_90.
- For multiple sequence columns, e.g. ppi aa_seq1/aa_seq2, search each chain separately.
  Merge sample flags with any semantics by default: True if either chain qualifies, matching
  ppi 10_recompute_all_ood.py (A<thr)|(B<thr) and (A==0)|(B==0).

m8 cache: if --m8-cache exists, replay flags without mmseqs2;
otherwise run mmseqs2 and write the m8 there for reuse.
m8 has three columns query,target,pident, with query=unique_id;
multiple-column keys use unique_id|<seq-col>. Only compatible generated caches can be reused.

Usage:
    # Single-chain task: fold definition.
    python seq_homology_ood.py --test cath_test.csv \
        --ref cath_train.csv cath_val.csv --out test_flagged.csv \
        --m8-cache ood_mmseqs_test_vs_trainval.m8
    # Paired-chain task: ppi definition.
    python seq_homology_ood.py --test ppi_test.csv \
        --ref ppi_train.csv ppi_val.csv --seq-col aa_seq1,aa_seq2 \
        --out ppi_flagged.csv
    # Replay existing m8 without mmseqs2.
    python seq_homology_ood.py --test cath_test.csv --ref cath_train.csv \
        --out replay.csv --m8-cache existing.m8
"""
import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

DEFAULT_THRESHOLDS = "90,80,70,60,50,40,30"


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--test", required=True, type=Path, help="Test CSV containing sequence columns")
    ap.add_argument("--ref", required=True, type=Path, nargs="+",
                    help="Reference CSVs (train/val); combine multiple files by union")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--seq-col", default="aa_seq",
                    help="Comma-separated sequence columns, e.g. aa_seq1,aa_seq2 for paired chains")
    ap.add_argument("--id-col", default="unique_id")
    ap.add_argument("--thresholds", default=DEFAULT_THRESHOLDS,
                    help="Comma-separated pident thresholds: maximum pident < threshold means True")
    ap.add_argument("--merge", choices=["any", "all"], default="any",
                    help="Sample-level merge for multiple chains: any means either chain qualifies, as in ppi")
    ap.add_argument("--sensitivity", default="7.5",
                    help="mmseqs2 easy-search sensitivity: fold uses 7.5, repaired lbs uses 7")
    ap.add_argument("--m8-cache", type=Path, default=None,
                    help="m8 cache: replay if present; otherwise search and save for reuse")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--tmp-dir", type=Path, default=None,
                    help="Temporary FASTA/mmseqs directory; created and cleaned automatically by default")
    ap.add_argument("--default-col", default="Default",
                    help="Default column; process Default=False output rows according to --nondefault")
    ap.add_argument("--nondefault", choices=["nan", "false"], default="nan",
                    help="Default=False values: nan follows the project-wide generation-time convention; "
                         "false is retained only for legacy compatibility")
    return ap.parse_args()


def read_csv(p):
    return pd.read_csv(p, dtype=str, keep_default_na=False)


def query_key(uid, col, multi):
    return f"{uid}|{col}" if multi else uid


def collect_records(df, seq_cols, id_col, multi):
    """Collect deduplicated (key, seq) records; keep the first per key."""
    records = {}
    cols = [id_col] + seq_cols
    for row in df[cols].itertuples(index=False):
        uid, seqs = row[0], row[1:]
        for col, seq in zip(seq_cols, seqs):
            if not seq:
                continue
            key = query_key(uid, col, multi)
            if key not in records:
                records[key] = seq
    return records


def write_fasta(records, path):
    with open(path, "w") as fh:
        for key, seq in records.items():
            fh.write(f">{key}\n{seq}\n")


def run_mmseqs(q_fa, t_fa, m8_path, tmp_dir, sensitivity, threads):
    if shutil.which("mmseqs") is None:
        sys.exit("Error: mmseqs (mmseqs2) not found on PATH. Install it and retry, "
                 "or use --m8-cache to replay an existing m8 without mmseqs2.")
    cmd = ["mmseqs", "easy-search", str(q_fa), str(t_fa), str(m8_path),
           str(tmp_dir / "mmseqs_tmp"), "-s", str(sensitivity),
           "--format-output", "query,target,pident",
           "--threads", str(threads)]
    print("mmseqs2:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def parse_m8_max_pident(m8_path):
    """Maximum pident per query; skip q==t self-alignments, as in repaired lbs."""
    max_pident = {}
    with open(m8_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            q, t, pident = parts[0], parts[1], float(parts[2])
            if q == t:
                continue
            if pident > max_pident.get(q, 0.0):
                max_pident[q] = pident
    return max_pident


def main():
    args = parse_args()
    seq_cols = [c.strip() for c in args.seq_col.split(",") if c.strip()]
    multi = len(seq_cols) > 1
    thresholds = [float(t) for t in args.thresholds.split(",")]
    thr_labels = [t.rstrip("0").rstrip(".") if "." in t else t
                  for t in args.thresholds.split(",")]

    test = read_csv(args.test)
    refs = pd.concat([read_csv(p) for p in args.ref], ignore_index=True)
    for col in seq_cols:
        for name, df in (("test", test), ("ref", refs)):
            if col not in df.columns:
                sys.exit(f"Error: {name} lacks sequence column {col!r}; available: {list(df.columns)}.")
    if args.id_col not in test.columns:
        sys.exit(f"Error: test lacks ID column {args.id_col!r}")
    print(f"test={len(test)}, ref={len(refs)} ({len(args.ref)} files)")

    tmp_ctx = None
    replay = args.m8_cache is not None and args.m8_cache.exists()
    if replay:
        m8_path = args.m8_cache
        print(f"Replaying existing m8: {m8_path} (without mmseqs2)")
    else:
        tmp_ctx = (tempfile.TemporaryDirectory(prefix="seq_homology_")
                   if args.tmp_dir is None else None)
        tmp_dir = args.tmp_dir or Path(tmp_ctx.name)
        tmp_dir.mkdir(parents=True, exist_ok=True)
        q_fa, t_fa = tmp_dir / "query.fasta", tmp_dir / "target.fasta"
        q_records = collect_records(test, seq_cols, args.id_col, multi)
        # Deduplicate targets by sequence only; keys do not affect classification.
        t_records = {}
        seen_seqs = set()
        for row in refs[seq_cols].itertuples(index=False):
            for seq in row:
                if seq and seq not in seen_seqs:
                    seen_seqs.add(seq)
                    t_records[f"t{len(t_records)}"] = seq
        write_fasta(q_records, q_fa)
        write_fasta(t_records, t_fa)
        print(f"Unique query sequences: {len(q_records)}; unique target sequences: {len(t_records)}")
        m8_path = args.m8_cache or (tmp_dir / "result.m8")
        m8_path.parent.mkdir(parents=True, exist_ok=True)
        run_mmseqs(q_fa, t_fa, m8_path, tmp_dir, args.sensitivity, args.threads)

    max_pident = parse_m8_max_pident(m8_path)
    hit_keys = set(max_pident)
    print(f"Queries with hits: {len(hit_keys)}")

    # Per-chain maximum pident (0.0 without a hit) and flags.
    chain_cols = []
    for col in seq_cols:
        mp_col = "max_pident" if not multi else f"max_pident_{col}"
        keys = test[args.id_col].map(lambda u: query_key(u, col, multi))
        test[mp_col] = keys.map(lambda k: max_pident.get(k, 0.0))
        chain_cols.append((col, mp_col, keys.isin(hit_keys)))

    def merge(flag_per_chain):
        out = flag_per_chain[0]
        for f in flag_per_chain[1:]:
            out = out | f if args.merge == "any" else out & f
        return out

    for thr, label in zip(thresholds, thr_labels):
        col = f"seq_Redundancy_{label}"
        test[col] = merge([(test[mp] < thr) for _, mp, _ in chain_cols])
        test[col] = test[col].map({True: "True", False: "False"})
    test["OOD_Orphan"] = merge([~has_hit for _, _, has_hit in chain_cols])
    test["OOD_Orphan"] = test["OOD_Orphan"].map({True: "True", False: "False"})

    # Check threshold monotonicity and subset relations.
    flags = {label: test[f"seq_Redundancy_{label}"] == "True" for label in thr_labels}
    asc = sorted(zip(thresholds, thr_labels))
    for (_, lo), (_, hi) in zip(asc, asc[1:]):
        assert bool((flags[lo] <= flags[hi]).all()), f"Threshold monotonicity violated: {lo} vs {hi}"
    lowest = asc[0][1]
    orphan = test["OOD_Orphan"] == "True"
    assert bool((orphan <= flags[lowest]).all()), "OOD_Orphan is not a subset of the lowest-threshold column!"
    print("Subset relations verified: "
          f"OOD_Orphan ⊆ {' ⊆ '.join(f'seq_Redundancy_{l}' for _, l in asc)}")

    # Default=False extreme-length rows: output columns receive NaN (not computed)
    # directly by project convention; finalize_ood_columns.py --mode nan is fallback validation.
    result_cols = [f"seq_Redundancy_{label}" for label in thr_labels] + ["OOD_Orphan"]
    if args.default_col in test.columns:
        nd = test[args.default_col] != "True"
        if nd.any():
            fill = "" if args.nondefault == "nan" else "False"
            test.loc[nd, result_cols] = fill
            print(f"Default=False: {int(nd.sum())} rows set to {'NaN' if fill == '' else 'False'}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    test.to_csv(args.out, index=False, lineterminator="\n")
    print(f"-> {args.out} ({len(test)} rows)")
    for label in thr_labels:
        print(f"  seq_Redundancy_{label}: {int(flags[label].sum())}")
    print(f"  OOD_Orphan: {int(orphan.sum())}")
    if tmp_ctx is not None:
        tmp_ctx.cleanup()


if __name__ == "__main__":
    main()
