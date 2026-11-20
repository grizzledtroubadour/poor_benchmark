#!/usr/bin/env python3
r"""Structural-homology OOD computation (homology-ood-annotation step 2).

Canonical definition based on fold_classification/scripts/annotate_test_ood.py:
- Search: query=test structures, target=train+val structures.
- Tool: foldseek easy-search with query,target,alntmscore output.
- Maximum reference alntmscore per test structure; no hit means 0.0.
  TM-score_<T> = (max alntmscore < T)
- Query IDs remove paths and structure extensions, matching fold normalization:
  `str.replace(r"\.pdb$","")` followed by the final path component.

--test-structs / --ref-structs each accept a PDB/mmCIF directory
(all *.pdb / *.cif files) or a manifest with one structure path per line.

m8 cache: replay existing --m8-cache results without foldseek;
otherwise run foldseek and save m8 to that path for reuse.

Output: one row per test structure: query_id, max_alntmscore, TM-score_<T>...

Usage:
    python struct_homology_ood.py --test-structs pdbs_test.txt \
        --ref-structs pdbs_trainval.txt --out tm_flags.csv \
        --m8-cache ood_foldseek_test_vs_trainval.m8
    # Replay existing m8 without foldseek.
    python struct_homology_ood.py --test-structs pdbs_test.txt \
        --ref-structs pdbs_trainval.txt --out replay.csv --m8-cache existing.m8
"""
import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

DEFAULT_THRESHOLDS = "0.9,0.8,0.7,0.6,0.5,0.4,0.3"
STRUCT_EXTS = (".pdb", ".cif", ".ent")


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--test-structs", required=True, type=Path,
                    help="Test structure directory or manifest with one path per line")
    ap.add_argument("--ref-structs", required=True, type=Path,
                    help="Train+val structure directory or manifest")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--thresholds", default=DEFAULT_THRESHOLDS,
                    help="Comma-separated alntmscore thresholds: maximum < threshold means True")
    ap.add_argument("--m8-cache", type=Path, default=None,
                    help="m8 cache: replay if present, otherwise search and save")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--tmp-dir", type=Path, default=None)
    return ap.parse_args()


def struct_list(spec):
    """Directory -> all structure files; manifest -> one path per line."""
    if spec.is_dir():
        files = sorted(p for p in spec.iterdir() if p.suffix.lower() in STRUCT_EXTS)
    else:
        files = []
        with open(spec) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    files.append(Path(line))
    missing = [p for p in files if not p.exists()]
    if missing:
        sys.exit(f"Error: {spec} lists {len(missing)} missing structure files, "
                 f"for example {missing[0]}")
    return files


def normalize_id(name):
    """Normalize query IDs: strip .pdb/.cif/.ent case-insensitively and take the basename.
    foldseek createdb removes any structure extension, including mmCIF inputs
    such as 8RRI_K.CIF, so strip those extensions too.
    Otherwise every .cif query would incorrectly appear to have no hit."""
    s = str(name)
    low = s.lower()
    for ext in STRUCT_EXTS:
        if low.endswith(ext):
            s = s[: -len(ext)]
            break
    return s.split("/")[-1]


def run_foldseek(q_dir, t_dir, m8_path, tmp_dir, threads):
    if shutil.which("foldseek") is None:
        sys.exit("Error: foldseek not found on PATH. Install it and retry, "
                 "or use --m8-cache to replay existing m8 without foldseek.")
    cmd = ["foldseek", "easy-search", str(q_dir), str(t_dir), str(m8_path),
           str(tmp_dir / "foldseek_tmp"),
           "--format-output", "query,target,alntmscore",
           "--threads", str(threads)]
    print("foldseek:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def parse_m8_max_tm(m8_path):
    """Maximum alntmscore per query; skip q==t self-alignments."""
    max_tm = {}
    with open(m8_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            q, t = normalize_id(parts[0]), normalize_id(parts[1])
            tm = float(parts[2])
            if q == t:
                continue
            if tm > max_tm.get(q, 0.0):
                max_tm[q] = tm
    return max_tm


def main():
    args = parse_args()
    thresholds = [float(t) for t in args.thresholds.split(",")]
    thr_labels = [t.rstrip("0").rstrip(".") if "." in t else t
                  for t in args.thresholds.split(",")]

    replay = args.m8_cache is not None and args.m8_cache.exists()
    tmp_ctx = None
    if replay:
        m8_path = args.m8_cache
        print(f"Replaying existing m8: {m8_path} (without foldseek)")
        # Replay still requires the test manifest to determine output rows.
        test_files = struct_list(args.test_structs)
        query_ids = [normalize_id(p.name) for p in test_files]
    else:
        test_files = struct_list(args.test_structs)
        ref_files = struct_list(args.ref_structs)
        query_ids = [normalize_id(p.name) for p in test_files]
        if len(set(query_ids)) != len(query_ids):
            sys.exit("Error: duplicate test structure IDs after removing extensions prevent one-to-one matching")
        print(f"Test structures: {len(test_files)}; reference structures: {len(ref_files)}")
        tmp_ctx = (tempfile.TemporaryDirectory(prefix="struct_homology_")
                   if args.tmp_dir is None else None)
        tmp_dir = args.tmp_dir or Path(tmp_ctx.name)
        tmp_dir.mkdir(parents=True, exist_ok=True)
        # foldseek requires directories; use symlinks instead of copying, as in fold.
        q_dir, t_dir = tmp_dir / "fs_query", tmp_dir / "fs_target"
        q_dir.mkdir(exist_ok=True)
        t_dir.mkdir(exist_ok=True)
        for p in test_files:
            dst = q_dir / p.name
            if not dst.exists():
                dst.symlink_to(p.resolve())
        for p in ref_files:
            dst = t_dir / p.name
            if not dst.exists():
                dst.symlink_to(p.resolve())
        m8_path = args.m8_cache or (tmp_dir / "result.m8")
        m8_path.parent.mkdir(parents=True, exist_ok=True)
        run_foldseek(q_dir, t_dir, m8_path, tmp_dir, args.threads)

    max_tm = parse_m8_max_tm(m8_path)
    print(f"Queries with hits: {len(max_tm)} / {len(query_ids)}")

    rows = []
    flags_by_label = {label: [] for label in thr_labels}
    for qid in query_ids:
        tm = max_tm.get(qid, 0.0)
        row = {"query_id": qid, "max_alntmscore": tm}
        for thr, label in zip(thresholds, thr_labels):
            v = tm < thr
            row[f"TM-score_{label}"] = str(v)
            flags_by_label[label].append(v)
        rows.append(row)

    # Check threshold monotonicity.
    asc = sorted(zip(thresholds, thr_labels))
    for (_, lo), (_, hi) in zip(asc, asc[1:]):
        for a, b in zip(flags_by_label[lo], flags_by_label[hi]):
            assert a <= b, f"Threshold monotonicity violated: {lo} vs {hi} (query row)"

    out = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False, lineterminator="\n")
    print(f"-> {args.out} ({len(out)} rows)")
    for label in thr_labels:
        print(f"  TM-score_{label}: {sum(flags_by_label[label])}")
    if tmp_ctx is not None:
        tmp_ctx.cleanup()


if __name__ == "__main__":
    main()
