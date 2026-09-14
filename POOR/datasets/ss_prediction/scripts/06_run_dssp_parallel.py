#!/usr/bin/env python3
"""
Parallel DSSP runner using shell script approach.
Splits work into N subsets and runs them in parallel.
"""

import csv
import os
import subprocess
import sys
from pathlib import Path

TASK_DIR = Path(__file__).resolve().parent.parent  # datasets/ss_prediction
PDB_DIR = TASK_DIR / "pdbs"
DSSP_DIR = TASK_DIR / "output" / "dssp"
OUTPUT_CSV = TASK_DIR / "output" / "dssp_labels.csv"
NUM_WORKERS = min(24, os.cpu_count() or 4)

DSSP_DIR.mkdir(parents=True, exist_ok=True)

MAP_3STATE = {
    "H": "H", "G": "H", "I": "H", "P": "H",
    "E": "E", "B": "E",
    "T": "C", "S": "C", " ": "C", "C": "C",
}


def parse_dssp(dssp_path):
    seq = []
    ss8 = []
    with open(dssp_path) as f:
        lines = f.readlines()
    start_idx = None
    for i, line in enumerate(lines):
        if line.strip().startswith("#  RESIDUE AA STRUCTURE"):
            start_idx = i + 1
            break
    if start_idx is None:
        return None, None, None
    for line in lines[start_idx:]:
        if len(line) < 17:
            continue
        aa = line[13] if len(line) > 13 else " "
        ss = line[16] if len(line) > 16 else " "
        if aa == "!":
            continue
        if aa == " " or aa == "X":
            seq.append("X")
            ss8.append("C")
            continue
        seq.append(aa)
        ss8.append(ss if ss != " " else "C")
    if not seq:
        return None, None, None
    ss3 = "".join(MAP_3STATE.get(s, "C") for s in ss8)
    return "".join(seq), "".join(ss8), ss3


def run_dssp_on_file(pdb_path):
    stem = pdb_path.stem
    dssp_path = DSSP_DIR / f"{stem}.dssp"
    if dssp_path.exists() and dssp_path.stat().st_size > 0:
        seq, ss8, ss3 = parse_dssp(dssp_path)
        if seq:
            return stem, seq, ss8, ss3, "ok"
    try:
        subprocess.run(
            ["mkdssp", str(pdb_path), str(dssp_path), "--output-format", "dssp"],
            capture_output=True, text=True, timeout=60, check=False
        )
        if not dssp_path.exists():
            return stem, None, None, None, "no_output"
        seq, ss8, ss3 = parse_dssp(dssp_path)
        if seq:
            return stem, seq, ss8, ss3, "ok"
        return stem, None, None, None, "parse_fail"
    except Exception as e:
        return stem, None, None, None, f"error:{e}"


def worker(subset_file, out_csv):
    pdb_files = [line.strip() for line in open(subset_file) if line.strip()]
    results = []
    ok = skip = fail = 0
    for i, path in enumerate(pdb_files, 1):
        result = run_dssp_on_file(Path(path))
        status = result[-1]
        if status == "ok":
            ok += 1
            results.append({"chain_key": result[0], "seq": result[1], "ss8": result[2], "ss3": result[3]})
        elif status == "skipped":
            skip += 1
            results.append({"chain_key": result[0], "seq": result[1], "ss8": result[2], "ss3": result[3]})
        else:
            fail += 1
        if i % 200 == 0:
            print(f"  {i}/{len(pdb_files)} ok={ok} fail={fail}", flush=True)
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["chain_key", "seq", "ss8", "ss3"])
        writer.writeheader()
        writer.writerows(results)
    print(f"Done: ok={ok} skip={skip} fail={fail} -> {out_csv}")


def main():
    pdb_files = sorted(PDB_DIR.glob("*.pdb"))
    print(f"Total PDB files: {len(pdb_files)}")

    # Split into subsets
    subset_dir = TASK_DIR / "output" / "dssp_subsets"
    subset_dir.mkdir(parents=True, exist_ok=True)

    per_worker = (len(pdb_files) + NUM_WORKERS - 1) // NUM_WORKERS
    subsets = []
    for i in range(NUM_WORKERS):
        start = i * per_worker
        end = min((i + 1) * per_worker, len(pdb_files))
        subset_file = subset_dir / f"subset_{i}.txt"
        with open(subset_file, "w") as f:
            for p in pdb_files[start:end]:
                f.write(str(p) + "\n")
        subsets.append(subset_file)

    # Launch workers
    pids = []
    for i, subset_file in enumerate(subsets):
        out_csv = subset_dir / f"results_{i}.csv"
        cmd = [sys.executable, __file__, "worker", str(subset_file), str(out_csv)]
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        pids.append((p, out_csv))
        print(f"Launched worker {i} (PID {p.pid})")

    # Wait and collect
    all_results = []
    all_errors = []
    for i, (p, out_csv) in enumerate(pids):
        stdout, _ = p.communicate()
        print(f"Worker {i} finished (returncode {p.returncode}):")
        print(stdout[-500:] if len(stdout) > 500 else stdout)
        if out_csv.exists():
            with open(out_csv) as f:
                reader = csv.DictReader(f)
                for row in reader:
                    all_results.append(row)
        else:
            all_errors.append(f"Missing {out_csv}")

    # Write combined output
    with open(OUTPUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["chain_key", "seq", "ss8", "ss3"])
        writer.writeheader()
        writer.writerows(all_results)

    print(f"\nWrote {len(all_results)} entries to {OUTPUT_CSV}")
    if all_errors:
        print(f"Errors: {all_errors[:5]}")


if __name__ == "__main__":
    if len(sys.argv) >= 4 and sys.argv[1] == "worker":
        worker(sys.argv[2], sys.argv[3])
    else:
        main()
