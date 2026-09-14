#!/usr/bin/env python3
"""Step 2: Extract single-chain structure files for all synthetic design chains.

Sources: POOR/data/pdb/{PDB}.pdb.gz (preferred) or {PDB}.cif.gz
- PDB format: filter ATOM/HETATM lines by chain column (col 22).
- CIF format: pure text-line filtering on _atom_site.auth_asym_id
  (BioPython MMCIFParser chokes on large CIFs; text filtering is proven).

Output (staging, pre-dedup): output/chains_all/{unique_id}.pdb|.cif
Log: output/extract_structures.log
"""

import gzip
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parent.parent
POOR = BASE.parent.parent
PDB_RAW_DIR = POOR / "data" / "pdb"
CHAINS_DIR = BASE / "output" / "chains_all"
CHAINS_DIR.mkdir(parents=True, exist_ok=True)


def find_raw(pdb_id):
    for stem in (pdb_id, pdb_id.upper()):
        for ext in (".pdb.gz", ".pdb", ".cif.gz", ".cif"):
            p = PDB_RAW_DIR / f"{stem}{ext}"
            if p.exists():
                return p
    return None


def extract_chain_pdb(raw_path, chain_id, out_path):
    opener = gzip.open if raw_path.suffix == ".gz" else open
    n = 0
    with opener(raw_path, "rt") as fin, open(out_path, "w") as fout:
        for line in fin:
            # Multi-model (NMR) entries: keep model 1 only (fixed post-hoc by
            # fix_multimodel_pdbs.py, now archived in output/scripts_archive/).
            if line.startswith("ENDMDL"):
                break
            if line.startswith(("ATOM  ", "HETATM")):
                if line[21:22] == chain_id:
                    fout.write(line)
                    n += 1
            elif line.startswith("TER") and n > 0 and line[21:22] == chain_id:
                fout.write(line)
    if n == 0:
        os.unlink(out_path)
        return False
    return True


def extract_chain_cif(raw_path, chain_id, out_path):
    opener = gzip.open if raw_path.suffix == ".gz" else open
    cols = {}
    in_atom = False
    chain_col = None
    n = 0
    with opener(raw_path, "rt") as fin, open(out_path, "w") as fout:
        for line in fin:
            s = line.rstrip("\n")
            if s.startswith("_atom_site."):
                in_atom = True
                field = s.split(".", 1)[1].strip()
                cols[field] = len(cols)
                if field == "auth_asym_id":
                    chain_col = cols[field]
                fout.write(line)
                continue
            if s.startswith(("#", "loop_")):
                in_atom = False
                fout.write(line)
                continue
            if not in_atom:
                if not s.startswith("_"):
                    fout.write(line)
                continue
            parts = s.split()
            if len(parts) < 5 or parts[0] not in ("ATOM", "HETATM"):
                continue
            if chain_col is not None and parts[chain_col] != chain_id:
                continue
            fout.write(line)
            n += 1
    if n == 0:
        os.unlink(out_path)
        return False
    return True


def process(row):
    uid, pdb_id, chain_id = row
    for ext in (".pdb", ".cif"):
        p = CHAINS_DIR / f"{uid}{ext}"
        if p.exists():
            return uid, True, ""
    raw = find_raw(pdb_id)
    if raw is None:
        return uid, False, "raw file missing"
    is_cif = ".cif" in raw.name
    out_path = CHAINS_DIR / f"{uid}{'.cif' if is_cif else '.pdb'}"
    try:
        ok = (extract_chain_cif if is_cif else extract_chain_pdb)(
            raw, chain_id, out_path)
    except Exception as e:  # noqa: BLE001
        if out_path.exists():
            os.unlink(out_path)
        return uid, False, f"{type(e).__name__}: {e}"
    return uid, ok, "" if ok else "no atoms matched"


def main():
    df = pd.read_csv(BASE / "output" / "design_chains_all.csv")
    tasks = list(df[["unique_id", "pdb_id", "chain_id"]].itertuples(
        index=False, name=None))
    print(f"chains to extract: {len(tasks)}")

    ok, failed = 0, []
    with ProcessPoolExecutor(max_workers=16) as ex:
        futs = {ex.submit(process, t): t for t in tasks}
        for i, fut in enumerate(as_completed(futs)):
            if i % 200 == 0:
                print(f"  {i}/{len(tasks)} ...")
            uid, success, err = fut.result()
            if success:
                ok += 1
            else:
                failed.append((uid, err))

    print(f"\nsuccess: {ok}, failed: {len(failed)}")
    if failed:
        with open(BASE / "output" / "extract_structures.log", "w") as f:
            for uid, err in failed:
                f.write(f"{uid}\t{err}\n")
        for uid, err in failed[:20]:
            print(f"  FAILED {uid}: {err}")


if __name__ == "__main__":
    main()
