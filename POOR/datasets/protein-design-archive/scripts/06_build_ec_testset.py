#!/usr/bin/env python3
"""Step 6 (v2): Build the EC-task extra test set base table.

Label source (strict, per-chain experimental annotation ONLY — no transfer):
  - chain-level EC parsed from the chain's OWN raw PDB entry
      .pdb: COMPND MOLECULE/CHAIN/EC sub-fields -> {chain: ECs}
      .cif: _entity.pdbx_ec (entity level) mapped to chains via
            _atom_site.(auth_asym_id, label_entity_id)
    UNION SIFTS data/sifts/sifts_ec_annotations.json for that chain.
  - v2: NO within-cluster transfer. Exact-duplicate groups (identical
    aa_seq, i.e. literally the same protein) require CONSISTENT L3 label
    sets across annotated members; conflicting groups are dropped.

Filter:
  - L4->L3 normalization; multi-label ';'-joined; all L3 in EC train labels
  - overlap removal vs EC train/val: pdb_id match OR mmseqs (-s 7)
    max pident >= 95

Outputs:
  output/design_ec_annotations.csv   (per-chain EC hits)
  output/ec_dropped_conflict.csv     (exact-dup groups with conflicting L3)
  output/design_ec_mmseqs_vs_trainval.m8
  output/ec_design_base.csv
  output/ec_removed_overlap.csv
"""

import gzip
import json
import re
import subprocess
from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parent.parent
POOR = BASE.parent.parent
OUT_DIR = BASE / "output"
PDB_RAW_DIR = POOR / "data" / "pdb"
SIFTS_EC = POOR / "data" / "sifts" / "sifts_ec_annotations.json"
FUNC = POOR / "datasets" / "func_prediction"

MMSEQS_M8 = OUT_DIR / "design_ec_mmseqs_vs_trainval.m8"
TMP = OUT_DIR / "s06_tmp"
TMP.mkdir(exist_ok=True)

PIDENT_CUTOFF = 95.0


def find_raw(pdb_id):
    for stem in (pdb_id, pdb_id.upper()):
        for ext in (".pdb.gz", ".cif.gz"):
            p = PDB_RAW_DIR / f"{stem}{ext}"
            if p.exists():
                return p
    return None


def parse_pdb_compnd(path):
    compnd = []
    with gzip.open(path, "rt", errors="ignore") as f:
        for line in f:
            if line.startswith("COMPND"):
                compnd.append(line[10:].rstrip())
            elif compnd and not line.startswith(("COMPND", "SOURCE")):
                break
    result = {}
    for mol in re.split(r"MOLECULE:", "".join(compnd)):
        m_ch = re.search(r"CHAIN:\s*([^;]+)", mol)
        m_ec = re.search(r"EC:\s*([0-9.\-, ]+)", mol)
        if m_ch and m_ec:
            ecs = {e.strip() for e in m_ec.group(1).split(",") if e.strip()}
            for c in [c.strip() for c in m_ch.group(1).split(",")]:
                result.setdefault(c, set()).update(ecs)
    return result


def parse_cif_ec(path):
    entity_ec, chain_entity = {}, {}
    with gzip.open(path, "rt", errors="ignore") as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        if lines[i].strip() == "loop_":
            j = i + 1
            cols = []
            while j < len(lines) and lines[j].strip().startswith("_"):
                cols.append(lines[j].strip())
                j += 1
            if cols and all(c.startswith("_entity.") for c in cols):
                idc = cols.index("_entity.id") if "_entity.id" in cols else None
                ecc = cols.index("_entity.pdbx_ec") if "_entity.pdbx_ec" in cols else None
                if idc is not None and ecc is not None:
                    k = j
                    while k < len(lines) and not lines[k].startswith(("#", "loop_", "_")):
                        parts = lines[k].split()
                        if len(parts) == len(cols) and parts[ecc] not in ("?", "."):
                            entity_ec[parts[idc]] = parts[ecc]
                        k += 1
            if cols and all(c.startswith("_atom_site.") for c in cols):
                ac = (cols.index("_atom_site.auth_asym_id")
                      if "_atom_site.auth_asym_id" in cols else None)
                ec_col = (cols.index("_atom_site.label_entity_id")
                          if "_atom_site.label_entity_id" in cols else None)
                if ac is not None and ec_col is not None:
                    k = j
                    while k < len(lines) and not lines[k].startswith(("#", "loop_", "_")):
                        parts = lines[k].split()
                        if len(parts) == len(cols) and parts[0] in ("ATOM", "HETATM"):
                            chain_entity[parts[ac]] = parts[ec_col]
                        k += 1
            i = j
        else:
            i += 1
    result = {}
    for ch, ent in chain_entity.items():
        if ent in entity_ec:
            result.setdefault(ch, set()).add(entity_ec[ent])
    return result


_raw_cache = {}


def chain_ec_pdbfile(pdb_id, chain_id):
    if pdb_id not in _raw_cache:
        raw = find_raw(pdb_id)
        if raw is None:
            _raw_cache[pdb_id] = {}
        elif ".cif" in raw.name:
            _raw_cache[pdb_id] = parse_cif_ec(raw)
        else:
            _raw_cache[pdb_id] = parse_pdb_compnd(raw)
    return _raw_cache[pdb_id].get(chain_id, set())


def to_l3(ec):
    parts = ec.split(".")
    if len(parts) >= 4 and parts[3] not in ("-", ""):
        return ".".join(parts[:3]) + ".-"
    return ec if ec.endswith("-") else ".".join(parts[:3]) + ".-"


def write_fasta(df, path):
    with open(path, "w") as f:
        for uid, seq in zip(df["unique_id"], df["aa_seq"]):
            f.write(f">{uid}\n{seq}\n")


def main():
    allc = pd.read_csv(OUT_DIR / "design_chains_all.csv")  # D+U only
    sifts = json.load(open(SIFTS_EC))

    # ---- 1. per-chain EC (own annotation only)
    rows = []
    for _, r in allc.iterrows():
        ecs = set(chain_ec_pdbfile(r["pdb_id"], r["chain_id"]))
        ann = sifts.get(f"{r['pdb_id'].upper()}_{r['chain_id']}")
        if ann:
            ecs |= {a["ec_number"] for a in ann.get("annotations", [])
                    if a.get("ec_number")}
        if ecs:
            rows.append({"unique_id": r["unique_id"], "pdb_id": r["pdb_id"],
                         "chain_id": r["chain_id"], "aa_seq": r["aa_seq"],
                         "ec_l4": ";".join(sorted(ecs))})
    ec_df = pd.DataFrame(rows)
    ec_df.to_csv(OUT_DIR / "design_ec_annotations.csv", index=False)
    print(f"per-chain EC hits (D+U): {len(ec_df)}")

    # ---- 2. exact-dup groups: require consistent L3 across members
    ec_df["l3_set"] = ec_df["ec_l4"].map(
        lambda s: frozenset(to_l3(e) for e in s.split(";")))
    grp = ec_df.groupby("aa_seq").agg(
        n=("unique_id", "count"),
        l3_sets=("l3_set", lambda s: set(s)),
        ec_l4_union=("ec_l4", lambda ss: ";".join(sorted(
            {e for s in ss for e in s.split(";")}))),
    )
    conflict = grp[grp["l3_sets"].map(len) > 1]
    if len(conflict):
        conflict.reset_index().to_csv(OUT_DIR / "ec_dropped_conflict.csv",
                                      index=False)
    good = grp[grp["l3_sets"].map(len) == 1]
    print(f"exact-dup groups with EC: {len(grp)}, conflicting (dropped): "
          f"{len(conflict)}, consistent: {len(good)}")

    # map to dedup'd representatives
    dedup = pd.read_csv(OUT_DIR / "design_chains_dedup_exact.csv")
    cand = dedup[dedup["aa_seq"].isin(good.index)].copy()
    cand["ec_l4_set"] = cand["aa_seq"].map(good["ec_l4_union"])
    cand["label"] = cand["aa_seq"].map(
        lambda s: ";".join(sorted(next(iter(good.loc[s, "l3_sets"])))))
    print(f"representatives with own-group EC: {len(cand)}")

    # ---- 3. label filter: all L3 in EC train labels
    train = pd.read_csv(FUNC / "splits" / "ec_train.csv")
    val = pd.read_csv(FUNC / "splits" / "ec_val.csv")
    train_l3 = set()
    for l in train["label"]:
        train_l3.update(str(l).split(";"))
    ok = cand["label"].map(lambda l: set(l.split(";")) <= train_l3)
    print(f"label filter: {len(cand)} -> {ok.sum()}")
    cand = cand[ok].copy()

    # ---- 4. overlap removal vs EC train/val
    tv_pdb = (set(train["unique_id"].str.split("_").str[0].str.lower())
              | set(val["unique_id"].str.split("_").str[0].str.lower()))
    cand["overlap_pdbid"] = cand["pdb_id"].isin(tv_pdb)

    if MMSEQS_M8.exists():
        print(f"mmseqs2: reusing {MMSEQS_M8.name}")
    else:
        q_fa = TMP / "design.fasta"
        t_fa = TMP / "trainval.fasta"
        write_fasta(cand, q_fa)
        write_fasta(pd.concat([train, val], ignore_index=True), t_fa)
        cmd = ["mmseqs", "easy-search", str(q_fa), str(t_fa), str(MMSEQS_M8),
               str(TMP / "mmseqs_tmp"), "-s", "7",
               "--format-output", "query,target,pident"]
        print("mmseqs2:", " ".join(cmd))
        subprocess.run(cmd, check=True)

    hits = pd.read_csv(MMSEQS_M8, sep="\t", header=None,
                       names=["query", "target", "pident"], usecols=[0, 1, 2])
    max_pident = hits.groupby("query")["pident"].max()
    cand["max_pident_tv"] = cand["unique_id"].map(max_pident)
    cand["overlap_seq95"] = cand["max_pident_tv"].fillna(0.0) >= PIDENT_CUTOFF

    removed = cand[cand["overlap_pdbid"] | cand["overlap_seq95"]]
    kept = cand[~(cand["overlap_pdbid"] | cand["overlap_seq95"])].copy()
    removed.to_csv(OUT_DIR / "ec_removed_overlap.csv", index=False)
    print(f"overlap removal: pdb_id {cand['overlap_pdbid'].sum()}, "
          f"seq>=95% {cand['overlap_seq95'].sum()}, removed {len(removed)}")

    kept.to_csv(OUT_DIR / "ec_design_base.csv", index=False)
    print(f"\nec_design_base.csv: {len(kept)} rows")
    print(f"  labels: {sorted(set(';'.join(kept['label']).split(';')))}")
    print(f"  seq_len < 60: {(kept['seq_len'] < 60).sum()}, "
          f"> 1000: {(kept['seq_len'] > 1000).sum()}")


if __name__ == "__main__":
    main()
