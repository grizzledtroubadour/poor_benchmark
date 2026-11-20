#!/usr/bin/env python3
"""Repair struct_label for structure-residue/aa_seq length mismatches using contiguous resseq anchors.

These structures omit residues; sequence alignment in add_struct_label.py can place gaps
arbitrarily in low-complexity regions and misalign interface labels (DATA_PROCESS.md §5.7).
Use retained PDB resseq numbering for deterministic mapping:
            - Enumerate first occurrences of (chain,resseq,icode); process contiguous segments separately.
            - Anchor each segment to a unique exact match in aa_seq (sequence identity is already 100%).
                        Map unique anchors in resseq order; ambiguous or absent matches receive -1.
            - Across segments, require monotonic non-overlapping aa_seq positions; otherwise reject the row.
This removes arbitrary alignment-gap placement.

Usage: python3 14_fix_struct_label_resseq.py [--apply]
"""
import os, sys
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
TD = os.path.join(ROOT, "datasets", "ppis_prediction")
SPLITS = ["ppis_train.csv", "ppis_val.csv", "ppis_test.csv"]

RES3 = {
    'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I',
    'LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S',
    'THR':'T','VAL':'V','TRP':'W','TYR':'Y','MSE':'M','SEC':'C','PYL':'K',
}

def parse_pdb(path):
    res=[]; seen=set()
    for l in open(path, errors="ignore"):
        if not l.startswith(("ATOM","HETATM")): continue
        key=(l[21], l[22:26], l[26])
        if key in seen: continue
        seen.add(key)
        res.append((int(l[22:26]), RES3.get(l[17:20].strip().upper(), 'X')))
    return res

def fix_row(residues, aa_seq, labels):
    """Return (new_struct_label, stats), or (None, reason)."""
    n = len(residues)
    # Split into contiguous resseq segments.
    runs=[]; cur=[0]
    for i in range(1, n):
        if residues[i][0] == residues[i-1][0] + 1:
            cur.append(i)
        else:
            runs.append(cur); cur=[i]
    runs.append(cur)
    pos_map = {}
    unanchored = 0
    for run in runs:
        sub = "".join(residues[i][1] for i in run)
        hits=[]; start=0
        while True:
            j = aa_seq.find(sub, start)
            if j==-1: break
            hits.append(j); start=j+1
        if len(hits)==1:
            j0=hits[0]
            for k,i in enumerate(run):
                pos_map[i]=j0+k
        else:
            unanchored += len(run)
    # Check monotonicity.
    mapped = sorted(pos_map.values())
    if len(mapped)!=len(set(mapped)):
        return None, "overlap"
    new=[labels[pos_map[i]] if i in pos_map else -1 for i in range(n)]
    return new, {"unanchored": unanchored, "runs": len(runs)}

def main():
    apply = "--apply" in sys.argv
    for f in SPLITS:
        df = pd.read_csv(os.path.join(TD,"splits",f), dtype=str, keep_default_na=False)
        fixed=0
        for idx, r in df.iterrows():
            residues = parse_pdb(os.path.join(TD,"pdbs",r["struct_file"]))
            n=len(residues)
            cur=[int(x) for x in r["struct_label"].strip("[]").split()]
            if len(cur)!=n or n==len(r["aa_seq"]):
                continue  # Direct/equal-length paths are safe.
            labels=[int(x) for x in r["label"].strip("[]").split()]
            new, st = fix_row(residues, r["aa_seq"], labels)
            if new is None:
                print(f"  [FAIL] {r['unique_id'][:12]} ({f}): {st}"); continue
            if new!=cur:
                d = np.array(new)!=np.array(cur)
                d1 = int((d & ((np.array(new)==1)|(np.array(cur)==1))).sum())
                print(f"  [FIX] {r['unique_id'][:12]} ({f}): {int(d.sum())} pos changed, {d1} at interface, unanchored={st['unanchored']}")
                if apply:
                    df.at[idx,"struct_label"]="[ "+" ".join(map(str,new))+" ]"
                fixed+=1
        if apply and fixed:
            df.to_csv(os.path.join(TD,"splits",f), index=False, lineterminator="\n")
        print(f"{f}: {fixed} rows fixed")
if __name__=="__main__":
    main()
