#!/usr/bin/env python3
"""Crop ppis structures to the representative chain matching aa_seq.

For each sample where the structure file does not already contain a chain
exactly equal to aa_seq, this script:
  1. Selects the chain best matching aa_seq.
  2. Globally aligns aa_seq to that chain and keeps exactly one residue per
     aa_seq position, skipping structural insertions.  This guarantees that
     after cropping: len(structure residues) == len(aa_seq) == len(labels).
  3. If the alignment identity is too low or aa_seq contains deletions with
     respect to the chain, the sample is flagged for removal because the
     structure is corrupted/placeholder.

With --apply:
  - Recoverable structures are overwritten in pdbs/; originals are moved to
    output/quarantine_structures/.
  - Unrecoverable rows are removed from the split CSVs; originals are backed
    up to output/quarantine_splits/.

Without --apply: report only.
"""
import os, sys, re, warnings, shutil
from difflib import SequenceMatcher
warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from Bio.Data.PDBData import protein_letters_3to1_extended as MAP3TO1
import pandas as pd


def one_letter(resname):
    return MAP3TO1.get(resname.strip().upper(), None)

TD = os.path.join(ROOT, "datasets", "ppis_prediction")
FILES = ["ppis_train.csv", "ppis_val.csv", "ppis_test.csv"]

# Minimum fraction of aa_seq residues that must match the chosen chain.
IDENTITY_THR = 0.95


def is_coord(t):
    if "." not in t:
        return False
    try:
        float(t)
        return True
    except ValueError:
        return False


def parse_residues(path):
    """Return ordered residues: list of {"lines": [...], "resname": ..., "ol": ...}."""
    residues = []
    cur_desc = None
    cur = None
    with open(path, errors="ignore") as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            p = line.split()
            resname = p[3]
            ci = None
            for i in range(4, len(p)):
                if is_coord(p[i]):
                    ci = i
                    break
            desc = tuple(p[4:ci]) if ci else (p[4],)
            if desc != cur_desc:
                cur_desc = desc
                cur = {"desc": desc, "resname": resname, "lines": []}
                residues.append(cur)
            cur["lines"].append(line)
    for r in residues:
        r["ol"] = one_letter(r["resname"])
    return residues


def chain_sequences(residues):
    """Map chain_id -> (sequence, list of all residue indices in residues)."""
    chain_idx = {}
    for i, r in enumerate(residues):
        desc = r["desc"]
        cid = desc[0][0] if desc and desc[0] else "_"
        chain_idx.setdefault(cid, []).append(i)
    out = {}
    for cid, idxs in chain_idx.items():
        seq = "".join(residues[i]["ol"] for i in idxs if residues[i]["ol"] is not None)
        out[cid] = (seq, idxs)
    return out


def map_query_to_target(target_seq, query_seq):
    """Global-align query (aa_seq) to target (chain). Return dict query_pos->target_pos, identity.

    query positions that align to a gap in target (structural deletions) are
    absent from the dict.
    """
    from Bio.Align import PairwiseAligner

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    if not target_seq or not query_seq:
        return {}, 0.0
    alignment = aligner.align(target_seq, query_seq)[0]
    coord = alignment.aligned  # shape (2, n_blocks, 2)
    q2t = {}
    matches = 0
    aligned_len = 0
    for k in range(coord.shape[1]):
        t0, t1 = int(coord[0][k][0]), int(coord[0][k][1])
        q0, q1 = int(coord[1][k][0]), int(coord[1][k][1])
        L = t1 - t0
        aligned_len += L
        for i in range(L):
            qp = q0 + i
            tp = t0 + i
            q2t[qp] = tp
            if target_seq[tp] == query_seq[qp]:
                matches += 1
    identity = matches / len(query_seq) if query_seq else 0.0
    return q2t, identity


def decide_crop(aa, residues):
    """Return (residue_indices_to_keep, how, best_chain_id, new_seq, identity) or None."""
    aa = aa.replace("|", "").strip()
    chains = chain_sequences(residues)
    if not chains:
        return None

    # 1) exact chain match
    for cid, (sseq, idxs) in chains.items():
        if sseq == aa:
            return idxs, "EXACT_CHAIN", cid, sseq, 1.0

    # 2) exact substring of a chain
    for cid, (sseq, idxs) in chains.items():
        pos = sseq.find(aa)
        if pos >= 0:
            keep = idxs[pos : pos + len(aa)]
            new_seq = "".join(residues[i]["ol"] for i in keep)
            return keep, "SUBSTRING", cid, new_seq, 1.0

    # 3) global alignment per chain, pick best identity
    best = None
    for cid, (sseq, idxs) in chains.items():
        q2t, identity = map_query_to_target(sseq, aa)
        if len(q2t) != len(aa):
            # structural deletions -> cannot represent every aa position
            continue
        if identity < IDENTITY_THR:
            continue
        keep = [idxs[q2t[q]] for q in range(len(aa))]
        new_seq = "".join(residues[i]["ol"] for i in keep)
        if best is None or identity > best[4]:
            best = (keep, "GLOBAL_ALIGN", cid, new_seq, identity)
    return best


def crop_and_write(path, keep_idxs, residues):
    """Overwrite path with cropped PDB containing kept residues."""
    kept_lines = []
    for i in keep_idxs:
        kept_lines.extend(residues[i]["lines"])
    with open(path, "w") as f:
        f.write("HEADER    CROPPED TO REPRESENTATIVE CHAIN\n")
        f.writelines(kept_lines)
        f.write("END\n")


def main():
    apply = "--apply" in sys.argv
    restrict = set(a for a in sys.argv[1:] if not a.startswith("--"))

    quarantine_struct = os.path.join(TD, "output", "quarantine_structures")
    quarantine_splits = os.path.join(TD, "output", "quarantine_splits")
    if apply:
        os.makedirs(quarantine_struct, exist_ok=True)
        os.makedirs(quarantine_splits, exist_ok=True)

    rows = []
    for fn in FILES:
        fp = os.path.join(TD, "splits", fn)
        df = pd.read_csv(fp)
        for _, r in df.iterrows():
            rows.append((fn, str(r["unique_id"]), str(r["aa_seq"])))
    if restrict:
        rows = [r for r in rows if r[1] in restrict]

    stats = {
        "checked": 0,
        "already_exact": 0,
        "cropped": 0,
        "unrecoverable": 0,
        "no_structure": 0,
    }
    recoverable = []   # (fn, uid, aa_len, old_struct_len, new_struct_len, how, ident, cid)
    unrecoverable = [] # (fn, uid, aa_len, best_struct_len, best_ident, best_cid, reason)

    for fn, uid, aa in rows:
        p = os.path.join(TD, "pdbs", f"{uid}.pdb")
        stats["checked"] += 1
        if not os.path.exists(p):
            stats["no_structure"] += 1
            unrecoverable.append((fn, uid, len(aa), 0, 0.0, None, "missing"))
            continue
        residues = parse_residues(p)
        result = decide_crop(aa, residues)
        if result is None:
            # gather diagnostic info
            chains = chain_sequences(residues)
            best_ident = 0.0
            best_cid = None
            best_len = 0
            reason = "no_chain"
            for cid, (sseq, idxs) in chains.items():
                q2t, ident = map_query_to_target(sseq, aa)
                if ident > best_ident:
                    best_ident = ident
                    best_cid = cid
                    best_len = len(sseq)
                    if len(q2t) != len(aa):
                        reason = "deletions_in_structure"
                    else:
                        reason = f"low_identity_{ident:.3f}"
            unrecoverable.append((fn, uid, len(aa), best_len, round(best_ident, 4), best_cid, reason))
            continue
        keep_idxs, how, cid, new_seq, identity = result
        if how == "EXACT_CHAIN":
            stats["already_exact"] += 1
            continue
        stats["cropped"] += 1
        recoverable.append((fn, uid, len(aa), len(residues), len(keep_idxs), how, round(identity, 4), cid))
        if apply:
            shutil.copy(p, os.path.join(quarantine_struct, f"{uid}.pdb"))
            crop_and_write(p, keep_idxs, residues)

    print("==== PPIS structure crop report ====")
    print(f"checked:          {stats['checked']}")
    print(f"already exact:    {stats['already_exact']}")
    print(f"cropped:          {stats['cropped']}")
    print(f"unrecoverable:    {len(unrecoverable)}")
    print(f"no structure:     {stats['no_structure']}")

    from collections import Counter
    print("\nPer-split unrecoverable:", Counter(u[0] for u in unrecoverable))
    print("Per-split cropped:", Counter(r[0] for r in recoverable))

    print("\nSample unrecoverable (lowest identity first):")
    for u in sorted(unrecoverable, key=lambda x: x[4])[:30]:
        print(f"  {u[1]:34} {u[0]:18} aa_len={u[2]:5} struct_len={u[3]:5} ident={u[4]:.3f} chain={u[5]} reason={u[6]}")
    if len(unrecoverable) > 30:
        print(f"  ... (+{len(unrecoverable)-30} more)")

    print("\nSample cropped:")
    for r in recoverable[:20]:
        print(f"  {r[1]:34} {r[0]:18} aa_len={r[2]:5} old_struct={r[3]:5} new_struct={r[4]:5} {r[5]:14} ident={r[6]:.3f} chain={r[7]}")
    if len(recoverable) > 20:
        print(f"  ... (+{len(recoverable)-20} more)")

    if apply:
        removed = 0
        for fn in FILES:
            fp = os.path.join(TD, "splits", fn)
            df = pd.read_csv(fp)
            bad_uids = {u[1] for u in unrecoverable if u[0] == fn}
            if bad_uids:
                shutil.copy(fp, os.path.join(quarantine_splits, fn))
                df = df[~df["unique_id"].isin(bad_uids)].reset_index(drop=True)
                df.to_csv(fp, index=False)
                removed += len(bad_uids)
        print(f"\n[apply] removed {removed} unrecoverable rows from splits")
        print(f"[apply] originals quarantined to {quarantine_struct} and {quarantine_splits}")
    else:
        print("\n[DRY-RUN] no files modified. Use --apply to execute.")


if __name__ == "__main__":
    main()
