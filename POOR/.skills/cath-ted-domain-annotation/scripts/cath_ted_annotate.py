#!/usr/bin/env python3
"""CATH+TED domain annotation for PDB chains OR UniProt accessions (task-agnostic).

Two input modes (auto-detected from CSV columns):
  chain mode     input has `pdb`,`chain` columns.
                 CATH = name match against cath-domain-list.txt; TED domains are
                 kept only if >= --overlap of the domain falls within the chain's
                 UniProt range (SIFTS segment mapping).
  uniprot mode   input has a `uniprot` column (sample = full-length protein,
                 e.g. AFDB-keyed datasets). CATH = union over all PDB chains mapped
                 to the accession in SIFTS (pdb_chain_cath + sifts_chain_uniprot,
                 domains resolved by name match against cath-domain-list.txt);
                 TED labels are taken as-is (no overlap filter).

Chain-mode output columns: cath_topos, cath_sfs, ted_topos, ted_sfs, topos, sfs,
uniprot, annotation_source, cath_domains (id:uni-ranges), ted_domains,
cath_chain_cov, aug_chain_cov.
Uniprot-mode output columns: cath_topos, cath_sfs, ted_topos, ted_sfs, topos, sfs,
annotation_source, cath_domains (bare ids), ted_domains.
topos/sfs = union (CATH priority + TED gap-filling) -- use these downstream.
annotation_source: cath | ted | both | none.

TED label levels: H-level -> C.A.T + C.A.T.H ; T-level -> C.A.T only.
Comma-separated double labels (foldseek,foldclass calls) are parsed as a union.

TED API results are cached project-wide (default: data/ted/ted_api_cache.json).

Usage:
  python cath_ted_annotate.py --chains my_chains.csv --out my_labels.csv
  python cath_ted_annotate.py --chains my_accs.csv --out my_labels.csv --offline
"""
import argparse
import csv
import gzip
import json
import threading
import time
import urllib.error
import urllib.request
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
API = "https://ted.cathdb.info/api/v1/uniprot/summary/{}"

rate_lock = threading.Lock()
last_req = [0.0]


def fetch_ted(acc):
    with rate_lock:
        wait = 0.09 - (time.time() - last_req[0])
        if wait > 0:
            time.sleep(wait)
        last_req[0] = time.time()
    for attempt in range(3):
        try:
            with urllib.request.urlopen(API.format(acc), timeout=30) as r:
                js = json.loads(r.read().decode())
            data = js.get("data", [])
            return {"status": "ok" if data else "not_found", "data": data}
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return {"status": "not_found", "data": []}
            time.sleep(2 ** attempt * 2)
        except Exception:
            time.sleep(2 ** attempt * 2)
    return {"status": "error", "data": []}


def parse_range_str(s, sep="_"):
    """'2-21_29-69' or '2-21,29-69' -> [(2,21),(29,69)]"""
    segs = []
    for part in s.replace(",", sep).split(sep):
        a, b = part.split("-")
        segs.append((int(a), int(b)))
    return segs


def merge_ranges(ranges):
    if not ranges:
        return []
    rs = sorted(ranges)
    merged = [list(rs[0])]
    for b, e in rs[1:]:
        if b <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([b, e])
    return [(b, e) for b, e in merged]


def overlap_len(segs_a, segs_b):
    return sum(max(0, min(e1, e2) - max(b1, b2) + 1)
               for b1, e1 in segs_a for b2, e2 in segs_b)


def total_len(segs):
    return sum(e - b + 1 for b, e in segs)


def fmt_ranges(segs):
    return ",".join(f"{b}-{e}" for b, e in segs)


def ted_labels_from_record(rec):
    """All labeled domains of a TED record -> (topos, sfs, [(ted_id, label, chopping)])."""
    tt, ts, doms = set(), set(), []
    if not rec or rec.get("status") != "ok":
        return tt, ts, doms
    for d in rec["data"]:
        label = d.get("cath_label", "-")
        if label == "-":
            continue
        doms.append((d["ted_id"], label, d.get("chopping", "")))
        for alt in label.split(","):
            parts = alt.split(".")
            if len(parts) < 3:
                continue
            tt.add(".".join(parts[:3]))
            if d.get("cath_assignment_level") == "H" and len(parts) >= 4:
                ts.add(".".join(parts[:4]))
    return tt, ts, doms


def query_ted(accs, cache, args):
    todo = [a for a in accs if cache.get(a, {}).get("status") in (None, "error")]
    print(f"unique UniProt: {len(accs)}, cached: {len(accs) - len(todo)}, to query: {len(todo)}")
    if todo and args.offline:
        print("offline mode: skipping API queries")
    elif todo:
        lock = threading.Lock()
        done = [0]
        t0 = time.time()

        def work(acc):
            res = fetch_ted(acc)
            with lock:
                cache[acc] = res
                done[0] += 1
                if done[0] % 500 == 0:
                    args.cache.parent.mkdir(parents=True, exist_ok=True)
                    json.dump(cache, open(args.cache, "w"))
                    dt = time.time() - t0
                    print(f"  {done[0]}/{len(todo)}  {done[0]/dt:.1f} req/s", flush=True)

        args.cache.parent.mkdir(parents=True, exist_ok=True)
        with ThreadPoolExecutor(max_workers=args.threads) as ex:
            list(ex.map(work, todo))
        json.dump(cache, open(args.cache, "w"))
    return cache


def load_cath_labels_by_domain(cath_list):
    """domain_id -> (topo, sf)"""
    dom_label = {}
    with open(cath_list) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.split()
            if len(p) < 5:
                continue
            topo = f"{p[1]}.{p[2]}.{p[3]}"
            dom_label[p[0].lower()] = (topo, f"{topo}.{p[4]}")
    return dom_label


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chains", type=Path, required=True,
                    help="CSV with pdb,chain columns OR a uniprot column "
                         "(extra columns preserved)")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--domains-out", type=Path, default=None,
                    help="optional per-domain CSV (chain mode only)")
    ap.add_argument("--cache", type=Path, default=ROOT / "data/ted/ted_api_cache.json")
    ap.add_argument("--overlap", type=float, default=0.5,
                    help="chain mode: min fraction of TED domain covered by the chain "
                         "UniProt range (default 0.5)")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--offline", action="store_true", help="use cache only, no API calls")
    ap.add_argument("--cath-list", type=Path,
                    default=ROOT / "data/CATHv44/cath-classification-data/cath-domain-list.txt")
    ap.add_argument("--cath-boundaries", type=Path,
                    default=ROOT / "data/CATHv44/cath-classification-data/cath-domain-boundaries-seqreschopping.txt")
    ap.add_argument("--sifts-dir", type=Path, default=ROOT / "data/sifts")
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.chains)))
    mode = "uniprot" if ("uniprot" in rows[0] and "pdb" not in rows[0]) else "chain"
    print(f"mode: {mode}; input rows: {len(rows)}")

    cache = json.load(open(args.cache)) if args.cache.exists() else {}

    if mode == "uniprot":
        out_rows = run_uniprot_mode(rows, cache, args)
        extra_cols = ["cath_topos", "cath_sfs", "ted_topos", "ted_sfs", "topos", "sfs",
                      "annotation_source", "cath_domains", "ted_domains"]
    else:
        out_rows, dom_rows = run_chain_mode(rows, cache, args)
        extra_cols = ["cath_topos", "cath_sfs", "ted_topos", "ted_sfs", "topos", "sfs",
                      "uniprot", "annotation_source", "cath_domains", "ted_domains",
                      "cath_chain_cov", "aug_chain_cov"]

    fieldnames = [c for c in list(rows[0].keys()) + extra_cols
                  if c in out_rows[0] or c in rows[0]]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(out_rows)

    if mode == "chain" and args.domains_out:
        args.domains_out.parent.mkdir(parents=True, exist_ok=True)
        with open(args.domains_out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["source", "domain_id", "pdb", "chain",
                                              "uniprot", "label", "uni_ranges",
                                              "mapped_frac", "in_chain_frac", "kept"])
            w.writeheader()
            w.writerows(dom_rows)
        print(f"domains written: {len(dom_rows)}  ->  {args.domains_out}")

    from collections import Counter
    stat = Counter(r["annotation_source"] for r in out_rows)
    cov = sum(1 for r in out_rows if r["topos"]) / len(out_rows) * 100
    print(f"annotation_source: {dict(stat)}")
    print(f"union topo coverage: {cov:.1f}%  ->  {args.out}")


# ---------------------------------------------------------------- chain mode
def run_chain_mode(rows, cache, args):
    for r in rows:
        r["pdb"] = r["pdb"].lower()
    needed = {(r["pdb"], r["chain"]) for r in rows}
    needed_pdbs = {p for p, _ in needed}
    print(f"unique (pdb,chain): {len(needed)}")

    # CATH: chain -> domain ids with labels
    dom_label = {}
    cath_doms = defaultdict(list)
    with open(args.cath_list) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.split()
            if len(p) < 5:
                continue
            key = (p[0][:4].lower(), p[0][4])
            if key not in needed:
                continue
            topo = f"{p[1]}.{p[2]}.{p[3]}"
            dom = p[0].lower()
            dom_label[dom] = (topo, f"{topo}.{p[4]}")
            cath_doms[key].append(dom)
    print(f"chains with CATH: {len(cath_doms)}, domains: {len(dom_label)}")

    # CATH domain boundaries (SEQRES numbering)
    dom_chop = {}
    with open(args.cath_boundaries) as f:
        for line in f:
            p = line.split()
            if len(p) == 2 and p[0].lower() in dom_label:
                try:
                    dom_chop[p[0].lower()] = parse_range_str(p[1], sep=",")
                except ValueError:
                    pass
    print(f"CATH domains with boundaries: {len(dom_chop)}/{len(dom_label)}")

    # SIFTS: uniprot accession + segment mapping (needed pdbs only)
    chain2uni = {}
    chain_segs = defaultdict(list)
    with gzip.open(args.sifts_dir / "pdb_chain_cath.tsv.gz", "rt") as f:
        for line in f:
            if line.startswith("#") or line.startswith("PDB"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) >= 4 and p[0] in needed_pdbs:
                chain2uni[(p[0], p[1])] = p[2]
    with gzip.open(args.sifts_dir / "sifts_chain_uniprot.tsv.gz", "rt") as f:
        for line in f:
            if line.startswith("#") or line.startswith("PDB"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[0] not in needed_pdbs:
                continue
            key = (p[0], p[1])
            chain2uni.setdefault(key, p[2])
            try:
                chain_segs[key].append((int(p[3]), int(p[4]), int(p[7]), int(p[8])))
            except ValueError:
                pass
    print(f"chains without UniProt mapping: {sum(1 for c in needed if c not in chain2uni)}")

    accs = sorted({chain2uni[c] for c in needed if c in chain2uni})
    cache = query_ted(accs, cache, args)

    out_rows, dom_rows = [], []
    for r in rows:
        key = (r["pdb"], r["chain"])
        uni = chain2uni.get(key)
        uni_ranges = merge_ranges([(sb, se) for _, _, sb, se in chain_segs.get(key, [])])
        chain_len = total_len(uni_ranges)

        ct, cs = set(), set()
        cath_dom_strs, cath_cov_ranges = [], []
        for dom in cath_doms.get(key, []):
            topo, sf = dom_label[dom]
            segs = dom_chop.get(dom)
            uni_segs = []
            if segs:
                for a, b in segs:
                    for rb, re, sb, se in chain_segs.get(key, []):
                        lo, hi = max(a, rb), min(b, re)
                        if lo <= hi:
                            uni_segs.append((sb + lo - rb, sb + hi - rb))
                uni_segs = merge_ranges(uni_segs)
            mapped_frac = (total_len(uni_segs) / total_len(segs)) if segs else None
            in_chain_frac = (overlap_len(uni_segs, uni_ranges) / total_len(segs)) if segs else None
            ct.add(topo)
            cs.add(sf)
            cath_dom_strs.append(f"{dom}:{fmt_ranges(uni_segs)}" if uni_segs else f"{dom}:")
            cath_cov_ranges.extend(uni_segs)
            if args.domains_out:
                dom_rows.append({"source": "cath", "domain_id": dom, "pdb": key[0],
                                 "chain": key[1], "uniprot": uni or "", "label": sf,
                                 "uni_ranges": fmt_ranges(uni_segs),
                                 "mapped_frac": f"{mapped_frac:.3f}" if mapped_frac is not None else "",
                                 "in_chain_frac": f"{in_chain_frac:.3f}" if in_chain_frac is not None else "",
                                 "kept": "True"})

        tt, ts = set(), set()
        ted_dom_strs, ted_cov_ranges = [], []
        rec = cache.get(uni or "", {})
        if uni and rec.get("status") == "ok" and uni_ranges:
            level_by_id = {d["ted_id"]: d.get("cath_assignment_level") for d in rec["data"]}
            _, _, tdoms = ted_labels_from_record(rec)
            for ted_id, label, chopping in tdoms:
                try:
                    segs = parse_range_str(chopping)
                except Exception:
                    continue
                dom_len = total_len(segs)
                ov_frac = overlap_len(segs, uni_ranges) / dom_len if dom_len else 0
                kept = ov_frac >= args.overlap
                if args.domains_out:
                    dom_rows.append({"source": "ted", "domain_id": ted_id,
                                     "pdb": key[0], "chain": key[1], "uniprot": uni,
                                     "label": label, "uni_ranges": chopping,
                                     "mapped_frac": "", "in_chain_frac": f"{ov_frac:.3f}",
                                     "kept": str(kept)})
                if not kept:
                    continue
                for alt in label.split(","):
                    parts = alt.split(".")
                    if len(parts) < 3:
                        continue
                    tt.add(".".join(parts[:3]))
                    if level_by_id.get(ted_id) == "H" and len(parts) >= 4:
                        ts.add(".".join(parts[:4]))
                ted_dom_strs.append(f"{ted_id}:{chopping}")
                ted_cov_ranges.extend(segs)

        topos, sfs = ct | tt, cs | ts
        source = ("both" if ct and tt else "cath" if ct else
                  "ted" if tt else "none")
        cath_cov = overlap_len(merge_ranges(cath_cov_ranges), uni_ranges) / chain_len if chain_len else 0
        aug_cov = overlap_len(merge_ranges(cath_cov_ranges + ted_cov_ranges),
                              uni_ranges) / chain_len if chain_len else 0
        out_rows.append({**r,
                         "cath_topos": ";".join(sorted(ct)), "cath_sfs": ";".join(sorted(cs)),
                         "ted_topos": ";".join(sorted(tt)), "ted_sfs": ";".join(sorted(ts)),
                         "topos": ";".join(sorted(topos)), "sfs": ";".join(sorted(sfs)),
                         "uniprot": uni or "", "annotation_source": source,
                         "cath_domains": ";".join(cath_dom_strs),
                         "ted_domains": ";".join(ted_dom_strs),
                         "cath_chain_cov": f"{cath_cov:.3f}",
                         "aug_chain_cov": f"{aug_cov:.3f}"})
    return out_rows, dom_rows


# -------------------------------------------------------------- uniprot mode
def run_uniprot_mode(rows, cache, args):
    accs_needed = sorted({r["uniprot"] for r in rows if r["uniprot"]})
    acc_set = set(accs_needed)

    # SIFTS reverse map: accession -> PDB chains (both files)
    acc2chains = defaultdict(set)
    with gzip.open(args.sifts_dir / "pdb_chain_cath.tsv.gz", "rt") as f:
        for line in f:
            if line.startswith("#") or line.startswith("PDB"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) >= 4 and p[2] in acc_set:
                acc2chains[p[2]].add((p[0], p[1]))
    with gzip.open(args.sifts_dir / "sifts_chain_uniprot.tsv.gz", "rt") as f:
        for line in f:
            if line.startswith("#") or line.startswith("PDB"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3 and p[2] in acc_set:
                acc2chains[p[2]].add((p[0], p[1]))
    rev_chains = set().union(*acc2chains.values()) if acc2chains else set()
    print(f"accessions: {len(accs_needed)}, with PDB chains: {len(acc2chains)}, "
          f"reverse chains: {len(rev_chains)}")

    # CATH labels for those chains (name match against domain list)
    chain_labels = defaultdict(lambda: [set(), set(), []])  # chain -> topos, sfs, dom ids
    with open(args.cath_list) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.split()
            if len(p) < 5:
                continue
            key = (p[0][:4].lower(), p[0][4])
            if key not in rev_chains:
                continue
            topo = f"{p[1]}.{p[2]}.{p[3]}"
            chain_labels[key][0].add(topo)
            chain_labels[key][1].add(f"{topo}.{p[4]}")
            chain_labels[key][2].append(p[0].lower())

    acc_cath = {}
    for acc, chains in acc2chains.items():
        tt, ss, doms = set(), set(), []
        for c in chains:
            lt, ls, ld = chain_labels.get(c, [set(), set(), []])
            tt |= lt
            ss |= ls
            doms += ld
        acc_cath[acc] = (tt, ss, doms)
    n_cath = sum(1 for a in accs_needed if acc_cath.get(a, (set(),))[0])
    print(f"accessions with CATH (reverse): {n_cath}")

    cache = query_ted(accs_needed, cache, args)

    out_rows = []
    for r in rows:
        acc = r["uniprot"]
        ct, cs, cdoms = acc_cath.get(acc, (set(), set(), []))
        tt, ts, tdoms = ted_labels_from_record(cache.get(acc, {}))
        topos, sfs = ct | tt, cs | ts
        source = ("both" if ct and tt else "cath" if ct else
                  "ted" if tt else "none")
        out_rows.append({**r,
                         "cath_topos": ";".join(sorted(ct)), "cath_sfs": ";".join(sorted(cs)),
                         "ted_topos": ";".join(sorted(tt)), "ted_sfs": ";".join(sorted(ts)),
                         "topos": ";".join(sorted(topos)), "sfs": ";".join(sorted(sfs)),
                         "annotation_source": source,
                         "cath_domains": ";".join(sorted(set(cdoms))),
                         "ted_domains": ";".join(f"{tid}:{ch}" for tid, _, ch in tdoms)})
    return out_rows


if __name__ == "__main__":
    main()
