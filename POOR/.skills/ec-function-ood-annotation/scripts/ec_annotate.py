#!/usr/bin/env python3
"""EC 标注：为样本生成统一的 EC 编号集合（ec-function-ood-annotation skill 第 1 步）。

两种输入模式（按列名自动识别）：
- 链模式：输入含 `pdb`、`chain` 列 → SIFTS 链级 EC（sifts_chain_ec.tsv.gz）
- UniProt 模式：输入含 `uniprot` 列 → ① SIFTS 反查（pdb_chain_uniprot 的
  SP_PRIMARY → 链 → 链级 EC 并集）；② 剩余 accession 走 UniProt REST 批量查询
  （fields=accession,ec），写入项目级共享缓存 data/uniprot_ec_cache.json；
  `--offline` 仅用缓存，不访问 API

输出：输入列之后追加 `ec_numbers`（`;` 分隔完整 EC 编号）、`has_ec`。
层级解析（L3/L4 截断）由 ec_ood.py 完成，本脚本只负责原始 EC 集合。

用法：
    python ec_annotate.py --input my_chains.csv --out my_ec_labels.csv
    python ec_annotate.py --input my_accs.csv --out my_ec_labels.csv --offline
"""
import argparse
import csv
import gzip
import json
import time
from collections import defaultdict
from pathlib import Path

import requests

POOR = Path(__file__).resolve().parents[3]
SIFTS_EC_TSV = POOR / "data" / "sifts" / "sifts_chain_ec.tsv.gz"
SIFTS_PDB_CHAIN_UNIPROT = POOR / "data" / "sifts" / "pdb_chain_uniprot.tsv.gz"
DEFAULT_CACHE = POOR / "data" / "uniprot_ec_cache.json"
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
BATCH_SIZE = 100


def load_sifts_chain_ec():
    m = defaultdict(set)
    with gzip.open(SIFTS_EC_TSV, "rt") as f:
        next(f)
        for line in f:
            p = line.strip().split("\t")
            if len(p) >= 4 and p[3] and p[3] != "-":
                m[(p[0].lower(), p[1])].add(p[3])
    return m


def load_uniprot_to_chains():
    m = defaultdict(list)
    with gzip.open(SIFTS_PDB_CHAIN_UNIPROT, "rt") as f:
        next(f)  # 注释行
        next(f)  # 表头
        for line in f:
            p = line.strip().split("\t")
            if len(p) >= 3:
                m[p[2]].append((p[0].lower(), p[1]))
    return m


def fetch_ec_batch(ids):
    """UniProt REST 批量查询 EC（同 ppi_prediction/scripts/fetch_uniprot_ec.py）。"""
    query = " OR ".join(f"accession:{x}" for x in ids)
    params = {"query": query, "format": "tsv", "fields": "accession,ec", "size": len(ids)}
    r = requests.get(UNIPROT_SEARCH_URL, params=params, timeout=60)
    r.raise_for_status()
    out = {}
    for line in r.text.strip().split("\n")[1:]:
        parts = line.split("\t")
        if len(parts) >= 2:
            out[parts[0]] = [x.strip() for x in parts[1].split(";") if x.strip()]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--offline", action="store_true", help="仅使用缓存，不访问 UniProt API")
    ap.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    args = ap.parse_args()

    with open(args.input) as f:
        rows = list(csv.DictReader(f))
    fieldnames = list(rows[0].keys())
    if {"pdb", "chain"} <= set(fieldnames):
        mode = "chain"
    elif "uniprot" in fieldnames:
        mode = "uniprot"
    else:
        raise SystemExit("无法识别输入模式：需要 (pdb, chain) 列或 uniprot 列")
    print(f"mode: {mode}; input rows: {len(rows)}")

    chain_ec = load_sifts_chain_ec()
    ec_of_row = {}
    if mode == "chain":
        for r in rows:
            ec_of_row[r["pdb"], r["chain"]] = chain_ec.get((r["pdb"].lower(), r["chain"]), set())
        print(f"chains with EC (SIFTS): {sum(bool(v) for v in ec_of_row.values())}/{len(rows)}")
    else:
        accs = sorted({r["uniprot"] for r in rows})
        u2c = load_uniprot_to_chains()
        for a in accs:
            s = set()
            for ck in u2c.get(a, []):
                s |= chain_ec.get(ck, set())
            ec_of_row[a] = s
        covered = sum(bool(ec_of_row[a]) for a in accs)
        print(f"accessions: {len(accs)}, with EC via SIFTS 反查: {covered}")

        # UniProt API 补缺
        cache = json.load(open(args.cache)) if args.cache.exists() else {}
        missing = [a for a in accs if not ec_of_row[a] and a not in cache]
        print(f"cached: {sum(a in cache for a in accs)}, to query: {len(missing)}")
        if missing and not args.offline:
            for i in range(0, len(missing), BATCH_SIZE):
                batch = missing[i:i + BATCH_SIZE]
                try:
                    res = fetch_ec_batch(batch)
                except Exception as e:
                    print(f"  batch {i} 查询失败: {e}")
                    res = {}
                for a in batch:
                    cache[a] = res.get(a, [])
                time.sleep(0.5)
                if (i // BATCH_SIZE) % 10 == 0:
                    print(f"  {i + len(batch)}/{len(missing)}")
            args.cache.parent.mkdir(parents=True, exist_ok=True)
            json.dump(cache, open(args.cache, "w"))
        for a in accs:
            if not ec_of_row[a]:
                ec_of_row[a] = set(cache.get(a, []))
        print(f"with EC after API 补缺: {sum(bool(ec_of_row[a]) for a in accs)}/{len(accs)}")

    if "ec_numbers" not in fieldnames:
        fieldnames += ["ec_numbers", "has_ec"]
    for r in rows:
        key = (r["pdb"], r["chain"]) if mode == "chain" else r["uniprot"]
        ecs = ec_of_row.get(key, set())
        r["ec_numbers"] = ";".join(sorted(ecs))
        r["has_ec"] = str(bool(ecs))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    print(f"-> {args.out}")


if __name__ == "__main__":
    main()
