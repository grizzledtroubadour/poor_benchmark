#!/usr/bin/env python3
"""
Stage 4: Global sequence clustering with MMseqs2 at 95% identity.

Single stage (v3, 2026-09-02): run easy-linclust once on all stage2 chains
(length >=31, including short 31–59-aa chains), supplying FASTA directly. mmseqs 18.8cc5c
easy-linclust accepts FASTA; prebuilt database inputs trigger createdb errors in this version.
Choose the longest representative per cluster, breaking ties lexicographically by chain_key.

Reproduction note: published data use the historical two-stage route: cluster >=60-aa and
<60-aa pools separately, cross-deduplicate, then remove <=30-aa representatives (see history
in DATA_PROCESS.md). Under longest-first selection, the single-stage route is equivalent
except for linclust boundary tie-breaking; small differences are allowed.

Inputs:
    ss_prediction/output/stage2_chain_sequences.fa
    ss_prediction/output/stage2_chain_metadata.csv

Outputs:
    ss_prediction/output/stage4_mmseqs_clusters.tsv
    ss_prediction/output/stage4_cluster_representatives.csv
    ss_prediction/output/stage4_nonredundant_sequences.fa
    ss_prediction/output/stage4_clustering_stats.json

Usage:
    python ss_prediction/scripts/04_run_mmseqs2_clustering.py
"""

import json
import logging
import subprocess
import pandas as pd
from pathlib import Path
from collections import defaultdict

TASK_DIR = Path(__file__).resolve().parent.parent  # datasets/ss_prediction
OUTPUT_DIR = TASK_DIR / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(OUTPUT_DIR / "mmseqs2_clustering.log", mode="w"),
        logging.StreamHandler(),
    ],
)
logger = logging.getLogger(__name__)

# Config
INPUT_FASTA = OUTPUT_DIR / "stage2_chain_sequences.fa"
METADATA_CSV = OUTPUT_DIR / "stage2_chain_metadata.csv"
TMP_DIR = OUTPUT_DIR / "mmseqs2_tmp"
CLUSTER_TSV = OUTPUT_DIR / "stage4_mmseqs_clusters.tsv"
REP_CSV = OUTPUT_DIR / "stage4_cluster_representatives.csv"
REP_FASTA = OUTPUT_DIR / "stage4_nonredundant_sequences.fa"
STATS_JSON = OUTPUT_DIR / "stage4_clustering_stats.json"

MIN_SEQ_ID = 0.95
COV = 0.8


def run_mmseqs2():
    """Run MMseqs2 easy-linclust on all stage2 chains (len>=31)."""
    TMP_DIR.mkdir(parents=True, exist_ok=True)
    out_prefix = TMP_DIR / "clust"

    logger.info("Running easy-linclust at 95% identity on all chains...")
    subprocess.run(
        [
            "mmseqs", "easy-linclust",
            str(INPUT_FASTA), str(out_prefix), str(TMP_DIR),
            "--min-seq-id", str(MIN_SEQ_ID),
            "-c", str(COV),
            "--cov-mode", "0",
            "--threads", "32",
        ],
        check=True,
    )

    import shutil
    shutil.copy(str(out_prefix) + "_cluster.tsv", CLUSTER_TSV)
    logger.info(f"Saved cluster TSV to {CLUSTER_TSV}")


def parse_clusters(tsv_path) -> dict:
    """Parse MMseqs2 cluster TSV: rep -> member."""
    clusters = defaultdict(list)
    with open(tsv_path) as f:
        for line in f:
            rep, member = line.strip().split("\t")
            clusters[rep].append(member)
    return clusters


def select_representatives(clusters: dict, meta_df: pd.DataFrame) -> pd.DataFrame:
    """Select the best representative from each cluster.

    Rule: longest sequence per cluster, then lexicographic chain_key."""
    meta_df = meta_df.copy()
    meta_df["chain_key"] = meta_df["pdb_id"] + "_" + meta_df["chain_id"]
    lookup = meta_df.set_index("chain_key")

    reps = []
    for rep_key, members in clusters.items():
        candidates = [lookup.loc[m] for m in members if m in lookup.index]
        if not candidates:
            continue
        best = max(candidates, key=lambda row: (row["num_residues"], row.name))
        reps.append(best)

    rep_df = pd.DataFrame(reps)
    rep_df = rep_df.reset_index().rename(columns={"index": "chain_key"})
    return rep_df


def main():
    if not INPUT_FASTA.exists():
        logger.error(f"Input FASTA not found: {INPUT_FASTA}")
        return

    logger.info("Loading metadata...")
    meta_df = pd.read_csv(METADATA_CSV)
    meta_df["num_residues"] = meta_df["num_residues"].astype(int)

    if not CLUSTER_TSV.exists():
        run_mmseqs2()
    else:
        logger.info(f"Reusing existing cluster file: {CLUSTER_TSV}")

    logger.info("Parsing clusters...")
    clusters = parse_clusters(CLUSTER_TSV)
    logger.info(f"Total clusters: {len(clusters)}")
    singletons = sum(1 for m in clusters.values() if len(m) == 1)
    multi = len(clusters) - singletons
    logger.info(f"Singleton clusters: {singletons}, Multi-member: {multi}")

    logger.info("Selecting cluster representatives...")
    rep_df = select_representatives(clusters, meta_df)
    logger.info(f"Selected {len(rep_df)} representative chains")

    rep_df.to_csv(REP_CSV, index=False)
    logger.info(f"Saved {len(rep_df)} representatives to {REP_CSV}")

    rep_keys = set(rep_df["chain_key"].tolist())
    with open(REP_FASTA, "w") as fout:
        with open(INPUT_FASTA) as fin:
            write_seq = False
            for line in fin:
                if line.startswith(">"):
                    header = line[1:].strip().split()[0]
                    write_seq = header in rep_keys
                    if write_seq:
                        fout.write(line)
                elif write_seq:
                    fout.write(line)
    logger.info(f"Saved representative FASTA to {REP_FASTA}")

    stats = {
        "total_chains": len(meta_df),
        "clusters": len(clusters),
        "singleton_clusters": int(singletons),
        "multi_member_clusters": int(multi),
        "representative_chains": len(rep_df),
        "min_seq_id": MIN_SEQ_ID,
        "coverage_threshold": COV,
        "reduction_rate": round((len(meta_df) - len(rep_df)) / len(meta_df) * 100, 2),
    }
    with open(STATS_JSON, "w") as f:
        json.dump(stats, f, indent=2)
    logger.info(f"Saved stats to {STATS_JSON}")

    print("\n=== Clustering Summary ===")
    for k, v in stats.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
