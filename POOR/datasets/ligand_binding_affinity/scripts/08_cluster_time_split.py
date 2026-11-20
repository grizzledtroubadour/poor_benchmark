#!/usr/bin/env python3
"""
Script: 08_cluster_time_split.py
Purpose:
  1. Partition samples into cliques with pairwise ECFP4 Tanimoto >= 0.6.
  2. Split each cluster chronologically by release_year in a 7:1:2 ratio.
  3. Move val/test samples with maximum training similarity <0.6 into train.
  4. Balance ratios (incorporating the former 12_balance_val_from_test.py / 13_adjust_val_to_train_9to1.py):
     - After fallback filtering, test is oversized and val undersized; move the earliest third of test to val.
     - Then move the earliest m val samples to train to obtain train:val ≈ 9:1.
       （m = round((9V - T) / 10)）。

Usage:
    cd ligand_binding_affinity/
    python scripts/08_cluster_time_split.py

Input:
    output/pdbbind_v2020_dataset.csv

Output:
    splits/train.csv
    splits/val.csv
    splits/test.csv
    output/cluster_split_stats.json
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import DataStructs

N_BITS = 2048
SIM_THRESHOLD = 0.6


def fp_from_on_bits(on_bits):
    fp = DataStructs.ExplicitBitVect(N_BITS)
    if on_bits:
        fp.SetBitsFromList([int(i) for i in on_bits])
    return fp


def clique_partition(fps, threshold=SIM_THRESHOLD):
    """
    Greedy clique partition: every pair within a cluster has similarity >= threshold.
    Return a list of clusters, each containing sample indices.
    """
    n = len(fps)
    neighbors = [set() for _ in range(n)]

    print(f"[Clustering] Computing {n} x {n} similarity matrix (threshold {threshold})...")
    for i in range(n):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[i + 1:])
        for offset, sim in enumerate(sims):
            if sim >= threshold:
                j = i + 1 + offset
                neighbors[i].add(j)
                neighbors[j].add(i)
        if (i + 1) % 1000 == 0:
            print(f"[Clustering] Processed {i + 1}/{n} samples")

    unassigned = set(range(n))
    clusters = []

    print("[Clustering] Partitioning into cliques...")
    while unassigned:
        seed = min(unassigned)
        clique = [seed]
        candidates = neighbors[seed] & unassigned
        unassigned.remove(seed)

        while candidates:
            v = min(candidates)
            clique.append(v)
            unassigned.remove(v)
            candidates = candidates & neighbors[v]

        clusters.append(clique)

    return clusters


def split_cluster(indices, df):
    """Split one cluster chronologically in a 7:1:2 ratio and return (index, split) pairs."""
    n = len(indices)
    if n == 0:
        return []

    # Sort by release_year, breaking ties stably by pdb_id.
    sub = df.iloc[indices].copy()
    sub = sub.sort_values(by=["release_year", "pdb_id"], kind="mergesort")
    sorted_indices = sub.index.tolist()

    n_train = int(n * 0.7)
    n_val = int(n * 0.1)
    n_test = n - n_train - n_val

    splits = ["train"] * n_train + ["val"] * n_val + ["test"] * n_test
    # splits correspond to sorted_indices in chronological order.
    return list(zip(sorted_indices, splits))


def compute_max_train_tanimoto(query_fps, train_fps):
    """Compute maximum Tanimoto similarity to training for each query."""
    max_sims = []
    for qfp in query_fps:
        sims = DataStructs.BulkTanimotoSimilarity(qfp, train_fps)
        max_sims.append(float(max(sims)) if sims else 0.0)
    return np.array(max_sims)


def main():
    parser = argparse.ArgumentParser(description="Cluster-time split with Tanimoto filter")
    parser.add_argument(
        "--input",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "output" / "pdbbind_v2020_dataset.csv",
        help="Input final dataset CSV",
    )
    parser.add_argument(
        "--splits-dir",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "splits",
        help="Output directory for split CSVs",
    )
    parser.add_argument(
        "--stats",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "output" / "cluster_split_stats.json",
        help="Output split statistics JSON",
    )
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    total = len(df)
    print(f"[Input] Total samples: {total}")

    # Build the fingerprint list.
    fps = []
    valid_mask = []
    for _, row in df.iterrows():
        ecfp4_str = row.get("ligand_ecfp4")
        if pd.isna(ecfp4_str) or not ecfp4_str:
            fps.append(None)
            valid_mask.append(False)
        else:
            on_bits = json.loads(ecfp4_str)
            fps.append(fp_from_on_bits(on_bits))
            valid_mask.append(True)

    valid_indices = np.where(valid_mask)[0].tolist()
    valid_fps = [fps[i] for i in valid_indices]
    print(f"[Input] Valid fingerprints: {len(valid_fps)}")

    # Partition valid fingerprints into cliques.
    clusters = clique_partition(valid_fps, threshold=SIM_THRESHOLD)
    print(f"[Clustering] Generated {len(clusters)} clusters")

    # Cluster size distribution
    sizes = [len(c) for c in clusters]
    size_counter = pd.Series(sizes).value_counts().sort_index()
    print(f"[Clustering] Cluster size distribution (first 10):")
    print(size_counter.head(10))

    # Map global indices to cluster IDs / sizes.
    cluster_id_map = {}
    cluster_size_map = {}
    for cid, c in enumerate(clusters):
        for local_idx in c:
            global_idx = valid_indices[local_idx]
            cluster_id_map[global_idx] = cid
            cluster_size_map[global_idx] = len(c)

    # Invalid fingerprints form individual clusters (marked -1).
    for i in range(total):
        if not valid_mask[i]:
            cluster_id_map[i] = -1
            cluster_size_map[i] = 1

    # Split each cluster chronologically.
    split_labels = np.empty(total, dtype=object)
    for cid, c in enumerate(clusters):
        global_c = [valid_indices[i] for i in c]
        assignments = split_cluster(global_c, df)
        for idx, label in assignments:
            split_labels[idx] = label

    # Assign samples with invalid fingerprints directly to train.
    for i in range(total):
        if not valid_mask[i]:
            split_labels[i] = "train"

    # Step 3: move val/test samples with maximum training similarity <0.6 to train.
    print("[Filtering] Checking maximum val/test Tanimoto similarity to train...")
    moved_total = 0
    iteration = 0
    while True:
        iteration += 1
        train_mask = split_labels == "train"
        val_mask = split_labels == "val"
        test_mask = split_labels == "test"

        train_fps = [fps[i] for i in range(total) if train_mask[i] and valid_mask[i]]
        if not train_fps:
            break

        moved_this_round = 0
        for mask, name in [(val_mask, "val"), (test_mask, "test")]:
            query_indices = np.where(mask)[0].tolist()
            query_fps = [fps[i] for i in query_indices if valid_mask[i]]
            if not query_fps:
                continue
            max_sims = compute_max_train_tanimoto(query_fps, train_fps)
            for local_pos, global_idx in enumerate(query_indices):
                if not valid_mask[global_idx]:
                    # Move invalid fingerprints directly to train.
                    split_labels[global_idx] = "train"
                    moved_this_round += 1
                elif max_sims[local_pos] < SIM_THRESHOLD:
                    split_labels[global_idx] = "train"
                    moved_this_round += 1

        print(f"[Filtering] Round  {iteration}  moved  {moved_this_round} samples to train")
        moved_total += moved_this_round
        if moved_this_round == 0:
            break

    # Step 4: balance ratios (formerly 12_balance_val_from_test.py).
    # After fallback filtering, move the earliest third of the oversized test set to val.
    print("[Balancing] Moving the earliest third of test to val...")
    test_idx = np.where(split_labels == "test")[0]
    test_sub = df.iloc[test_idx].sort_values(by=["release_year", "pdb_id"], kind="mergesort")
    n_move_test_to_val = len(test_sub) // 3
    move_ids = set(test_sub.index[:n_move_test_to_val])
    for i in move_ids:
        split_labels[i] = "val"
    print(f"[Balancing] Moving  {n_move_test_to_val} samples to val "
          f"(train: {(split_labels == 'train').sum()}, val: {(split_labels == 'val').sum()}, "
          f"test: {(split_labels == 'test').sum()})")

    # Step 5: train:val ≈ 9:1 (formerly 13_adjust_val_to_train_9to1.py).
    # (T + m) / (V - m) = 9  =>  m = (9V - T) / 10
    print("[Balancing] Moving the earliest val samples to train for train:val ≈ 9:1...")
    T = int((split_labels == "train").sum())
    V = int((split_labels == "val").sum())
    m = int(round((9.0 * V - T) / 10.0))
    m = max(0, min(m, V))
    val_idx = np.where(split_labels == "val")[0]
    val_sub = df.iloc[val_idx].sort_values(by=["release_year", "pdb_id"], kind="mergesort")
    move_ids = set(val_sub.index[:m])
    for i in move_ids:
        split_labels[i] = "train"
    n_move_val_to_train = m
    print(f"[Balancing] Moving from val  {n_move_val_to_train} samples to train "
          f"(train:val = {(split_labels == 'train').sum() / (split_labels == 'val').sum():.3f} : 1)")

    # Final statistics
    train_df = df[split_labels == "train"].copy()
    val_df = df[split_labels == "val"].copy()
    test_df = df[split_labels == "test"].copy()

    train_df["split"] = "train"
    val_df["split"] = "val"
    test_df["split"] = "test"

    train_df["cluster_id"] = train_df.index.map(cluster_id_map)
    val_df["cluster_id"] = val_df.index.map(cluster_id_map)
    test_df["cluster_id"] = test_df.index.map(cluster_id_map)

    train_df["cluster_size"] = train_df.index.map(cluster_size_map)
    val_df["cluster_size"] = val_df.index.map(cluster_size_map)
    test_df["cluster_size"] = test_df.index.map(cluster_size_map)

    # Compute final max_train_tanimoto
    print("[Fingerprints] Computing final val/test max_train_tanimoto...")
    final_train_fps = [fps[i] for i in train_df.index if valid_mask[i]]
    for subset_df in [val_df, test_df]:
        sims_list = []
        for idx in subset_df.index:
            if valid_mask[idx]:
                sims = DataStructs.BulkTanimotoSimilarity(fps[idx], final_train_fps)
                sims_list.append(float(max(sims)) if sims else 0.0)
            else:
                sims_list.append(np.nan)
        subset_df["max_train_tanimoto"] = sims_list

    train_df["max_train_tanimoto"] = np.nan

    # Save
    args.splits_dir.mkdir(parents=True, exist_ok=True)
    train_path = args.splits_dir / "train.csv"
    val_path = args.splits_dir / "val.csv"
    test_path = args.splits_dir / "test.csv"

    train_df.to_csv(train_path, index=False)
    val_df.to_csv(val_path, index=False)
    test_df.to_csv(test_path, index=False)

    def affinity_counts(d):
        return d["affinity_type"].value_counts().to_dict() if len(d) > 0 else {}

    def describe_series(s):
        s_clean = s.dropna()
        if len(s_clean) == 0:
            return {}
        return {
            "count": int(len(s_clean)),
            "mean": float(s_clean.mean()),
            "std": float(s_clean.std()),
            "min": float(s_clean.min()),
            "25%": float(s_clean.quantile(0.25)),
            "50%": float(s_clean.quantile(0.50)),
            "75%": float(s_clean.quantile(0.75)),
            "max": float(s_clean.max()),
        }

    stats = {
        "total": total,
        "similarity_threshold": SIM_THRESHOLD,
        "clustering": {
            "num_clusters": len(clusters),
            "cluster_size_distribution": {
                "min": int(min(sizes)),
                "max": int(max(sizes)),
                "mean": round(sum(sizes) / len(sizes), 2),
                "singleton_clusters": int(sum(1 for s in sizes if s == 1)),
                "top_sizes": size_counter.head(10).to_dict(),
            },
        },
        "final_counts": {
            "train": len(train_df),
            "val": len(val_df),
            "test": len(test_df),
        },
        "final_ratios": {
            "train": round(len(train_df) / total, 4),
            "val": round(len(val_df) / total, 4),
            "test": round(len(test_df) / total, 4),
        },
        "moved_to_train": moved_total,
        "moved_from_test_to_val": n_move_test_to_val,
        "moved_from_val_to_train": n_move_val_to_train,
        "train_val_ratio": round(len(train_df) / len(val_df), 3),
        "year_range": {
            "train": {"min": int(train_df["release_year"].min()), "max": int(train_df["release_year"].max())},
            "val": {"min": int(val_df["release_year"].min()), "max": int(val_df["release_year"].max())},
            "test": {"min": int(test_df["release_year"].min()), "max": int(test_df["release_year"].max())},
        },
        "affinity_type_counts": {
            "train": affinity_counts(train_df),
            "val": affinity_counts(val_df),
            "test": affinity_counts(test_df),
        },
        "max_train_tanimoto": {
            "val": describe_series(val_df["max_train_tanimoto"]),
            "test": describe_series(test_df["max_train_tanimoto"]),
        },
        "tanimoto_threshold_counts": {
            "val": {
                "<0.6": int((val_df["max_train_tanimoto"] < 0.6).sum()),
                "0.6-0.7": int(((val_df["max_train_tanimoto"] >= 0.6) & (val_df["max_train_tanimoto"] < 0.7)).sum()),
                "0.7-0.8": int(((val_df["max_train_tanimoto"] >= 0.7) & (val_df["max_train_tanimoto"] < 0.8)).sum()),
                "0.8-0.9": int(((val_df["max_train_tanimoto"] >= 0.8) & (val_df["max_train_tanimoto"] < 0.9)).sum()),
                ">=0.9": int((val_df["max_train_tanimoto"] >= 0.9).sum()),
            },
            "test": {
                "<0.6": int((test_df["max_train_tanimoto"] < 0.6).sum()),
                "0.6-0.7": int(((test_df["max_train_tanimoto"] >= 0.6) & (test_df["max_train_tanimoto"] < 0.7)).sum()),
                "0.7-0.8": int(((test_df["max_train_tanimoto"] >= 0.7) & (test_df["max_train_tanimoto"] < 0.8)).sum()),
                "0.8-0.9": int(((test_df["max_train_tanimoto"] >= 0.8) & (test_df["max_train_tanimoto"] < 0.9)).sum()),
                ">=0.9": int((test_df["max_train_tanimoto"] >= 0.9).sum()),
            },
        },
        "split_files": {
            "train": str(train_path),
            "val": str(val_path),
            "test": str(test_path),
        },
    }

    with open(args.stats, "w", encoding="utf-8") as f:
        json.dump(stats, f, indent=2)

    print(f"\n[Output] Training CSV: {train_path} ({len(train_df)})")
    print(f"[Output] Validation CSV: {val_path} ({len(val_df)})")
    print(f"[Output] Test CSV: {test_path} ({len(test_df)})")
    print(f"[Output] Statistics JSON: {args.stats}")
    print(f"\n[Summary] Final ratios: train {stats['final_ratios']['train']}, "
          f"val {stats['final_ratios']['val']}, test {stats['final_ratios']['test']}")
    print(f"[Summary] Moved  {moved_total} samples to train (similarity fallback)")
    print(f"[Summary] Ratio balancing: test→val {n_move_test_to_val}  samples, val→train {n_move_val_to_train}  samples")
    print("\n=== Final max_train_tanimoto distribution ===")
    print("val:")
    print(val_df["max_train_tanimoto"].describe())
    print("test:")
    print(test_df["max_train_tanimoto"].describe())


if __name__ == "__main__":
    main()
