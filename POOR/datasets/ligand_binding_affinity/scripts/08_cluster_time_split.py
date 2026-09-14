#!/usr/bin/env python3
"""
脚本：08_cluster_time_split.py
用途：
  1. 按 ECFP4 Tanimoto >= 0.6 对样本进行团划分聚类（同一簇内两两相似度均 >= 0.6）
  2. 每个簇按 release_year 进行 7:1:2 时间划分
  3. 检查 val/test 中是否存在与训练集最大相似度 < 0.6 的样本，将其移入 train
  4. 比例平衡（已并入原 12_balance_val_from_test.py / 13_adjust_val_to_train_9to1.py）：
     - 兜底过滤后 test 占比偏高、val 偏低，先将 test 中按时间最早的 1/3 样本移入 val；
     - 再将 val 中按时间最早的 m 个样本移入 train，使 train:val ≈ 9:1
       （m = round((9V - T) / 10)）。

使用：
    cd ligand_binding_affinity/
    python scripts/08_cluster_time_split.py

输入：
    output/pdbbind_v2020_dataset.csv

输出：
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
    贪心团划分：将样本划分成若干簇，保证每个簇内任意两个样本相似度 >= threshold。
    返回簇列表，每个簇是样本索引列表。
    """
    n = len(fps)
    neighbors = [set() for _ in range(n)]

    print(f"[聚类] 计算 {n} x {n} 相似度矩阵（阈值 {threshold}）...")
    for i in range(n):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[i + 1:])
        for offset, sim in enumerate(sims):
            if sim >= threshold:
                j = i + 1 + offset
                neighbors[i].add(j)
                neighbors[j].add(i)
        if (i + 1) % 1000 == 0:
            print(f"[聚类] 已处理 {i + 1}/{n} 个样本")

    unassigned = set(range(n))
    clusters = []

    print("[聚类] 进行团划分...")
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
    """对单个簇按时间 7:1:2 切分，返回每个索引对应的 split 标签。"""
    n = len(indices)
    if n == 0:
        return []

    # 按 release_year 升序，同年份按 pdb_id 稳定排序
    sub = df.iloc[indices].copy()
    sub = sub.sort_values(by=["release_year", "pdb_id"], kind="mergesort")
    sorted_indices = sub.index.tolist()

    n_train = int(n * 0.7)
    n_val = int(n * 0.1)
    n_test = n - n_train - n_val

    splits = ["train"] * n_train + ["val"] * n_val + ["test"] * n_test
    # splits 按时间顺序对应 sorted_indices
    return list(zip(sorted_indices, splits))


def compute_max_train_tanimoto(query_fps, train_fps):
    """计算 query 列表与训练集的最大 Tanimoto 相似度数组。"""
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
    print(f"[输入] 总样本数: {total}")

    # 构建指纹列表
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
    print(f"[输入] 有效指纹数: {len(valid_fps)}")

    # 团划分聚类（仅在有效指纹上进行）
    clusters = clique_partition(valid_fps, threshold=SIM_THRESHOLD)
    print(f"[聚类] 共生成 {len(clusters)} 个簇")

    # 簇大小分布
    sizes = [len(c) for c in clusters]
    size_counter = pd.Series(sizes).value_counts().sort_index()
    print(f"[聚类] 簇大小分布（前 10）:")
    print(size_counter.head(10))

    # 建立全局索引 -> 簇ID / 簇大小
    cluster_id_map = {}
    cluster_size_map = {}
    for cid, c in enumerate(clusters):
        for local_idx in c:
            global_idx = valid_indices[local_idx]
            cluster_id_map[global_idx] = cid
            cluster_size_map[global_idx] = len(c)

    # 无效指纹样本单独成簇（标记为 -1）
    for i in range(total):
        if not valid_mask[i]:
            cluster_id_map[i] = -1
            cluster_size_map[i] = 1

    # 对每个簇进行时间划分
    split_labels = np.empty(total, dtype=object)
    for cid, c in enumerate(clusters):
        global_c = [valid_indices[i] for i in c]
        assignments = split_cluster(global_c, df)
        for idx, label in assignments:
            split_labels[idx] = label

    # 无效指纹样本直接划入 train
    for i in range(total):
        if not valid_mask[i]:
            split_labels[i] = "train"

    # 第三步：移动 val/test 中与训练集最大相似度 < 0.6 的样本到 train
    print("[过滤] 检查 val/test 与 train 的最大 Tanimoto...")
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
                    # 无效指纹直接移入 train
                    split_labels[global_idx] = "train"
                    moved_this_round += 1
                elif max_sims[local_pos] < SIM_THRESHOLD:
                    split_labels[global_idx] = "train"
                    moved_this_round += 1

        print(f"[过滤] 第 {iteration} 轮移动 {moved_this_round} 个样本到 train")
        moved_total += moved_this_round
        if moved_this_round == 0:
            break

    # 第四步：比例平衡（原 12_balance_val_from_test.py）
    # 兜底过滤后 test 占比偏高、val 偏低，将 test 按时间最早的 1/3 移到 val
    print("[平衡] 将 test 中按时间最早的 1/3 样本移入 val...")
    test_idx = np.where(split_labels == "test")[0]
    test_sub = df.iloc[test_idx].sort_values(by=["release_year", "pdb_id"], kind="mergesort")
    n_move_test_to_val = len(test_sub) // 3
    move_ids = set(test_sub.index[:n_move_test_to_val])
    for i in move_ids:
        split_labels[i] = "val"
    print(f"[平衡] 从 test 移动 {n_move_test_to_val} 个样本到 val "
          f"(train: {(split_labels == 'train').sum()}, val: {(split_labels == 'val').sum()}, "
          f"test: {(split_labels == 'test').sum()})")

    # 第五步：train:val ≈ 9:1（原 13_adjust_val_to_train_9to1.py）
    # (T + m) / (V - m) = 9  =>  m = (9V - T) / 10
    print("[平衡] 将 val 中按时间最早的样本移入 train，使 train:val ≈ 9:1...")
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
    print(f"[平衡] 从 val 移动 {n_move_val_to_train} 个样本到 train "
          f"(train:val = {(split_labels == 'train').sum() / (split_labels == 'val').sum():.3f} : 1)")

    # 最终统计
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

    # 计算最终 max_train_tanimoto
    print("[指纹] 计算最终 val/test 的 max_train_tanimoto...")
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

    # 保存
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

    print(f"\n[输出] 训练集 CSV: {train_path} ({len(train_df)})")
    print(f"[输出] 验证集 CSV: {val_path} ({len(val_df)})")
    print(f"[输出] 测试集 CSV: {test_path} ({len(test_df)})")
    print(f"[输出] 统计 JSON: {args.stats}")
    print(f"\n[汇总] 最终比例: train {stats['final_ratios']['train']}, "
          f"val {stats['final_ratios']['val']}, test {stats['final_ratios']['test']}")
    print(f"[汇总] 共移动 {moved_total} 个样本到 train（相似度兜底）")
    print(f"[汇总] 比例平衡：test→val {n_move_test_to_val} 个，val→train {n_move_val_to_train} 个")
    print("\n=== 最终 max_train_tanimoto 分布 ===")
    print("val:")
    print(val_df["max_train_tanimoto"].describe())
    print("test:")
    print(test_df["max_train_tanimoto"].describe())


if __name__ == "__main__":
    main()
