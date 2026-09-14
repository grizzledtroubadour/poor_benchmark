#!/usr/bin/env python3
"""
脚本：03_clean_metadata.py
用途：对 PDBbind v2020.R1 元数据进行清洗过滤。

使用：
    cd binding_affinity/
    python scripts/03_clean_metadata.py

默认过滤规则（可通过命令行参数开关）：
    --exclude-ic50                排除 IC50 样本
    --exclude-inequality          排除亲和力带不等号（< / <= / ~）的样本
    --exclude-incomplete-ligand   排除 [Incomplete ligand] 标记样本
    --exclude-covalent            排除 [Covalent complex] 标记样本
    --exclude-different-protein   排除 [Different protein in assay] 标记样本
    --exclude-different-ligand    排除 [Different ligand in assay] 标记样本
    --exclude-uncommon-element    排除 [Uncommon element] 标记样本（默认保留）

输入：
    output/pdbbind_v2020_metadata.csv

输出：
    output/pdbbind_v2020_metadata_cleaned.csv
    output/pdbbind_v2020_metadata_cleaned_stats.json
"""

import argparse
import json
import pandas as pd
from pathlib import Path


def clean_metadata(
    input_csv: Path,
    output_csv: Path,
    stats_path: Path,
    exclude_ic50: bool,
    exclude_inequality: bool,
    exclude_incomplete_ligand: bool,
    exclude_covalent: bool,
    exclude_different_protein: bool,
    exclude_different_ligand: bool,
    exclude_uncommon_element: bool,
):
    df = pd.read_csv(input_csv)
    total = len(df)

    stats = {
        "total_before": total,
        "filters": {},
    }

    # 依次应用过滤规则，记录每步移除数量
    # 注意：使用副本，避免修改原始 df
    cleaned = df.copy()

    if exclude_ic50:
        mask = cleaned["affinity_type"] == "IC50"
        removed = mask.sum()
        cleaned = cleaned[~mask]
        stats["filters"]["exclude_ic50"] = {
            "removed": int(removed),
            "remaining": int(len(cleaned)),
        }

    if exclude_inequality:
        # 只保留精确值 '='
        mask = cleaned["affinity_operator"] != "="
        removed = mask.sum()
        cleaned = cleaned[~mask]
        stats["filters"]["exclude_inequality"] = {
            "removed": int(removed),
            "remaining": int(len(cleaned)),
        }

    if exclude_incomplete_ligand:
        mask = cleaned["incomplete_ligand"].notna() & (cleaned["incomplete_ligand"] != "")
        removed = mask.sum()
        cleaned = cleaned[~mask]
        stats["filters"]["exclude_incomplete_ligand"] = {
            "removed": int(removed),
            "remaining": int(len(cleaned)),
        }

    if exclude_covalent:
        mask = cleaned["covalent_complex"].notna() & (cleaned["covalent_complex"] != "")
        removed = mask.sum()
        cleaned = cleaned[~mask]
        stats["filters"]["exclude_covalent"] = {
            "removed": int(removed),
            "remaining": int(len(cleaned)),
        }

    if exclude_different_protein:
        mask = cleaned["different_protein"].notna() & (cleaned["different_protein"] != "")
        removed = mask.sum()
        cleaned = cleaned[~mask]
        stats["filters"]["exclude_different_protein"] = {
            "removed": int(removed),
            "remaining": int(len(cleaned)),
        }

    if exclude_different_ligand:
        mask = cleaned["different_ligand"].notna() & (cleaned["different_ligand"] != "")
        removed = mask.sum()
        cleaned = cleaned[~mask]
        stats["filters"]["exclude_different_ligand"] = {
            "removed": int(removed),
            "remaining": int(len(cleaned)),
        }

    if exclude_uncommon_element:
        mask = cleaned["uncommon_element"].notna() & (cleaned["uncommon_element"] != "")
        removed = mask.sum()
        cleaned = cleaned[~mask]
        stats["filters"]["exclude_uncommon_element"] = {
            "removed": int(removed),
            "remaining": int(len(cleaned)),
        }

    # 最终统计
    stats["total_after"] = int(len(cleaned))
    stats["removed_total"] = int(total - len(cleaned))
    stats["retention_rate"] = round(len(cleaned) / total, 4) if total > 0 else 0.0

    # 亲和力类型分布
    stats["affinity_type_counts_after"] = cleaned["affinity_type"].value_counts().to_dict()

    # 年份范围
    stats["release_year_range"] = {
        "min": int(cleaned["release_year"].min()),
        "max": int(cleaned["release_year"].max()),
    }

    # pK 统计
    stats["pK_stats"] = cleaned["pK"].describe().to_dict()

    # 保存
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    cleaned.to_csv(output_csv, index=False)

    with open(stats_path, "w", encoding="utf-8") as f:
        json.dump(stats, f, indent=2)

    print(f"[输入] 清洗前样本数: {total}")
    print(f"[输出] 清洗后样本数: {len(cleaned)}")
    print(f"[输出] 保留比例: {stats['retention_rate']*100:.2f}%")
    print(f"[输出] 清洗后 CSV: {output_csv}")
    print(f"[输出] 统计 JSON: {stats_path}")
    print()
    print("=== 各过滤步骤 ===")
    for name, s in stats["filters"].items():
        print(f"  {name}: 移除 {s['removed']}，剩余 {s['remaining']}")
    print()
    print(f"=== 最终亲和力类型分布 ===")
    for k, v in stats["affinity_type_counts_after"].items():
        print(f"  {k}: {v}")


def main():
    parser = argparse.ArgumentParser(description="Clean PDBbind v2020.R1 metadata")
    parser.add_argument(
        "--input",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "output" / "pdbbind_v2020_metadata.csv",
        help="Input metadata CSV",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "output" / "pdbbind_v2020_metadata_cleaned.csv",
        help="Output cleaned metadata CSV",
    )
    parser.add_argument(
        "--stats",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "output" / "pdbbind_v2020_metadata_cleaned_stats.json",
        help="Output stats JSON",
    )

    # 过滤开关
    parser.add_argument("--exclude-ic50", action="store_true", default=True, help="Exclude IC50 samples")
    parser.add_argument("--exclude-inequality", action="store_true", default=True, help="Exclude samples with inequality operators")
    parser.add_argument("--exclude-incomplete-ligand", action="store_true", default=True, help="Exclude [Incomplete ligand] samples")
    parser.add_argument("--exclude-covalent", action="store_true", default=True, help="Exclude [Covalent complex] samples")
    parser.add_argument("--exclude-different-protein", action="store_true", default=True, help="Exclude [Different protein in assay] samples")
    parser.add_argument("--exclude-different-ligand", action="store_true", default=True, help="Exclude [Different ligand in assay] samples")
    parser.add_argument("--exclude-uncommon-element", action="store_true", default=False, help="Exclude [Uncommon element] samples")

    # 允许用 --no-exclude-xxx 关闭默认过滤
    parser.add_argument("--no-exclude-ic50", dest="exclude_ic50", action="store_false")
    parser.add_argument("--no-exclude-inequality", dest="exclude_inequality", action="store_false")
    parser.add_argument("--no-exclude-incomplete-ligand", dest="exclude_incomplete_ligand", action="store_false")
    parser.add_argument("--no-exclude-covalent", dest="exclude_covalent", action="store_false")
    parser.add_argument("--no-exclude-different-protein", dest="exclude_different_protein", action="store_false")
    parser.add_argument("--no-exclude-different-ligand", dest="exclude_different_ligand", action="store_false")
    parser.add_argument("--no-exclude-uncommon-element", dest="exclude_uncommon_element", action="store_false")

    args = parser.parse_args()

    clean_metadata(
        args.input,
        args.output,
        args.stats,
        args.exclude_ic50,
        args.exclude_inequality,
        args.exclude_incomplete_ligand,
        args.exclude_covalent,
        args.exclude_different_protein,
        args.exclude_different_ligand,
        args.exclude_uncommon_element,
    )


if __name__ == "__main__":
    main()
