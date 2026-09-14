#!/usr/bin/env python3
"""Stage selected monomer structures from data/pdb_structures_all/ into pdbs/.

补缺口脚本（2026-08 新增）：pdbs/ 目录历史上由 splits 引用的结构文件从
data/pdb_structures_all/ 复制而来，当时无独立脚本。本脚本固化该步骤：

需要的文件清单 =
    splits/ppi_{train,val,test}.csv 的 struct_file1 / struct_file2 并集
    （struct_file 列由 09_rebuild_with_xtal_negatives.py 按
      output/structure_selection.csv 的选定结构生成）
  ∪ output/splits_pre_xtal_backup/（XTAL 重构前历史 splits）引用文件
    （默认包含，可用 --no-pre-xtal 关闭；现有 pdbs/ 共 46,402 个文件
      = 当前 splits 引用 44,178 ∪ pre-XTAL 历史引用 2,224）

复制规则（幂等）：
  - 目标已存在 → 跳过；
  - 源文件存在 → 复制；
  - 源文件缺失 → 记入 missing 清单（当前 data/pdb_structures_all/ 仅保留
    10,690 个文件，多数历史文件已移出；全新环境重放需先运行
    06_download_structures_and_filter_test.py 重新下载全部 50,306 个选定结构）。

用法：
  python 07_stage_pdbs.py                      # 幂等 staging 到 pdbs/
  python 07_stage_pdbs.py --dry-run            # 只统计期望清单与缺失，不写文件
  python 07_stage_pdbs.py --dest DIR [--limit N]  # 重定向目标（验证用）
"""
import argparse
import json
import logging
import shutil
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

TASK_DIR = Path(__file__).resolve().parent.parent  # datasets/ppi_prediction
SPLITS_DIR = TASK_DIR / "splits"
PRE_XTAL_DIR = TASK_DIR / "output" / "splits_pre_xtal_backup"
SRC_DIR = TASK_DIR / "data" / "pdb_structures_all"
DEST_DIR = TASK_DIR / "pdbs"
SUMMARY_JSON = TASK_DIR / "output" / "stage_pdbs_summary.json"

STRUCT_COLS = ["struct_file1", "struct_file2"]


def collect_split_files(csv_path):
    df = pd.read_csv(csv_path, usecols=lambda c: c in STRUCT_COLS)
    files = set()
    for col in df.columns:
        files.update(f for f in df[col].dropna().astype(str) if f and f != "nan")
    return files


def collect_needed_files(include_pre_xtal):
    needed = set()
    for split in ["train", "val", "test"]:
        path = SPLITS_DIR / f"ppi_{split}.csv"
        files = collect_split_files(path)
        logger.info(f"{path.name}: {len(files)} unique structure files")
        needed.update(files)
    n_current = len(needed)
    if include_pre_xtal and PRE_XTAL_DIR.exists():
        for path in sorted(PRE_XTAL_DIR.glob("*.csv")):
            files = collect_split_files(path)
            logger.info(f"pre-XTAL {path.name}: {len(files)} unique structure files")
            needed.update(files)
        logger.info(
            f"needed total: {len(needed)} (current splits {n_current}, "
            f"pre-XTAL extra {len(needed) - n_current})"
        )
    return needed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--src", type=str, default=str(SRC_DIR), help="源结构目录")
    parser.add_argument("--dest", type=str, default=str(DEST_DIR), help="目标目录（默认 pdbs/）")
    parser.add_argument("--no-pre-xtal", action="store_true",
                        help="不包含 pre-XTAL 历史 splits 引用的文件")
    parser.add_argument("--dry-run", action="store_true", help="只统计，不复制")
    parser.add_argument("--limit", type=int, default=0, help="最多复制 N 个文件（验证用，0=不限）")
    args = parser.parse_args()

    src = Path(args.src)
    dest = Path(args.dest)
    needed = collect_needed_files(include_pre_xtal=not args.no_pre_xtal)

    existing = {f.name for f in dest.glob("*.pdb")} if dest.exists() else set()
    to_copy = sorted(needed - existing)
    missing_src = [f for f in to_copy if not (src / f).exists()]
    copyable = [f for f in to_copy if (src / f).exists()]
    if args.limit:
        copyable = copyable[: args.limit]

    logger.info(
        f"needed={len(needed)}, already in dest={len(needed & existing)}, "
        f"to copy={len(to_copy)} (src available {len(to_copy) - len(missing_src)}, "
        f"src missing {len(missing_src)})"
    )

    n_copied = 0
    if not args.dry_run:
        dest.mkdir(parents=True, exist_ok=True)
        for f in copyable:
            shutil.copy2(src / f, dest / f)
            n_copied += 1
        logger.info(f"copied {n_copied} files -> {dest}")
    else:
        logger.info("dry-run: no files written")

    summary = {
        "needed_files": len(needed),
        "already_in_dest": len(needed & existing),
        "copied": n_copied,
        "src_missing": len(missing_src),
        "src_missing_files": missing_src,
    }
    out = SUMMARY_JSON if dest == DEST_DIR else dest / "stage_pdbs_summary.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as fh:
        json.dump(summary, fh, indent=2)
    logger.info(f"summary -> {out}")


if __name__ == "__main__":
    main()
