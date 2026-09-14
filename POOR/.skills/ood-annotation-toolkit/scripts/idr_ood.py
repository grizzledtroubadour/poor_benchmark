#!/usr/bin/env python3
"""IDR（内在无序区）OOD 计算（ood-annotation-toolkit）。

统一工具，但**口径存在任务间分裂**，必须显式选 mode：
- `--mode residue-ratio`：disorder>0.5 的残基占比（ss / func / fold 口径，阈值 0.3）
- `--mode region-ratio`：最长连续无序区长度 / 序列长（ppi 口径，阈值 0.3）
阈值用 `--threshold` 调整（affinity / optimal_ph 历史上用 0.1，重算这些任务时对齐）。

序列清洗同样存在分裂，用 `--clean` 选择：
- `strip`（默认）：直接删除非 20 种标准氨基酸字符（fold / ss 口径）
- `map`：B→N, Z→Q, J→L, U→C, O→K，其余非标准字符→A（ppi 口径）

缓存：按清洗后序列增量缓存（--cache，csv：
sequence, residue_ratio, max_region_ratio），两个 mode 的指标同时缓存，
换 mode / 阈值重跑不需重算；缓存逐 batch 落盘，中断可续。

多序列列任务（--seq-col aa_seq1,aa_seq2）：每条链各自计算，
OOD_IDR = 任一链超过阈值（ppi 的 any 语义）。

用法：
    python idr_ood.py --test ssp_test.csv --out test_idr.csv \
        --mode residue-ratio --threshold 0.3 --cache output/idr_cache.csv
    python idr_ood.py --test ppi_test.csv --seq-col aa_seq1,aa_seq2 \
        --mode region-ratio --clean map --out ppi_idr.csv --cache output/idr_cache.csv
"""
import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

STD_AA = "ACDEFGHIKLMNPQRSTVWY"


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--test", required=True, type=Path, help="test csv（含序列列）")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--seq-col", default="aa_seq",
                    help="序列列，逗号分隔多列（双链任务如 aa_seq1,aa_seq2）")
    ap.add_argument("--mode", required=True, choices=["residue-ratio", "region-ratio"],
                    help="residue-ratio=无序残基占比（ss/func/fold）；"
                         "region-ratio=最长连续无序区占比（ppi）")
    ap.add_argument("--threshold", type=float, default=0.3,
                    help="OOD_IDR 判定阈值（默认 0.3；affinity/pH 历史用 0.1）")
    ap.add_argument("--clean", choices=["strip", "map"], default="strip",
                    help="非标准残基清洗方式（strip=fold/ss，map=ppi）")
    ap.add_argument("--cache", type=Path, default=None,
                    help="序列级增量缓存 csv（存在则复用，计算后追加写回）")
    ap.add_argument("--batch-size", type=int, default=512)
    ap.add_argument("--default-col", default="Default",
                    help="Default 列名；存在时 Default=False 行的 OOD_IDR 按 --nondefault 处理")
    ap.add_argument("--nondefault", choices=["nan", "false"], default="nan",
                    help="Default=False 行的取值（nan=全项目统一惯例，流水线位置直接写入；"
                         "false=旧惯例，仅为兼容保留）")
    return ap.parse_args()


def clean_sequence(seq, mode):
    seq = seq.upper()
    if mode == "strip":
        return re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", seq)
    seq = seq.replace("B", "N").replace("Z", "Q").replace("J", "L")
    seq = seq.replace("U", "C").replace("O", "K")
    return re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "A", seq)


def metrics_from_scores(scores):
    """从 per-residue disorder 分数计算两种口径的指标。"""
    arr = np.asarray(scores, dtype=float)
    if len(arr) == 0:
        return np.nan, np.nan
    binary = arr > 0.5
    residue_ratio = float(binary.mean())
    # 最长连续 True 区段长度 / 序列长
    max_run = run = 0
    for b in binary:
        run = run + 1 if b else 0
        if run > max_run:
            max_run = run
    return residue_ratio, max_run / len(arr)


def main():
    args = parse_args()
    seq_cols = [c.strip() for c in args.seq_col.split(",") if c.strip()]
    multi = len(seq_cols) > 1

    test = pd.read_csv(args.test, dtype=str, keep_default_na=False)
    for col in seq_cols:
        if col not in test.columns:
            sys.exit(f"错误：test 缺少序列列 {col!r}")
    print(f"test={len(test)}，seq_cols={seq_cols}，mode={args.mode}，"
          f"threshold={args.threshold}，clean={args.clean}")

    # 收集去重后的清洗序列
    clean_seqs = set()
    for col in seq_cols:
        clean_seqs.update(clean_sequence(s, args.clean) for s in test[col] if s)
    clean_seqs = sorted(clean_seqs)
    print(f"去重清洗序列：{len(clean_seqs)}")

    cache = {}
    if args.cache is not None and args.cache.exists():
        cdf = pd.read_csv(args.cache, dtype=str, keep_default_na=False)
        for r in cdf.itertuples(index=False):
            cache[r.sequence] = (float(r.residue_ratio), float(r.max_region_ratio))
        print(f"缓存命中：{len(cache)}（{args.cache}）")

    missing = [s for s in clean_seqs if s not in cache]
    print(f"待计算：{len(missing)}")
    if missing:
        try:
            import metapredict as meta
        except ImportError:
            sys.exit("错误：当前环境未安装 metapredict（需要 v3）。"
                     "请切换到含 metapredict 的 python（如 /root/miniconda3/bin/python），"
                     "或提供已覆盖这些序列的 --cache 仅做判定。")
        print(f"metapredict version: {meta.__version__}")
        empty = [s for s in missing if len(s) == 0]
        for s in empty:
            cache[s] = (np.nan, np.nan)
        todo = [s for s in missing if len(s) > 0]
        for i in range(0, len(todo), args.batch_size):
            chunk = todo[i:i + args.batch_size]
            scores = meta.predict_disorder_batch(chunk, show_progress_bar=False)
            for seq, sc in zip(chunk, scores):
                cache[seq] = metrics_from_scores(sc[1])
            print(f"  {min(i + args.batch_size, len(todo))}/{len(todo)}")
            if args.cache is not None:  # 逐 batch 落盘，中断可续
                save_cache(cache, args.cache)
        if args.cache is not None:
            save_cache(cache, args.cache)

    metric_idx = 0 if args.mode == "residue-ratio" else 1

    def metric_of(seq):
        if not seq:
            return np.nan
        return cache.get(clean_sequence(seq, args.clean), (np.nan, np.nan))[metric_idx]

    flags = []
    ratio_cols = []
    for col in seq_cols:
        ratio_col = "idr_ratio" if not multi else f"idr_ratio_{col}"
        test[ratio_col] = test[col].map(metric_of)
        ratio_cols.append(ratio_col)
        flags.append(test[ratio_col] > args.threshold)
    ood = flags[0]
    for f in flags[1:]:
        ood = ood | f
    test["OOD_IDR"] = ood.map({True: "True", False: "False"})

    # Default=False（极端长度预挑出）行：OOD_IDR 直接置 NaN（"未计算"语义，
    # 全项目统一惯例，流水线位置写入；finalize_ood_columns.py --mode nan 仅作兜底校验）
    if args.default_col in test.columns:
        nd = test[args.default_col] != "True"
        if nd.any():
            fill = "" if args.nondefault == "nan" else "False"
            test.loc[nd, "OOD_IDR"] = fill
            print(f"  Default=False {int(nd.sum())} 行 OOD_IDR 置 {'NaN' if fill == '' else 'False'}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    test.to_csv(args.out, index=False, lineterminator="\n")
    n_nan = int(sum(test[c].isna().sum() for c in ratio_cols))
    print(f"-> {args.out}（{len(test)} 行）")
    print(f"  OOD_IDR=True: {int(ood.sum())}，指标为 NaN 的链（判 False）: {n_nan}")


def save_cache(cache, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(
        [(s, v[0], v[1]) for s, v in cache.items()],
        columns=["sequence", "residue_ratio", "max_region_ratio"])
    df.to_csv(path, index=False, lineterminator="\n")


if __name__ == "__main__":
    main()
