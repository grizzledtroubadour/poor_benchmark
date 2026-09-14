#!/usr/bin/env python3
"""聚合 SSP 全量数据表（ssp_full_data.csv）。

从历史断链环节补建：划分脚本（08_create_ssp_datasets.py）的输入
``output/ssp_full_data.csv`` 原由未留存脚本生成，本脚本按现行产物重建：

1. 行集与 ``aa_seq``/``file_type`` 取自 ``output/stage4_cluster_representatives.csv``
   （40,067 条代表链；``aa_seq`` 为链全长序列，``file_type`` 为源结构格式）；
2. 标签源为 ``output/dssp_labels.csv`` 的 ``ss8``（逐已解析残基的权威 DSSP
   8-state），映射 ``H0 G1 I2 E3 B4 T5 S6 P7 C8``；
3. ``aa_seq`` 比 DSSP 序列长（DSSP 会丢弃缺骨架原子的残基）：序列一致时直接
   映射；不一致时用 ``difflib.SequenceMatcher``（默认 autojunk=True）对齐——
   equal 块直接映射，replace 块按位置顺序映射，未覆盖的残基置 ``-1``
   （ignore index）。该口径逐行复现了历史文件的全部 40,067 行。

输出：``output/ssp_full_data.csv``，列 ``unique_id, aa_seq, file_type, labels``
（labels 为 ``[8 8 3 ...]`` 风格空格分隔字符串，长度 == len(aa_seq)）。

已知历史偏差：现存文件的 ``file_type`` 有 787 行误记为 ``pdb``（实际源仅
``.cif.gz``，见 stage2_chain_metadata.csv）；该列下游（08 划分脚本）不使用，
本脚本按代表链元数据产出正确值，差异记录在 DATA_PROCESS.md。

用法：
    python scripts/07_aggregate_ssp_full_data.py                    # 写入 output/
    python scripts/07_aggregate_ssp_full_data.py --out output/_verify_ssp/ssp_full_data.csv
"""
import argparse
import difflib
from pathlib import Path

import pandas as pd

TASK_DIR = Path(__file__).resolve().parent.parent
OUTPUT_DIR = TASK_DIR / "output"

DSSP_LABELS = OUTPUT_DIR / "dssp_labels.csv"
REPRESENTATIVES = OUTPUT_DIR / "stage4_cluster_representatives.csv"

SS8_TO_INT = {"H": 0, "G": 1, "I": 2, "E": 3, "B": 4, "T": 5, "S": 6, "P": 7,
              "C": 8, " ": 8}


def build_labels(aa_seq: str, dssp_seq: str, ss8: str) -> list:
    """把 DSSP ss8 映射到全长 aa_seq；DSSP 缺失残基置 -1。"""
    if aa_seq == dssp_seq:
        return [SS8_TO_INT[c] for c in ss8]
    lab = [-1] * len(aa_seq)
    sm = difflib.SequenceMatcher(None, aa_seq, dssp_seq)  # autojunk=True（默认）
    for tag, a0, a1, d0, d1 in sm.get_opcodes():
        if tag == "equal":
            for k in range(a1 - a0):
                lab[a0 + k] = SS8_TO_INT[ss8[d0 + k]]
        elif tag == "replace":
            for k in range(min(a1 - a0, d1 - d0)):
                lab[a0 + k] = SS8_TO_INT[ss8[d0 + k]]
        # delete -> 保持 -1；insert -> DSSP 多余残基忽略
    return lab


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path,
                    default=OUTPUT_DIR / "ssp_full_data.csv",
                    help="输出路径（默认 output/ssp_full_data.csv）")
    args = ap.parse_args()

    dssp = pd.read_csv(DSSP_LABELS)
    rep = pd.read_csv(REPRESENTATIVES)
    # 行序与 dssp_labels.csv（按 chain_key 排序）一致，与历史文件相同
    df = dssp[["chain_key", "seq", "ss8"]].merge(
        rep[["chain_key", "seq", "file_type"]],
        on="chain_key", suffixes=("_dssp", ""),
    )
    assert len(df) == len(rep) == len(dssp)

    n_aligned = 0
    rows = []
    for r in df.itertuples():
        lab = build_labels(r.seq, r.seq_dssp, r.ss8)
        assert len(lab) == len(r.seq)
        n_aligned += r.seq != r.seq_dssp
        rows.append({
            "unique_id": r.chain_key,
            "aa_seq": r.seq,
            "file_type": r.file_type,
            "labels": "[" + " ".join(map(str, lab)) + "]",
        })
    print(f"chains: {len(rows)}, aligned via SequenceMatcher: {n_aligned}")

    out = pd.DataFrame(rows, columns=["unique_id", "aa_seq", "file_type", "labels"])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    print(f"saved: {args.out} ({len(out)} rows)")


if __name__ == "__main__":
    main()
