#!/usr/bin/env python3
"""构建 fold_classification 划分样本池（domain pool）。

从历史断链环节补建：v2 划分脚本（04_build_topology_split.py）的样本池原直接
取自 v1 C.A 粒度划分备份（output/splits_ca_backup/），本脚本从一手产物重建
等价的样本池：

1. 输入 ``output/cath_s95_summary.csv``（62,901 个成功提取的 S95 域）；
2. ``aa_seq`` 取自 CATH 参考序列 ``cath-domain-seqs-S95.fa``（域级权威序列，
   与 v1 池完全一致；PDB ATOM 记录常缺失 loop 残基，不作序列源）；
3. 剔除 ``length_expected > 1000`` 的超长域（5 个：1ej6B00、1u6gC00、
   3a6pA00、3egwA02、3gjxD00；口径与 v1 池一致——按 CATH 列表长度判定，
   池中保留 3w3uA01 这类参考序列略超 1000 但 length_expected=999 的域）；
4. ``struct_file = {domain_id}.pdb``，并断言对应文件存在于 ``pdbs/``。

输出：``output/domain_pool.csv``（62,896 行，列 unique_id, aa_seq, struct_file）。
与 ``output/splits_ca_backup/`` 三个 v1 划分文件的并集做行级 diff 验证等价。

用法：
    python scripts/03_build_domain_pool.py                      # 写入 output/domain_pool.csv
    python scripts/03_build_domain_pool.py --out output/_verify_pool/domain_pool.csv
"""
import argparse
from pathlib import Path

import pandas as pd

TASK_DIR = Path(__file__).resolve().parent.parent
ROOT = TASK_DIR.parent.parent                       # 仓库根（POOR）
SUMMARY_CSV = TASK_DIR / "output" / "cath_s95_summary.csv"
PDB_DIR = TASK_DIR / "pdbs"
CATH_SEQ_S95 = (
    ROOT / "data" / "CATHv44" / "sequence-data" / "cath-domain-seqs-S95.fa"
)

MAX_LEN = 1000  # length_expected > 1000 的域剔除


def load_fasta(fasta_file: Path) -> dict:
    """加载 CATH 域参考序列，键为 domain_id（如 12asA00）"""
    seqs = {}
    current_id = None
    current_seq = []
    with open(fasta_file) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_id is not None:
                    seqs[current_id] = "".join(current_seq)
                header = line[1:].strip()
                current_id = header.split("|")[-1].split("/")[0].strip()
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_id is not None:
            seqs[current_id] = "".join(current_seq)
    return seqs


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path,
                    default=TASK_DIR / "output" / "domain_pool.csv",
                    help="输出路径（默认 output/domain_pool.csv）")
    args = ap.parse_args()

    summary = pd.read_csv(SUMMARY_CSV)
    fasta = load_fasta(CATH_SEQ_S95)
    print(f"summary 行数: {len(summary)}, FASTA 条目: {len(fasta)}")

    dropped = summary[summary["length_expected"] > MAX_LEN]
    pool = summary[summary["length_expected"] <= MAX_LEN].copy()
    print(f"剔除 length_expected>{MAX_LEN}: {len(dropped)} "
          f"({', '.join(dropped['domain_id'])})")

    pool["unique_id"] = pool["domain_id"]
    pool["aa_seq"] = pool["unique_id"].map(fasta)
    assert pool["aa_seq"].notna().all(), "域缺失 CATH 参考序列"
    pool["struct_file"] = pool["unique_id"] + ".pdb"

    missing_pdb = [u for u in pool["unique_id"]
                   if not (PDB_DIR / f"{u}.pdb").exists()]
    assert not missing_pdb, f"pdbs/ 缺文件: {missing_pdb[:5]}..."

    out = pool[["unique_id", "aa_seq", "struct_file"]].sort_values("unique_id")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    print(f"样本池已保存: {args.out} ({len(out)} 行)")


if __name__ == "__main__":
    main()
