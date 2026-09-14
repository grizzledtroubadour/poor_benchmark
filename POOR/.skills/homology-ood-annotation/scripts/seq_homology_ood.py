#!/usr/bin/env python3
"""序列同源 OOD 计算（homology-ood-annotation skill 第 1 步）。

统一口径（以 fold_classification/scripts/annotate_test_ood.py 为基准）：
- 搜索方向：query = test 序列，target = train+val 序列
- 工具：mmseqs2 `easy-search -s 7.5`（默认 e-value <= 1e-3 门槛，
  lbs 修复版用 -s 7，差异见 SKILL.md 历史实现表）
- 判定：每个 test 样本取对 ref 的 max pident；无 hit 按 0.0 计
  seq_Redundancy_<T> = (max pident < T)；OOD_Orphan = 无任何 hit
  → 保证 OOD_Orphan ⊆ seq_Redundancy_30 ⊆ ... ⊆ seq_Redundancy_90
- 多序列列任务（如 ppi 的 aa_seq1/aa_seq2）：每条链各自搜索，
  样本级判定默认按 any 语义合并（任一链满足即 True，与 ppi
  10_recompute_all_ood.py 的 (A<thr)|(B<thr)、(A==0)|(B==0) 一致）

m8 缓存：--m8-cache 指向的 m8 已存在时仅重放判定（无需 mmseqs2）；
不存在时运行 mmseqs2 并把 m8 写到该路径供下次复用。
m8 格式固定为 `query,target,pident` 三列（query 名 = unique_id；
多列模式为 `unique_id|<seq-col>`，自行生产的缓存才能复用）。

用法：
    # 单链任务（fold 口径）
    python seq_homology_ood.py --test cath_test.csv \
        --ref cath_train.csv cath_val.csv --out test_flagged.csv \
        --m8-cache ood_mmseqs_test_vs_trainval.m8
    # 双链任务（ppi 口径）
    python seq_homology_ood.py --test ppi_test.csv \
        --ref ppi_train.csv ppi_val.csv --seq-col aa_seq1,aa_seq2 \
        --out ppi_flagged.csv
    # 仅重放已有 m8（不需要 mmseqs2）
    python seq_homology_ood.py --test cath_test.csv --ref cath_train.csv \
        --out replay.csv --m8-cache existing.m8
"""
import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

DEFAULT_THRESHOLDS = "90,80,70,60,50,40,30"


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--test", required=True, type=Path, help="test csv（含序列列）")
    ap.add_argument("--ref", required=True, type=Path, nargs="+",
                    help="参考集 csv（train、val，可多个，取并集）")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--seq-col", default="aa_seq",
                    help="序列列，逗号分隔多列（双链任务如 aa_seq1,aa_seq2）")
    ap.add_argument("--id-col", default="unique_id")
    ap.add_argument("--thresholds", default=DEFAULT_THRESHOLDS,
                    help="逗号分隔的 pident 阈值（max pident < 阈值 → True）")
    ap.add_argument("--merge", choices=["any", "all"], default="any",
                    help="多序列列任务的样本级合并语义（any=任一链满足，ppi 口径）")
    ap.add_argument("--sensitivity", default="7.5",
                    help="mmseqs2 easy-search -s 参数（fold 用 7.5，lbs 修复版用 7）")
    ap.add_argument("--m8-cache", type=Path, default=None,
                    help="m8 缓存路径：存在则仅重放；不存在则搜索后写入复用")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--tmp-dir", type=Path, default=None,
                    help="临时目录（fasta、mmseqs tmp）；默认自动创建并清理")
    ap.add_argument("--default-col", default="Default",
                    help="Default 列名；存在时 Default=False 行的结果列按 --nondefault 处理")
    ap.add_argument("--nondefault", choices=["nan", "false"], default="nan",
                    help="Default=False 行的取值（nan=全项目统一惯例，流水线位置直接写入；"
                         "false=旧惯例，仅为兼容保留）")
    return ap.parse_args()


def read_csv(p):
    return pd.read_csv(p, dtype=str, keep_default_na=False)


def query_key(uid, col, multi):
    return f"{uid}|{col}" if multi else uid


def collect_records(df, seq_cols, id_col, multi):
    """去重收集 (key, seq)；同 key 取首条。"""
    records = {}
    cols = [id_col] + seq_cols
    for row in df[cols].itertuples(index=False):
        uid, seqs = row[0], row[1:]
        for col, seq in zip(seq_cols, seqs):
            if not seq:
                continue
            key = query_key(uid, col, multi)
            if key not in records:
                records[key] = seq
    return records


def write_fasta(records, path):
    with open(path, "w") as fh:
        for key, seq in records.items():
            fh.write(f">{key}\n{seq}\n")


def run_mmseqs(q_fa, t_fa, m8_path, tmp_dir, sensitivity, threads):
    if shutil.which("mmseqs") is None:
        sys.exit("错误：未在 PATH 找到 mmseqs（mmseqs2）。请安装 mmseqs2 后重试，"
                 "或用 --m8-cache 指向已有 m8 仅重放判定（重放模式不需要 mmseqs2）。")
    cmd = ["mmseqs", "easy-search", str(q_fa), str(t_fa), str(m8_path),
           str(tmp_dir / "mmseqs_tmp"), "-s", str(sensitivity),
           "--format-output", "query,target,pident",
           "--threads", str(threads)]
    print("mmseqs2:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def parse_m8_max_pident(m8_path):
    """每个 query 的 max pident；跳过 q==t 自比对（与 lbs 修复版一致）。"""
    max_pident = {}
    with open(m8_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            q, t, pident = parts[0], parts[1], float(parts[2])
            if q == t:
                continue
            if pident > max_pident.get(q, 0.0):
                max_pident[q] = pident
    return max_pident


def main():
    args = parse_args()
    seq_cols = [c.strip() for c in args.seq_col.split(",") if c.strip()]
    multi = len(seq_cols) > 1
    thresholds = [float(t) for t in args.thresholds.split(",")]
    thr_labels = [t.rstrip("0").rstrip(".") if "." in t else t
                  for t in args.thresholds.split(",")]

    test = read_csv(args.test)
    refs = pd.concat([read_csv(p) for p in args.ref], ignore_index=True)
    for col in seq_cols:
        for name, df in (("test", test), ("ref", refs)):
            if col not in df.columns:
                sys.exit(f"错误：{name} 缺少序列列 {col!r}（现有列：{list(df.columns)}）")
    if args.id_col not in test.columns:
        sys.exit(f"错误：test 缺少 id 列 {args.id_col!r}")
    print(f"test={len(test)}, ref={len(refs)} ({len(args.ref)} 个文件)")

    tmp_ctx = None
    replay = args.m8_cache is not None and args.m8_cache.exists()
    if replay:
        m8_path = args.m8_cache
        print(f"重放已有 m8：{m8_path}（不调用 mmseqs2）")
    else:
        tmp_ctx = (tempfile.TemporaryDirectory(prefix="seq_homology_")
                   if args.tmp_dir is None else None)
        tmp_dir = args.tmp_dir or Path(tmp_ctx.name)
        tmp_dir.mkdir(parents=True, exist_ok=True)
        q_fa, t_fa = tmp_dir / "query.fasta", tmp_dir / "target.fasta"
        q_records = collect_records(test, seq_cols, args.id_col, multi)
        # target 只需按序列去重（键名不参与判定）
        t_records = {}
        seen_seqs = set()
        for row in refs[seq_cols].itertuples(index=False):
            for seq in row:
                if seq and seq not in seen_seqs:
                    seen_seqs.add(seq)
                    t_records[f"t{len(t_records)}"] = seq
        write_fasta(q_records, q_fa)
        write_fasta(t_records, t_fa)
        print(f"query 去重序列：{len(q_records)}，target 去重序列：{len(t_records)}")
        m8_path = args.m8_cache or (tmp_dir / "result.m8")
        m8_path.parent.mkdir(parents=True, exist_ok=True)
        run_mmseqs(q_fa, t_fa, m8_path, tmp_dir, args.sensitivity, args.threads)

    max_pident = parse_m8_max_pident(m8_path)
    hit_keys = set(max_pident)
    print(f"有 hit 的 query：{len(hit_keys)}")

    # 每链 max pident（无 hit → 0.0）与 per-chain 判定
    chain_cols = []
    for col in seq_cols:
        mp_col = "max_pident" if not multi else f"max_pident_{col}"
        keys = test[args.id_col].map(lambda u: query_key(u, col, multi))
        test[mp_col] = keys.map(lambda k: max_pident.get(k, 0.0))
        chain_cols.append((col, mp_col, keys.isin(hit_keys)))

    def merge(flag_per_chain):
        out = flag_per_chain[0]
        for f in flag_per_chain[1:]:
            out = out | f if args.merge == "any" else out & f
        return out

    for thr, label in zip(thresholds, thr_labels):
        col = f"seq_Redundancy_{label}"
        test[col] = merge([(test[mp] < thr) for _, mp, _ in chain_cols])
        test[col] = test[col].map({True: "True", False: "False"})
    test["OOD_Orphan"] = merge([~has_hit for _, _, has_hit in chain_cols])
    test["OOD_Orphan"] = test["OOD_Orphan"].map({True: "True", False: "False"})

    # 阶梯单调性与子集关系自检
    flags = {label: test[f"seq_Redundancy_{label}"] == "True" for label in thr_labels}
    asc = sorted(zip(thresholds, thr_labels))
    for (_, lo), (_, hi) in zip(asc, asc[1:]):
        assert bool((flags[lo] <= flags[hi]).all()), f"阶梯单调性破坏: {lo} vs {hi}"
    lowest = asc[0][1]
    orphan = test["OOD_Orphan"] == "True"
    assert bool((orphan <= flags[lowest]).all()), "OOD_Orphan 不是最低阈值列的子集！"
    print("子集关系验证通过："
          f"OOD_Orphan ⊆ {' ⊆ '.join(f'seq_Redundancy_{l}' for _, l in asc)}")

    # Default=False（极端长度预挑出）行：结果列直接置 NaN（"未计算"语义，
    # 全项目统一惯例，流水线位置写入；finalize_ood_columns.py --mode nan 仅作兜底校验）
    result_cols = [f"seq_Redundancy_{label}" for label in thr_labels] + ["OOD_Orphan"]
    if args.default_col in test.columns:
        nd = test[args.default_col] != "True"
        if nd.any():
            fill = "" if args.nondefault == "nan" else "False"
            test.loc[nd, result_cols] = fill
            print(f"Default=False {int(nd.sum())} 行结果列置 {'NaN' if fill == '' else 'False'}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    test.to_csv(args.out, index=False, lineterminator="\n")
    print(f"-> {args.out}（{len(test)} 行）")
    for label in thr_labels:
        print(f"  seq_Redundancy_{label}: {int(flags[label].sum())}")
    print(f"  OOD_Orphan: {int(orphan.sum())}")
    if tmp_ctx is not None:
        tmp_ctx.cleanup()


if __name__ == "__main__":
    main()
