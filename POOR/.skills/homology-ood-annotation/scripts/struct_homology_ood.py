#!/usr/bin/env python3
r"""结构同源 OOD 计算（homology-ood-annotation skill 第 2 步）。

统一口径（以 fold_classification/scripts/annotate_test_ood.py 为基准）：
- 搜索方向：query = test 结构，target = train+val 结构
- 工具：foldseek `easy-search`，输出 `query,target,alntmscore`
- 判定：每个 test 结构取对 ref 的 max alntmscore；无 hit 按 0.0 计
  TM-score_<T> = (max alntmscore < T)
- query id 取结构文件名去路径与扩展名（与 fold 脚本对 foldseek 输出的
  归一化一致：`str.replace(r"\.pdb$","")` + 取路径末段）

输入结构：--test-structs / --ref-structs 各接受一个 PDB/mmCIF 目录
（取其中全部 *.pdb / *.cif）或一个清单文件（每行一个结构文件路径）。

m8 缓存：--m8-cache 指向的 m8 已存在时仅重放判定（无需 foldseek）；
不存在时运行 foldseek 并把 m8 写到该路径供下次复用。

输出：每个 test 结构一行：query_id, max_alntmscore, TM-score_<T>...

用法：
    python struct_homology_ood.py --test-structs pdbs_test.txt \
        --ref-structs pdbs_trainval.txt --out tm_flags.csv \
        --m8-cache ood_foldseek_test_vs_trainval.m8
    # 仅重放已有 m8（不需要 foldseek）
    python struct_homology_ood.py --test-structs pdbs_test.txt \
        --ref-structs pdbs_trainval.txt --out replay.csv --m8-cache existing.m8
"""
import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

DEFAULT_THRESHOLDS = "0.9,0.8,0.7,0.6,0.5,0.4,0.3"
STRUCT_EXTS = (".pdb", ".cif", ".ent")


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--test-structs", required=True, type=Path,
                    help="test 结构目录或清单文件（每行一个结构路径）")
    ap.add_argument("--ref-structs", required=True, type=Path,
                    help="train+val 结构目录或清单文件")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--thresholds", default=DEFAULT_THRESHOLDS,
                    help="逗号分隔的 alntmscore 阈值（max < 阈值 → True）")
    ap.add_argument("--m8-cache", type=Path, default=None,
                    help="m8 缓存路径：存在则仅重放；不存在则搜索后写入复用")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--tmp-dir", type=Path, default=None)
    return ap.parse_args()


def struct_list(spec):
    """目录 → 全部结构文件；文件 → 每行一个路径。"""
    if spec.is_dir():
        files = sorted(p for p in spec.iterdir() if p.suffix.lower() in STRUCT_EXTS)
    else:
        files = []
        with open(spec) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    files.append(Path(line))
    missing = [p for p in files if not p.exists()]
    if missing:
        sys.exit(f"错误：{spec} 中 {len(missing)} 个结构文件不存在，"
                 f"如 {missing[0]}")
    return files


def normalize_id(name):
    """fold 脚本的 query 归一化：去结构扩展名（.pdb/.cif/.ent，大小写不敏感）、
    取路径末段。foldseek createdb 对任何结构扩展名都会剥掉，mmCIF 输入
    （如 8RRI_K.CIF）在 m8 中同样以去扩展名形式出现，必须一并剥除，
    否则 .cif query 会全部误判为无 hit。"""
    s = str(name)
    low = s.lower()
    for ext in STRUCT_EXTS:
        if low.endswith(ext):
            s = s[: -len(ext)]
            break
    return s.split("/")[-1]


def run_foldseek(q_dir, t_dir, m8_path, tmp_dir, threads):
    if shutil.which("foldseek") is None:
        sys.exit("错误：未在 PATH 找到 foldseek。请安装 foldseek 后重试，"
                 "或用 --m8-cache 指向已有 m8 仅重放判定（重放模式不需要 foldseek）。")
    cmd = ["foldseek", "easy-search", str(q_dir), str(t_dir), str(m8_path),
           str(tmp_dir / "foldseek_tmp"),
           "--format-output", "query,target,alntmscore",
           "--threads", str(threads)]
    print("foldseek:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def parse_m8_max_tm(m8_path):
    """每个 query 的 max alntmscore；跳过 q==t 自比对。"""
    max_tm = {}
    with open(m8_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            q, t = normalize_id(parts[0]), normalize_id(parts[1])
            tm = float(parts[2])
            if q == t:
                continue
            if tm > max_tm.get(q, 0.0):
                max_tm[q] = tm
    return max_tm


def main():
    args = parse_args()
    thresholds = [float(t) for t in args.thresholds.split(",")]
    thr_labels = [t.rstrip("0").rstrip(".") if "." in t else t
                  for t in args.thresholds.split(",")]

    replay = args.m8_cache is not None and args.m8_cache.exists()
    tmp_ctx = None
    if replay:
        m8_path = args.m8_cache
        print(f"重放已有 m8：{m8_path}（不调用 foldseek）")
        # 重放模式仍需 test 清单以确定输出行
        test_files = struct_list(args.test_structs)
        query_ids = [normalize_id(p.name) for p in test_files]
    else:
        test_files = struct_list(args.test_structs)
        ref_files = struct_list(args.ref_structs)
        query_ids = [normalize_id(p.name) for p in test_files]
        if len(set(query_ids)) != len(query_ids):
            sys.exit("错误：test 结构 id（文件名去扩展名）存在重复，无法一一对应")
        print(f"test 结构：{len(test_files)}，ref 结构：{len(ref_files)}")
        tmp_ctx = (tempfile.TemporaryDirectory(prefix="struct_homology_")
                   if args.tmp_dir is None else None)
        tmp_dir = args.tmp_dir or Path(tmp_ctx.name)
        tmp_dir.mkdir(parents=True, exist_ok=True)
        # foldseek 需要结构文件目录，用符号链接避免拷贝（fold 脚本同款做法）
        q_dir, t_dir = tmp_dir / "fs_query", tmp_dir / "fs_target"
        q_dir.mkdir(exist_ok=True)
        t_dir.mkdir(exist_ok=True)
        for p in test_files:
            dst = q_dir / p.name
            if not dst.exists():
                dst.symlink_to(p.resolve())
        for p in ref_files:
            dst = t_dir / p.name
            if not dst.exists():
                dst.symlink_to(p.resolve())
        m8_path = args.m8_cache or (tmp_dir / "result.m8")
        m8_path.parent.mkdir(parents=True, exist_ok=True)
        run_foldseek(q_dir, t_dir, m8_path, tmp_dir, args.threads)

    max_tm = parse_m8_max_tm(m8_path)
    print(f"有 hit 的 query：{len(max_tm)} / {len(query_ids)}")

    rows = []
    flags_by_label = {label: [] for label in thr_labels}
    for qid in query_ids:
        tm = max_tm.get(qid, 0.0)
        row = {"query_id": qid, "max_alntmscore": tm}
        for thr, label in zip(thresholds, thr_labels):
            v = tm < thr
            row[f"TM-score_{label}"] = str(v)
            flags_by_label[label].append(v)
        rows.append(row)

    # 阶梯单调性自检
    asc = sorted(zip(thresholds, thr_labels))
    for (_, lo), (_, hi) in zip(asc, asc[1:]):
        for a, b in zip(flags_by_label[lo], flags_by_label[hi]):
            assert a <= b, f"阶梯单调性破坏: {lo} vs {hi}（query 行）"

    out = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False, lineterminator="\n")
    print(f"-> {args.out}（{len(out)} 行）")
    for label in thr_labels:
        print(f"  TM-score_{label}: {sum(flags_by_label[label])}")
    if tmp_ctx is not None:
        tmp_ctx.cleanup()


if __name__ == "__main__":
    main()
