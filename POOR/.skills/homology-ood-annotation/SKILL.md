---
name: homology-ood-annotation
description: 为任意任务按 POOD-Benchmark 统一口径计算搜索型同源 OOD 列：mmseqs2 序列同源（seq_Redundancy_90..30 + OOD_Orphan）与 foldseek 结构同源（TM-score_0.9..0.3）。搜索方向固定为 test→train+val，支持 m8 缓存复用与双链任务（aa_seq1/aa_seq2）。当需要为数据集补充或重算序列/结构同源性 OOD 列时使用。
---

# 搜索型同源 OOD 标注（mmseqs2 序列 + foldseek 结构）

## 描述

本 skill 封装了 POOD-Benchmark 统一的**搜索型**同源 OOD 流程：`seq_homology_ood.py`（mmseqs2 easy-search，序列）与 `struct_homology_ood.py`（foldseek easy-search，结构）。两者均为"test 样本对 train+val 参考集搜索，取最大相似度按阈值阶梯判定"，与 `cath-ted-domain-annotation`（分类标签留一法）互补。

## 何时使用

- 任意任务需要补充或重算 `seq_Redundancy_90` ~ `seq_Redundancy_30` 与 `OOD_Orphan`
- 任意任务需要补充或重算 `TM-score_0.9` ~ `TM-score_0.3`
- 已有 mmseqs2 / foldseek 的 m8 结果文件，只需重放阈值判定（复跑、阈值调整）

## 核心口径（各任务必须一致）

| 决策点 | 约定 | 理由 |
|:---|:---|:---|
| 搜索方向 | **query = test，target = train+val** | 与 fold `annotate_test_ood.py` 一致；反向搜索会漏掉"test 有多个同源"的语义 |
| 序列搜索 | mmseqs2 `easy-search -s 7.5 --format-output query,target,pident`（easy-search 默认 e-value ≤ 1e-3 门槛） | fold 现行实现；easy-search 显著性门槛保证 hit 均为显著同源。**-s 7（lbs）vs 7.5 的实测差异可忽略**（lbs 8,864 链两档对比：seq_Redundancy ≥50 各列零差异、seq40 Δ=3、seq30 Δ=9、Orphan Δ=18/932，仅影响 pident 22.7–42.0 的检测极限弱 hit，见 lbs `output/verify/sens_s{7,75}.m8`） |
| 序列判定 | 每个 test 样本取 max pident（无 hit 按 0.0 计）；`seq_Redundancy_<T>` = max pident < T（T ∈ 90,80,...,30） | 阶梯阈值，步长 10% |
| Orphan 语义 | **无任何显著 hit → `OOD_Orphan=True`** | easy-search 默认 e-value ≤ 1e-3，"无 hit" 等价于旧 lbs 修复版的 "best e-value > 1e-3"；保证 **OOD_Orphan ⊆ seq_Redundancy_30 ⊆ ... ⊆ seq_Redundancy_90**（脚本内置断言自检） |
| 结构搜索 | foldseek `easy-search --format-output query,target,alntmscore`；query id 取文件名去路径与 `.pdb` 后缀 | fold 现行实现的归一化方式 |
| 结构判定 | 每个 test 结构取 max alntmscore（无 hit 按 0.0 计）；`TM-score_<T>` = max alntmscore < T（T ∈ 0.9,...,0.3） | 阶梯阈值，步长 0.1 |
| 双链任务合并 | 两条链各自搜索，样本级判定按 **any 语义**合并：`seq_Redundancy_<T>` = (A<T) \| (B<T)，`OOD_Orphan` = (A 无 hit) \| (B 无 hit)；TM-score 同理 | 与 ppi `10_recompute_all_ood.py` 的 `(A_max_id<thr)\|(B_max_id<thr)`、`(A==0)\|(B==0)` 一致（**不是取两链 max pident 再判定**——任一侧链无同源即将整个 pair 判为 OOD）；`--merge all` 可选 |
| 自比对 | 解析 m8 时跳过 q==t 行 | 与 lbs 修复版一致；防止 test/ref 有重叠 id 时 max pident 恒 100 |
| m8 缓存 | `--m8-cache` 已存在 → 仅重放判定（不需要外部工具）；不存在 → 搜索后写入该路径 | 搜索是重计算，判定是轻量重放 |

## 数据来源（各任务历史实现与已知差异）

| 任务 | 实现 | 与本 skill 口径的差异 |
|:---|:---|:---|
| fold_classification | `scripts/05_annotate_test_ood.py`（**本 skill 基准**） | 无差异；m8 缓存 `output/ood_mmseqs_test_vs_trainval.m8` / `ood_foldseek_test_vs_trainval.m8` 可直接 `--m8-cache` 重放 |
| ligand_binding_site | `output/scripts_archive/24_fix_seq_homology_easysearch.py`（easy-search 修复终版，已归档） | `-s 7`（非 7.5）；format-output 多一列 evalue，Orphan 用 best e-value > 1e-3 判定（与"无 hit"等价）；query 按 chain_key（`{pdb}_{chain}`）去重而非 unique_id |
| ppi_prediction | `scripts/10_recompute_all_ood.py` | 2026-08-27 起 mmseqs2 与 foldseek 均统一为 **easy-search**（skill 脚本重算写回，缓存 `output/ood_mmseqs_test_vs_trainval_easysearch.m8`、`ood_foldseek_test_vs_trainval_easysearch.m8`，均可 `--m8-cache` 重放）；此前 mmseqs 为旧式 `search --max-seqs 10 -e 100 -c 0.0`（seq30 515→3,653、Orphan 211→3,017），foldseek 为 `search -a -e inf` 保留全部 hit（旧 max_tmscore 分布在 [0.3,0.4) 出现 0 条的断层伪影、TM-score_0.3 与 0.4 两列退化相等，重算后 0.9: 4,827→8,172、0.3: 1,610→943）；旧缓存仅供溯源。双链 any 合并即本 skill 默认 |
| func_prediction | `scripts/orphan_mark.py`（**已归档**，现行口径见任务 DATA_PROCESS.md §4.3 easy-search -s 7） | 历史参考集是 UniRef50 而非 train+val，且 Orphan 还要求 fident ≥ 0.30——语义不同，不可用本 skill 重放其旧 OOD_Orphan |
| lbs 旧流程 | `scripts/11_build_ood_splits.py`（搜索调用已同步为 easy-search -s 7） | 旧式 `mmseqs search -e 100 -c 0.0` 弱比对顶满 max pident，曾导致 Orphan(1,285) ≫ SeqID<30%(134) 倒挂；由 24_fix（已归档）修复 |

> 重算既有任务时优先沿用该任务的 m8 缓存重放；重新搜索则用本 skill 默认参数，并在 DATA_PROCESS.md 注明与历史实现的参数差异。

## 用法

### 1. 序列同源（seq_homology_ood.py）

```bash
# 单链任务（fold 口径）
python .skills/homology-ood-annotation/scripts/seq_homology_ood.py \
    --test cath_test.csv --ref cath_train.csv cath_val.csv \
    --out test_flagged.csv --m8-cache output/ood_mmseqs_test_vs_trainval.m8 \
    --threads 16

# 双链任务（ppi 口径，any 合并）
python .skills/homology-ood-annotation/scripts/seq_homology_ood.py \
    --test ppi_test.csv --ref ppi_train.csv ppi_val.csv \
    --seq-col aa_seq1,aa_seq2 --out ppi_flagged.csv --threads 16

# 仅重放已有 m8（不需要 mmseqs2；注意 m8 query 键需与 unique_id 对应）
python .skills/homology-ood-annotation/scripts/seq_homology_ood.py \
    --test cath_test.csv --ref cath_train.csv \
    --out replay.csv --m8-cache existing.m8
```

输出 = test 全行 + `max_pident`（双链为 `max_pident_<seq-col>` 每链一列）
+ 7 个 `seq_Redundancy_*`（True/False 字符串）+ `OOD_Orphan`。

关键参数：`--seq-col`（默认 `aa_seq`）、`--thresholds`（默认 `90,80,70,60,50,40,30`）、
`--merge {any,all}`、`--sensitivity`（默认 7.5）、`--id-col`、`--tmp-dir`。

### 2. 结构同源（struct_homology_ood.py）

```bash
# 目录模式（取目录下全部 *.pdb/*.cif）
python .skills/homology-ood-annotation/scripts/struct_homology_ood.py \
    --test-structs test_pdbs/ --ref-structs trainval_pdbs/ \
    --out tm_flags.csv --m8-cache output/ood_foldseek_test_vs_trainval.m8

# 清单模式（每行一个结构文件路径）
python .skills/homology-ood-annotation/scripts/struct_homology_ood.py \
    --test-structs test_pdbs.txt --ref-structs trainval_pdbs.txt --out tm_flags.csv
```

输出 = 每个 test 结构一行：`query_id, max_alntmscore, TM-score_0.9 ... TM-score_0.3`。

### 3. 写回 splits（调用方职责）

与 ec skill 相同：先备份、整表 `dtype=str, keep_default_na=False` 读入、仅覆写目标列、
Default=False 行按任务惯例置 NaN（可用 `ood-annotation-toolkit` 的
`finalize_ood_columns.py --mode nan` 统一处理）、写回后逐单元格 diff 验证。

## 验证要求

- 脚本内置断言：`OOD_Orphan ⊆ seq_Redundancy_30 ⊆ ... ⊆ seq_Redundancy_90` 阶梯单调性；
  TM-score 列同理（破坏即报错退出）
- 回归基准（fold，12,911 行 test，`--m8-cache` 重放，2026-08-27 实测）：与
  `splits/cath_test.csv` 的 7 个 seq_Redundancy 列 + OOD_Orphan、7 个 TM-score 列逐行一致
- 新任务抽查 ≥5 条 OOD_Orphan=True 样本：确认其序列/结构对 train+val 确无显著 hit
- 有 hit query 数应与 m8 行数同量级；异常偏少时先检查 query 键是否匹配（双链缓存的
  键为 `unique_id|<seq-col>`，历史任务缓存键可能是 protein_id / chain_key，不可混用）

## 注意事项

1. **缓存键兼容**：本 skill 自产 m8 的 query 键 = `unique_id`（双链为 `unique_id|<seq-col>`）；
   fold 的缓存键同为 unique_id 可直接重放；ppi（protein_id）/ lbs（chain_key）的历史
   缓存键不同，重放前需确认键空间一致，否则全部判为 Orphan
2. **阈值列名**：`--thresholds` 写入列名时去掉尾零（90→`seq_Redundancy_90`，
   0.9→`TM-score_0.9`），自定义阈值注意与既有列名对齐
3. foldseek 需要结构文件目录，脚本用符号链接组建 query/target 目录（fold 同款做法，
   避免拷贝 ~10 GB）；`--tmp-dir` 指向的目录会保留 fasta/m8 等中间产物
4. mmseqs2 / foldseek 不在 PATH 时给出明确报错并提示改用 `--m8-cache` 重放
5. Default=False（极端长度）行的 NaN 惯例不在本 skill 处理——本 skill 输出全行计算值，
   NaN 化由 `ood-annotation-toolkit/finalize_ood_columns.py` 统一执行
