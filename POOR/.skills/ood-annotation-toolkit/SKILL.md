---
name: ood-annotation-toolkit
description: POOD-Benchmark 计算型 OOD 与收尾通用工具箱：metapredict IDR OOD（idr_ood.py，两种历史口径可选）、回归任务值域分 bin LongTail OOD（value_bin_longtail_ood.py）、Default=False 行 NaN 化与 InD 列重算（finalize_ood_columns.py）、划分文件一致性与泄露检查（check_splits.py，只读）。当需要为数据集计算 IDR/LongTail-bin OOD 列、统一执行 NaN/InD 收尾惯例、或体检 splits 时使用。
---

# OOD 标注工具箱（计算型 OOD 与收尾）

## 描述

本 skill 收编了 POOD-Benchmark 四个横切工具：IDR 无序 OOD、值域分 bin LongTail OOD、OOD 列收尾（NaN 惯例 + InD 重算）、splits 一致性体检。与 `homology-ood-annotation`（搜索型同源）、`cath-ted-domain-annotation`（结构分类）、`ec-function-ood-annotation`（EC 功能）互补。

## 何时使用

- 任意任务需要计算 `idr_ratio` + `OOD_IDR`（metapredict v3）
- 回归任务（kcat / affinity / optimal_ph 类）需要 `OOD_LongTail_*Bin_le*` 列
- 任意 test csv 需要统一执行 Default=False→NaN 惯例或重算 `InD` 列
- 需要检查 train/val/test 泄露、列结构、Default/InD 统计（只读体检）

## 核心口径

### idr_ood.py（IDR 无序 OOD）——口径存在任务间分裂，必须显式选 mode

| 决策点 | 约定 |
|:---|:---|
| 预测工具 | metapredict v3（实测 3.0.2），残基 disorder > 0.5 判无序 |
| `--mode residue-ratio` | 无序残基占比（ss / func / fold / kcat / affinity / optimal_ph 口径；affinity 的"区域长度求和/序列长"数学上与此等价） |
| `--mode region-ratio` | **最长连续无序区**长度 / 序列长（ppi 口径，源自已归档的 `18_compute_ood_idr.py`） |
| `--threshold` | 默认 0.3（ss/func/fold/ppi）；**affinity / optimal_ph 历史用 0.1** |
| `--clean strip`（默认） | 删除非 20 种标准氨基酸字符（fold / ss 口径） |
| `--clean map` | B→N, Z→Q, J→L, U→C, O→K，其余→A（ppi 口径） |
| 双链任务 | `--seq-col aa_seq1,aa_seq2`，OOD_IDR = 任一链超阈值（any 语义，ppi 口径） |
| 缓存 | `--cache` 按清洗后序列增量缓存（sequence, residue_ratio, max_region_ratio 双指标同时存），逐 batch 落盘，换 mode/阈值重跑零重算 |

各任务历史口径差异表（重算时必须对齐）：

| 任务 | mode | 阈值 | clean | 来源 |
|:---|:---|:---|:---|:---|
| ss_prediction | residue-ratio | 0.3 | strip（仅删 XBZJU，与 strip 近似） | `scripts/09_run_idr_prediction.py` |
| func_prediction | residue-ratio | 0.3 | strip | DATA_PROCESS.md §IDR |
| fold_classification | residue-ratio | 0.3 | strip | `scripts/05_annotate_test_ood.py`（缓存 `output/idr_ratio_full.csv`） |
| ppi_prediction | **region-ratio** | 0.3 | **map** | `output/scripts_archive/18_compute_ood_idr.py` / `scripts/10_recompute_all_ood.py` |
| ligand_binding_affinity | residue-ratio | **0.1** | strip（另去 `\|` 链分隔符） | `scripts/10_compute_structure_idr_ood.py` |
| enzyme_optimal_ph | residue-ratio | **0.1** | strip | `scripts/04_compute_ood_markers.py`（原 compute_ood_idr.py 已归档并入此脚本） |
| enzyme_kinetics | residue-ratio（仅 `idr_ratio` 连续列） | — | strip | DATA_PROCESS.md |

### value_bin_longtail_ood.py（回归值域分 bin LongTail）

| 决策点 | 约定 |
|:---|:---|
| 判定 | test Default=True 行所属 bin 的 **train** 频次 <= 阈值 → True；label 缺失/超出 bin 范围 → False |
| clip 单独成箱 | `--clip lo,hi`：==lo 与 ==hi 的值各自独立成 bin（kcat 删失值口径，防止边界堆积遮蔽尾部） |
| Default=False 行 | 默认 NaN（`--nondefault nan`，全项目统一惯例，流水线位置直接写入）；`--nondefault false` 仅为兼容保留 |
| 分箱 | `pd.cut` 默认左开右闭 `(a,b]`，边 = arange(start, end, width) |

各任务历史参数（重算时对齐）：

| 任务 | label | bin | clip | 阈值 | 列前缀 |
|:---|:---|:---|:---|:---|:---|
| kcat | log10(kcat) | [-6.5, 7.0) 宽 0.5 | -6,6 | 50,100 | `OOD_LongTail_KcatBin` |
| ligand_binding_affinity | -log10(K) | [-0.5, 16.5) 宽 0.5 | 无 | 50,100 | `OOD_LongTail_AffBin` |
| enzyme_optimal_ph | pH | [1.5, 13.5) 宽 0.5 | 无 | 50,100 | `OOD_LongTail_pHbin` |

### finalize_ood_columns.py（NaN 惯例 + InD 收尾）

| 决策点 | 约定 |
|:---|:---|
| `--mode nan` | Default==False 行：除 `OOD_ExtremeShort`/`OOD_ExtremeLong` 外所有 **OOD_ 前缀**列置 NaN（空单元格，表"未计算"）；**只动 OOD_ 前缀列**，seq_Redundancy_*/TM-score_*/idr_ratio 不在范围（与原 `nondefault_ood_nan.py` 一致） |
| 豁免开关 | `--exempt-nondefault` 跳过 nan（预留给整表 Default=False 属刻意设计的特殊文件；现有任务均不需要——PDA 两个 design 集合为整表 Default=True） |
| `--mode ind` | `InD` = 所有 `OOD_*` + `seq_Redundancy_*` + `TM-score_*` 均为 False（NaN 视为 False）；已有 InD → 原位覆写；有 `Pure_ID`（fold/PDA 旧名）→ 原位更名重算；否则末尾追加 |
| 字节保护 | 整表 `dtype=str, keep_default_na=False` 读入，写回前断言非目标列逐单元格不变 |

### check_splits.py（只读体检）

泄露（train∩val / train∩test / val∩test unique_id）、核心列完整性、废弃列、label 空值、train/val 不应携带 OOD 列、列集对比（train==val；test 仅多 OOD 块与 `idr_ratio` 等合法标注列）、Default/InD 统计、Default=False 行 NaN 惯例违例、`OOD_Orphan ⊆ seq_Redundancy_30` 子集关系。任何问题 exit 1。

## 用法

```bash
S=.skills/ood-annotation-toolkit/scripts

# 1. IDR（ss/func/fold 口径）
python $S/idr_ood.py --test ssp_test.csv --out ssp_idr.csv \
    --mode residue-ratio --threshold 0.3 --cache output/idr_cache.csv
# IDR（ppi 口径，双链）
python $S/idr_ood.py --test ppi_test.csv --seq-col aa_seq1,aa_seq2 \
    --mode region-ratio --clean map --out ppi_idr.csv --cache output/idr_cache.csv

# 2. 值域分 bin LongTail（kcat 口径）
python $S/value_bin_longtail_ood.py --train kcat_train.csv --test kcat_test.csv \
    --out kcat_test.csv --bin-start -6.5 --bin-end 7.0 --bin-width 0.5 \
    --clip -6,6 --thresholds 50,100 --col-prefix OOD_LongTail_KcatBin

# 3. 收尾（先 NaN 化再重算 InD；--dry-run 只报告不写回）
python $S/finalize_ood_columns.py --test-csv cath_test.csv --mode both --dry-run
python $S/finalize_ood_columns.py --test-csv cath_test.csv --mode both

# 4. 体检（只读；ppi 需指定双链核心列）
python $S/check_splits.py --train ppi_train.csv --val ppi_val.csv --test ppi_test.csv \
    --seq-cols aa_seq1,aa_seq2 --struct-cols struct_file1,struct_file2
```

## 验证要求

- `idr_ood.py` 回归基准：fold 口径（residue-ratio @0.3, strip）应复现 `splits/cath_test.csv` 的 `idr_ratio`/`OOD_IDR`；换缓存重跑秒级完成
- `finalize_ood_columns.py --dry-run` 对 fold / ppi 现状应报告 InD 零差异、NaN 无清除（2026-08-27 实测两者均已合规）
- `check_splits.py` 对 fold 与 ppi 全部通过（2026-08-27 实测）
- 写回类脚本（value_bin / finalize）遵循惯例：先备份、只覆写目标列、写后 diff 验证

## 注意事项

1. **IDR 不可混口径**：ppi 的 region-ratio 与其他任务的 residue-ratio 不可互换；同一 `--cache` 文件双指标并存，可安全跨任务复用，但判定列必须按本任务历史 mode 出
2. **metapredict 环境**：不在所有 python 环境中可用（本机 `/root/miniconda3/bin/python` 有 3.0.2）；缺库时报错并提示改用缓存
3. **重复 unique_id**（kcat 同酶多底物）：本工具箱按行写值，不涉及 uid 去重映射问题；EC 类任务的 groupby 注意事项见 ec skill
4. **未收编的 analysis 脚本**：`migrate_column_schema.py` 是一次性列名迁移（9 项已完成的 rename，无复用价值），`survey_columns.py` 的全库多任务列频率清点功能未吸收（check_splits.py 只做单任务三文件列清单）——两者逻辑不在本 toolkit 内，原件已删除
5. optimal_ph 的 LongTail bin 用 `--nondefault false` 对齐历史；新任务一律用默认 `nan`
