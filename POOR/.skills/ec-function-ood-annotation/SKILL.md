---
name: ec-function-ood-annotation
description: 为任意任务的样本生成统一的 EC 编号标注，并按 POOD-Benchmark 统一口径计算 OOD_NewEC_L3/L4（all-not-in、train+val 参考）与 OOD_LongTail_EC_L3/L4_le5/le10（all-not-in、train 频次、含 0）。当需要为数据集补充或重算 EC 功能新颖性 / 长尾功能 OOD 列时使用。
---

# EC 功能标注与 NewEC / LongTailEC OOD 计算

## 描述

本 skill 封装了 POOD-Benchmark 统一的 EC 功能 OOD 流程：先为样本生成 EC 编号集合（SIFTS 链级 / UniProt 反查 + API 补缺），再按统一口径计算 6 个 OOD 列。与 `cath-ted-domain-annotation`（结构分类 OOD）互补，覆盖功能维度。

## 何时使用

- 任意任务需要补充或重算 `OOD_NewEC_L3` / `OOD_NewEC_L4`
- 任意任务需要补充或重算 `OOD_LongTail_EC_{L3,L4}_le{5,10}`
- 需要为样本生成统一来源的 EC 编号集合

## 核心口径（各任务必须一致）

| 决策点 | 约定 | 理由 |
|:---|:---|:---|
| 层级解析 `levels()` | EC 段数 ≥3 → L3=前三段；段数 ==4 → L4=完整编号；**允许部分 "-"**（如 `3.4.22.-` → L3=`3.4.22`） | 与统一驱动脚本一致；丢弃部分编号会漏标 |
| **NewEC** 参考集 | **train+val** 出现的 EC 集合 | 2026-08-27 定稿口径 |
| **NewEC** 判定 | **all-not-in**：有 EC 注释且所有 EC（对应层级）都不在参考集 → True；无注释 → False | 多 EC 样本只有"全新功能"才算 NewEC；any 语义过松（ppis 曾 142→162） |
| **LongTailEC** 参考集 | **train-only**，逐样本计数（携带该类别的 train 样本数） | 与 func/optimal_ph/kcat 的 LongTail 口径一致 |
| **LongTailEC** 判定 | **all-not-in**：有 EC 注释且所有 EC（对应层级）train 频次都 < 5 / < 10（**含 0**）→ True；无注释 → False | 含 0 使 NewEC ⊆ LongTailEC，两场景语义嵌套 |
| Default/极端行 | **Default=False（极端长度预挑出）行一律置 NaN**（2026-08-27 起全项目统一，表示"未计算"，与 False 区分；2026-08-28 起由本 skill 脚本与各任务流水线脚本在生成位置直接写入，`.skills/ood-annotation-toolkit/scripts/finalize_ood_columns.py --mode nan` 仅作兜底校验） | PDA 两个 design 集合为整表 Default=True（额外测试集、无预挑出样本），不触发该惯例 |

> NewEC 与 LongTailEC 均为 all-not-in（2026-08-27 定稿）：NewEC 衡量"完全陌生的功能"，
> LongTailEC 衡量"所有功能均为长尾"。初版 LongTailEC 曾为 any 语义，后统一定稿为 all-not-in；
> func_ec 的 `OOD_LongTail` 是预测标签级（L3）口径，与 LongTailEC（L4）层级不同，不构成语义冲突。

## 数据来源

| 文件/接口 | 用途 |
|:---|:---|
| `data/sifts/sifts_chain_ec.tsv.gz` | PDB+链 → EC（链模式唯一来源） |
| `data/sifts/pdb_chain_uniprot.tsv.gz` | accession → PDB 链反查（UniProt 模式来源 ①） |
| `data/uniprot_ec_cache.json` | **项目级共享 UniProt EC 缓存**（UniProt 模式来源 ②，跨任务复用） |
| UniProt REST `rest.uniprot.org/uniprotkb/search?fields=accession,ec` | 缓存未命中时批量查询（100/批，0.5s 间隔） |

> **EC 源漂移警告**：UniProt 当前 EC 注释与数据集原始标注可能不同（实测：kcat 数据集中
> Q71SQ4/Q9BP39 的 `1.11.1.15` 在 UniProt 已删除，P91883 由 `1.11.1.15` 更新为 `1.11.1.24`）。
> 重算既有任务时必须沿用该任务的历史 EC 源（见下表），新任务用本 skill 默认源即可。

### 各任务历史 EC 源对照（重算时保持一致）

| 任务 | EC 源 |
|:---|:---|
| kcat | `datasets/enzyme_kinetics_prediction/data/kcat_with_reactant_ecfp4.csv`（uniprot→ec，数据集原始标注；个别行为多 EC 字符串，需按 `;` 拆分） |
| optimal_ph | `datasets/enzyme_optimal_ph/data/metadata.csv`（ec_id；**782 个 UniProt 为 `; ` 分隔的多 EC 字符串，必须拆分**，2026-08-27 修复） |
| ppi | UniProt API 缓存（`output/uniprot_ec_mapping.csv`）∪ SIFTS 链级，bio/xtal 行双路并集，random 负样本仅 UniProt |
| ligand_affinity | SIFTS 链级，按 pdb 全链并集 |
| ss / lbs | SIFTS 链级 |
| ppis | SIFTS 链级，seq_hash→pair_meta→(pdb,chain) 映射链 |
| fold / PDA | `data/sifts/sifts_ec_annotations.json`（NewEC 为 EC7 L1 prior，**不在本 skill 口径内**） |
| func_ec | SIFTS 链级（NewEC/LongTailEC 均**仅 L4 列**——L3 是预测目标层级，其 `OOD_LongTail` 已覆盖 L3 级长尾；LongTailEC 不得用 label 列计算，label 实为 `x.x.x.-` 形式的 L3） |
| func go_mf/go_bp/go_cc | SIFTS 链级（NewEC L3+L4、LongTailEC 四列共六列，2026-08-27 落地，脚本 `add_go_ec_ood.py`（本 skill scripts/ 下）；GO 版 OOD 不可行——标签词表封闭致 NewGO 恒 0、全量 GO 注释 all-not-in 也恒 0） |

## 用法

### 1. EC 标注（ec_annotate.py）

```bash
# 链模式（输入含 pdb, chain 列）
python .skills/ec-function-ood-annotation/scripts/ec_annotate.py \
    --input my_chain_list.csv --out my_ec_labels.csv

# UniProt 模式（输入含 uniprot 列；--offline 仅用缓存）
python .skills/ec-function-ood-annotation/scripts/ec_annotate.py \
    --input my_acc_list.csv --out my_ec_labels.csv --offline
```

输出在输入列后追加：`ec_numbers`（`;` 分隔完整 EC）、`has_ec`。

### 2. OOD 计算（ec_ood.py）

```bash
# 单文件含 train/val/test（split 列标识，';' 分隔多归属）：
python .skills/ec-function-ood-annotation/scripts/ec_ood.py \
    --labels all_ec_labels.csv --split-col split --out test_flagged.csv

# test 与参考集分文件（ref 文件也需 split 列以区分 train/val）：
python .skills/ec-function-ood-annotation/scripts/ec_ood.py \
    --labels test_ec.csv --ref-labels trainval_ec.csv --out test_flagged.csv
```

输出仅 test 行，追加 6 列：`OOD_NewEC_L3/L4`、`OOD_LongTail_EC_L3_le5/le10`、`OOD_LongTail_EC_L4_le5/le10`。

### 3. 写入 splits（调用方职责）

1. **先备份** test csv（各任务 `output/backup_*` 目录）
2. 屏蔽 Default=False（极端长度）行为 NaN（全项目统一惯例）
3. L4 列就地重算、L3/LongTail 列追加末尾（与 2026-08-27 统一时一致）
4. 验证行数/列序/unique_id 不变，其余列值不变

### 4. 批量重跑既有任务（batch 驱动）

本 skill `scripts/` 下收编了 2026-08-27 全项目统一口径时的三个 batch 驱动脚本
（自 `analysis/scripts/` 移入，仓库根按 `Path(__file__).resolve().parents[3]` 解析）：

| 脚本 | 作用 | 覆盖任务 |
|:---|:---|:---|
| `add_unified_newec.py` | NewEC L3/L4 统一为 all-not-in + train+val 参考并落盘 | kcat、optimal_ph、ligand_affinity、ppi、ss、lbs、ppis、func_ec（仅 L4） |
| `add_unified_longtail_ec.py` | LongTailEC 四列统一为 all-not-in + train 频次（含 0）并落盘 | 上述 + fold、PDA cath_design（非 ES 行；fold/PDA 用 `sifts_ec_annotations.json` 源） |
| `add_go_ec_ood.py` | 把 NewEC + LongTailEC 六列扩展到 GO 任务 | func go_mf/go_bp/go_cc |

三个脚本内置各任务的 EC 源适配表（见"各任务历史 EC 源对照"），直接写回
`splits/` 并输出 `output/*_summary.json`。**警告**：它们不做备份、不执行
NaN 收尾——重跑前自行备份，重跑后用 `ood-annotation-toolkit` 的
`finalize_ood_columns.py --mode both` 收尾、`check_splits.py` 体检。
口径细节与单任务新标注优先用 `ec_annotate.py` + `ec_ood.py` 两阶段流程。

## 验证要求

- 链模式回归基准（ss，Default 子集，2026-08-27 实测）：NewEC_L3=1、NewEC_L4=246、
  LongTail_L3=37/77、LongTail_L4=739/1113（all-not-in）—— 本 skill 复算逐行零差异
- 抽查 ≥5 条 NewEC=True 样本：其所有 EC（对应层级）确认不在 train+val 集合
- 计数按**逐行**统计并注明；kcat 同一 accession 多行（多底物），逐行与唯一 accession 计数差异巨大（如 L4<10：2,355 vs 788）

## 注意事项

1. **重复 unique_id**：逐行写回时不能用 `flags.loc[uids]`（重复标签会扩展行数），用 dict 映射（`groupby(level=0).first().to_dict()`）
2. **pandas 浮点重写**：读写含 `idr_ratio` 等 float 列的 csv 会损失末位精度（17 位→16 位），写回后与备份 diff 并恢复原始文本
3. **fold 特殊**：`OOD_NewEC` 是 EC7（L1）划分前预留 prior，**不适用**本 skill 口径；func_ec 只保留 L4（L3 是其分类目标）
4. **optimal_ph 警告**：train 仅 ~8.6k 样本且按 UniProt 计频次，L4 级 LongTail 占比 64–78%，区分度有限，L3 级更有意义
5. UniProt API 查询失败的批次按空集合写入缓存，重跑前需人工剔除对应缓存键
6. ref 文件缺 `split` 列时全部按 train 计频次（会高估 LongTail 频次），双文件模式务必带 split 列
