---
name: cath-ted-domain-annotation
description: 为任意任务的 PDB 链样本生成统一的 CATH+TED Topology/Superfamily 域标注，并按 POOD-Benchmark 口径计算 OOD_FoldHoldout / OOD_SuperfamilyHoldout。当需要为数据集（ppis、ligand_binding_affinity、func_prediction、ss_prediction 等）补充或重算结构分类标注与 CATH 留一法 OOD 列时使用。
---

# CATH+TED 域标注与 Holdout OOD 计算

## 描述

本 skill 封装了 POOD-Benchmark 统一的结构域标注流程：**CATH v4.4 实验标注优先，TED putative 标注补缺（并集语义）**，输出链级 Topology（C.A.T）与 Superfamily（C.A.T.H）标签，并可进一步按标准口径计算 `OOD_FoldHoldout` / `OOD_SuperfamilyHoldout`。

## 何时使用

- 任意任务需要为 PDB 链样本补充 Topology / Superfamily 标注
- 需要计算或重算 `OOD_FoldHoldout` / `OOD_SuperfamilyHoldout`（CATH 留一法）列
- 数据集中存在未被 CATH 收录的链（较新 PDB 结构），仅 CATH 口径覆盖率不足时

## 核心口径（各任务必须一致）

| 决策点 | 约定 | 理由 |
|:---|:---|:---|
| 标签来源优先级 | CATH v4.4 实验标注为准，TED 仅补缺（并集） | TED 为 foldseek/foldclass 计算预测，且会漏标 ~6-7% CATH 已有标签，不能替代 |
| TED 域过滤 | 仅保留 chopping 与链 UniProt 区间重叠 **≥50% 域长**的域 | TED 域基于 UniProt 全长，防止链外区域的域混入；25%/75% 敏感性已验证稳健 |
| TED 标签级别 | `H` 级 → C.A.T + C.A.T.H；`T` 级 → 仅 C.A.T | T 级无超家族信息，不可用于 SuperfamilyHoldout |
| Holdout 参考集 | **train-only**，train 自身标注同样用 CATH+TED 并集（对称增强） | 只增强 test 会导致 TED-only 标签假性"未见"，虚增 Holdout |
| Holdout 判定 | any 语义：任一标签未在参考集出现 → True；**无注释 → False** | 与既有任务脚本一致 |
| 样本粒度 | 由输入链清单决定：链级任务一链一行；复合物级任务（如 ligand）先按 PDB 聚合全部链取并集 | 与旧流程 16_cath_ec_ood.py（已归档至 `ligand_binding_affinity/output/scripts_archive/`）/ add_cath_holdout_ood.py 口径一致 |

## 数据来源

| 文件 | 用途 |
|:---|:---|
| `data/CATHv44/cath-classification-data/cath-domain-list.txt` | CATH v4.4 官方域分类（实验标注基线） |
| `data/CATHv44/cath-classification-data/cath-domain-boundaries-seqreschopping.txt` | CATH 域边界（SEQRES 编号，换算为 UniProt 坐标后使用） |
| `data/sifts/pdb_chain_cath.tsv.gz` | PDB+链 → UniProt（SP_PRIMARY，优先来源） |
| `data/sifts/sifts_chain_uniprot.tsv.gz` | PDB+链 → UniProt 及 SP_BEG/SP_END 区间（TED 域重叠过滤依据） |
| `data/ted/ted_api_cache.json` | **项目级共享 TED API 缓存**（跨任务复用，已含 9,889 个 accession） |
| TED REST API `https://ted.cathdb.info/api/v1/uniprot/summary/<acc>` | 缓存未命中时在线查询（~10 req/s 并发安全） |

> Zenodo 的 TED 批量文件（record 13908086）在本机网络常被限流（403）；API+缓存是已验证的主路径。

## 用法

### 1. 标注（cath_ted_annotate.py）

支持两种输入模式（按列名自动识别）：

- **链模式**：输入含 `pdb`、`chain` 列。CATH 名称匹配 + TED 50% 域重叠过滤
- **UniProt 模式**：输入含 `uniprot` 列（样本为全长蛋白，如 AFDB 键数据集）。
  CATH 经 SIFTS 反查（`pdb_chain_cath` ∪ `sifts_chain_uniprot` 的 SP_PRIMARY →
  链 → `cath-domain-list` 名称匹配，取该 accession 全部 PDB 链的并集）；
  TED 按 accession 直接取全部标注域，**不做重叠过滤**（样本即全长序列）

```bash
python .skills/cath-ted-domain-annotation/scripts/cath_ted_annotate.py \
    --chains my_chain_list.csv --out output/my_labels.csv          # 链模式
python .skills/cath-ted-domain-annotation/scripts/cath_ted_annotate.py \
    --chains my_acc_list.csv --out output/my_labels.csv --offline  # UniProt 模式
```

链模式输出在输入列之后追加：

| 列 | 说明 |
|:---|:---|
| `cath_topos` / `cath_sfs` | CATH 实验标注（`;` 分隔集合） |
| `ted_topos` / `ted_sfs` | TED 补缺标注（已经过重叠过滤） |
| `topos` / `sfs` | 并集，**下游 Holdout 计算使用这两列** |
| `uniprot` | 链对应的 UniProt accession |
| `annotation_source` | `cath` / `ted` / `both` / `none` |
| `cath_domains` | CATH 域及其 UniProt 坐标区间，如 `1qqfA00:5-120,200-260`（无映射时域 ID 后为空） |
| `ted_domains` | 保留的 TED 域及其 chopping（UniProt 坐标） |
| `cath_chain_cov` / `aug_chain_cov` | 链 UniProt 区间被 CATH / CATH+TED 标注域覆盖的比例 |

**CATH 域边界**：来自 `data/CATHv44/cath-classification-data/cath-domain-boundaries-seqreschopping.txt`
（SEQRES 编号），经 SIFTS 段映射（`RES_BEG/RES_END ↔ SP_BEG/SP_END`）换算为 UniProt 坐标。
CATH 标签本身仍按名称匹配保留（不受边界影响）；边界用于区间分析与覆盖率统计。
用 `--domains-out domains.csv` 可额外输出域级明细表（含 `mapped_frac`、`in_chain_frac`、
`kept` 列），其中 `in_chain_frac < 0.5` 的 CATH 域基本对应 SIFTS 无该链 UniProt 映射的情形。

关键参数：`--overlap`（默认 0.5）、`--threads`（默认 8）、`--offline`（仅用缓存，不访问 API）、`--cache`（默认项目共享缓存）。

### 2. Holdout 计算（cath_ted_holdout.py）

```bash
# 单文件含 train/val/test 全量（split 列标识，';' 分隔多归属）：
python .skills/cath-ted-domain-annotation/scripts/cath_ted_holdout.py \
    --labels output/my_labels.csv --split-col split --out output/test_flagged.csv

# 或 test 与参考集分文件：
python .skills/cath-ted-domain-annotation/scripts/cath_ted_holdout.py \
    --labels test_labels.csv --ref-labels trainval_labels.csv --out test_flagged.csv
```

输出追加 `OOD_FoldHoldout` / `OOD_SuperfamilyHoldout` 两列（`True`/`False`），仅含 test 行。

> **硬性要求**：各任务 wrapper 写回 splits 时的 holdout 重算**必须调用
> `cath_ted_holdout.py`**，不得自行重实现判定逻辑（参考集范围、any 语义、
> 无注释 → False 等口径以此脚本为唯一实现）。

### 3. 复合物级任务（PDB 级并集）

先自行把链级标注按样本聚合：`topos`/`sfs` 对样本内所有链取并集，再跑 holdout 脚本
（参考 `datasets/ligand_binding_affinity/output/scripts_archive/16_cath_ec_ood.py` 的 groupby 写法）。

## 验证要求

- 抽查 ≥20 条链：标签表中 TED 标签必须是该 accession API 原始返回的子集（重叠过滤后）
- 覆盖率与 Holdout 计数应与下列基准同量级（2026-08-26，ppis / ligand 实测）：

| 任务 | 仅 CATH topo 覆盖 | CATH+TED | FoldHoldout 变化 | SuperfamilyHoldout 变化 |
|:---|:---|:---|:---|:---|
| ppis test (Default) | 92.8% | 97.4% | 64 → 65 | 231 → 235 |
| ligand test | 70.8% | 90.9% | 24 → 38 | 45 → 73 |

> ligand 为 2026-08-26 正式 splits 更新后的终值（含链清单播种修复 +1/+1；另修复了 TED 逗号分隔双标签的解析，全库 82 个域，标签取两方法并集）。

- 双标注样本 CATH vs TED 一致性基准：完全一致 88-91%、有交集 >99%；明显偏离需排查映射错误

## 注意事项

1. **对称增强**：train/val/test 必须在同一次 annotate 运行中处理，参考集才有 TED 补缺
2. TED 基于 AFDB **v4** 模型；本地 v6 结构按 accession 对应即可，无需处理版本差异
3. TED 标签为 putative：接入正式 splits 时在 `DATA_PROCESS.md` 注明"CATH 实验标注优先，TED 补缺"
4. 约 8% accession 不在 TED 中（`not_found`，多为新 UniProt 条目），这些链仍按无注释处理
5. TED 域与 CATH 域边界口径不同，TED 为已有 CATH 标注链"多出"的标签（ppis 333 链）含边界差异成分，解读互补性时需注意
6. TED API 的 `cath_label` 极少数情况下为逗号分隔的双标签（`foldseek,foldclass` 两方法各给一个，全库 82 个域），脚本按并集处理——**直接 `split(".")` 会产生 `3.30.160,3` 之类的伪标签**，曾导致假性 Holdout 翻转，已实现防护
7. **UniProt 模式 caveat**：同一 accession 的不同 PDB 构建体可能覆盖不同区域/突变体，反查并集可能略微过覆盖；pair 级任务（如 ppi_prediction）按双方标签并集判定，进一步放大 any 语义的敏感度
