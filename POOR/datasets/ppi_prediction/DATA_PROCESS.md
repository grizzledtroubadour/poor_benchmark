# PPI Prediction 数据集处理文档

> 记录从原始 PPI 数据到最终蛋白质-蛋白质相互作用预测数据集的完整处理流程，可按本文档从原始数据一次性复现。
> 处理位置：`datasets/ppi_prediction/`（蛋白质-蛋白质相互作用预测任务主目录）
> 数据源：[PINDER](https://github.com/pinder-org/pinder) (release `2024-02`)
> 参考论文：Kovtun et al., *PINDER: The Protein INteraction Dataset and Evaluation Resource*, 2024

## 任务定义与当前规模

本数据集用于**蛋白质-蛋白质相互作用二分类预测**。每个样本是一个蛋白对（蛋白 A / 蛋白 B），标签为 `1`（相互作用）或 `0`（无相互作用）。

当前 splits 规模（以 `splits/` 实际文件统计为准）：

| 文件 | 行数 | 组成（正 bio / XTAL 负 / 随机负） |
|------|------|------|
| `ppi_train.csv` | **58,968** | 29,484 / 10,301 / 19,183 |
| `ppi_val.csv` | **8,394** | 4,197 / 2,104 / 2,093 |
| `ppi_test.csv` | **16,264** | 8,132 / 3,947 / 4,185（另含 OOD 标记列，共 37 列） |

Test 内部：`Default=True` 14,939 / `Default=False`（`OOD_ExtremeLong`，任一链 > 1,000）1,325；正负样本严格 1:1。

> 所有脚本均从任务目录（`datasets/ppi_prediction/`）下运行。被取代的历史脚本保留在 `output/scripts_archive/`（见归档清单）。

---

## 阶段一：源数据获取与解析（PINDER 索引 + 序列获取）

> 目标：获取 PPI 原始注释表与全部蛋白序列。无需下载 PINDER 完整 700 GB 结构数据。

### 1.1 数据来源

位置：`datasets/ppi_prediction/data/`

| 文件 | 用途 | 大小 | 说明 |
|------|------|------|------|
| `data/index.parquet` | 核心入口 | ~119 MB | 每个 dimer 的元数据与 split 信息 |
| `data/metadata.parquet` | 详细注释 | ~100 MB | 结构方法、分辨率、接口统计、标签等 |

**核心数据规模（原始 PINDER 2024-02）**：总样本数 2,319,564（index.parquet）；唯一 PDB ID ~62,706；唯一 UniProt ID 47,098；正样本 (BIO) 1,180,362；负样本候选 (XTAL) 873,090。

**关键字段**：

| 字段 | 类型 | 说明 |
|------|------|------|
| `id` | string | PINDER system ID，如 `7rzb__A1_A0A229LVN5--7rzb__A2_A0A229LVN5` |
| `pdb_id` | string | PDB 四位代码 |
| `uniprot_R` / `uniprot_L` | string | 受体 / 配体的 UniProt ID |
| `cluster_id` | string | 接口聚类 ID，用于去冗余 |
| `split` | string | `train` / `val` / `test` / `invalid`（PINDER 原始划分，本数据集不沿用，见阶段四时间切分） |
| `label` | string | `BIO`=生物学相互作用（正样本），`XTAL`=晶体接触（负样本候选） |
| `complex_type` | string | `homomer` / `heteromer` |
| `length1` / `length2` | int | 受体 / 配体序列长度 |

**手动下载命令**：

```bash
aria2c -x 8 -s 8 -o data/index.parquet     "https://storage.googleapis.com/pinder/2024-02/index.parquet"
aria2c -x 8 -s 8 -o data/metadata.parquet  "https://storage.googleapis.com/pinder/2024-02/metadata.parquet"
```

### 1.2 UniProt 序列获取

**脚本**：`scripts/01_fetch_uniprot_sequences.py`

- 可选地从项目级共享目录 `<repo>/data/uniprot/uniprot_sprot.fasta` 预加载 Swiss-Prot 序列（文件不存在时自动跳过，当前未提供该文件）
- 通过 UniProt REST API (`rest.uniprot.org/uniprotkb/stream`) 批量下载缺失序列，支持断点续传
- 输出：`data/uniprot_sequences_all.fasta`
- 覆盖情况：UniProt API 初始覆盖 PINDER 中 96.33% 的 UniProt ID（45,368 / 47,098），剩余 1,730 个 ID 由下一步回补

### 1.3 缺失序列从 PDB 单体结构回补

**脚本**：`scripts/05_fetch_missing_sequences_from_pdb.py`

- 基于 `output/ppi_pairs_corrected_seq95.csv` 找出缺失 UniProt 序列的蛋白
- 通过 `aria2c` 批量下载对应 PINDER holo 单体 PDB 文件到 `data/pdb_structures_missing/`（3,566 个）
- 使用 Bio.PDB 解析 ATOM 记录并恢复氨基酸序列，追加到 `data/uniprot_sequences_all.fasta`
- 缺失判断为 `p not in seqs or not seqs[p]`（空序列记录视为缺失）
- 自动级联重跑 `04_build_pairs_with_corrected_ranking.py` 与 `09_rebuild_with_xtal_negatives.py`，更新后续 pipeline

**回补效果**：成功回补 1,691 个蛋白的序列；最终数据集中仅 1 条样本因序列缺失被过滤，序列覆盖率达 **99.999%**。

> 修复注记（2026-08-28，已落数据）：train 样本 `4p5x__A1_O75746--4p5x__A2_O75746` 的 `aa_seq1/aa_seq2` 曾为空，根源是 `>O75746` 空 fasta 记录（UniProt stream API 空响应产物，全库唯一）且旧缺失判断只看键是否存在。修复：两个 fasta 中的空记录补入 UniProt 当前 678aa 序列（SLC25A12）；`splits/ppi_train.csv` 该行同步补齐（仅 2 个单元格，备份 `output/backup_before_o75746_fix/`）；05 号脚本缺失判断改为同时检查空序列。修复后三个 split 空链样本均为 0。

---

## 阶段二：标签聚合与计算（蛋白对标签 + representative 选择 + EC 缓存）

> 目标：按无序蛋白对去重并固化正/负标签来源，为每个蛋白对选择最佳 representative 结构，缓存 UniProt EC 注释。
> 线性执行顺序为 02 → 03（阶段三去冗余）→ 04；此处按主题归组。

### 2.1 正/负样本标签来源

PINDER 的 `metadata.label` 字段已提供生物学相关性注释：

| 标签 | 含义 | 用途 |
|------|------|------|
| `BIO` | PRODIGY-cryst 预测为生物学相互作用 | **正样本** (`label=1`) |
| `XTAL` | 预测为晶体接触（非生理相互作用） | **难负样本** (`label=0`) |

负样本构造与 1:1 平衡在阶段四划分时完成（见阶段四 Step 1）。

### 2.2 蛋白对去重（整合所有 split）

**脚本**：`scripts/02_build_unique_protein_pairs.py`

- 整合 train / val / test / invalid 四类数据
- 按无序蛋白对 `{protein_A_id, protein_B_id}` 去重
- 若同一蛋白对同时有 BIO 和 XTAL，保留 BIO
- 输出：`output/ppi_unique_pairs.csv`（74,779 行）

### 2.3 Representative 选择（修正后口径）

**脚本**：`scripts/04_build_pairs_with_corrected_ranking.py`

选择优先级（按顺序）：

1. **是否有 apo / predicted 单体结构**：`high` 质量 apo > `low` 质量 apo > predicted > 仅 holo
2. **有效分辨率更高**：PINDER `metadata.parquet` 仅提供 holo 复合物分辨率，未直接给出 apo/predicted 的分辨率或置信度，暂以 holo 分辨率作为代理（见注意事项）
3. **总长度更长**
4. **日期最新**

输出：`output/ppi_pairs_corrected_uniprot.csv`（74,779 行）、`output/ppi_pairs_corrected_seq95.csv`（65,214 行，**下游主输入**）

### 2.4 UniProt EC 编号缓存

**脚本**：`scripts/08_fetch_uniprot_ec.py`

- 通过 UniProt REST API 批量获取 EC 编号（分号分隔）
- 输入：`data/uniprot_sequences_all.fasta`；输出：`output/uniprot_ec_mapping.csv`（47,059 条）
- 供阶段四 EC 功能 OOD（UniProt + SIFTS 双路并集）使用

---

## 阶段三：去冗余与数据精化（序列去冗余 + 结构选择/staging）

### 3.1 序列同一性去重

**脚本**：`scripts/03_dedup_by_sequence_identity.py`

- 使用 mmseqs2 对 UniProt 序列做 95% identity 聚类
- 若两对蛋白的每条链都落在同一序列 cluster 中，则视为冗余
- 输出：`output/ppi_unique_pairs_seq95_dedup.csv`（63,189 行）

### 3.2 单体结构选择与下载

**脚本**：`scripts/06_download_structures_and_filter_test.py`

- 为每个 representative pair（BIO + XTAL）的两条链选择最佳单体结构，优先级：**apo 单体 > predicted 单体 (AlphaFold) > holo 单体**
- 通过 `aria2c` 批量下载选中的单体 PDB 文件到 `data/pdb_structures_all/`
- 输出：
  - `output/structure_selection.csv`（65,214 行，representative → 选中结构文件映射，含 `both_chains_apo_or_predicted` 标记）
  - `output/ppi_test_structure_filtered.csv`（诊断性产物：按结构可得性过滤后的 test 集；**最终 splits 由 09 号脚本生成**，该过滤逻辑已并入 09）

**结构源分布**（基于 `structure_selection.csv`，65,214 条 representative）：

| R 链源 | L 链源 | 数量 |
|--------|--------|------|
| apo    | apo    | 5,648 |
| apo    | predicted | 3,864 |
| apo    | holo   | 405 |
| predicted | apo | 2,975 |
| predicted | predicted | 48,130 |
| predicted | holo | 642 |
| holo   | apo    | 216 |
| holo   | predicted | 695 |
| holo   | holo   | 2,639 |

R 链：predicted 51,747 / apo 9,917 / holo 3,550；L 链：predicted 52,689 / apo 8,839 / holo 3,686；两条链均为 apo/predicted 的 pair：60,617 / 65,214。选定结构文件共 50,306 个，下载全部成功（`output/structure_download_summary.json`）。

### 3.3 结构 staging 到 pdbs/

**脚本**：`scripts/07_stage_pdbs.py`（2026-08 补缺口新增；历史上该步骤无独立脚本）

- 需要的文件清单 = `splits/ppi_{train,val,test}.csv` 的 `struct_file1`/`struct_file2` 并集（44,178 个；struct_file 列由 09 号脚本按 `structure_selection.csv` 生成）∪ `output/splits_pre_xtal_backup/`（XTAL 重构前历史 splits）引用文件（2,224 个，`--no-pre-xtal` 可关闭）
- 复制规则（幂等）：目标已存在 → 跳过；源存在 → 复制；源缺失 → 记入 missing 清单
- **验证结果（2026-08）**：`--dry-run` 生成的期望清单（46,402 个）与现有 `pdbs/` 文件清单双向 diff 均为 0；小批量真实复制与 `pdbs/` 对应文件逐字节一致
- 注意：当前 `data/pdb_structures_all/` 仅保留 10,690 个文件（多数历史文件已移出），全新环境重放需先运行 06 重新下载全部 50,306 个选定结构；structure_selection 选定但从未被任何 splits 版本引用的 3,904 个文件不属于 staging 范围

### 3.4 结构完整性校验（诊断）

**脚本**：`scripts/diag_validate_structures.py`（纯诊断，不进主流程）

- 使用 Bio.PDB 解析 `data/pdb_structures_all/` 下所有 PDB 文件，检查空文件、有效 ATOM 残基、可解析性
- 输出：`output/structure_validation_summary.json`、`output/structure_validation_failed.csv`
- 校验结果：全部 **50,306** 个下载的 PDB 文件通过完整性校验（empty=0, no_atoms=0, parse_error=0）

---

## 阶段四：数据划分与 OOD 标注

### Step 1: XTAL 难负样本 + 时间切分重建 splits

**脚本**：`scripts/09_rebuild_with_xtal_negatives.py`

负样本由两部分组成，在 train/val/test 划分**之后**按每个 split 独立构造：

1. **XTAL 难负样本**：PINDER 标注的 XTAL 对全部保留为负样本。这些是"长得像界面但非生物相互作用"的难负样本，来自真实复合物结构，蛋白池与 BIO 高度重叠（65.2%），消除了"随机配对 vs 真实复合物"的浅层捷径。
2. **随机配对补齐**：当 XTAL 数量不足时，从同 split 的**全部蛋白池**（BIO + XTAL 蛋白 ID 合并去重）中随机抽取两个蛋白组成负样本对。**两个蛋白允许相同**（模拟同源二聚体负样本），消除 `seq1==seq2 → label=1` 的 homodimer shortcut。排除所有在 `output/ppi_unique_pairs.csv` 中标记为 BIO 或 XTAL 的已知蛋白对。
3. **1:1 比例**：每个 split 中负样本总量（XTAL + 随机）等于正样本数量。
4. **`pair_type` 列**：每行标注样本类型——`bio`（正样本）、`xtal`（XTAL 负样本）、`random`（随机配对负样本），train/val/test 三个文件均含此列。

**划分原则**：

| 原则 | 说明 |
|------|------|
| **时间切分** | 按 BIO 样本的 `release_date` 排序，累积 70% / 80% 处切分（**70:10:20**），XTAL 用相同时间节点划分 |
| **全局长度硬剔除** | 划分前剔除任一链长度 > 2,000 的样本（不进入任何 split） |
| **train/val 长度限制** | 仅保留单链 ≤ 1,000 的样本（> 1,000 的归入 test 作为 `OOD_ExtremeLong`） |
| **test 结构过滤** | 仅保留双链均有 apo 或 predicted 结构的样本（inline 完成） |
| **正负样本平衡** | 每 split 中 BIO 正样本与 XTAL+随机负样本 1:1 |

- 输入：`output/ppi_pairs_corrected_seq95.csv`、`output/structure_selection.csv`、`output/ppi_unique_pairs.csv`、`data/uniprot_sequences_all.fasta`
- 随机配对 seed = 42/142/242；test 集同时写入 `Default` 列（任一链 > 1,000 → `False`）
- 输出：`splits/ppi_{train,val,test}.csv`、`output/xtal_rebuild_summary.json`

> 历史注记：旧版由 `21_resample_negative_pairs.py` 纯随机配对生成全部负样本（不允许自配对），导致 homodimer shortcut 与"随机配对 vs 真实复合物"捷径，模型 AUROC 虚高；09 号脚本（21 → 23 → 09 迭代）引入 XTAL 难负样本 + 允许自配对后解决。旧脚本已归档，XTAL 重构前的旧 splits 备份在 `output/splits_pre_xtal_backup/`。

### Step 2: OOD 主体重算（极端长度 / Default / NaN 收尾为本地独有逻辑）

**脚本**：`scripts/10_recompute_all_ood.py`（基于最终 `splits/ppi_{train,val,test}.csv` 一次性计算）

- **OOD_ExtremeShort**：蛋白对中**任意一条链**长度 < 60 → `True`（Default 子集 715 条；另有 16 条"一短一长"混合对属 Default=False 行，此列随收尾置 NaN）
- **OOD_ExtremeLong**：任意一条链长度 > 1,000 → `True`（1,325 条）；`OOD_ExtremeLong=True` 的行 `Default=False`（14,939 条 Default=True）
- PPI 任务对极端长度样本**只做标记、不剔除**（train/val 在 Step 1 已限制 ≤ 1,000）
- **Default=False 行收尾**：仅保留 `OOD_ExtremeLong`（预挑出原因列），其余全部 OOD 列（含 `OOD_ExtremeShort`、`seq_Redundancy_*`、`TM-score_*`）置 NaN，表示"未计算"——与其他任务 Default=False 行惯例统一。备份：`output/backup_before_nondefault_full_nan/`
- 脚本内 mmseqs2 / Foldseek / IDR / EC 各段已标注 **[superseded]**（见 Step 3–6，重算以 skill 链为准）

### Step 3: 序列同源性（seq_Redundancy_* / OOD_Orphan，仅 Default 子集）

**正式口径**：`.skills/homology-ood-annotation/scripts/seq_homology_ood.py`

- 参考库 = train+val 全部唯一蛋白序列，query = test 蛋白；mmseqs2 `easy-search -s 7.5`（默认 e-value ≤ 1e-3 显著性门槛）
- 缓存：`output/ood_mmseqs_test_vs_trainval_easysearch.m8`（skill 格式：query 键 `unique_id|aa_seq{1,2}`，可 `--m8-cache` 重放）
- **seq_Redundancy_XX**（XX ∈ 90..30，步长 10）：任意一条链 max identity < XX → `True`
- **OOD_Orphan**：任意一条链在 train+val 中无显著 hit → `True`

> 历史注记：旧式 `mmseqs search --max-seqs 10 -e 100 -c 0.0` 参数曾严重低估低阈值 OOD（seq_Redundancy_30 仅 515、Orphan 仅 211），2026-08-27 已改用 easy-search 经 skill 脚本重算写回（备份 `output/backup_before_easysearch_unify/`）；旧缓存 `output/ood_mmseqs_test_vs_trainval.m8` 仅供溯源。本地 10 号脚本的 mmseqs 段已标注 superseded。

### Step 4: 结构相似性（TM-score_*，仅 Default 子集）

**正式口径**：`.skills/homology-ood-annotation/scripts/struct_homology_ood.py`

- query = test 蛋白单体结构（`pdbs/`），target = train+val 蛋白单体结构库
- Foldseek `easy-search`（默认 e-value ≤ 10 门槛）`--format-output query,target,alntmscore`
- 缓存：`output/ood_foldseek_test_vs_trainval_easysearch.m8`（skill 格式，可 `--m8-cache` 重放）
- **TM-score_X.X**（0.9..0.3，步长 0.1）：任意一条链 max alntmscore < 阈值 → `True`

> 历史注记：旧式 `foldseek search -a -e inf` 保留全部弱 hit，曾抬高 max alntmscore、压低低阈值 OOD（TM-score_0.9: 4,827→8,172、TM-score_0.3: 1,610→943），2026-08-27 已统一为 easy-search 重算（备份 `output/backup_before_foldseek_easysearch_unify/`）；旧缓存 `output/ood_foldseek_tmscore.csv` 仅供溯源。本地 10 号脚本的 foldseek 段已标注 superseded。

### Step 5: 内在无序区域（OOD_IDR，仅 Default 子集）

**统一工具**：`.skills/ood-annotation-toolkit/scripts/idr_ood.py`（metapredict v3；本任务口径 = `--mode region-ratio --clean map --threshold 0.3`）

- 对 test 集所有唯一蛋白序列做残基级 disorder 预测，缓存于 `output/idr_predictions.csv`（7,418 条序列）
- 非标准氨基酸映射后预测（B→N, Z→Q, J→L, U→C, O→K，其余→A）
- 单链判定：最长连续 disorder > 0.5 区域长度占序列全长 > 0.3 → IDR-containing（**region-ratio 口径**）
- **OOD_IDR**：任意一条链为 IDR-containing → `True`
- 本地 10 号脚本的 IDR 段已标注 superseded

### Step 6: EC 功能 OOD（OOD_NewEC_* / OOD_LongTail_EC_*，仅 Default 子集）

**统一口径**：`.skills/ec-function-ood-annotation` 的批量驱动 `add_unified_newec.py` / `add_unified_longtail_ec.py`

**EC 收集（双路并集）**：

- **UniProt REST API** 缓存：`scripts/08_fetch_uniprot_ec.py` 产物 `output/uniprot_ec_mapping.csv`（47,059 条）
- **SIFTS** 链级 EC：`<repo>/data/sifts/sifts_chain_ec.tsv.gz`，按 PDB ID + chain 补充
- bio/xtal 行：解析 `unique_id` 得 (pdb, chain, uniprot)，取 UniProt 缓存 ∪ SIFTS 链级并集；random 负样本：仅 UniProt 缓存
- 层级解析：段数 ≥ 3 → L3 取前三段；段数 = 4 → L4 完整编号；允许部分 `-`（如 `3.4.22.-` → L3=`3.4.22`）

**统一口径（all-not-in）**：

| 列 | 参考集 | 判定 |
|------|--------|------|
| `OOD_NewEC_L4` / `OOD_NewEC_L3` | train+val 出现的 EC 集合 | 样本有 EC 注释且其**所有** EC（对应层级）均未在参考集出现 → `True`；无注释 → `False` |
| `OOD_LongTail_EC_L3_le5` / `OOD_LongTail_EC_L3_le10` / `OOD_LongTail_EC_L4_le5` / `OOD_LongTail_EC_L4_le10` | **train** 逐样本频次 | 样本有 EC 注释且其**所有** EC（对应层级）train 频次都 < 5 / < 10（**含 0**）→ `True`；无注释 → `False` |

> 历史注记：旧列 `OOD_LongTail_EC_L4`（train+val 频次 < 10、any 语义、655 条）已被四列取代并删除；`OOD_NewEC_L4` 由 any 口径（245 条）就地重算为 all-not-in（203 条）。备份：`output/backup_before_newec_unify/`、`output/backup_before_longtail_unify/`、`output/backup_before_longtail_allnotin/`。本地 10 号脚本的 EC 段已标注 superseded。

### Step 7: CATH+TED 域标注（OOD_FoldHoldout / OOD_SuperfamilyHoldout，仅 Default 子集）

**脚本**：`scripts/11_ted_augment_cath_ood.py`（wrapper：标注 → 对级聚合 → 调 skill 的 `cath_ted_holdout.py` 判定 → 写回 splits）

**口径**（UniProt accession 模式，**对级**：`protein_A_id` / `protein_B_id` 双方标签取并集）：

- 标签 = CATH v4.4 实验标注（经 SIFTS UniProt→PDB 链反查并集）∪ TED putative 标注（REST API，缓存于项目级 `<repo>/data/ted/ted_api_cache.json`，按全长不过滤；H 级提供 C.A.T.H，T 级仅 C.A.T）
- 参考集 = **Train** 的 CATH+TED 并集（train/val/test 对称增强；当前参考集 1,262 个 Topology / 4,711 个 Superfamily）
- any 语义：蛋白对双方并集中任一 Topology / Superfamily 未在参考集出现 → `True`；双方均无注释 → `False`；`Default=False` 行置 NaN

**流程**：

1. `cath_ted_annotate.py --chains output/ppi_ted_aug_acc_list.csv --out output/ppi_acc_cath_ted_labels.csv --offline`（UniProt 模式，39,455 个 accession 全部缓存命中）
2. 聚合为对级标签表 `output/ppi_pair_cath_ted_labels.csv`（`unique_id, topos, sfs, split`，样本级标签 = 双方 accession 并集）
3. `cath_ted_holdout.py --labels <对级标签表> --split-col split --out output/ppi_pair_holdout_flags.csv`（skill 统一判定）
4. 写回 `splits/ppi_test.csv`（`Default=False` 行置 NaN）

**标注来源分布**（39,455 个 accession）：both 19,756 / ted-only 14,066 / cath-only 1,365 / none 4,268。Topology 覆盖率：仅 CATH 53.5% → CATH+TED **89.2%**。

**Holdout 结果**（test，Default=True n=14,939）：`OOD_FoldHoldout` 44 条；`OOD_SuperfamilyHoldout` 265 条；双方均无注释样本（→ False）950 条。

**输出文件**：`output/ppi_acc_cath_ted_labels.csv`（accession 级标签）、`output/ppi_ted_aug_acc_list.csv`（键清单，历史一次性生成）、`output/ppi_pair_cath_ted_labels.csv` / `output/ppi_pair_holdout_flags.csv`（中间产物）、`output/backup_before_ted_aug/`（更新前备份）。

### Step 8: NaN 惯例与 InD 收尾

**统一工具**：`.skills/ood-annotation-toolkit/scripts/finalize_ood_columns.py`

- **NaN 惯例**：`Default=False` 行（极端长度预挑出样本）除 `OOD_ExtremeLong`（预挑出原因列）外，**全部** OOD 列均为 NaN（未计算，与 `False` 区分）——包括 `OOD_ExtremeShort` 与搜索型同源列 `seq_Redundancy_*` / `TM-score_*`。备份：`output/backup_before_nondefault_nan/`、`output/backup_before_nondefault_full_nan/`
- **`InD`**：纯 in-distribution 标识——所有 `OOD_*` + `seq_Redundancy_*` + `TM-score_*` 标记均为 `False`（NaN 视为 `False`）时为 `True`，由 `finalize_ood_columns.py --mode ind` 在所有 OOD 列齐备后重算。当前 **1,944 / 16,264（12.0%）**。备份：`output/backup_before_ind/`
- 体检：`.skills/ood-annotation-toolkit/scripts/check_splits.py`

### OOD 统计摘要（以当前 `splits/ppi_test.csv` 实测为准）

test 共 **16,264** 行（正 8,132 / 负 8,132），`Default=True` 14,939 / `Default=False` 1,325。

| OOD 列 | True 数 | 备注 |
|--------|--------:|------|
| `seq_Redundancy_90` | 10,655 | 仅 Default 子集，Default=False 行 NaN |
| `seq_Redundancy_80` | 9,856 | 仅 Default 子集，Default=False 行 NaN |
| `seq_Redundancy_70` | 9,027 | 仅 Default 子集，Default=False 行 NaN |
| `seq_Redundancy_60` | 8,098 | 仅 Default 子集，Default=False 行 NaN |
| `seq_Redundancy_50` | 6,982 | 仅 Default 子集，Default=False 行 NaN |
| `seq_Redundancy_40` | 5,496 | 仅 Default 子集，Default=False 行 NaN |
| `seq_Redundancy_30` | 3,653 | 仅 Default 子集，Default=False 行 NaN |
| `OOD_Orphan` | 3,017 | 仅 Default 子集，Default=False 行 NaN |
| `TM-score_0.9` | 8,172 | 仅 Default 子集，Default=False 行 NaN |
| `TM-score_0.8` | 5,345 | 仅 Default 子集，Default=False 行 NaN |
| `TM-score_0.7` | 3,771 | 仅 Default 子集，Default=False 行 NaN |
| `TM-score_0.6` | 2,843 | 仅 Default 子集，Default=False 行 NaN |
| `TM-score_0.5` | 2,186 | 仅 Default 子集，Default=False 行 NaN |
| `TM-score_0.4` | 1,510 | 仅 Default 子集，Default=False 行 NaN |
| `TM-score_0.3` | 943 | 仅 Default 子集，Default=False 行 NaN |
| `OOD_ExtremeShort` | 715 | 仅 Default 子集，Default=False 行 NaN（16 条混合对随之置空） |
| `OOD_ExtremeLong` | 1,325 | = Default=False 行数（唯一保留的预挑出原因列） |
| `OOD_IDR` | 3,743 | 仅 Default 子集，Default=False 行 NaN |
| `OOD_NewEC_L4` | 203 | 仅 Default 子集，Default=False 行 NaN |
| `OOD_NewEC_L3` | 10 | 仅 Default 子集，Default=False 行 NaN |
| `OOD_LongTail_EC_L3_le5` | 17 | 仅 Default 子集，Default=False 行 NaN |
| `OOD_LongTail_EC_L3_le10` | 72 | 仅 Default 子集，Default=False 行 NaN |
| `OOD_LongTail_EC_L4_le5` | 444 | 仅 Default 子集，Default=False 行 NaN |
| `OOD_LongTail_EC_L4_le10` | 698 | 仅 Default 子集，Default=False 行 NaN |
| `OOD_FoldHoldout` | 44 | 仅 Default 子集 |
| `OOD_SuperfamilyHoldout` | 265 | 仅 Default 子集 |
| `InD` | 1,944 | 12.0% |

---

## 输出文件汇总

```
datasets/ppi_prediction/
├── DATA_PROCESS.md                      # 本文档
├── scripts/
│   ├── 01_fetch_uniprot_sequences.py    # 阶段一：从 UniProt API 获取序列
│   ├── 02_build_unique_protein_pairs.py # 阶段二：整合所有 split 按蛋白对去重
│   ├── 03_dedup_by_sequence_identity.py # 阶段三：mmseqs2 95% 序列同一性去重
│   ├── 04_build_pairs_with_corrected_ranking.py  # 阶段二：修正 representative 选择
│   ├── 05_fetch_missing_sequences_from_pdb.py    # 阶段一：从 PINDER PDB 回补缺失 UniProt 序列（级联重跑 04/09）
│   ├── 06_download_structures_and_filter_test.py # 阶段三：单体结构选择、下载与测试集结构过滤
│   ├── 07_stage_pdbs.py                 # 阶段三：选定结构 staging 到 pdbs/（2026-08 补缺口新增）
│   ├── 08_fetch_uniprot_ec.py           # 阶段二：UniProt EC 编号缓存
│   ├── 09_rebuild_with_xtal_negatives.py  # 阶段四：XTAL 难负样本 + 随机补齐 + 时间切分重建 splits
│   ├── 10_recompute_all_ood.py          # 阶段四：OOD 主体（mmseqs/foldseek/IDR/EC 段 superseded；长度列/Default/NaN 收尾保留）
│   ├── 11_ted_augment_cath_ood.py       # 阶段四：CATH+TED 域标注（OOD_FoldHoldout / OOD_SuperfamilyHoldout）
│   └── diag_validate_structures.py      # 诊断：结构完整性校验（不进主流程）
├── data/                                # 原始数据与任务级缓存
│   ├── index.parquet / metadata.parquet # PINDER 2024-02 注释表
│   ├── uniprot_sequences_all.fasta      # 全部蛋白序列
│   ├── pdb_structures_all/              # 下载的单体 PDB 结构（历史 50,306 个，当前保留 10,690 个，见 3.3）
│   └── pdb_structures_missing/          # 序列回补用 holo 单体 PDB（3,566 个）
├── output/                              # 中间结果、缓存与备份
│   ├── structure_selection.csv          # representative → 单体结构选择（65,214 行）
│   ├── uniprot_ec_mapping.csv           # UniProt → EC 编号缓存（47,059 条）
│   ├── idr_predictions.csv              # metapredict IDR 预测缓存（7,418 条序列）
│   ├── ood_mmseqs_test_vs_trainval_easysearch.m8     # mmseqs2 缓存（skill 格式，可重放）
│   ├── ood_foldseek_test_vs_trainval_easysearch.m8   # Foldseek 缓存（skill 格式，可重放）
│   ├── ood_mmseqs_test_vs_trainval.m8 / ood_foldseek_tmscore.csv  # 旧参数缓存（仅供溯源）
│   ├── ppi_acc_cath_ted_labels.csv      # accession 级 CATH+TED 标签
│   ├── scripts_archive/                 # 被取代的历史脚本（保留原文件名）
│   ├── backup_before_*/                 # 各次 splits 写回前备份
│   └── splits_pre_xtal_backup/          # XTAL 重构前的旧 splits
├── pdbs/                                # 最终数据集对应的结构文件（46,402 个）
└── splits/                              # 最终数据集划分文件
    ├── ppi_train.csv                    # 58,968 行，9 列（无 OOD 列）
    ├── ppi_val.csv                      # 8,394 行，9 列（无 OOD 列）
    └── ppi_test.csv                     # 16,264 行，37 列（含全部 OOD 标记）
```

**列说明**（train/val/test 共有）：

| 列 | 说明 |
|------|------|
| `unique_id` | 蛋白对唯一 ID（bio/xtal 为 PINDER representative_id，随机负样本为 `neg_xxxx`） |
| `protein_A_id` / `protein_B_id` | 两条链的 UniProt ID |
| `aa_seq1` / `aa_seq2` | 两条链的氨基酸序列 |
| `struct_file1` / `struct_file2` | 两条链对应的单体 PDB 文件名（仅文件名，位于 `pdbs/`） |
| `label` | 样本级标签，`1`=相互作用，`0`=无相互作用 |
| `pair_type` | 样本类型（`bio` / `xtal` / `random`） |

**`splits/ppi_test.csv` 完整列清单**（37 列，以实际表头为准）：

```
unique_id, protein_A_id, protein_B_id, aa_seq1, aa_seq2,
struct_file1, struct_file2, label, Default, pair_type,
OOD_IDR, OOD_NewEC_L4, OOD_ExtremeShort, OOD_ExtremeLong,
seq_Redundancy_90, seq_Redundancy_80, seq_Redundancy_70, seq_Redundancy_60,
seq_Redundancy_50, seq_Redundancy_40, seq_Redundancy_30, OOD_Orphan,
TM-score_0.9, TM-score_0.8, TM-score_0.7, TM-score_0.6, TM-score_0.5,
TM-score_0.4, TM-score_0.3,
OOD_FoldHoldout, OOD_SuperfamilyHoldout,
OOD_LongTail_EC_L3_le5, OOD_LongTail_EC_L3_le10,
OOD_LongTail_EC_L4_le5, OOD_LongTail_EC_L4_le10,
OOD_NewEC_L3, InD
```

---

## 数据质量检查

1. **泄露检查**：Train ∩ Val、Train ∩ Test、Val ∩ Test 的 `unique_id` 交集均为空（已核对，0 重叠）；随机负样本排除了 `output/ppi_unique_pairs.csv` 中全部已知 BIO/XTAL 蛋白对。
2. **分布一致性**：三个 split 正负样本严格 1:1；负样本中 XTAL 占比 train 34.9% / val 50.1% / test 48.5%（XTAL 数量受时间段内可得性限制，不足部分由随机配对补齐）。
3. **长度约束**：划分前全局剔除任一链 > 2,000 的样本；train/val 单链 ≤ 1,000；test 保留 > 1,000 样本并标记 `OOD_ExtremeLong` / `Default=False`。
4. **结构完整性**：`data/pdb_structures_all/` 下载的 50,306 个 PDB 全部通过 Bio.PDB 解析校验（诊断脚本 `diag_validate_structures.py`）。
5. **序列覆盖率**：99.999%（仅 1 条样本因序列缺失被过滤）。
6. **splits 写回验证**：历次 OOD 列重算/追加均核对行数、列序、`unique_id` 与备份一致（备份见 `output/backup_before_*/`）。
7. **pdbs/ 清单验证（2026-08）**：`07_stage_pdbs.py --dry-run` 期望清单（当前 splits 44,178 ∪ pre-XTAL 2,224 = 46,402）与 `pdbs/` 实际文件清单双向 diff 为 0。

---

## 诊断与质控工具

| 脚本 | 用途 |
|------|------|
| `scripts/diag_validate_structures.py` | Bio.PDB 解析校验 `data/pdb_structures_all/` 全部 PDB 文件（空文件 / 无 ATOM / 解析错误），产出 `output/structure_validation_summary.json` 与失败清单 |

---

## 归档清单（`output/scripts_archive/`）

| 脚本 | 归档原因 |
|------|----------|
| `02_build_ppi_pairs.py` / `03_dedup_ppi_pairs.py` | 初版配对/去冗余，已被 02/03 取代 |
| `06_build_pairs_with_custom_ranking.py` | 自定义 ranking 试验版，已被 04（修正口径）取代 |
| `07_build_negative_pairs.py` / `08_time_split_dataset.py` | 旧负样本/时间划分链路，已被 09 取代（05 的级联重跑链已同步改为 04 → 09） |
| `15_mark_ood.py` / `17_build_ood_splits.py` / `18_compute_ood_idr.py` / `19_compute_ec_ood.py` / `22_recompute_homology_newfold_ood.py` | 分散的旧 OOD 脚本，已被 10 一站式重算 + skill 统一口径取代 |
| `16_fix_negative_structure_labels.py` / `20_add_pdb_file_columns.py` | 一次性列修复，修复已落数据 |
| `21_resample_negative_pairs.py` | 纯随机负样本旧方案（homodimer shortcut），已被 09 的 XTAL+随机补齐取代 |

---

## 注意事项

1. **与 `ppis_prediction/` 的区别**：`ppis_prediction/` 面向**结合界面预测**（残基级标签），`ppi_prediction/` 面向**蛋白对相互作用预测**（样本级二分类）。
2. **分辨率代理**：PINDER 元数据仅含 holo 复合物分辨率，representative 选择暂以 holo 分辨率作为 apo/predicted 质量的代理；解析 apo 真实分辨率与 AlphaFold pLDDT 是后续改进项。
3. **`Default` 口径**：本任务 `Default=False` 仅由 `OOD_ExtremeLong` 触发（1,325 行），极端长度样本只做标记、不从任何 split 剔除。`Default=False` 行仅保留 `OOD_ExtremeLong`，其余全部 OOD 列（含 `seq_Redundancy_*` / `TM-score_*` / `OOD_ExtremeShort`）按 NaN 惯例置空——与其他任务 Default=False 行口径统一。
4. **结构搜索/序列搜索口径**：2026-08-27 起 mmseqs2 与 foldseek 均统一为 easy-search（skill 脚本重算写回）；本地 10 号脚本对应段落已标注 superseded，旧缓存仅供溯源。
5. **EC 口径**：`OOD_NewEC_*` / `OOD_LongTail_EC_*` 统一为 all-not-in（NewEC 参考 train+val，LongTailEC 参考 train 频次且含 0），旧 `OOD_LongTail_EC_L4`（655 条）已删除；统一口径由 `.skills/ec-function-ood-annotation` 的批量驱动脚本落盘。
6. **TED API 逗号分隔双标签**（foldseek/foldclass）已按并集解析，无逗号伪标签残留；TED 缓存全部命中，标注过程离线可复现。
7. **`data/pdb_structures_all/` 现状（2026-08 核实）**：目录当前仅保留 10,690 个 PDB（含 6,786 个 splits 引用文件与 3,904 个选定但从未被任何 splits 版本引用的文件），多数历史文件已移出；`pdbs/`（46,402 个）为最终权威结构目录。全新环境重放 staging 前需先重跑 06 重新下载。
8. **历史 splits 依赖**：`07_stage_pdbs.py` 的完整 46,402 文件清单包含 `output/splits_pre_xtal_backup/` 引用的 2,224 个历史文件；删除该备份目录后 staging 只能复现当前 splits 引用的 44,178 个文件。
