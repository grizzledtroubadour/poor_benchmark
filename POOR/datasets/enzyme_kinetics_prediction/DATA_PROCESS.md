# 酶动力学预测（kcat）数据集处理文档

> 本文档记录 `datasets/enzyme_kinetics_prediction/` 从 CatPred-DB 原始数据到最终划分的一次性复现流程。所有统计数字以当前 `splits/*.csv` 实际内容重算为准。

## 任务定义与当前规模

数据来源为 **CatPred-DB** 的 **kcat** 回归任务（label = `log10(kcat)`，结构为 AlphaFold v6 预测模型，样本键为 UniProt accession，同一 accession 可多底物多行）：

- Zenodo：`https://zenodo.org/records/14775076`（约 1.1 GB，含处理后数据集、PDB 记录、原始 BRENDA/SABIO CSV）
- GitHub：`https://github.com/maranasgroup/CatPred-DB`（约 397 MB）

当前 `splits/` 规模（实测）：

| 文件 | 行数 | UniProt 数 | 列数 | 说明 |
|------|-----:|-----------:|-----:|------|
| `splits/kcat_train.csv` | 15,479 | 4,888 | 6 | 训练集（6 核心列） |
| `splits/kcat_val.csv` | 2,150 | 662 | 6 | 验证集（同 train） |
| `splits/kcat_test.csv` | 5,002 | 1,486 | **37** | 测试集（6 核心列 + 31 OOD/标注列；Default=True 4,412 行，极端长度 590 行） |

处理流程分四个阶段：

- **阶段一**：源数据获取与解析（CatPred-DB 下载、AlphaFold 结构、PDB 过滤、Q9Y233 截断）
- **阶段二**：标签聚合与计算（反应物过滤 + ECFP4 指纹）
- **阶段三**：去冗余与数据精化（本任务无序列聚类去冗余；精化措施见该节说明）
- **阶段四**：数据划分与 OOD 标注（长度 OOD 划分 → 标准列 → 全部 OOD 维度 → NaN/InD 收尾）

## 目录结构

```
enzyme_kinetics_prediction/
├── data/                                   # 任务级原始数据、中间缓存与处理日志（*.log）
│   ├── catpred-db-github/                  # GitHub 仓库内容
│   ├── zenodo/                             # Zenodo 解压内容
│   ├── kcat_filtered_with_pdb.csv          # PDB 过滤后数据（22,736 行）
│   └── kcat_with_reactant_ecfp4.csv        # 反应物过滤 + ECFP4 指纹（22,690 行）
├── scripts/                                # 数据处理脚本（01-08，按执行顺序编号）
├── splits/                                 # 最终划分文件（仅 3 个）
├── pdbs/                                   # AlphaFold v6 结构文件（解压后 .pdb，7,058 个）
├── output/                                 # 中间结果、统计摘要、OOD 子集、备份
│   ├── scripts_archive/                    # 已归档的历史脚本
│   └── backup_*/                           # 各次 splits 写回前的备份
└── DATA_PROCESS.md                         # 本文档
```

> **注意**：本任务处理日志以 `data/*.log` 形式存放，未单设 `logs/`；`pdbs/` 与 `splits/` 存放可直接用于模型训练与评测的最终产物。

---

# 阶段一：源数据获取与解析

> 目标：从 CatPred-DB 原始发布得到「有结构」的 kcat 数据集（**22,736** 行，7,057 个 UniProt）。

## 1.1 数据下载（无脚本环节）

- GitHub 仓库克隆至 `data/catpred-db-github/`（日志 `data/catpred-db-github_clone.log`）
- Zenodo 解压至 `data/zenodo/`

## 1.2 AlphaFold 结构获取

- **脚本**：`scripts/01_download_alphafold_missing.py`
- **输入**：`data/kcat_missing_alphafold_uniprots.txt`
- **输出**：项目级 `data/alphafold_structures/AF-*-F1-model_v6.pdb.gz`；日志 `data/alphafold_download_progress.log`
- **说明**：kcat 数据集共涉及 **7,232** 个唯一 UniProt。项目级结构库初始覆盖 4,522 个；经 AlphaFold DB API 补下载（8 线程、断点续传）后最终覆盖 **7,057** 个（97.6%）。

## 1.3 PDB 过滤与准备

- **脚本**：`scripts/02_process_kcat_with_pdbs.py`
- **输入**：`data/catpred-db-github/datasets/processed/kcat_max_wt_singleSeqs_wpdbs.csv` + 项目级 `data/alphafold_structures/`
- **输出**：`data/kcat_filtered_with_pdb.csv`（22,736 行）、`pdbs/*.pdb.gz`、`data/kcat_uniprots_removed_no_pdb.txt`、`data/kcat_pdb_sequence_consistency.csv`、`data/kcat_substrate_stats.json`
- **说明**：
  - 保留有本地 PDB 文件的样本：**22,736** 行、**7,057** 个 UniProt（无结构剔除 461 个 UniProt）
  - 序列一致性检查：7,045 / 7,057（99.83%）完全匹配
  - 底物统计：有效 `reactant_smiles` 覆盖 99.8%；唯一 `reaction_smiles` 12,114；主要 EC 大类为 Hydrolases(3)、Oxidoreductases(1)、Transferases(2)、Lyases(4)；酶+完整 reactant_smiles 完全重复 0 行

- **脚本**：`scripts/03_decompress_pdbs_and_truncate_q9y233.py`
- **输入**：`pdbs/*.pdb.gz`、`data/kcat_filtered_with_pdb.csv`
- **输出**：`pdbs/*.pdb`（解压）、Q9Y233 截断结构、`data/q9y233_truncation_mapping.json`
- **说明**：解压全部 `.pdb.gz` 为 `.pdb`；Q9Y233 的 AlphaFold 模型与数据集序列不一致（外部数据源补救），在 PDB 序列中滑窗定位数据集序列并截断保留对应 **779** 个残基。当前 `pdbs/` 中截断版本已替换为 `AF-Q9Y233-F1-model_v6.pdb`，原始全长保留为 `AF-Q9Y233-F1-model_v6_full.pdb`（`pdbs/` 共 7,058 个文件 = 7,057 个 UniProt + 1 个全长备份）。

---

# 阶段二：标签聚合与计算

> 目标：得到「有反应物 SMILES、有 ECFP4 指纹」的干净 kcat 数据集（**22,690** 行）。

## 2.1 反应物过滤与 ECFP4 指纹

- **脚本**：`scripts/04_compute_reactant_ecfp4.py`
- **输入**：`data/kcat_filtered_with_pdb.csv`
- **输出**：`data/kcat_with_reactant_ecfp4.csv`（22,690 行，含 `ecfp4_bit_0`~`ecfp4_bit_2047` 与 `ecfp4` 01 字符串）、`output/kcat_ecfp4_summary.json`；日志 `data/compute_reactant_ecfp4.log`
- **说明**：
  - 仅使用 `reactant_smiles`，排除 46 行无反应物数据，最终保留 **22,690** 行（7,051 个 UniProt）
  - RDKit Morgan 指纹（半径=2，2048 位）逐反应物计算，多分子按位 OR 合并；平均每个合并指纹激活 **57.6** 个 bit
  - 历史上的 Step 3（更新 `kcat_random_*` 划分）依赖已归档的 random splits，属死代码，已从脚本中删除

## 2.2 目标值与长度分布（22,690 行，实测）

- label = `log10_value`（log10(kcat)），上游 CatPred-DB 已将其截断（删失）至 **[-6, 6]**（BRENDA 原始范围 2.24e-11 ~ 9.3e8 s⁻¹）
- 均值 **0.89**，标准差 **1.65**，最小 -6.00，最大 6.00；该删失对值域 LongTail 的影响见阶段四 4.6

| 长度范围 | 数量 | 占比 |
|----------|-----:|-----:|
| < 60 | 60 | 0.3% |
| 60 – 1000 | 22,041 | 97.1% |
| > 1000 | 589 | 2.6% |

---

# 阶段三：去冗余与数据精化

> 本任务**不做序列聚类去冗余**（同一 UniProt 的多底物行是任务设计的一部分）；样本级防泄漏由阶段四的 **UniProt 级划分**保证。本阶段的精化措施分布在其上下游，集中说明如下：

| 精化措施 | 位置 | 数量 |
|----------|------|------|
| 无 PDB 结构剔除 | 阶段一 1.3 | 461 个 UniProt |
| 无反应物 SMILES 剔除 | 阶段二 2.1 | 46 行 |
| 正常长度 val/test 反应物相似度 < 0.6 的 UniProt 整体移入 train | 阶段四 4.1 | 69 个 UniProt |
| 极端长度且与 train 反应物相似度 < 0.6 的 UniProt 整体剔除 | 阶段四 4.1 | 21 行（5 个 UniProt） |
| 单链长度 > 2000 剔除（全项目统一规则，无独立脚本） | 阶段四 4.7 | 38 行（10 个 UniProt） |
| 列名 schema 迁移（`pdb_file`→`struct_file`、`labels`→`label`、删 `file_type`，全项目统一步骤，无本任务脚本） | 首轮 OOD 之后 | — |

---

# 阶段四：数据划分与 OOD 标注

> 目标：UniProt 级别无重叠的 7:1:2 划分；极端长度样本整体进入 test 作为 OOD 候选；再按统一口径补齐全部 OOD 维度。
>
> **统一口径**：各 OOD 维度的唯一实现沉淀在四个 skill 中，重算时必须调用对应 skill 脚本，不得另起口径：
> - `.skills/homology-ood-annotation`：`seq_Redundancy_*` / `OOD_Orphan`（mmseqs2 easy-search）/ `TM-score_*`（Foldseek，支持 `--m8-cache` 重放）
> - `.skills/ood-annotation-toolkit`：`idr_ratio` / `OOD_IDR`（metapredict）、值域 bin LongTail、NaN/InD 收尾、splits 体检
> - `.skills/cath-ted-domain-annotation`：`OOD_FoldHoldout` / `OOD_SuperfamilyHoldout`（CATH+TED 域标注，经 `scripts/07_ted_augment_cath_ood.py` 调用）
> - `.skills/ec-function-ood-annotation`：`OOD_NewEC_*` / `OOD_LongTail_EC_*`（all-not-in 统一口径）

## 4.1 基于长度的 OOD 训练/验证/测试划分

- **脚本**：`scripts/05_create_length_ood_split.py`
- **输入**：`data/kcat_with_reactant_ecfp4.csv`
- **输出**：`splits/kcat_{train,val,test}.csv`、`output/kcat_split_summary.json`；日志 `data/create_length_ood_split.log`
- **划分策略**：
  1. 分离极端长度样本（< 60 或 > 1000，649 行）
  2. 正常长度样本按 **UniProt** 随机 7:1:2 划分（`seed=42`），同一蛋白不出现在多个 split
  3. 正常 val/test 与 train 的反应物最大 Tanimoto（ECFP4，阈值 **0.6**）< 0.6 的 UniProt 整体移入 train（val 26 + test 43，共 **69** 个）
  4. 极端长度样本与扩充后 train 的最大相似度 < 0.6 的 UniProt 整体剔除（**5** 个 UniProt、**21** 行）；其余极端长度行全部进入 test 作为 OOD 候选

**划分结果（脚本输出时点）**：

| Split | 行数 | UniProt 数 | 长度范围 | 极端长度数 |
|-------|-----:|-----------:|----------|-----------:|
| Train | 15,479 | 4,888 | 60–1000 | 0 |
| Val | 2,150 | 662 | 60–1000 | 0 |
| Test | 5,040 | 1,496 | 混合 | 628 |
| 已剔除 | 21 | 5 | 极端 | — |

> 后续 4.7 按全项目统一规则从 test 剔除 38 条超长（>2000）样本，test 最终为 **5,002** 行 / **1,486** 个 UniProt（极端长度 590）。

## 4.2 POOD 标准列与长度 OOD 列

- **脚本**：`scripts/06_create_pood_compatible_splits.py`（首轮一次性流水线；**仅"加标准列 + 长度 OOD"部分为现行有效功能**，其 mmseqs2/Foldseek/metapredict 段已 superseded，见 4.3–4.4）
- **输入/输出**：`splits/kcat_{train,val,test}.csv`（就地改写）
- **说明**：加标准列（`unique_id`、`file_type`、`aa_seq`、`labels`——后经全项目统一 schema 迁移为现行的 `struct_file`/`label` 并删除 `file_type`，见注意事项）；计算 `Default` / `OOD_ExtremeLength` / `OOD_ExtremeLong` / `OOD_ExtremeShort`；摘要 `output/kcat_pood_summary.json`（旧口径计数，仅供参考）。
- **unique_id 生成规则**：`unique_id = AF-<acc>-F1-model_v6`（蛋白级标识，与 `struct_file` 同名）。kcat 为（酶, 反应物）对级数据，同一 UniProt 多底物多行；为保证 split 内唯一（README "unique within a task"），对有重复的 unique_id 追加反应物内容哈希后缀 `-<md5(reactant_smiles)[:8]>`（`(uniprot, reactant_smiles)` 组合已验证无重复）；未重复的保持不变。2026-09 修复前仅为蛋白级 id，split 内大量重复（train 10,591 / val 1,488 / test 3,516 行），修复在源头（`add_standard_columns`）完成，备份 `output/backup_before_unique_id_fix/`

**最终文件格式**：Train/Val 为 6 列（`unique_id, aa_seq, struct_file, ligand_smiles, ligand_ecfp4, label`）；Test 为 6 列 + 31 个 OOD/标注列（`label` 强制保存为字符串；`ligand_smiles` 为反应物 SMILES，多分子以 `.` 分隔；`ligand_ecfp4` 为激活 bit 索引列表字符串）。

## 4.3 序列同源与结构同源 OOD

- **skill 脚本**（重算唯一入口）：
  - `.skills/homology-ood-annotation/scripts/seq_homology_ood.py`（mmseqs2 easy-search，query=test、target=train+val，取最大 identity：`seq_Redundancy_<T>` = max identity < T%，T=90…30；无显著 hit 或最佳 e-value > 1e-3 → `OOD_Orphan=True`）
  - `.skills/homology-ood-annotation/scripts/struct_homology_ood.py`（Foldseek，取最大 alntmscore：`TM-score_<T>` = max TM-score < T，T=0.9…0.3，支持 `--m8-cache` 重放）
- **说明**：序列搜索按 `aa_seq` 去重、结构搜索按 UniProt 去重后进行，再映射回逐行（同一 UniProt 多底物多行共享标记）
- **[superseded]**：`scripts/06_create_pood_compatible_splits.py` 内联的 mmseqs2/Foldseek 段为首轮历史实现，重算一律走上述 skill 脚本
- **2026-08-30 统一验证**：序列侧按 skill easy-search（`-s 7.5`）重跑写回（m8 缓存 `analysis/output/seq_ood_unify_assessment/kcat/mmseqs_easysearch.m8`，备份 `output/backup_before_seq_easysearch_unify/`）；与旧值仅 34/4,412 个 Default 样本有边界变化（seq50/40/30 与 Orphan 共 83 个单元格，全部 T→F），`InD` 366 不变——旧列与 easy-search 口径实质等价。

## 4.4 IDR 无序 OOD

- **skill 脚本**：`.skills/ood-annotation-toolkit/scripts/idr_ood.py`（metapredict v3.0.2；本任务口径：`idr_ratio` = disorder score > 0.5 的残基占比（连续标注列，全 test 行有值）；`OOD_IDR` = 连续无序区域占比 > 0.3）
- **[superseded]**：`scripts/06_create_pood_compatible_splits.py` 与已归档 `reformat_splits_to_pood_format.py` 的内联 metapredict 计算

## 4.5 CATH+TED 域标注：OOD_FoldHoldout / OOD_SuperfamilyHoldout

- **脚本**：`scripts/07_ted_augment_cath_ood.py`（wrapper）
- **输入**：`splits/kcat_{train,val,test}.csv`；TED 缓存 `data/ted/ted_api_cache.json`（项目级共享）
- **输出**：`output/kcat_ted_aug_acc_list.csv`（7,036 个 accession 键清单）、`output/kcat_acc_cath_ted_labels.csv`、`output/kcat_cath_ted_holdout.csv`、`output/backup_before_ted_aug/`
- **流程**：
  0. 从 splits 生成 accession 清单（脚本内置，原为无脚本的历史一次性产物）
  1. `.skills/cath-ted-domain-annotation/scripts/cath_ted_annotate.py --chains output/kcat_ted_aug_acc_list.csv --out output/kcat_acc_cath_ted_labels.csv --offline`
  2. `.skills/cath-ted-domain-annotation/scripts/cath_ted_holdout.py --labels ... --split-col splits --out output/kcat_cath_ted_holdout.csv`
  3. wrapper 按 accession 把 holdout 结果映射回 `splits/kcat_test.csv` 逐行写回（同一 accession 多底物多行），`Default=False` 行置 NaN。accession 提取用正则 `^AF-(.+?)-F1-model_v6(?:-[0-9a-f]{8})?$`，兼容 unique_id 的反应物哈希后缀形态（见 4.2）
- **口径**（UniProt accession 模式）：标签 = CATH v4.4 实验标注（经 SIFTS UniProt→PDB 链反查并集）∪ TED putative 标注（按全长不过滤；H 级提供 C.A.T.H，T 级仅 C.A.T）；参考集 = **train** 的 CATH+TED 并集（train/val/test 对称增强）；any 语义；无注释 → False
- **覆盖率**：标注来源分布 both 2,700 / ted-only 4,103 / cath-only 42 / none 191；Topology 覆盖率仅 CATH 39.0% → CATH+TED **97.3%**；TED 缓存全部命中（离线运行，0 次新查询）
- **结果**（当前 splits，Default=True n=4,412）：`OOD_FoldHoldout` True **43**、`OOD_SuperfamilyHoldout` True **160**、无注释（→ False）**109** 行

## 4.6 EC 功能 OOD 与值域 bin LongTail

**NewEC（`OOD_NewEC_L4` / `OOD_NewEC_L3`）**

- **skill 脚本**：`.skills/ec-function-ood-annotation/scripts/add_unified_newec.py`（摘要 `output/newec_unified_summary.json`，备份 `output/backup_before_newec_unify/`）
- **口径**：参考集 = train+val 出现的 EC 集合；**all-not-in**——样本有 EC 注释且其所有 EC（对应层级）都不在 train+val → True；无注释 → False；层级解析允许部分 `-`
- **EC 源**：`data/kcat_with_reactant_ecfp4.csv`（uniprot → ec，数据集原始标注；个别 UniProt 为 `;` 分隔多 EC，需拆分）。重算必须沿用该历史源，不可用 UniProt 当前注释（存在源漂移）
- **结果**（Default n=4,412，逐行计）：`OOD_NewEC_L4` True **853**；`OOD_NewEC_L3` True **3**
- **取代关系**：旧列 `OOD_NewFunction`（EC L4 any 语义，1,224 条）已删除，由 `OOD_NewEC_L4` 取代（旧脚本 `add_additional_ood_labels.py` 已归档）

**LongTailEC（四列）**

- **skill 脚本**：`.skills/ec-function-ood-annotation/scripts/add_unified_longtail_ec.py`（摘要 `output/longtail_ec_unified_summary.json`，备份 `output/backup_before_longtail_unify/`、`output/backup_before_longtail_allnotin/`）
- **口径**：参考集 = **train-only** 逐行频次；**all-not-in**——样本有 EC 注释且其所有 EC（对应层级）在 train 中频次都 < 5 / < 10（**含 0**）→ True；无 EC 注释 → False；`Default=False` 行置 NaN
- **结果**（Default n=4,412，逐行计）：`OOD_LongTail_EC_L3_le5` **28**、`L3_le10` **75**、`L4_le5` **1,524**、`L4_le10` **2,163**
- **取代关系**：旧列 `OOD_LongTailFunction`（217 条）与过渡列 `OOD_LongTail_EC_L3`（85 条）均已删除

**值域分 bin LongTail（`OOD_LongTail_KcatBin_le50` / `le100`，回归任务专属）**

- **脚本**：`scripts/08_add_value_bin_longtail_ood.py`；逻辑已沉淀至 `.skills/ood-annotation-toolkit/scripts/value_bin_longtail_ood.py`（重算入口）
- **口径**：bin 宽 0.5，边 [-6.5, 7.0)；**clip 边界值单独成箱**（label==-6.0、==6.0 各作为独立 bin 统计 train 频次 3 / 92，避免 ==6.0 的删失堆积遮蔽右尾）；test `Default=True` 行所属 bin 的 train 计数 ≤ 50 / 100 → True；label 缺失/无法分箱 → False；`Default=False` 行置 NaN
- **结果**（Default n=4,412）：le50 True **7**、le100 True **133**（摘要 `output/kcat_longtail_kcatbin_summary.json`，备份 `output/backup_before_kcatbin_clipfix/`）
- **分布形态**：train label（n=15,479）近似正态（mean=0.88，std=1.67，skew=0.11），≤50 bin 为 (-6.5,-3.5]，≤100 加 (-4.0,-3.0] 与 (5.0,5.5]，两端对称

## 4.7 超长剔除与 NaN/InD 收尾

- **超长剔除**（无脚本环节）：按全项目统一规则（单链长度 > 2000 一律剔除），从 test 移除 **38** 行（**10** 个 UniProt）：test 由 5,040 → **5,002** 行，极端长度 628 → 590；备份 `output/len2000_removal_backup/kcat_test.before_len2000.csv`。当前 test 长度范围 16–1928
- **NaN 惯例与兜底校验**：`Default=False` 行除 `OOD_ExtremeShort/Long` 外所有 `OOD_*`、`seq_Redundancy_*`、`TM-score_*` 列为 NaN（空单元格，"未计算"语义）。2026-08-28 起该 NaN 由各生成脚本（`06` 及 skill 标注脚本）在流水线位置直接写入；`.skills/ood-annotation-toolkit/scripts/finalize_ood_columns.py --mode nan` 仅作兜底校验（当前数据上为 no-op；历史备份 `output/backup_before_nondefault_nan/`）；`idr_ratio` 为连续标注列不受影响
- **InD 重算**：`.skills/ood-annotation-toolkit/scripts/finalize_ood_columns.py --mode ind`——所有 `OOD_*` + `seq_Redundancy_*` + `TM-score_*` 均为 False（NaN 视为 False）→ `InD=True`（备份 `output/backup_before_ind/`）；当前 `InD` True **366**（占 test 7.3%）

## 4.8 splits 体检

以 `.skills/ood-annotation-toolkit/scripts/check_splits.py`（只读）检查：各 split 内 `unique_id` 唯一；train∩val、train∩test、val∩test 的 `unique_id` / UniProt accession 交集均为 **0**；train/val 为 6 核心列、test 为 37 列；Default=False 行非长度 OOD 列全部为空。

```bash
python ../../.skills/ood-annotation-toolkit/scripts/check_splits.py \
    --train splits/kcat_train.csv --val splits/kcat_val.csv --test splits/kcat_test.csv
```

## 4.9 OOD 统计摘要（当前 splits 实测）

测试集共 **5,002** 行：`Default=True` **4,412** 行，`Default=False`（极端长度）**590** 行。逐行统计；除长度 OOD 外其余维度的 True 全部落在 Default=True 子集内。

| OOD 列名 | True 数 | 说明 |
|----------|--------:|------|
| `Default` | 4,412 | 正常长度测试样本 |
| `OOD_ExtremeShort` | 59 | 长度 < 60 |
| `OOD_ExtremeLong` | 531 | 长度 > 1000（剔除 >2000 后） |
| `seq_Redundancy_90` | 3,936 | 与 train+val 最大一致性 < 90% |
| `seq_Redundancy_80` | 3,415 | < 80% |
| `seq_Redundancy_70` | 3,063 | < 70% |
| `seq_Redundancy_60` | 2,645 | < 60% |
| `seq_Redundancy_50` | 2,145 | < 50% |
| `seq_Redundancy_40` | 1,218 | < 40% |
| `seq_Redundancy_30` | 386 | < 30% |
| `OOD_Orphan` | 173 | 无显著同源 hit |
| `TM-score_0.9` | 956 | max TM-score < 0.9 |
| `TM-score_0.8` | 339 | < 0.8 |
| `TM-score_0.7` | 106 | < 0.7 |
| `TM-score_0.6` | 13 | < 0.6 |
| `TM-score_0.5` | 6 | < 0.5 |
| `TM-score_0.4` | 6 | < 0.4 |
| `TM-score_0.3` | 1 | < 0.3 |
| `idr_ratio` | — | 无序残基比例（连续值，全 test 行） |
| `OOD_IDR` | 48 | 连续无序区域 > 30% |
| `OOD_FoldHoldout` | 43 | CATH+TED Topology 留一法 |
| `OOD_SuperfamilyHoldout` | 160 | CATH+TED Superfamily 留一法 |
| `OOD_NewEC_L4` | 853 | 所有 EC L4 均不在 train+val |
| `OOD_NewEC_L3` | 3 | 所有 EC L3 均不在 train+val |
| `OOD_LongTail_EC_L3_le5` | 28 | 所有 EC L3 train 频次 < 5 |
| `OOD_LongTail_EC_L3_le10` | 75 | < 10 |
| `OOD_LongTail_EC_L4_le5` | 1,524 | 所有 EC L4 train 频次 < 5 |
| `OOD_LongTail_EC_L4_le10` | 2,163 | < 10 |
| `OOD_LongTail_KcatBin_le50` | 7 | 所属 0.5-log bin 在 train ≤ 50 样本 |
| `OOD_LongTail_KcatBin_le100` | 133 | ≤ 100 样本 |
| `InD` | 366 | 纯 ID 样本（占 test 7.3%） |

**分布一致性（当前 splits 实测）**：

| 维度 | Train | Val | Test |
|------|-------|-----|------|
| 样本数 | 15,479 | 2,150 | 5,002 |
| UniProt 数 | 4,888 | 662 | 1,486 |
| 平均长度 | 410.9 | 421.5 | 489.7 |
| 长度范围 | 61–990 | 72–965 | 16–1928 |
| 长度 < 60 | 0 | 0 | 59 |
| 长度 > 1000 | 0 | 0 | 531 |

## 4.10 工具与参数

| 工具 | 版本 | 用途 | 核心参数 |
|------|------|------|----------|
| mmseqs2 | v18.8cc5c | 序列同源性搜索 | `easy-search`，`--format-output query,target,fident,evalue` |
| Foldseek | commit 718d421 | 结构相似性搜索 | `search -a`，`--format-output query,target,alntmscore` |
| metapredict | v3.0.2 | 内在无序区域预测 | 默认 disorder 预测 |
| RDKit | — | ECFP4 指纹 | Morgan 半径=2，2048 bit |

---

# 复现路径（一次性）

```bash
cd datasets/enzyme_kinetics_prediction

# 阶段一：源数据获取与解析（1.1 下载为无脚本环节）
python scripts/01_download_alphafold_missing.py        # 1.2 补下载 AlphaFold 结构
python scripts/02_process_kcat_with_pdbs.py            # 1.3 PDB 过滤 + 序列一致性
python scripts/03_decompress_pdbs_and_truncate_q9y233.py  # 1.3 解压 + Q9Y233 截断

# 阶段二：标签聚合与计算
python scripts/04_compute_reactant_ecfp4.py            # 2.1 反应物过滤 + ECFP4

# 阶段四：划分与 OOD 标注
python scripts/05_create_length_ood_split.py           # 4.1 长度 OOD 划分
python scripts/06_create_pood_compatible_splits.py     # 4.2 标准列 + 长度 OOD（mmseqs/foldseek/IDR 段已 superseded）
# 4.3 序列/结构同源：.skills/homology-ood-annotation/scripts/seq_homology_ood.py、struct_homology_ood.py
# 4.4 IDR：.skills/ood-annotation-toolkit/scripts/idr_ood.py
python scripts/07_ted_augment_cath_ood.py              # 4.5 CATH+TED（调 skill 两脚本）
# 4.6 EC 统一口径：.skills/ec-function-ood-annotation/scripts/add_unified_newec.py、add_unified_longtail_ec.py
python scripts/08_add_value_bin_longtail_ood.py        # 4.6 值域 bin LongTail（或 toolkit value_bin_longtail_ood.py）
# 4.7 超长剔除（>2000，无脚本环节）+ finalize_ood_columns.py --mode nan / --mode ind

# 4.8 体检（只读）
python ../../.skills/ood-annotation-toolkit/scripts/check_splits.py \
    --train splits/kcat_train.csv --val splits/kcat_val.csv --test splits/kcat_test.csv
```

---

# 输出文件汇总

| 文件/目录 | 说明 |
|-----------|------|
| `splits/kcat_{train,val,test}.csv` | 最终划分（15,479 / 2,150 / 5,002 行） |
| `pdbs/` | AlphaFold v6 结构（7,058 个 `.pdb`，含 Q9Y233 全长备份） |
| `data/kcat_filtered_with_pdb.csv` / `data/kcat_with_reactant_ecfp4.csv` | 阶段一/二产物（22,736 / 22,690 行） |
| `output/kcat_split_summary.json` | 4.1 划分摘要（划分时点，test=5,040） |
| `output/kcat_pood_summary.json` | 首轮 OOD 标注摘要（旧口径计数，仅供参考） |
| `output/kcat_ecfp4_summary.json` | 2.1 ECFP4 摘要 |
| `output/kcat_longtail_kcatbin_summary.json` | 值域 bin LongTail 摘要 |
| `output/longtail_ec_unified_summary.json` / `newec_unified_summary.json` | EC 统一口径摘要 |
| `output/kcat_acc_cath_ted_labels.csv` / `kcat_ted_aug_acc_list.csv` / `kcat_cath_ted_holdout.csv` | CATH+TED 标注与 holdout 结果 |
| `output/kcat_test_ood_*.csv` | 各 OOD 场景子集（首轮生成，新旧两套命名残留，见注意事项） |
| `output/backup_*/`、`output/len2000_removal_backup/` | 各次 splits 写回前的备份 |
| `output/scripts_archive/` | 已归档的历史脚本 |

---

# 归档清单（`output/scripts_archive/`，不再使用）

| 脚本 | 归档原因 |
|------|----------|
| `add_additional_ood_labels.py` | 旧 EC 口径（`OOD_NewFunction` / `OOD_LongTailFunction`），已被 ec-function skill 的 all-not-in 统一口径（`OOD_NewEC_*` / `OOD_LongTail_EC_*`）取代，旧列已从 splits 删除 |
| `reformat_splits_to_pood_format.py` | 6 核心列重排 + `idr_ratio` + 极端行 NaN 的功能已被全项目统一 schema 迁移与 `finalize_ood_columns.py` 覆盖；其产出的旧列名（`pdb_file`/`labels`/`file_type`）已失效 |
| `analyze_catpred_data.py`、`analyze_distributions.py`、`compare_zenodo_github.py` | 探索性统计，不进主流程 |
| `split_kcat_random.py`、`create_ecfp4_splits.py`、`move_low_similarity_to_train.py`、`create_ood_test_labels.py`、`finalize_splits_and_replace_sequences.py` | 已被取代的划分/写回脚本（random splits 路线） |

# 注意事项

1. **旧 OOD 列已废弃**：`OOD_NewFunction`、`OOD_LongTailFunction`、`OOD_LongTail_EC_L3` 已被统一口径列取代并从 splits 删除；`output/kcat_pood_summary.json` 与部分 `kcat_test_ood_*.csv` 子集仍是旧口径产物，仅供参考，不要当作当前 splits 统计来源。
2. **OOD 子集 csv 新旧两套命名残留**：`output/` 下同时存在 `kcat_test_ood_ExtremeShort.csv`（新）与 `kcat_test_ood_OOD_ExtremeShort.csv`（旧）等两套命名，为历次脚本不同命名约定所致，未做清理。
3. **列名 schema 迁移**：首轮 `reformat_splits_to_pood_format.py`（已归档）产出的核心列名为 `pdb_file`/`labels`/`file_type`，后经全项目统一 schema 迁移改为现行的 `struct_file`/`label`（并删除 `file_type` 列）；以当前 splits 表头为准。
4. **无脚本环节**：>2000 剔除 38 行（4.7）与 schema 迁移（注意事项 3）均为全项目统一规则执行时的一次性操作，仅保留备份佐证；复现时按相应规则执行即可。
5. **历史统计口径变化**：2026-08-27 前的 OOD 统计（如旧文档 `seq_Redundancy_90`=4,558）把极端长度行计入分母/分子；现行 NaN 惯例下所有非长度 OOD 列仅对 Default=True 子集统计（4,412 行），全文数字以 4.9 为准。
6. **逐行 vs 逐 accession 计数**：kcat 同一 UniProt 多底物多行，EC 类 OOD 的逐行计数远大于按唯一 accession 计数；文档统一采用逐行口径（与 splits 逐行写回一致）。
7. **TED 为 putative 标注**：CATH 实验标注优先、TED 补缺（并集）；TED API 极少数逗号分隔双标签已按并集解析，无逗号伪标签残留。
8. **EC 源漂移**：EC 功能 OOD 的 EC 源必须沿用 `data/kcat_with_reactant_ecfp4.csv`（数据集原始标注），不可用 UniProt 当前注释重查。
9. **Q9Y233 截断结构**：`pdbs/AF-Q9Y233-F1-model_v6.pdb` 为截断后版本（779 残基），全长备份为 `AF-Q9Y233-F1-model_v6_full.pdb`；splits 的 `struct_file` 指向前者。
