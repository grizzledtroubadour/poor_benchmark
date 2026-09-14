# PPIS Prediction 数据集处理文档

> 记录从 DIPS-Plus v1.3.0（Zenodo）原始数据到最终蛋白质-蛋白质结合界面预测（Protein-Protein Interaction Site Prediction, PPIS）数据集的完整处理流程，可按本文档从原始数据一次性复现。
> 处理位置：`datasets/ppis_prediction/`
> 数据源：https://zenodo.org/records/8140981
> 参考论文：Morehead et al., *DIPS-Plus: The enhanced database of interacting protein structures for interface prediction*, Scientific Data 2023

**任务定义**：单链残基级二分类——给定一条蛋白质链的序列与结构，预测每个残基是否处于蛋白质-蛋白质结合界面。最终交付为 **Chain-Centric**（唯一链序列为中心）数据集：每个样本唯一对应一条链、一套聚合界面标签。早期 Pair-Centric（复合物对级别）流程仅作历史参考，压缩记录于文末"历史方案"一节。

**当前 splits 规模**（`splits/`，实测）：

| 文件 | 行数 | 说明 |
|------|------|------|
| `ppis_train.csv` | **6,043** | 训练集 |
| `ppis_val.csv` | **672** | 验证集 |
| `ppis_test.csv` | **1,980** | 测试集（1,685 Default + 295 ExtremeLength）+ OOD 列 |

核心字段：`unique_id`（序列 MD5，即 `seq_hash`）、`aa_seq`（从 PDB 提取）、`struct_file`（`{unique_id}.pdb`，位于 `pdbs/`）、`label`（逐残基界面标签，`1` 界面 / `0` 非界面 / `-1` 标注缺失，与 `aa_seq` 严格等长）、`struct_label`（按 PDB 结构残基顺序对齐的标注，见 Stage 4.7）。

---

## Stage 1: 源数据获取与解析

> 目标：获取 DIPS-Plus v1.3.0 原始数据，解压、索引并解析为元数据表与 pair 级界面标签。

### 1.1 数据下载 — `scripts/01_download_dips_plus.py`

| 项目 | 内容 |
|------|------|
| Zenodo Record ID | 8140981（DOI `10.5281/zenodo.8140981`），版本 1.3.0 |
| 文件 | `final_raw_dips.tar.gz`，15.7 GB，MD5 `a4af4e14162aa88a59c05fef5d562088` ✅ |
| 下载工具 | `aria2c`（4 连接，断点续传）；`--msa` 可选下载 MSA 归档（未下载，见注意事项） |

### 1.2 解压与结构探查 — `scripts/02_extract_and_inspect.py`

解压至 `data/raw/`（约 88 GB）：`pairs-postprocessed*.txt` 官方划分列表、`*_pdb_metadata.csv`、以及按 PDB ID 前两位分桶（`aa/`…`zz/`）的 pair 文件。Pair 文件命名 `{pdb_id}.pdb{model}_{pair_index}.{ext}`（dill 或 hdf5）。

### 1.3 建立文件索引 — `scripts/03_build_pair_index.py`

扫描全部 pair 文件：唯一复合物 **42,112** 条，dill 与 HDF5 并存，优先使用 HDF5（经 `hickle` 读取）。

输出：`output/dips_plus_pair_index.csv`（42,112 行）、`output/dips_plus_split_index.csv`、`output/extraction_summary.json`。

### 1.4 解析元数据 — `scripts/04_parse_pair_files.py`

逐 pair 读取 HDF5（hickle 序列化的 8 元素结构：`complex_name` / `df0` / `df1` / `pos_idx` / `neg_idx` / `srcs` / `id` / `sequences`），输出 `output/dips_plus_metadata.csv`（42,112 行）：`unique_id`、`pdb_id`、`file_path`、`file_format`、`receptor_seq`/`ligand_seq` 及长度、链 ID、原子/残基数、`pos_pairs`/`neg_pairs`、`has_idr_*` 等。

### 1.5 Pair 级界面标签 — `scripts/05_generate_interface_labels.py`

由 pair 文件中预采样的 `pos_idx`（6 Å 原子距离阈值的正例界面原子对）推导残基级标签：残基的任意原子出现在 `pos_idx` 中即标 1。

输出：`output/interface_labels/{unique_id}_labels.npz`（42,112 个，`receptor_labels`/`ligand_labels`）、`output/interface_label_distribution.csv`、`output/interface_label_summary.json`（总体 pair 级界面残基比例 14.8%）。

> 本脚本一度随 Pair-Centric 旧线归档，2026-08 恢复为正式步骤（pair 级标签是 Chain-Centric 聚合的上游环节，也是理解 `pos_idx` 标签语义的参照实现）。

---

## Stage 2: 标签聚合与计算

> 目标：以唯一链序列为中心，聚合该序列在所有 pair 中的界面残基（并集），生成残基级界面标签。

### 2.1 Chain-Centric 标签聚合 — `scripts/06_build_chain_centric_dataset.py`

1. 以完整氨基酸序列作为链的唯一标识（MD5，记 `seq_hash`）。
2. 遍历全部 42,112 条 pair，提取受体和配体的界面标签（界面定义同 §1.5 的 `pos_idx` 口径）。
3. 对每条唯一序列，将其在所有 pair 中出现时的界面残基位置取**并集（union）**。

| 指标 | 数值 |
|------|------|
| 原始 pair 数 | 42,112 |
| 唯一链序列数 | **22,223** |
| 出现多次的唯一链 | 12,422（55.9%） |
| 平均界面残基比例 | 28.3%（并集覆盖多个结合区域，高于 pair 级 14.8%） |

输出：`output/chain_centric_metadata.csv`（22,223 行）、`output/chain_interface_labels/{seq_hash}_labels.npz`、`output/chain_centric_summary.json`。

---

## Stage 3: 去冗余与数据精化

> 目标：95% identity 去冗余得到代表链，并为每条代表链提取坐标最完整、链归属正确的 PDB 结构。

### 3.1 序列去冗余 — `scripts/07_dedup_chain_centric.py`

mmseqs2 `cluster`（95% identity / 80% coverage / cov-mode 0）对 22,223 条唯一链聚类，每簇保留 mmseqs2 自动选择的代表序列。

| 指标 | 数值 |
|------|------|
| 代表链数 | **8,850** |
| 压缩比例 | 60.2% |
| 平均/最大聚类大小 | 2.51 / 64 |

输出：`output/chain_centric_dedup_metadata.csv`、`output/chain_centric_dedup_summary.json`。

> 历史版本（`dedup_chain_centric_and_link_structures.py`）曾从 Pair-Centric 旧线产物为代表链建结构符号链接（`output/chain_representative_structures/`），该死代码已裁掉——结构来源字段由 §3.2 直接回填。

### 3.2 代表结构提取 — `scripts/08_extract_representative_structures.py`

1. 扫描全部 42,112 条 pair 元数据，为每条唯一序列选**原子数最多**的最佳出现。
2. 加载该 pair 的 `df0`/`df1`，**按序列匹配**选择与目标序列对应的 df（完全匹配 → 唯一长度匹配 → 相似度最高）；DIPS-Plus pair 文件中 `df0`/`df1` 与元数据 receptor/ligand 可能互换，不能按 role 直取。
3. 提取坐标写为单链 PDB（统一链 ID `A`），并回填 `chain_centric_dedup_metadata.csv` 的 `structure_source_unique_id` / `structure_source_role`（`{role}->{matched_df}`）。

> 序列匹配选 df 的逻辑原属一次性补丁 `fix_representative_structures.py`（初版按 role 直取导致约 6,786 条链提取错误），已融入本脚本；其 PDB 列格式 bug（resname 后多空格致 chain/resseq 右移）不带入，本脚本直接产出标准列布局（坐标位于 cols 31-38，见脚本内注释）。

| 指标 | 数值 |
|------|------|
| 有结构的代表链 | 8,850 / 8,850（100%） |
| 修正后 PDB 与 metadata 序列一致 | 8,741 |
| 长度一致但残基类型差异（多为修饰残基） | 63 |
| 序列长度仍不一致 | 46（经 §4.3 对齐映射或 §4.4 清理处理） |

输出：`pdbs/{seq_hash}.pdb`（8,850 个）。

---

## Stage 4: 数据划分与 OOD 标注

> 目标：72:8:20 划分 + 多维 OOD 标记 + 统一口径收尾，生成 `splits/` 三个文件。

### 4.1 数据划分 — `scripts/09_split_chain_centric.py`

1. 链长不在 **[60, 1000]** 的链独立为 `test_extremelength`。
2. 其余按 **72 : 8 : 20** 随机划分 Train / Val / Test（seed=42）。代表链已 95% 去冗余，随机划分无序列泄露。

划分时点（8,850 条）：Train 6,153 / Val 683 / Test 1,710 / ExtremeLength 304。

输出：`output/chain_centric_splits/chain_centric_{train,val,test,test_extremelength}.csv`、`output/chain_centric_split_summary.json`。

### 4.2 序列 / 结构 / IDR OOD 初版标记 — `scripts/10_mark_ood_chain_centric.py`

参考集 = Train + Val，输出 `output/chain_centric_splits/chain_centric_test_ood.csv`、`output/ood_chain_centric_summary.json`。

- **序列同源**（mmseqs2 search，`OOD_LowHomology_30..90`、`OOD_Orphan`）与**结构相似**（foldseek，`OOD_NewFold_0.3..0.9`）两段：**superseded**——当前口径以 `.skills/homology-ood-annotation/`（`seq_homology_ood.py` / `struct_homology_ood.py`，easy-search）为准，且**两段现值均已按 skill 口径重算**：`TM-score_*` 于 2026-08-28（§4.4 裁剪后结构，见注意事项 4），`seq_Redundancy_*`/`OOD_Orphan` 于 2026-08-29（旧 `mmseqs search` 默认灵敏度系统性漏检 hit、高估 OOD，如 Orphan 522 → 458、seq_Redundancy_90 1,570 → 1,457，InD 64 → 134；备份 `output/backup_before_seq_easysearch_unify/`，重跑产物 `analysis/output/ppis_seq_easysearch_audit/`，m8 缓存 `output/ood_mmseqs_test_vs_trainval.m8`）。
- **OOD_IDR**：保留有效——使用 DIPS-Plus 自带 **flDPnn** 无序倾向（界面残基中 propensity > 0.5 的比例 > 0.3），与 ood-annotation-toolkit 的 metapredict 口径不同，不可混用。

### 4.3 最终三文件生成与标签对齐 — `scripts/11_build_final_ppis_splits.py`

1. 合并四个中间划分文件；对每个 `seq_hash` 从 `pdbs/{seq_hash}.pdb` 提取 `aa_seq`。
2. 加载 `output/chain_interface_labels/{seq_hash}_labels.npz` 标签；PDB 序列与标签等长则直接使用，否则用 **Biopython 全局对齐**映射（identity ≥ 0.85，或 ≥ 0.75 且 coverage ≥ 0.75；PDB 插入残基标 `-1`）；对齐失败跳过。
3. 测试集合并正常长度与极端长度样本，生成 `Default` / `OOD_ExtremeLong` / `OOD_ExtremeShort`，常规 OOD 列仅对 `Default=True` 判定。
4. 本脚本另内置 CATH 名称匹配 Holdout、SIFTS EC NewFunction、界面稀疏 LongTail（有效界面标注中正样本比例 < 5% → `OOD_LongTail`）的初版计算；CATH 与 EC 两轴后被 skill 统一口径重算（§4.5、§4.6），界面稀疏 `OOD_LongTail` 仍为当前值。

对齐统计：合计 8,849 条中 8,807 直接匹配、42 条对齐映射、跳过 1 条（`8b00203a…`，PDB 279 vs metadata 488 无法对齐）。

### 4.4 结构裁剪清理 — `scripts/12_crop_ppis_structures.py`

对每条样本检查 `pdbs/{unique_id}.pdb` 是否含与 `aa_seq` 完全一致的链；否则选最佳匹配链全局对齐并**裁剪**为一一对应（跳过结构插入），保证 `len(结构残基) == len(aa_seq) == len(label)`；对齐 identity < 0.95 或结构有缺失的样本判不可恢复，整行从 splits 删除。

- 删除不可恢复样本 **154**（Train −109 / Val −11 / Test −34），裁剪结构 **61**，清理后总规模 **8,695**。
- 隔离产物：`output/quarantine_structures/`（裁剪前原始 PDB）、`output/quarantine_splits/`（删除前 splits 备份）。
- 氨基酸三字母→单字母映射使用 Biopython `protein_letters_3to1_extended`（脚本内联，原依赖 `analysis/scripts/verify_aa_seq_vs_structure.py` 已移除）。

### 4.5 CATH+TED 域标注与 Holdout — `scripts/13_ted_augment_cath_ood.py`（任务 wrapper）

方法学：`.skills/cath-ted-domain-annotation/`（标注 `cath_ted_annotate.py` + Holdout `cath_ted_holdout.py`，判定口径以 skill 脚本为唯一实现）。

1. 链清单：seq_hash → （来源 pair, role) → (pdb, chain)（`src_chain` 映射，复用已归档 `add_cath_holdout_ood.py` 的逻辑）→ `output/ppis_ted_aug_chain_list.csv`（8,695 链）。
2. `cath_ted_annotate.py --offline`：CATH v4.4 名称匹配 ∪ TED putative 标注（REST API 缓存 `data/ted/ted_api_cache.json`，域重叠 ≥50% 域长）→ `output/ppis_chain_cath_ted_labels.csv`、`output/ppis_chain_domains.csv`。
3. `cath_ted_holdout.py`：参考集 = Train 的 CATH+TED 并集；any 语义（链任一 Topology/Superfamily 未在参考集出现 → True；无注释 → False）→ 就地更新 `splits/ppis_test.csv` 的 `OOD_FoldHoldout` / `OOD_SuperfamilyHoldout`（Default=False 行直接写 NaN，见 §4.7）。

当前值（test Default n=1,685）：`OOD_FoldHoldout` **65**、`OOD_SuperfamilyHoldout` **235**；标注来源分布 both 7,125 / cath-only 985 / ted-only 323 / none 262。仅 CATH 旧口径备份在 `output/backup_before_ted_aug/`。

### 4.6 EC 功能 OOD（统一口径）— `.skills/ec-function-ood-annotation/`

- **NewEC**（`add_unified_newec.py`）：参考集 = train+val EC 集合，all-not-in 语义；`OOD_NewEC_L4`（**144**）、`OOD_NewEC_L3`（**9**）。
- **LongTailEC**（`add_unified_longtail_ec.py`）：参考集 = train-only 频次，all-not-in（所有 EC train 频次 < 5 / < 10，含 0）→ 四列：`OOD_LongTail_EC_L3_le5` **62** / `L3_le10` **137** / `L4_le5` **515** / `L4_le10` **659**（Default n=1,685）。
- EC 数据源：`data/sifts/sifts_chain_ec.tsv.gz`（SIFTS 链级），映射链同 §4.5 的 `src_chain`；test Default 中有 EC 注释 873 条（51.8%）。
- 与界面稀疏 `OOD_LongTail`（§4.3）语义不同，二者并存。

### 4.7 收尾与体检（skill 链）

按序执行（脚本均在 `.skills/` 下）：

1. **homology-ood-annotation**（`seq_homology_ood.py` / `struct_homology_ood.py`）：seq_Redundancy/Orphan/TM-score 的当前口径实现，**两类列现值均为本步口径**——`TM-score_*` 2026-08-28 重算（query=test 1,980 结构 / target=train+val 6,715 结构，m8 缓存 `output/foldseek_easysearch_test_vs_trainval.m8`）；`seq_Redundancy_*`/`OOD_Orphan` 2026-08-29 重算（easy-search `-s 7.5`，test 1,980 / ref=train+val 去重序列，m8 缓存 `output/ood_mmseqs_test_vs_trainval.m8`）。
2. **ood-annotation-toolkit `idr_ood.py`**：**不适用**于本任务——`OOD_IDR` 用 DIPS-Plus flDPnn 口径（§4.2），不替换。
3. **cath-ted-domain-annotation**：即 §4.5。
4. **ec-function-ood-annotation**：即 §4.6（回归任务的 `value_bin_longtail_ood.py` 不适用，本任务为残基级分类）。
5. **finalize**（`finalize_ood_columns.py`）：
   - `--mode nan`：兜底校验——`Default=False`（295 行极端长度）行的所有 OOD 标记列（`OOD_*` 除 `OOD_ExtremeShort/Long`，含 `seq_Redundancy_*` / `TM-score_*`）应为 NaN。2026-08-28 起 NaN 由各生成脚本（`11`/`13` 及 skill 脚本）在流水线位置直接写入，本步在当前数据上为 no-op；同日已将历史遗留的 seq/TM 列 False 值统一为 NaN（备份 `output/backup_nondefault_seqtm_nan/`）。
   - `--mode ind`：`InD` = 全部 OOD/seq/TM 列均为 False（NaN 视为 False）→ 当前 **134 / 1,980（6.8%）**。
6. **struct_label**（`struct-label-alignment/add_struct_label.py`）：对 `struct_file` 每个结构残基（ATOM 按文件顺序、`(chain, resseq, icode)` 首次出现）给出与 PDB 残基序列全局比对对齐的标注，无法映射记 `-1`；三个 splits 文件均已追加该列（备份 `output/backup_before_struct_label/`）。
   - **resseq 锚定复核与修复（2026-09-02，`scripts/14_fix_struct_label_resseq.py`）**：13 条"结构残基数 ≠ aa_seq 长度"的行（结构缺失残基，走比对路径）经 resseq 编号锚定复核——其中 11 行确认现有 struct_label 与编号映射完全一致（比对结果本就正确）；2 行（train 的 `fc54381e1a5a…`、`16d04ef06e42…`）末端 C 端 tag 残基超出 aa_seq 覆盖范围，原比对误将其标为 0 / 1（界面），已按 resseq 映射修正为 `-1`（备份 `output/backup_20260902_structlabel_fix_ppis_train.csv`）。修复后本任务含 `-1` 的为这 2 行（各 1 个残基）。
   - 背景教训：比对式继承在"结构缺失残基落在重复/低复杂度区"时存在缺口放置歧义（SSP 148 行错帧事件，见 `ss_prediction/DATA_PROCESS.md` §5.5 与 skill 注意事项）；LBS 经同类全量复核无此问题（96,077 行直通 + 179 行等长点替换）。
7. **check_splits**（`check_splits.py`）：划分体检——Train∩Val / Train∩Test / Val∩Test 的 `unique_id` 交集均为 **0**（实测）；`len(label)==len(aa_seq)` 全量成立；`struct_file` 在 `pdbs/` 零缺失。

此外列名经项目统一 schema 迁移规范为现名（`file_type→struct_file`、`labels→label`、`OOD_LowHomology_*→seq_Redundancy_*`、`OOD_NewFold_*→TM-score_*`、`OOD_CATH_FoldHoldout→OOD_FoldHoldout`、`OOD_NewFunction→OOD_NewEC_L4`、`ExtremeLong/Short→OOD_ExtremeLong/Short`；无独立脚本，全项目统一步骤）。

### 4.8 最终划分统计（实测）

| 集合 | 样本数 | 占比（相对 normal 8,400） | 平均长度 | 平均界面比例 |
|------|--------|----------------------|---------|-------------|
| Train | 6,043 | 71.9% | 252.8 | 26.4% |
| Val | 672 | 8.0% | 262.3 | 26.3% |
| Test（Default） | 1,685 | 20.1% | 256.9 | 26.4% |
| Test_ExtremeLength | 295 | — | 173.1 | 35.9%（34 条极长 1,012–1,484 + 261 条极短 7–59） |

测试集 OOD True 数（1,980 行；除 Default/Extreme/InD 外均为 Default 子集内判定）：

| OOD 列 | True | | OOD 列 | True |
|--------|-----:|---|--------|-----:|
| `Default` | 1,685 | | `TM-score_0.9` / `0.5` / `0.3` | 770 / 47 / 20 |
| `OOD_ExtremeLong` / `Short` | 34 / 261 | | `OOD_IDR` | 28 |
| `seq_Redundancy_90` / `70` / `50` / `30` | 1,457 / 1,298 / 1,066 / 602 | | `OOD_FoldHoldout` / `SuperfamilyHoldout` | 65 / 235 |
| `OOD_Orphan` | 458 | | `OOD_NewEC_L4` / `L3` | 144 / 9 |
| `OOD_LongTail`（界面稀疏 <5%） | 48 | | LongTailEC L3_le5/le10、L4_le5/le10 | 62 / 137 / 515 / 659 |
| `InD` | 134 | | | |

---

## 输出文件汇总

```
datasets/ppis_prediction/
├── DATA_PROCESS.md                        # 本文档
├── splits/                                # 最终划分（仅三个文件）
│   ├── ppis_train.csv / ppis_val.csv / ppis_test.csv
├── pdbs/                                  # 代表链 PDB（{seq_hash}.pdb，8,833 个，见注意事项 7）
├── data/raw/                              # DIPS-Plus pair 文件（约 88 GB）
├── logs/                                  # 运行日志
├── scripts/                               # 现役脚本（01–13，见各阶段）
└── output/
    ├── dips_plus_pair_index.csv / dips_plus_metadata.csv
    ├── interface_labels/                  # pair 级标签（42,112 NPZ，§1.5）
    ├── chain_centric_metadata.csv / chain_interface_labels/  # 唯一链聚合标签（22,223 NPZ）
    ├── chain_centric_dedup_metadata.csv / chain_centric_splits/
    ├── ppis_chain_cath_ted_labels.csv / ppis_chain_domains.csv / ppis_ted_aug_chain_list.csv
    ├── newec_unified_summary.json / longtail_ec_unified_summary.json
    ├── backup_before_*/                   # 各次写回前的 splits 备份
    ├── quarantine_splits/ / quarantine_structures/  # §4.4 隔离产物
    └── scripts_archive/                   # 归档脚本（见归档清单）
```

## 数据质量检查

1. **泄露**：Train/Val/Test `unique_id` 两两交集均为 0（实测）；代表链已 95% 去冗余。
2. **一致性**：全部 8,695 样本满足 `len(label)==len(aa_seq)` 且结构残基与 `aa_seq` 一一对应；`-1` 掩码残基共 140 个（均在测试集）；`struct_file` 零缺失。
3. **分布**：Train/Val/Test(Default) 平均界面比例均约 26.4%，长度分布接近。

## 归档清单（`output/scripts_archive/`，含原因）

| 脚本 | 归档原因 |
|------|---------|
| `fix_representative_structures.py` | 序列匹配选 df 逻辑已融入 `08_extract_representative_structures.py`（其 resname 列格式 bug 不带入） |
| `18_fix_pdb_column_layout.py` | 一次性 PDB 列布局/拼接链修复，修复已落数据（备份 `output/backup_before_pdb_layout_fix/`），根因已修入 08 |
| `add_function_longtail_ood.py` | EC/界面 LongTail 初版补充；EC 列被 §4.6 统一口径取代、界面 LongTail 已在 11 内置 |
| `evaluate_cath_coverage.py` | 纯诊断打印（CATH 覆盖率），不进主流程 |
| `18_add_longtail_ec_ood.py` | 初版单列 LongTailEC（any 语义），已被统一 all-not-in 四列口径取代 |
| `add_cath_holdout_ood.py` | 仅 CATH 名称匹配 Holdout，被 §4.5 CATH+TED 统一标注取代 |
| `generate_interface_labels.py` | 2026-08 已取消归档，恢复为 `scripts/05_generate_interface_labels.py`（归档副本与现役一致，已去重） |
| `split_dataset.py`、`filter_by_pairwise_identity*.py`、`filter_by_sequence_identity.py`、`split_82_and_mark_lowhomology.py`、`mark_ood.py`、`mark_ood_82split.py`、`extract_chain_structures.py`、`fix_mismatched_structures.py` | Pair-Centric 旧线脚本（见"历史方案"），被 Chain-Centric 主线取代 |

## 注意事项

1. **磁盘空间**：核心归档 15.7 GB，解压后 `data/raw/` 约 **88 GB**；mmseqs2/foldseek 中间文件另需空间。
2. **MSA 归档未下载**：`interim_external_feats_dips_msas.tar.gz`（12.0 GB）因空间限制未下载。
3. **IDR 口径**：`OOD_IDR` 用 DIPS-Plus 自带 flDPnn（界面残基 propensity>0.5 比例 >0.3），与 toolkit metapredict 口径不同，重算不可混用。
4. **seq/TM 同源列的口径状态**：两段均已统一为 skill easy-search 口径。`TM-score_*` 于 2026-08-28 在 §4.4 裁剪后结构上重算（旧历史实现 `search -a 1` 系统性高估 OOD，如 0.9 档 1,600 → 770、0.3 档 1,365 → 20；备份 `output/backup_before_foldseek_easysearch_unify/`，重跑产物 `analysis/output/foldseek_easysearch_audit/ppis/`）。`seq_Redundancy_*`/`OOD_Orphan` 于 2026-08-29 重算（旧 `mmseqs search` 默认灵敏度漏检 hit、高估 OOD：seq_Redundancy_90 1,570 → 1,457、seq_Redundancy_30 686 → 602、Orphan 522 → 458，InD 64 → 134；备份 `output/backup_before_seq_easysearch_unify/`，重跑产物 `analysis/output/ppis_seq_easysearch_audit/`，m8 缓存 `output/ood_mmseqs_test_vs_trainval.m8`）。
5. **PDB 坐标列格式（本次整理发现并已修）**：`pdbs/` 数据本身为标准布局（坐标 cols 31-38，可正常解析）；但整理前 08 的模板在 resseq 后只留 3 空格，重跑产物坐标会整体左移一列（负坐标百位可丢负号）。已在 08 中修正为 4 列间隔，修正后 60 链抽样重提取与现有 pdbs **逐字节一致**。
6. **隔离产物**：`output/quarantine_structures/`（62 个裁剪前原始 PDB）与 `output/quarantine_splits/` 仅供溯源。
7. **`pdbs/` 文件数**：8,833 个，比 splits 引用的 8,695 多 138 个，为 §4.4 删除行的残留结构（按项目约束未删除）；以 splits 的 `struct_file` 为准。
8. **Test 集来源**：DIPS-Plus 无标准 Test 集，本项目对代表链随机 72:8:20 划分；如需与文献对比可改用官方 val 作 Test。
9. **无脚本环节**：列名 schema 迁移（§4.7）为全项目统一步骤，本目录无独立脚本。
10. **历史演进**：(a) Pair-Centric → Chain-Centric（见下节）；(b) CATH-only → CATH+TED 并集（§4.5，TED 为 putative 标注，"CATH 实验标注优先，TED 补缺"）；(c) EC 初版列 → 统一 all-not-in 口径（§4.6）。各次写回前备份在 `output/backup_before_*/`。
11. **数据引用**：使用 DIPS-Plus 请引用 `Morehead et al., Scientific Data 2023`。

---

## 历史方案：Pair-Centric 旧线（已归档）

> 早期按复合物对（pair）粒度构建，仅作对比参考；脚本已归档至 `output/scripts_archive/`，中间产物保留在 `output/`（`interface_labels/`、`structures/`、`dips_plus_metadata_filtered_global_95id.csv`、`ood_82split_summary.json` 等）。

1. **界面标签**：`generate_interface_labels.py`（已恢复为现役 05）由 `pos_idx` 推导 pair 级受体/配体残基标签。
2. **去冗余**：`filter_by_pairwise_identity_global.py` 要求受体-受体与配体-配体同时 ≥95% identity / 80% coverage，mmseqs2 自搜索 + 贪心最大独立集，42,112 → 10,775 条。
3. **划分**：`split_82_and_mark_lowhomology.py` 按 PDB ID 8:2 切 Test，池内 9:1 切 Train/Val（Train 7,821 / Val 870 / Test 2,084）。
4. **OOD**：`mark_ood_82split.py` 双链严格条件计算 LowHomology/NewFold/Orphan/ExtremeLength/IDR。
5. **结构**：`extract_chain_structures.py` 提取 21,550 个链级 PDB（`output/structures/`）。

**为何改为 Chain-Centric**：最终任务是单链残基级界面预测；pair 粒度下同一链序列在不同复合物中携带不同界面标注（同序列多标签冲突），且与单链模型输入不匹配。以唯一序列为中心聚合标签后，每个样本唯一对应一条链、一套标签。
