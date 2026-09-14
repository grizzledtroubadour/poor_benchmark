# 蛋白质-配体结合亲和力（Ligand Binding Affinity）数据集处理文档

> 任务名称：ligand_binding_affinity（回归，预测 pK = -log10(Kd/Ki)）
> 数据源：PDBbind+ v2020.R1 免费数据包（https://www.pdbbind-plus.org.cn/download）
> 处理位置：`datasets/ligand_binding_affinity/`
>
> 本文档描述**从原始数据一次性复现**当前 `splits/` 与 `pdbs/` 的完整线性流程（四阶段结构）。
> 统计数字与列名均已按当前 `splits/*.csv` 实际表头重算核对（只读）。

## 任务定义与当前规模

回归任务：由蛋白序列/结构 + 配体 SMILES/ECFP4 预测结合亲和力 pK（`label`，字符串保存，统一 6 位小数）。

| 划分 | 文件 | 样本数 | Kd / Ki | 说明 |
|:---|:---|---:|:---|:---|
| 训练集 | `splits/ligand_binding_affinity_train.csv` | 8,366 | 4,638 / 3,728 | 6 列 |
| 验证集 | `splits/ligand_binding_affinity_val.csv` | 930 | 516 / 414 | 6 列 |
| 测试集 | `splits/ligand_binding_affinity_test.csv` | 2,367 | 1,684 / 683 | 38 列（含完整 OOD 块） |
| **合计** | — | **11,663** | — | `pdbs/` 下 11,663 个 `{pdb_id}.pdb` |

整体比例约 **7.2 : 0.8 : 2.0**，train:val ≈ **9 : 1**；test `Default` 全为 True（LBA/PPI 惯例，见 §9 注意事项）。

数据包说明：**PDBbind v2020.R1**（PDBbind+ 免费 Demo 数据包），结构文件使用 PDBbind v2024 流程重新处理，质量优于旧版 PDBbind-CN v2020。索引文件 `INDEX_general_PL.2020R1.lst`，复合物 19,037（蛋白-配体），无预计算 pK（解析时自行计算）。

脚本总览（`scripts/`，按执行顺序编号）：

| 编号 | 脚本 | 阶段 |
|:---|:---|:---|
| 01 | `01_download_pdbbind2020.sh` | ① 数据包下载与解压 |
| 02 | `02_parse_index.py` | ① 索引解析 |
| 03 | `03_clean_metadata.py` | ② 标签清洗（Kd/Ki 筛选） |
| 04 | `04_parse_ligands.py` | ② 配体解析（SMILES + ECFP4） |
| 05 | `05_extract_protein_sequences.py` | ③ 按链蛋白序列提取（`\|` 分隔） |
| 06 | `06_build_dataset_csv.py` | ③ 样本级数据集落盘 |
| 07 | `07_copy_whole_protein_pdbs.py` | ③ 整蛋白结构落盘 |
| 08 | `08_cluster_time_split.py` | ④ 配体聚类 + 簇内时间切分 + 比例平衡 |
| 09 | `09_format_splits_to_skill.py` | ④ 统一 schema 格式化 + 序列同源 OOD |
| 10 | `10_compute_structure_idr_ood.py` | ④ 结构相似 + IDR OOD |
| 11 | `11_ted_augment_cath_ood.py` | ④ CATH+TED Fold/Superfamily Holdout |
| 12 | `12_add_value_bin_longtail_ood.py` | ④ 值域 AffBin LongTail OOD |

---

## 阶段一：源数据获取与解析

### 1.1 数据包下载与解压

**脚本**：`scripts/01_download_pdbbind2020.sh`
**输入**：`data/index.tar.gz`、`data/P-L.tar.gz`（手动放入）
**输出**：`data/index/INDEX_general_PL.2020R1.lst`、`data/P-L/{1981-2000,2001-2010,2011-2019}/{pdb_id}/`

```bash
cd ligand_binding_affinity
bash scripts/01_download_pdbbind2020.sh
```

解压后每个样本目录含 `{pdb_id}_protein.pdb`、`{pdb_id}_ligand.mol2`、`{pdb_id}_ligand.sdf`、`{pdb_id}_pocket.pdb`。

### 1.2 索引解析

**脚本**：`scripts/02_parse_index.py`
**输入**：`data/index/INDEX_general_PL.2020R1.lst`
**输出**：`output/pdbbind_v2020_metadata.csv`（19,037 条，解析失败 0；亲和力类型 Kd / Ki / IC50，pK = -log10(affinity_value) 在此计算）

---

## 阶段二：标签聚合与计算

### 2.1 标签清洗（Kd/Ki 筛选）

**脚本**：`scripts/03_clean_metadata.py`
**输入**：`output/pdbbind_v2020_metadata.csv`
**输出**：`output/pdbbind_v2020_metadata_cleaned.csv`（+ `_stats.json`）

排除规则：

1. `affinity_type == IC50`
2. `affinity_operator` 为 `<`、`<=`、`~` 等不等号
3. 注释标记 `[Incomplete ligand]`
4. 注释标记 `[Covalent complex]`
5. 注释标记 `[Different protein in assay]`
6. 注释标记 `[Different ligand in assay]`

> 关于 `[Uncommon element]`：指配体中含有不常见元素（如硼、金属等）。本次未排除，共保留 96 条，后续若模型对元素类型敏感可再处理。

| 项目 | 数值 |
|------|------|
| 清洗前 | 19,037 |
| 清洗后 | **11,664**（Kd 6,839 / Ki 4,825） |
| 保留率 | 61.3% |

### 2.2 配体解析（SMILES + ECFP4）

**脚本**：`scripts/04_parse_ligands.py`
**输入**：`output/pdbbind_v2020_metadata_cleaned.csv`
**输出**：`output/pdbbind_v2020_metadata_ligand.csv`（+ `_stats.json`、`ligand_parse_failures.csv`）

- 优先读取 `{pdb_id}_ligand.sdf`，失败则回退到 `{pdb_id}_ligand.mol2`；
- RDKit 生成：`ligand_smiles`（isomeric canonical SMILES）、`ligand_ecfp4`（ECFP4/Morgan 指纹，radius=2，nBits=2048，保存为 ON-bits 列表）、`ligand_num_heavy_atoms`。

| 项目 | 数值 |
|------|------|
| 输入 | 11,664 |
| 解析成功 | **11,663**（`2pll` 失败，剔除） |
| 解析失败 | 1 |

---

## 阶段三：去冗余与数据精化

> 本任务评估过序列/配体去冗余路线（95% 序列同一性、配体 95%、SMILES 95% 三版），
> 最终**未采用**，样本量保持清洗后的 11,663；相应脚本已归档（见 §8），
> 中间产物 `output/pdbbind_v2020_metadata_dedup95*.csv` 仅作留痕。
> 数据精化指序列-结构对齐与整蛋白结构落盘，如下三步。

### 3.1 按链蛋白序列提取

**脚本**：`scripts/05_extract_protein_sequences.py`
**输入**：`output/pdbbind_v2020_metadata_ligand.csv` + `data/P-L/**/{pdb_id}_protein.pdb`
**输出**：`output/pdbbind_v2020_metadata_seq.csv`（+ `_stats.json`、`sequence_extraction_failures.csv`）

- 从原始蛋白 PDB 中**按链**提取标准氨基酸与常见修饰氨基酸序列；
- 多链蛋白的各链序列用 `|` 分隔（单链不加分隔符），同步产出 `protein_length`（各链求和，不含 `|`）、`protein_chains`（`;` 分隔）、`num_protein_chains`；
- 按 (chain, res_seq, res_name) 去重且只保留蛋白残基，因此从原始文件提取与从 07 过滤后的 `pdbs/` 提取结果一致（历史注记：该按链口径原由已归档的 `17_reextract_sequences_from_pdbs.py` 事后重提取，现已并入本步一次做对）；
- 提取成功 **11,663** / 失败 0；多链样本 **6,393**，单链 5,270；平均长度 532.8 aa，范围 24–5,280 aa。

### 3.2 样本级数据集落盘

**脚本**：`scripts/06_build_dataset_csv.py`
**输入**：`output/pdbbind_v2020_metadata_seq.csv`
**输出**：`output/pdbbind_v2020_dataset.csv`、`output/protein_chain_extraction_stats.json`

- 保留序列提取成功的 11,663 条样本，补充 `protein_pdb_file`（`pdbs/{pdb_id}.pdb`）列；
- 替代已归档 `08_extract_protein_chains.py` 的链级结构路线（不再逐链写 PDB）；`chain_extraction_success` / `chain_extraction_error` 两列为兼容历史下游保留，取值与序列提取成功标记一致。

`pdbbind_v2020_dataset.csv` 核心字段：

| 字段 | 说明 |
|------|------|
| `pdb_id` | PDB 四位代码（即后续 `unique_id`） |
| `resolution` / `release_year` | 分辨率（Å）/ 释放年份（时间划分依据） |
| `affinity_type` / `affinity_value` / `pK` | Kd 或 Ki / 亲和力数值（M）/ 回归目标 |
| `ligand_smiles` / `ligand_ecfp4` / `ligand_num_heavy_atoms` | 配体特征 |
| `protein_sequence` / `protein_length` | 蛋白序列（多链 `\|` 分隔）/ 氨基酸总数 |
| `protein_chains` / `num_protein_chains` | 链 ID 列表 / 链数 |
| `protein_pdb_file` | 整蛋白 PDB 路径：`pdbs/{pdb_id}.pdb` |
| `protein_file` / `ligand_sdf_file` / `pocket_file` | 原始文件相对路径 |

### 3.3 整蛋白结构落盘

**脚本**：`scripts/07_copy_whole_protein_pdbs.py`
**输入**：`output/pdbbind_v2020_dataset.csv` + `data/P-L/**/{pdb_id}_protein.pdb`
**输出**：`pdbs/{pdb_id}.pdb`（11,663 个）、`output/whole_protein_pdb_stats.json`

- 文件名**仅使用 PDB ID**；过滤非蛋白残基的 `HETATM`（水、配体、辅因子等），仅保留标准与常见修饰氨基酸；
- 遇到第一个 `ENDMDL` 停止，避免 NMR 多模型重复。

---

## 阶段四：数据划分与 OOD 标注

### 4.1 配体聚类 + 簇内时间切分 + 比例平衡

**脚本**：`scripts/08_cluster_time_split.py`
**输入**：`output/pdbbind_v2020_dataset.csv`
**输出**：`splits/{train,val,test}.csv`（中间产物，09 消费后删除）、`output/cluster_split_stats.json`

1. **团划分聚类**：基于 ECFP4 指纹贪心团划分，保证同一簇内任意两样本 Tanimoto ≥ 0.6；孤立样本成单点簇（5,394 簇，单点 3,205，最大簇 271，平均 2.16）；
2. **簇内时间切分**：每簇按 `release_year` 升序（同年按 `pdb_id`）按 7:1:2 切 train/val/test；
3. **相似度兜底**：val/test 中与 train 最大 Tanimoto < 0.6 的样本迭代移入 train（共移动 **2,866** 个）；
4. **比例平衡**（已并入原 `12_balance_val_from_test.py` / `13_adjust_val_to_train_9to1.py`）：兜底后 test 偏高、val 偏低，先将 test 按时间最早的 1/3（1,183 个）移入 val，再将 val 按时间最早的 442 个移入 train，使 train:val ≈ 9:1（实际 8.996）。

最终划分：train 8,366（71.7%）/ val 930（8.0%）/ test 2,367（20.3%）。

### 4.2 统一 schema 格式化与序列同源性 OOD

**脚本**：`scripts/09_format_splits_to_skill.py`
**输入**：`output/pdbbind_v2020_dataset.csv` + `splits/{train,val,test}.csv`
**输出**：`splits/ligand_binding_affinity_{train,val,test}.csv`（删除中间产物）、`output/{test,train_val}_for_mmseqs.fasta`、`output/mmseqs_test_vs_trainval.tsv`

- 输出统一 schema：`unique_id, aa_seq, struct_file, ligand_smiles, ligand_ecfp4, label`；`label`（pK）强制字符串、统一 6 位小数；`struct_file` 为文件名（`{pdb_id}.pdb`）；
- test 追加 `Default`（本任务全 True）、`OOD_ExtremeShort` / `OOD_ExtremeLong`（蛋白总长 < 60 / > 1,000，仅标记不移除）、`seq_Redundancy_90..30`、`OOD_Orphan`；
- 序列同源性用 mmseqs2 计算 test 相对 train+val 的最大 identity（mmseqs2 只接受标准氨基酸，FASTA 写出时去除 `|`）；支持 `--mmseqs-cache` 重放已有 convertalis 结果。
- **superseded**：本脚本的 mmseqs2 段为历史首遍实现，统一口径与重放工具以 `.skills/homology-ood-annotation` 的 `scripts/seq_homology_ood.py` 为准（easy-search / `--m8-cache`）。

字段说明：

| 字段 | 说明 |
|------|------|
| `unique_id` | 样本唯一标识，即 `pdb_id` |
| `aa_seq` | 蛋白质氨基酸序列；多链样本用 `\|` 分隔各链 |
| `struct_file` | 蛋白结构文件名（仅文件名）：`{pdb_id}.pdb`，位于 `pdbs/` |
| `ligand_smiles` / `ligand_ecfp4` | 配体 canonical SMILES / ECFP4 ON-bits 列表 |
| `label` | 亲和力标签 `pK`，**强制保存为字符串，统一保留 6 位小数** |

### 4.3 结构相似性与 IDR OOD

**脚本**：`scripts/10_compute_structure_idr_ood.py`
**输入**：`splits/ligand_binding_affinity_{train,val,test}.csv`、`pdbs/`
**输出**：就地更新 test.csv（追加 `TM-score_0.9..0.3`、`idr_ratio`、`OOD_IDR`）、更新 `output/cluster_split_stats.json`

- **Foldseek-Multimer 复合物比对（2026-08-29 起）**：`easy-multimersearch`（query=test 2,367 复合物结构，target=train+val 9,296 结构），取每对复合物 `complexqtmscore`/`complexttmscore` 对称均值的最大值；foldseek createdb 会将多链 PDB 拆链（`<file>_<chain>`），聚合时归一回复合物 stem（单链文件为裸名）。复合物级 TM 最贴合"复合物-配体亲和力"语义。旧口径为单链 `search -a` + `qtmscore/ttmscore` 对称均值（注意与其他任务的 `alntmscore` 均不一致），三方对比（旧对称均值 / 单链 alntm max / 复合物 multimer）见 `analysis/output/foldseek_easysearch_audit/lba/multimer/three_way_compare.txt`，旧 test 备份于 `output/backup_before_multimer_tm_unify/`，per-stem 缓存 `output/foldseek_multimer_max_tm_per_stem.csv`；
- metapredict 预测内在无序区域占比，`idr_ratio > 0.1` → `OOD_IDR`；统一工具见 `.skills/ood-annotation-toolkit`（`scripts/idr_ood.py`）。

### 4.4 CATH Fold / Superfamily Holdout（CATH+TED 并集终值）

**脚本**：`scripts/11_ted_augment_cath_ood.py`
**输入**：`splits/ligand_binding_affinity_{train,val,test}.csv`、项目级 `data/sifts/`、`data/CATHv44/`
**输出**：`output/ligand_ted_aug_chain_list.csv`、`ligand_chain_cath_ted_labels.csv`、`ligand_chain_domains.csv`、`ligand_pdb_cath_ted_labels.csv`、`ligand_pdb_cath_ted_holdout.csv`；就地更新 test.csv 两列

终值口径为 **CATH+TED 并集**（CATH v4.4 实验标注优先，TED putative 补缺），统一流程见 `.skills/cath-ted-domain-annotation`；holdout 判定（train-only 参考、any 语义、无注释 → False）以该 skill 的 `scripts/cath_ted_holdout.py` 为唯一实现。流程：

1. **链清单**：11,663 个 PDB ×（SIFTS 两文件 ∪ `cath-domain-list.txt` 中出现的链）——仅 SIFTS 播种会漏 **944 条 CATH 已标注链**（342 个 PDB，含 87 个 test PDB 的 278 链）；
2. 调用 skill 的 `cath_ted_annotate.py`（`--offline`，TED 经 50% 域重叠过滤，H 级提供 C.A.T.H、T 级仅 C.A.T；API 缓存于项目级 `data/ted/ted_api_cache.json`）；
3. 链级标签按 PDB 聚合取并集（PDB 级口径）；**对称增强**——train 自身同样用 CATH+TED 并集作参考集，避免 test 的 TED-only 标签假性"未见"；
4. holdout 判定调用 skill 的 `cath_ted_holdout.py`，仅就地更新 test 两列（其余列不动）。

**数据规模**：27,352 链（11,663 PDB）；标注来源 both 18,990 / cath-only 3,219 / ted-only 2,115 / none 3,028；TED 缓存命中 3,811 个 accession（0 次新查询）。

**结果（test，n=2,367）**：`OOD_FoldHoldout` **38** 条、`OOD_SuperfamilyHoldout` **73** 条、无注释样本（静默 False）149 条。更新前备份：`output/backup_before_ted_aug/`。

> superseded：仅 CATH 的本地基线实现 `16_cath_ec_ood.py` 已归档（该口径下 27.5% 的测试 PDB 无任何标注，Fold/Superfamily 仅 24/45 条）；其 EC 部分同时被 §4.5 的 skill 统一口径取代。

### 4.5 EC 功能 OOD：`OOD_NewEC_*` / `OOD_LongTail_EC_*`（skill 落盘）

> 统一口径、层级解析与驱动脚本见 `.skills/ec-function-ood-annotation`（`scripts/add_unified_newec.py`、`scripts/add_unified_longtail_ec.py`）。本任务 EC 来源：`data/sifts/sifts_chain_ec.tsv.gz`，按 PDB 全链并集。

- `OOD_NewEC_L4` / `OOD_NewEC_L3`：参考集 = **train+val** 出现的 EC 集合；**all-not-in** 语义——样本有 EC 注释且其所有 EC（对应层级）都不在参考集 → True；无 EC 注释 → False；
- `OOD_LongTail_EC_{L3,L4}_le{5,10}`：参考集 = **train-only** 逐样本频次（**含 0**）；all-not-in——样本所有 EC（对应层级）train 频次都 < 5 / < 10 → True；无注释 → False；
- 旧列 `OOD_NewFunction`（any 口径，54 条）与 `OOD_LongTailFunction`（181 条）已由统一口径列**删除替换**；备份：`output/backup_before_newec_unify/`、`output/backup_before_longtail_allnotin/`。

### 4.6 值域 AffBin LongTail

**脚本**：`scripts/12_add_value_bin_longtail_ood.py`
**输入**：`splits/ligand_binding_affinity_{train,test}.csv`
**输出**：就地更新 test.csv（追加 `OOD_LongTail_AffBin_le50` / `le100`）、`output/ligand_longtail_affbin_summary.json`

- `label`（pK）按 **0.5 宽 bin**（范围 [-0.5, 16.5)）分箱，统计 **train** 各 bin 样本数；
- test 样本所属 bin 的 train 计数 ≤ 50 / ≤ 100 → True（本任务 test 全部 Default=True）；
- train pK 分布（n=8,366）近似正态（mean=6.26，std=1.97），LongTail 两端对称；test 阳性样本 label 横跨 (0.66, 12.30)，确为分布两端的极端亲和力样本；
- 本任务 label 无 clip 删失（对照：kcat 任务存在上游 CatPred-DB [-6, 6] 截断）；
- 通用工具见 `.skills/ood-annotation-toolkit`（`scripts/value_bin_longtail_ood.py`）。

### 4.7 InD 收尾与体检（skill 落盘）

> 统一收尾工具：`.skills/ood-annotation-toolkit` 的 `scripts/finalize_ood_columns.py`；体检用同 skill 的 `scripts/check_splits.py`。

- `InD`：所有 `OOD_*` + `seq_Redundancy_*` + `TM-score_*` 标记均为 False（NaN 视为 False）时为 True，由 `finalize_ood_columns.py --mode ind` 在所有 OOD 列齐备后统一重算；当前 **1,091 / 2,367（46.1%）**；备份见 `output/backup_before_ind/`；
- NaN 惯例（`--mode nan`）：本任务 test `Default` 全为 True（LBA/PPI 惯例），无 Default=False 行，故无需 NaN 化；
- `check_splits.py` 体检全部通过（泄露、子集关系等，见 §5）。

### 4.8 Test OOD 统计（终值，当前 splits 实测）

| OOD 维度 | 标记列 | True 数量 |
|----------|--------|----------|
| ID 测试集 | `Default` | 2,367 |
| 极端短链 / 极端长链 | `OOD_ExtremeShort` / `OOD_ExtremeLong` | 15 / 303 |
| 序列冗余 < 90%~30% | `seq_Redundancy_90` ~ `seq_Redundancy_30` | 575 / 498 / 454 / 409 / 368 / 294 / 241 |
| 无同源 hit | `OOD_Orphan` | 181 |
| 结构新颖（复合物 TM < 0.9~0.3） | `TM-score_0.9` ~ `TM-score_0.3` | 767 / 662 / 447 / 306 / 201 / 135 / 79 |
| 高内在无序区（idr > 0.1） | `OOD_IDR` | 93 |
| CATH+TED 超家族 / 折叠留一法 | `OOD_SuperfamilyHoldout` / `OOD_FoldHoldout` | 73 / 38 |
| EC 新功能 L4 / L3（统一口径） | `OOD_NewEC_L4` / `OOD_NewEC_L3` | 46 / 8 |
| EC 长尾 L3（le5 / le10） | `OOD_LongTail_EC_L3_le5` / `le10` | 24 / 37 |
| EC 长尾 L4（le5 / le10） | `OOD_LongTail_EC_L4_le5` / `le10` | 189 / 344 |
| 亲和力 bin 长尾（le50 / le100） | `OOD_LongTail_AffBin_le50` / `le100` | 19 / 65 |
| 纯 ID 样本 | `InD` | 1,247（52.7%） |

当前 test 完整 38 列表头：

```
unique_id, aa_seq, struct_file, ligand_smiles, ligand_ecfp4, label,
Default, OOD_ExtremeShort, OOD_ExtremeLong,
seq_Redundancy_90, ..., seq_Redundancy_30,
TM-score_0.9, ..., TM-score_0.3,
idr_ratio, OOD_IDR, OOD_Orphan,
OOD_SuperfamilyHoldout, OOD_FoldHoldout,
OOD_NewEC_L4, OOD_LongTail_AffBin_le50, OOD_LongTail_AffBin_le100,
OOD_LongTail_EC_L3_le5, OOD_LongTail_EC_L3_le10,
OOD_LongTail_EC_L4_le5, OOD_LongTail_EC_L4_le10,
OOD_NewEC_L3, InD
```

### 4.9 配体化学相似度（Tanimoto）

val/test 每个样本与训练集所有样本的最大 Tanimoto：

| 集合 | 均值 | 中位数 | 最小值 | 最大值 |
|------|------|--------|--------|--------|
| Val | 0.858 | 0.873 | 0.600 | 1.000 |
| Test | 0.846 | 0.835 | 0.600 | 1.000 |

| 相似度区间 | Val | Test |
|------------|-----|------|
| 0.6–0.7 | 179 | 495 |
| 0.7–0.8 | 198 | 551 |
| 0.8–0.9 | 103 | 290 |
| ≥ 0.9 | 450 | 1,031 |

> 所有 val/test 样本与训练集的最大 Tanimoto 均 ≥ 0.6（划分设计的下限）。

---

## 5. 数据质量检查

1. **泄露检查**：Train ∩ Val、Train ∩ Test、Val ∩ Test 的 `unique_id` 交集均为空（`ood-annotation-toolkit` 的 `check_splits.py` 体检全部通过，`OOD_Orphan ⊆ seq_Redundancy_30` 子集关系成立）；
2. **划分设计防泄露**：同一配体相似度簇（pairwise Tanimoto ≥ 0.6）不跨 train/val/test；val/test 与 train 最大 Tanimoto ≥ 0.6 为刻意保留的化学相似性下限；
3. **分布一致性**：长度、Kd/Ki 比例在各子集间一致；test pK 范围 (0.66, 12.30) 落在 train 分布内；
4. **序列-结构一致性**：`aa_seq` 与 `pdbs/` 结构逐样本对齐（按链提取口径，见 §3.1；历史上重提取验证不一致 0 条）；
5. **异常值处理**：`2pll` 配体解析失败被移除（§2.2）；IC50 与不等号/异常注释记录已在清洗阶段排除（§2.1）。

---

## 6. 输出文件汇总

```
ligand_binding_affinity/
├── DATA_PROCESS.md              # 本文档
├── data/                        # 原始数据
│   ├── index/
│   └── P-L/
├── pdbs/                        # 整蛋白 PDB 文件（11,663 个，{pdb_id}.pdb）
├── scripts/                     # 处理脚本（编号见文首总览）
├── output/
│   ├── pdbbind_v2020_metadata.csv
│   ├── pdbbind_v2020_metadata_cleaned.csv
│   ├── pdbbind_v2020_metadata_ligand.csv
│   ├── pdbbind_v2020_metadata_seq.csv
│   ├── pdbbind_v2020_dataset.csv          # 样本级数据集（阶段四输入）
│   ├── pdbbind_v2020_metadata_dedup95*.csv  # 旧去冗余路线留痕（未采用）
│   ├── ligand_parse_failures.csv
│   ├── whole_protein_pdb_stats.json
│   ├── cluster_split_stats.json
│   ├── {test,train_val}_for_mmseqs.fasta / mmseqs_test_vs_trainval.tsv
│   ├── ligand_ted_aug_chain_list.csv      # PDB×链清单（含 splits 归属）
│   ├── ligand_chain_cath_ted_labels.csv   # 链级 CATH+TED 标签
│   ├── ligand_chain_domains.csv           # 域级明细（81,029 行）
│   ├── ligand_pdb_cath_ted_labels.csv / ligand_pdb_cath_ted_holdout.csv
│   ├── ligand_longtail_affbin_summary.json
│   ├── scripts_archive/                   # 归档脚本（见 §8）
│   └── backup_before_*/                   # 各次 splits 更新前备份
└── splits/
    ├── ligand_binding_affinity_train.csv  # 8,366
    ├── ligand_binding_affinity_val.csv    # 930
    └── ligand_binding_affinity_test.csv   # 2,367，含完整 OOD 标记列
```

---

## 7. 诊断与质控工具

本任务无独立 `diag_*` 脚本；质控通过以下途径完成：

- `.skills/ood-annotation-toolkit/scripts/check_splits.py`：splits 体检（泄露、列完整性、子集关系）；
- `08_cluster_time_split.py` 输出的 `output/cluster_split_stats.json`：聚类与划分统计、max_train_tanimoto 分布；
- 各步骤的 `*_stats.json` / 失败清单 CSV（见 §6 目录树）。

---

## 8. 归档脚本清单（`output/scripts_archive/`）

| 脚本 | 归档原因 |
|:---|:---|
| `05_dedup_sequences_95.py` | 旧去冗余路线（95% 序列同一性），评估后未采用 |
| `06_dedup_by_ligand_95.py` | 旧去冗余路线（配体 95%），评估后未采用 |
| `07_dedup_by_smiles_95.py` | 旧去冗余路线（SMILES 95%），评估后未采用 |
| `08_extract_protein_chains.py` | 链级结构路线已被整蛋白路线取代；`pdbbind_v2020_dataset.csv` 现由 `06_build_dataset_csv.py` 生成 |
| `10_time_split_and_tanimoto.py` | 旧划分路线，已被 `08_cluster_time_split.py` 取代 |
| `12_balance_val_from_test.py` | 比例平衡逻辑已融入 `08_cluster_time_split.py` |
| `13_adjust_val_to_train_9to1.py` | 同上（train:val 9:1 调整已融入 08） |
| `16_cath_ec_ood.py` | 仅 CATH 基线与 EC 旧口径，已被 `11_ted_augment_cath_ood.py` + ec-function-ood-annotation skill 双重取代 |
| `17_reextract_sequences_from_pdbs.py` | 按链 `\|` 序列提取逻辑已融入 `05_extract_protein_sequences.py` |

---

## 9. 注意事项（历史演进、无脚本环节与已知偏差）

1. **LBA/PPI 惯例——Default 全 True**：本任务 train/val 保留极端长度样本（train 含 67 条 < 60 aa、1,008 条 > 1,000 aa 总长的样本），test 不做极端长度剔除，`Default=True` 覆盖全部 2,367 个测试样本；`OOD_ExtremeShort` / `OOD_ExtremeLong` 仅标记不移除，无 Default=False 行，故无 NaN 化需求。
2. **多链复合物**：`aa_seq` 以 `|` 分隔各链序列（test 1,264 条、train 4,591 条多链样本），下游工具使用前应先去除分隔符；长度统计为各链求和。
3. **超长单链集中在边界**：> 1,000 aa 的**单链**样本长度集中在 1,014–1,016 边界（train 25 条单链超长中 16 条落在 1,014–1,016，val 3 条为 1,015–1,016；test 仅 4 条，为 1,029–1,195）；其余超长样本均为多链复合物总长度超限；无单链 > 2,000 aa 样本。
4. **≤3aa 短链片段为原始沉积特征（2026-08-28 核查）**：train 7 条 / val 2 条 / test 2 条样本含 ≤3aa 的链片段，已逐样本核对 `data/P-L/` PDBbind 原始文件与 RCSB 镜像，全部来自 PDB 沉积本身而非提取 artifact：肽类抑制剂链（1eb1 仅存 PRO 2）、UNK 占位/短肽 stub（2j9n GLN+UNK、3e4a/6byz/4hva/3t6y/3t64/3t70 poly-Ala/His）、交替构象去重塌缩（5jjm Ser-Asp ×12 构象）、HETATM 沉积的修饰残基（4aze SEP→S，注意其记录名为 `HETATM17267` 无空格溢出格式）。按链提取口径（ATOM+HETATM、修饰残基映射）忠实保留了这些链，无需修复。
5. **浮点重序列化漂移（2026-08-28 验证发现）**：`output/pdbbind_v2020_dataset.csv` 经历史多次中间文件读写，`3o56` 的 `affinity_value` 存在 1 ulp 级漂移（当前 `7.939999999999993e-10`，由 05+06 重新生成为 `7.939999999999995e-10`，原始 `metadata.csv` 为 `7.939999999999999e-10`）；`pK` 字符串不受影响（逐位一致），对标签与划分无任何影响。属历史口径偏差，留痕不修复。
6. **整理验证记录（2026-08-28）**：新 `05`+`06` 链从 `pdbbind_v2020_metadata_ligand.csv` 重生成 `pdbbind_v2020_dataset.csv`，与现有文件 11,663 行 × 36 列逐值 diff 一致（仅上述 3o56 单格 1 ulp 偏差）；合并后的 `08_cluster_time_split.py` 全量重跑（5,394 簇、移动 2,866 + 1,183 + 442），train/val/test 样本归属与当前 splits **完全一致**。
7. **历史演进**（终值已融入正文，此处仅留痕）：FoldHoldout/SuperfamilyHoldout 由仅 CATH 基线（24/45）经 TED 增强重算为终值（38/73）；NewEC/LongTailEC 统一为 all-not-in 口径并删除旧列 `OOD_NewFunction`/`OOD_LongTailFunction`；AffBin LongTail 列后补；`InD` 由 toolkit 统一重算。各次更新的 splits 备份在 `output/backup_before_*/`。
