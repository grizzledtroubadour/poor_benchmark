# 酶最适 pH 预测（optimal_ph_prediction）数据集处理文档

> 任务名称：optimal_ph_prediction（回归，预测酶最适 pH）
> 数据来源：Zenodo 记录 [18405148](https://zenodo.org/records/18405148) —— *Dataset EnzyBase12k: A Curated Dataset for Predicting Enzyme pH Optima Using pHoptNN*
> 处理位置：`datasets/enzyme_optimal_ph/`
>
> 本文档描述**从原始数据一次性复现**当前 `splits/` 与 `pdbs/` 的完整线性流程（四阶段结构）。
> 统计数字与列名均已按当前 `splits/*.csv` 实际表头重算核对（只读）。

## 任务定义与当前规模

回归任务：由酶序列/结构预测其最适 pH（`label` = pH optimum，字符串保存）。

| 划分 | 文件 | 样本数 | 说明 |
|:---|:---|---:|:---|
| 训练集 | `splits/optimal_ph_prediction_train.csv` | 7,794 | AF 结构，4 列 |
| 验证集 | `splits/optimal_ph_prediction_val.csv` | 866 | AF 结构，4 列 |
| 测试集 | `splits/optimal_ph_prediction_test.csv` | 2,955 | PDB + 极端长度，34 列（含完整 OOD 块） |
| **合计** | — | **11,615** | `pdbs/` 下 11,615 个 `.cif` |

测试集内部：`Default=True` 2,717（原 PDB 样本）；`OOD_ExtremeLong` 211（`seq_length > 1000`）；`OOD_ExtremeShort` 27（`seq_length < 60` 的 24 条 + 2026-09-01 从 train 补移的 3 条，见 §4.2.1），后两者 `Default=False`。

脚本总览（`scripts/`，按执行顺序编号）：

| 编号 | 脚本 | 阶段 |
|:---|:---|:---|
| 01 | `01_verify_structures.py` | ① 结构解压与校验 |
| 02 | `02_split_dataset.py` | ④ 初分（极端长度留出 / PDB 作测试 / AF 9:1） |
| 03 | `03_prepare_splits_and_pdbs.py` | ④ 统一 schema 标准化 + pdbs 落盘 |
| 04 | `04_compute_ood_markers.py` | ④ 序列/结构同源 + Orphan + IDR OOD |
| 05 | `05_add_phbin_longtail_ood.py` | ④ pH 值域长尾 OOD |
| 06 | `06_ted_augment_cath_ood.py` | ④ CATH+TED Fold/Superfamily Holdout |
| 07 | `07_move_short_train_rows_to_test.py` | ④ 极端长度口径修正：aa_seq<60 的 train 行移入 test（§4.2.1） |
| 诊断 | `diag_analyze_metadata.py`、`diag_validate_sequence_structure.py`、`diag_cluster_sequences_mmseqs.py` | 诊断与质控（见 §7） |

---

## 阶段一：源数据获取与解析

### 1.1 Zenodo 文件列表与下载

| 文件名 | 大小（API 标注） | 说明 |
|:---|---:|:---|
| `metadata.csv` | ~6.4 MB | 元数据，11,615 行 × 13 列 |
| `DATASET.zip` | ~1.87 GB | 包含 `.pdb` 和 `.pqr` 格式结构文件 |
| `Dataset_mmcifs.zip` | ~559 MB | 包含 `.cif` 格式结构文件 |

```bash
# 元数据
curl -L -o data/metadata.csv "https://zenodo.org/api/records/18405148/files/metadata.csv/content"

# 结构文件（大文件，建议使用断点续传）
wget -c -t 0 --timeout=300 -O data/DATASET.zip \
  "https://zenodo.org/api/records/18405148/files/DATASET.zip/content"

wget -c -t 0 --timeout=300 -O data/Dataset_mmcifs.zip \
  "https://zenodo.org/api/records/18405148/files/Dataset_mmcifs.zip/content"
```

下载后的文件存放在 `data/` 目录，下载日志见 `logs/`。

### 1.2 结构解压与校验

**脚本**：`scripts/01_verify_structures.py`
**输入**：`data/DATASET.zip`、`data/Dataset_mmcifs.zip`、`data/metadata.csv`
**输出**：`data/dataset/`（`DATASET/EnzyBase12k_pbs|pqrs`）、`data/dataset_mmcifs/`（`EnzyBase12k_mmcifs`）、`output/structure_verification.json`

- 解压两个 zip（已解压则跳过），统计文件数量并将文件名与 `uniprot_id` / `pdb_id_final` 匹配；
- 校验结果：`.pdb` / `.pqr` / `.cif` 各 **11,615** 个，与 metadata 行数一致，文件名均与 `uniprot_id` 匹配。

### 1.3 元数据字段

| 字段 | 类型 | 说明 |
|:---|:---|:---|
| `uniprot_id` | str | UniProt 登录号（即后续的 `unique_id`） |
| `ph_optimum` | float | 酶最适 pH（回归标签） |
| `structure_source` | str | 结构来源：`PDB` / `AF` |
| `structure_mode` | str | 结构模式：`AF2_downloaded` / `AF3` / `AF3_template` |
| `pdb_id_final` | str | 最终 PDB ID（仅 PDB 来源） |
| `resolution` | float | 分辨率（Å，仅 PDB 来源） |
| `pdb_chain` | str | PDB 链标识（仅 PDB 来源） |
| `ec_id` | str | EC 编号（**782 个 UniProt 为 `; ` 分隔的多 EC 字符串**，见 §9） |
| `cut_to_uniprot_domain` | int | 是否截断至 UniProt 结构域 |
| `systematic_name` | str | 系统命名 |
| `organism` | str | 来源物种 |
| `seq_length` | int | 序列长度 |
| `uniprot_seq_cut` | str | 截断后的氨基酸序列 |

元数据核心统计（诊断脚本 `diag_analyze_metadata.py` 产出）：样本 11,615（UniProt 无重复）；pH 均值 7.20、中位 7.50、范围 1.5–12.5；AF 占 76.3%、PDB 占 23.7%；PDB 分辨率中位 1.95 Å（n=2,735）；序列长度中位 386、范围 32–1500；EC 注释覆盖率 98.6%。

---

## 阶段二：标签聚合与计算

本任务标签直接来自源数据：`label` = `metadata.csv` 的 `ph_optimum`，无多源聚合步骤。EC 注释（`ec_id`）作为阶段四 EC 功能 OOD 的输入，同样直接取自 `metadata.csv`（本任务历史 EC 源，见 §9 注意事项）。

---

## 阶段三：去冗余与数据精化

> 本任务最终**未按任何质量/冗余条件过滤样本**；本节步骤均为诊断性评估，结论仅作记录。

### 3.1 序列-结构一致性验证（诊断）

**脚本**：`scripts/diag_validate_sequence_structure.py`
**输入**：`data/metadata.csv` + 解压后的结构文件（优先 `.cif`）
**输出**：`output/sequence_structure_validation.csv` 及 `_summary.json`

将结构中提取的氨基酸序列（第一条链）与 `uniprot_seq_cut` 做 edlib 全局比对：验证 11,615 样本全部解析成功；identity 均值/中位 0.938/0.959；identity ≥ 0.90 占 77.6%，< 0.50 为 0。22.4% 样本因结构截断 identity 在 0.70–0.90 之间，无严重不一致；**未按 identity 过滤**。

### 3.2 MMseqs2 95% 聚类（诊断，未应用）

**脚本**：`scripts/diag_cluster_sequences_mmseqs.py`
**输入**：`data/metadata.csv` 的 `uniprot_seq_cut`
**输出**：`output/sequences.fa`、`output/mmseqs_cluster_95/`、`output/mmseqs_cluster_95_summary.json`

参数：min-seq-id 0.95、coverage 0.8、cov-mode 0。结果：11,615 序列 → 11,058 簇（单例 10,622，最大簇 9），冗余 557 条（4.8%）。**仅用于说明源数据冗余度低，未据此过滤**，最终划分保留全部 11,615 条。

### 3.3 候选过滤条件汇总

| 条件 | 建议阈值 | 是否应用 | 说明 |
|:---|:---|:---|:---|
| 序列长度 | 60–1000 | ✅（划分层面） | 界外样本移入测试集作极端长度 OOD，见 §4.1 |
| PDB 分辨率 | ≤ 3.0 Å | ❌ | 候选条件，未过滤 |
| pH 异常值 | 1.0–13.0 | ❌ | 实际范围 1.5–12.5，无需过滤 |
| 序列-结构一致性 | identity ≥ 0.90 | ❌ | 候选条件，未过滤 |
| 序列冗余 | MMseqs2 95% identity | ❌ | 仅诊断（§3.2），未过滤 |

---

## 阶段四：数据划分与 OOD 标注

### 4.1 初分

**脚本**：`scripts/02_split_dataset.py`
**输入**：`data/metadata.csv`
**输出**：`splits/{train,val,test,ood_extremelength}.csv`（中间产物，下一步消费后删除）、`output/split_summary.json`

1. **OOD-ExtremeLength**：`seq_length < 60` 或 `> 1000` 的样本（235 条）单独留出；
2. **Test**：剩余样本中所有 `structure_source == "PDB"` 的样本（2,717 条）；
3. **Train / Val**：剩余 AF 样本按 9:1 随机划分（`random_seed=42`）。

### 4.2 统一 schema 标准化与结构落盘

**脚本**：`scripts/03_prepare_splits_and_pdbs.py`
**输入**：`splits/{train,val,test,ood_extremelength}.csv` + 解压结构文件
**输出**：`splits/optimal_ph_prediction_{train,val,test}.csv`、`pdbs/*.cif`（11,615 个）、`output/split_mapping.json`

- 将 `ood_extremelength.csv` 并入 test 并删除初分中间文件；
- 原 PDB 测试样本 `Default=True`，极端长度样本 `Default=False`；
- 新增 `OOD_ExtremeLong`（`seq_length > 1000`）与 `OOD_ExtremeShort`（`seq_length < 60`）；
- `aa_seq` 从结构文件（优先 `.cif`）中提取；`label`（pH optimum）强制保存为字符串；
- 全部结构文件复制到 `pdbs/`，文件名 = `unique_id`（统一 `.cif`，无则回退 `.pdb`）；
- **直接产出统一 schema**：train/val 为 `unique_id, aa_seq, struct_file, label`（4 列），test 追加 `Default` / `OOD_ExtremeLong` / `OOD_ExtremeShort`（历史注记：初版产出为 `unique_id, file_type, aa_seq, labels`，经全项目 schema 迁移后定型，迁移记录见 §9）。

最终划分统计（按当前 splits 重算）：

| 数据集 | 样本数 | 结构来源 | pH 均值 | aa_seq 长度中位数 |
|:---|---:|:---|---:|---:|
| 训练集 | 7,794 | AF | 7.129 | 367 |
| 验证集 | 866 | AF | 7.139 | 363 |
| 测试集 | 2,955 | PDB + 极端长度 | 7.390 | 347 |

> pH 均值与长度中位数按当前 `splits/*.csv` 的 `label` / `aa_seq` 列重算（`aa_seq` 为结构提取序列，与 metadata 的 `seq_length` 口径不同）。

### 4.2.1 极端长度口径修正（2026-09-01）

**脚本**：`scripts/07_move_short_train_rows_to_test.py`

§4.1 的极端长度预留按 metadata `seq_length`（UniProt 切域长度）判定，而 `aa_seq` 为结构提取序列，两者口径不同。这导致 3 条 metadata 长度 ≥60 但结构序列 <60 的样本（`P05959` 55aa、`B1KRG2` 59aa、`P0C7A9` 46aa，label 均 8.0）误入 train。本步将这 3 条从 train 移入 test：`Default=False`、`OOD_ExtremeShort=True`、`OOD_ExtremeLong=False`、其余 OOD 列 NaN、`InD=False`（结构文件本就在 `pdbs/`，无需移动）。结果：train 7,797 → 7,794，test 2,952 → 2,955，`OOD_ExtremeShort` 24 → 27。变更前备份：`output/backup_20260901/`。

> 已知挂起：长端存在同源的口径错位（约 90 行 `OOD_ExtremeLong=True` 但 aa_seq 为 747–999aa），本次未统一，EL 仍按 `seq_length` 口径；如需统一为 `len(aa_seq)` 口径需重判 EL 标注并同步本节与统计。

### 4.3 OOD 标注总约定

除 `Default` / `OOD_ExtremeLong` / `OOD_ExtremeShort` / `InD` 外，所有 OOD 列（含 `seq_Redundancy_*`、`TM-score_*`）**仅对 `Default=True` 样本计算**；`Default=False`（极端长度）行这些列一律为 **NaN**（空单元格，表"未计算"），由 §4.10 的收尾步骤统一落实。

OOD 方法学统一口径见 `.skills/` 下的 `homology-ood-annotation`、`ood-annotation-toolkit`、`cath-ted-domain-annotation`、`ec-function-ood-annotation` 四个 skill；下面各步先列本地脚本，再列对应 skill 口径/复算工具。

### 4.4 序列同源性 OOD：`seq_Redundancy_*` 与 `OOD_Orphan`

**脚本**：`scripts/04_compute_ood_markers.py`
**输入**：`splits/optimal_ph_prediction_{train,val,test}.csv`
**输出**：就地更新 test.csv、`output/ood_markers_summary.json`

- 工具：MMseqs2 `easy-search`（`-s 7`），Default Test vs Train+Val，取每个 query 的最大 pident；
- `seq_Redundancy_{90..30}`：max pident **<** 阈值（步长 10%）→ True。**2026-08-29 起由 `≤` 统一为 `<`**（与全项目一致）：pident 一位小数精度下恰等于阈值的 35 行曾多标 True（各档 3/1/3/4/9/8/7 行），已重算修正（备份 `output/backup_before_threshold_lt_unify/`；重跑验证产物 `analysis/output/ph_threshold_audit/`，m8 缓存 `output/ood_mmseqs_test_vs_trainval.m8`）；P14618（pident 恰为 90.0）因此失去唯一 OOD 标记转入 `InD`（49 → 50）；
- `OOD_Orphan`：在 Train+Val 中无显著 hit → True。
- 统一口径与复算工具：`.skills/homology-ood-annotation`（`scripts/seq_homology_ood.py`，支持 easy-search 与 `--m8-cache` 重放）。

### 4.5 结构相似性 OOD：`TM-score_*`

**脚本**：`scripts/04_compute_ood_markers.py`（同一步内完成）

- 工具：Foldseek `easy-search`，Default Test vs Train+Val 结构（`pdbs/` 下 `.cif`），取最大 alntmscore；
- `TM-score_{0.9..0.3}`：max TM-score **<** 阈值（步长 0.1）→ True（2026-08-29 与 seq 列一并由 `≤` 统一为 `<`；实测无恰好等于阈值的样本，数值不变。m8 缓存 `output/ood_foldseek_test_vs_trainval.m8`）。
- 统一口径与复算工具：`.skills/homology-ood-annotation`（`scripts/struct_homology_ood.py`）。

### 4.6 无序区域 OOD：`OOD_IDR`

**脚本**：`scripts/04_compute_ood_markers.py`（IDR 口径已并入本脚本一次算对；原独立重算脚本 `compute_ood_idr.py` 已归档）

- 工具：metapredict v3（`predict_disorder_fasta` 批量预测），残基 disorder propensity > 0.5 判无序；
- 口径：**residue-ratio**（无序残基数 / 序列长），阈值 **> 0.1** → True，仅 Default 样本（Default=False 行在保存前直接置 NaN，见 §4.10）；
- 本任务历史口径为 0.1，与 ss/func/fold/ppi 的 0.3 不同；跨任务口径对照见 `.skills/ood-annotation-toolkit`（`scripts/idr_ood.py` 口径表）；
- 中间产物：`output/test_sequences_for_idr.fasta`；阳性 **39** 条（Default n=2,717），数据集以结构化酶为主。

### 4.7 折叠 / 超家族留一法 OOD：`OOD_FoldHoldout` / `OOD_SuperfamilyHoldout`

**脚本**：`scripts/06_ted_augment_cath_ood.py`（wrapper）
**输入**：`splits/optimal_ph_prediction_test.csv`、`output/optimal_ph_ted_aug_acc_list.csv`（键清单）
**输出**：`output/optimal_ph_acc_cath_ted_labels.csv`，就地更新 test.csv 两列

标注与 holdout 判定统一调用 `.skills/cath-ted-domain-annotation` 的 `scripts/cath_ted_annotate.py` 与 `scripts/cath_ted_holdout.py`（全项目唯一实现）。

**口径**（UniProt accession 模式，`unique_id` 即 accession）：

- 标签 = CATH v4.4 实验标注（经 SIFTS UniProt→PDB 链反查并集）∪ TED putative 标注（REST API，项目级共享缓存 `data/ted/ted_api_cache.json`，按全长不过滤；H 级提供 C.A.T.H，T 级仅 C.A.T）；
- 参考集 = **Train** 的 CATH+TED 并集（train/val/test 对称增强）；
- any 语义：任一样本 Topology / Superfamily 未在参考集出现 → True；无注释 → False。

**数据规模与覆盖率**：

- 键清单：11,615 个 UniProt accession（`output/optimal_ph_ted_aug_acc_list.csv`），TED 缓存全部命中（离线运行，0 次新查询）；
- 标注来源分布：ted-only 8,167 / both 2,787 / none 595 / cath-only 66；

| 口径 | Topology 覆盖率（全部 11,615 个 accession） |
|------|----------------|
| 仅 CATH | 24.6% |
| CATH+TED | **94.9%** |

**Holdout 结果**（test Default n=2,717）：`OOD_FoldHoldout` 91 条、`OOD_SuperfamilyHoldout` 311 条、无注释样本（→ False）65 条。更新前备份见 `output/backup_before_ted_aug/`。

### 4.8 EC 功能 OOD：`OOD_NewEC_*` 与 `OOD_LongTail_EC_*`（skill 落盘）

**统一口径**（方法学见 `.skills/ec-function-ood-annotation`；本任务由该 skill 的批量驱动 `scripts/add_unified_newec.py`、`scripts/add_unified_longtail_ec.py` 落盘）：

- EC 来源：`data/metadata.csv` 的 `ec_id`（本任务历史 EC 源）；层级解析允许部分 `-`（如 `3.4.22.-` → L3=`3.4.22`）；
- **多 EC 字符串拆分**：`ec_id` 含 `; ` 分隔的多 EC 字符串（782 个 UniProt），统一口径按拆分后计算；
- `OOD_NewEC_L3` / `OOD_NewEC_L4`：**all-not-in**——样本有 EC 注释且其所有 EC（对应层级）都不在 **train+val** 参考集 → True；无注释 → False；
- `OOD_LongTail_EC_{L3,L4}_le{5,10}`：**all-not-in**——样本有 EC 注释且其所有 EC（对应层级）在 **train** 中频次都 < 5 / < 10（**含 0**）→ True；无注释 → False；
- 落盘摘要：`output/newec_unified_summary.json`、`output/longtail_ec_unified_summary.json`；更新前备份：`output/backup_before_newec_unify/`、`output/backup_before_longtail_unify/`。

| 列 | True 数（Default n=2,717） |
|------|------:|
| `OOD_NewEC_L3` | 20 |
| `OOD_NewEC_L4` | 647 |
| `OOD_LongTail_EC_L3_le5` | 118 |
| `OOD_LongTail_EC_L3_le10` | 221 |
| `OOD_LongTail_EC_L4_le5` | 1,763 |
| `OOD_LongTail_EC_L4_le10` | 2,208 |

> 注意：本任务 train 仅 7,794 样本且按 UniProt 计 EC 频次，L4 级长尾占比高达 64.9%–81.3%，区分度有限，**L3 级更有意义**。

> superseded：EC 类 OOD 的本地首遍实现 `compute_additional_ood.py`（单一 EC 字符串、≤ 阈值口径的旧列 `OOD_NewFunction` / `OOD_LongTailFunction_L3_le{5,10}`）已归档，被上述 skill 统一口径取代。

### 4.9 pH 值域长尾 OOD：`OOD_LongTail_pHbin_*`

**脚本**：`scripts/05_add_phbin_longtail_ood.py`（抽取自已归档 `compute_additional_ood.py` 的 pHbin 部分，列名已适配统一 schema）
**输入**：`splits/optimal_ph_prediction_{train,test}.csv`
**输出**：就地更新 test.csv、更新 `output/ood_markers_summary.json`

- 分箱：pH 按 **[1.5, 13.5) 宽 0.5** 区间分箱（`pd.cut`）；
- 判定：Default 测试样本所属 bin 的 **train** 频次 ≤ 阈值 → True；
- 统一复算工具：`.skills/ood-annotation-toolkit`（`scripts/value_bin_longtail_ood.py`，其参数表含本任务口径）。

| 列 | 阈值 | True 数（Default n=2,717） |
|------|------|------:|
| `OOD_LongTail_pHbin_le50` | train bin 样本数 ≤ 50 | 25 |
| `OOD_LongTail_pHbin_le100` | train bin 样本数 ≤ 100 | 58 |

### 4.10 非 Default 行 NaN 校验与 `InD`（skill 落盘）

**工具**：`.skills/ood-annotation-toolkit/scripts/finalize_ood_columns.py`（先 `--mode nan` 再 `--mode ind`）；体检用同 skill 的 `scripts/check_splits.py`。

- **NaN 惯例**：`Default=False` 行除 `OOD_ExtremeShort` / `OOD_ExtremeLong` 外的所有 OOD 标记列（含 `seq_Redundancy_*`、`TM-score_*` 与全部 `OOD_*` 计算列）为 **NaN**（空单元格，表"未计算"），与 `False` 区分。2026-08-28 起该 NaN 由各流水线脚本（`04`/`05`/`06` 及 skill 标注脚本）在生成位置直接写入，`finalize_ood_columns.py --mode nan` 仅作**兜底校验**（当前数据上为 no-op）；历史备份见 `output/backup_before_nondefault_nan/`；
- **`InD`**：所有 `OOD_*` + `seq_Redundancy_*` + `TM-score_*` 标记均为 False（NaN 视为 False）→ True；结果 **50** 条（占测试集 1.7%）；备份见 `output/backup_before_ind/`；
- 读写惯例：`pd.read_csv(dtype=str, keep_default_na=False)` + `to_csv(index=False, lineterminator="\n")`，写回前备份、写后验证非目标列不变。

### 4.11 OOD 统计总表（当前 splits 实测，Default n=2,717）

| 列 | 定义摘要 | True 数 |
|:---|:---|---:|
| `seq_Redundancy_90` ~ `30` | max pident < 90/80/70/60/50/40/30 | 2,540 / 2,396 / 2,227 / 2,025 / 1,674 / 1,197 / 644 |
| `TM-score_0.9` ~ `0.3` | max TM-score < 0.9/0.8/0.7/0.6/0.5/0.4/0.3 | 926 / 411 / 171 / 40 / 8 / 4 / 3 |
| `OOD_Orphan` | Train+Val 无显著 hit | 422 |
| `OOD_IDR` | IDR 残基占比 > 0.1 | 39 |
| `OOD_FoldHoldout` | CATH+TED Topology 未在 Train 出现 | 91 |
| `OOD_SuperfamilyHoldout` | CATH+TED Superfamily 未在 Train 出现 | 311 |
| `OOD_NewEC_L3` / `L4` | EC L3/L4 all-not-in（参考 train+val） | 20 / 647 |
| `OOD_LongTail_EC_L3_le5` / `le10` | EC L3 train 频次 < 5 / < 10（含 0，all-not-in） | 118 / 221 |
| `OOD_LongTail_EC_L4_le5` / `le10` | EC L4 train 频次 < 5 / < 10（含 0，all-not-in） | 1,763 / 2,208 |
| `OOD_LongTail_pHbin_le50` / `le100` | pH bin train 频次 ≤ 50 / ≤ 100 | 25 / 58 |
| `InD` | 全部 OOD 标记为 False（NaN 视为 False） | 50 |

当前 test 完整 34 列表头：

```
unique_id, aa_seq, struct_file, label, Default, OOD_ExtremeLong, OOD_ExtremeShort,
seq_Redundancy_90, seq_Redundancy_80, seq_Redundancy_70, seq_Redundancy_60,
seq_Redundancy_50, seq_Redundancy_40, seq_Redundancy_30,
TM-score_0.9, TM-score_0.8, TM-score_0.7, TM-score_0.6, TM-score_0.5, TM-score_0.4, TM-score_0.3,
OOD_Orphan, OOD_IDR, OOD_NewEC_L3, OOD_LongTail_pHbin_le50, OOD_LongTail_pHbin_le100,
OOD_FoldHoldout, OOD_SuperfamilyHoldout,
OOD_LongTail_EC_L3_le5, OOD_LongTail_EC_L3_le10, OOD_LongTail_EC_L4_le5, OOD_LongTail_EC_L4_le10,
OOD_NewEC_L4, InD
```

---

## 5. 数据质量检查

### 5.1 泄露检查

| 检查项 | 结果 |
|:---|---:|
| Train ∩ Val | 0 |
| Train ∩ Test | 0 |
| Val ∩ Test | 0 |
| Train / Val / Test 内部重复 unique_id | 0 |

### 5.2 文件完整性

| 检查项 | 结果 |
|:---|---:|
| `splits/` 下最终文件数 | 3（train / val / test） |
| `pdbs/` 下结构文件数 | 11,615（全部 `.cif`） |
| 每个 split 样本均有对应结构文件 | ✅ |
| `label` 列保存为字符串 | ✅ |
| train / val 不含 OOD 列；test 含完整 OOD 块（34 列） | ✅ |
| `Default=False` 行计算型 OOD 列为 NaN | ✅ |
| `OOD_Orphan` ⊆ `seq_Redundancy_30`（422 ≤ 644） | ✅ |

### 5.3 异常值说明

- `OOD_ExtremeLong`（211）与 `OOD_ExtremeShort`（24）样本 `Default=False`，单独作为极端长度 OOD 子集评估，不参与其他 OOD 维度计算（对应列为 NaN）；
- `OOD_IDR` 阳性 39 条（1.4% Default），数据集以结构化酶为主，高无序样本稀少；
- `OOD_NewEC_L3` 阳性 20 条（0.7% Default），EC L3 功能类别在 Train+Val 与 Test 之间高度重叠；
- `OOD_LongTail_EC_L3_le10` 覆盖 221 条（8.1% Default），EC L3 类别存在明显长尾分布；L4 级长尾占比 64.9%–81.3%，区分度有限（见 §4.8）。

---

## 6. 输出文件汇总与目录结构

```
enzyme_optimal_ph/
├── data/                       # 原始与中间数据
│   ├── metadata.csv
│   ├── DATASET.zip / Dataset_mmcifs.zip
│   ├── dataset/                # 解压的 pdb/pqr（DATASET/EnzyBase12k_pbs|pqrs）
│   └── dataset_mmcifs/         # 解压的 cif（EnzyBase12k_mmcifs）
├── output/                     # 分析报告、中间结果与备份
│   ├── metadata_analysis.json / METADATA_ANALYSIS_REPORT.md / figures/
│   ├── structure_verification.json
│   ├── sequence_structure_validation.csv (+ _summary.json)
│   ├── sequences.fa / mmseqs_cluster_95/ (+ _summary.json)
│   ├── split_summary.json / split_mapping.json
│   ├── ood_markers_summary.json
│   ├── newec_unified_summary.json / longtail_ec_unified_summary.json
│   ├── optimal_ph_ted_aug_acc_list.csv / optimal_ph_acc_cath_ted_labels.csv
│   ├── test_sequences_for_idr.fasta
│   ├── backup_before_*/        # 各次 splits 更新前备份
│   └── scripts_archive/        # 归档的一次性/旧口径脚本（见 §8）
├── pdbs/                       # 最终结构文件（11,615 个 .cif）
├── splits/                     # 最终划分（仅 train / val / test 三个 csv）
├── scripts/                    # 处理脚本（编号见文首总览）
└── logs/                       # 下载日志
```

| 输出 | 位置 | 说明 |
|:---|:---|:---|
| 最终划分 | `splits/optimal_ph_prediction_{train,val,test}.csv` | 7,794 / 866 / 2,955 行；test 34 列 |
| 结构文件 | `pdbs/*.cif` | 11,615 个，文件名 = `unique_id` |
| CATH+TED 标签表 | `output/optimal_ph_acc_cath_ted_labels.csv` | accession 级并集标签（含域区间、标注来源） |
| 各阶段摘要 | `output/*.json` | 见目录树 |

---

## 7. 诊断与质控工具

以下脚本不进主流程，重命名为 `diag_*` 保留在 `scripts/`：

| 脚本 | 用途 | 输出 |
|:---|:---|:---|
| `diag_analyze_metadata.py` | 元数据全面统计与可视化 | `output/metadata_analysis.json`、`METADATA_ANALYSIS_REPORT.md`、`figures/` |
| `diag_validate_sequence_structure.py` | 序列-结构一致性 edlib 比对（§3.1） | `output/sequence_structure_validation.csv` |
| `diag_cluster_sequences_mmseqs.py` | MMseqs2 95% 聚类冗余度评估（§3.2） | `output/mmseqs_cluster_95*/` |

---

## 8. 归档脚本清单（`output/scripts_archive/`）

| 脚本 | 归档原因 |
|:---|:---|
| `compute_ood_idr.py` | IDR@0.1 口径已融入 `04_compute_ood_markers.py`（一次算对，不再两遍） |
| `compute_additional_ood.py` | pHbin 部分已抽出为 `05_add_phbin_longtail_ood.py`；EC 旧口径已被 ec-function-ood-annotation skill 取代 |
| `fix_labels_dtype.py` | 一次性修复（labels 列字符串化）；统一 schema 后 `03` 直接产出字符串 label |
| `restrict_ood_to_default.py` | 早期将非 Default 行 OOD 列置 False 的一次性脚本；2026-08-27 起全项目统一为 NaN 惯例（§4.10） |
| `update_ood_summary.py` | 上述修复的辅助摘要脚本，随之一并归档 |

---

## 9. 注意事项（历史演进、无脚本环节与已知偏差）

1. **列 schema 演进**：`03_prepare_splits_and_pdbs.py` 原始产出为 `unique_id, file_type, aa_seq, labels`，后经全项目统一列 schema 迁移为 `unique_id, aa_seq, struct_file, label`；当前脚本已直接产出统一 schema（2026-08-28 整理时验证：对 320 行子集重跑，全部列与当前 splits 逐值一致）。
2. **非 Default 行口径演进**：早期用 `restrict_ood_to_default.py` 将非 Default 行的 `seq_Redundancy_*` / `TM-score_*` / `OOD_Orphan` 统一设为 **False**；2026-08-27 起全项目统一为 **NaN 惯例**（§4.10），该脚本已归档。2026-08-28 起 `04`/`05`/`06` 及 skill 标注脚本在流水线位置直接写入 NaN，`finalize_ood_columns.py --mode nan` 仅作兜底校验（当前数据上为 no-op）。
3. **EC 列名演进**：旧列 `OOD_NewFunction` / `OOD_LongTailFunction_L3_le5` / `OOD_LongTailFunction_L3_le10` 已被统一口径的 `OOD_NewEC_L3/L4` 与 `OOD_LongTail_EC_{L3,L4}_le{5,10}` 替换（§4.8），旧列名不再存在于 splits。
4. **IDR 阈值口径**：本任务 IDR 阈值为 **0.1**（其他任务多为 0.3）；历史上由 `compute_ood_idr.py` 第二遍重算修正，现已融入 `04_compute_ood_markers.py` 一次算对。跨任务 IDR 口径对照见 `ood-annotation-toolkit` 的 `idr_ood.py` 口径表。
5. **MMseqs2 95% 聚类未应用**：§3.2 的 557 条冗余仅作诊断，最终划分保留全部 11,615 条样本。
6. **EC 多字符串**：`metadata.csv` 的 `ec_id` 有 782 个 UniProt 为 `; ` 分隔多 EC 字符串，所有 EC 统计必须先拆分（§4.8）。
7. **EC 源漂移**：重算 EC 类 OOD 必须沿用本任务历史 EC 源（`data/metadata.csv` 的 `ec_id`），不要改用 UniProt 实时注释（见 `ec-function-ood-annotation` skill 的 EC 源对照表）。
8. **TED 标注为 putative**：CATH 实验标注优先，TED 仅补缺；约 8% accession 不在 TED 中（`not_found`），按无注释处理（见 `cath-ted-domain-annotation` skill 注意事项）。
