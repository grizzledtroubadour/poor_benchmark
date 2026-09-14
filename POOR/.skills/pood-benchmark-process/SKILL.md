# POOD-Benchmark 数据处理规范

## 描述

本 skill 汇总了 POOD-Benchmark 项目（`raw_data/`）原始数据处理的全部规范，包括：

- 各任务目录内部文件夹组织约定
- 数据集划分文件格式与输出路径
- OOD（分布外）测试集定义与命名
- `DATA_PROCESS.md` 文档写作标准
- 数据质量检查与可复现性要求

## 何时使用

在以下场景调用本 skill：

- 在任意 POOD-Benchmark 任务目录（`fold_classification/`、`ppi_prediction/`、`func_prediction/`、`ss_prediction/`、`ppis_prediction/`、`ligand_binding_site/`、`antibody_antigen/`）中处理数据
- 编写、审阅或修改 `DATA_PROCESS.md`
- 生成 `{task_dir}/splits/` 下的划分文件
- 设计新的 OOD 测试集或调整 OOD 阈值
- 检查数据泄露、分布一致性、异常值

## 1. 任务目录内部文件夹组织

每个任务目录必须按以下约定组织子文件夹：

| 文件夹 | 用途 | 内容要求 |
|:---|:---|:---|
| `data/` | 任务专属数据来源 | 存放本任务所需的原始数据、外部下载数据及任务级中间缓存；不同任务的数据相互隔离，避免跨任务混用 |
| `output/` | 中间输出文件 | 存放处理过程中生成的所有中间结果、汇总表、统计摘要、缓存文件、备份文件等 |
| `logs/` | 处理日志 | 存放正在运行或已完成的脚本日志，便于后台监控任务进度、排查问题与复现运行过程 |
| `pdbs/` | 最终数据集对应的结构文件 | 存放与 `splits/` 中样本对应的蛋白质链结构文件（PDB / mmCIF 格式），作为模型输入的结构源 |
| `splits/` | 最终数据集划分文件 | **仅包含** `*_train.csv`、`*_val.csv`、`*_test.csv` 三个最终划分文件 |
| `scripts/` | 处理脚本 | 存放完成该任务数据处理流程的所有脚本文件，脚本命名应能反映执行阶段与功能 |

> **注意**：`data/`、`output/`、`logs/` 存放的是过程数据与日志；`pdbs/` 与 `splits/` 存放的是可直接用于模型训练与评测的最终产物。

## 2. 数据集划分输出路径

所有任务的最终数据集划分文件统一输出到任务目录下的 `splits/` 文件夹：

```
{task_dir}/splits/
```

例如：

- `fold_classification/splits/`
- `func_prediction/splits/`
- `ppi_prediction/splits/`

> `func_prediction/splits/` 需同时放置 `ec_*.csv` 与 `go_{bp,cc,mf}_*.csv`。

### 2.1 文件命名

| 集合 | 文件名 |
|:---|:---|
| 训练集 | `{task_name}_train.csv` |
| 验证集 | `{task_name}_val.csv` |
| 测试集 | `{task_name}_test.csv` |

`splits/` 目录下**仅包含**上述三个最终划分文件，不得存放中间文件。

## 3. 数据集划分文件格式

### 3.1 通用列结构

所有划分文件必须包含以下核心列：

| 字段 | 类型 | 说明 |
|------|------|------|
| `unique_id` | string | 样本唯一标识。单链任务通常为 `{pdb_id}_{chain}` 或 CATH 域 ID；PPI 任务为蛋白对组合 ID |
| `aa_seq` | string | 蛋白链氨基酸序列（单链任务）或对应蛋白序列。多链蛋白可用 `|` 分隔各链序列，下游工具使用时应先去除该分隔符。PPI 等 pair-level 任务改用 `aa_seq1`/`aa_seq2` |
| `struct_file` | string | 结构文件名，**仅文件名、不含目录**（如 `7E3T_A.pdb`），文件位于本任务的 `pdbs/`；扩展名（`.pdb` / `.cif`）即结构文件类型。PPI 等 pair-level 任务改用 `struct_file1`/`struct_file2` |
| `label` | string | **强制保存为字符串**，任务相关标签。多标签用 `;` 分隔；per-residue 标签以数组字符串形式保存 |

### 3.2 按任务类型的扩展列

#### 单链单标签任务（如 fold_classification）

```
unique_id, aa_seq, struct_file, label
```

#### 单链多标签任务（如 EC、GO）

```
unique_id, aa_seq, struct_file, label
```

- `label` 为多个标签，用英文分号 `;` 分隔
- 示例：
  - EC：`2.7.6.-;3.1.7.-`
  - GO：`GO:0006164;GO:0006189`

#### Pair-level 任务（如 PPI）

| 字段 | 类型 | 说明 |
|------|------|------|
| `unique_id` | string | 蛋白对唯一标识 |
| `protein_A_id` / `protein_B_id` | string | 蛋白 A / B 的标识符（如 UniProt ID） |
| `aa_seq1` / `aa_seq2` | string | 蛋白 A / B 的氨基酸序列 |
| `struct_file1` / `struct_file2` | string | 蛋白 A / B 的结构文件名（仅文件名，位于 `pdbs/`） |
| `label` | int/string | 相互作用标签，`1` 表示有相互作用，`0` 表示无 |

#### Per-residue 任务（如 secondary structure prediction）

```
unique_id, aa_seq, struct_file, label
```

- `label` 为与 `aa_seq` **逐残基对齐**的标签数组，以字符串形式保存
- **数组序列化统一为方括号包裹、单个空格分隔的整数**：`[v1 v2 v3 ... vn]`（如 `[0 0 1 0 -1]`），**不使用逗号分隔**；未解析/缺失残基以 `-1` 掩码
- 必须保证 `len(label) == len(aa_seq)`

#### 配体相关任务（如 kcat / LBA / LBS）

在核心列之外携带配体专属特征列，`ligand_smiles`、`ligand_ecfp4`

### 3.3 OOD 标记列（仅测试集）

测试集 `*_test.csv` 必须在通用列之后追加 OOD 标记列，取值为布尔值 `True` / `False`（极端长度样本的非长度维度 OOD 列可为 `NaN`，见 3.4）。

| 列名 | 说明 |
|------|------|
| `Default` | 没有通过 OOD-ExtremeLength 提前去除的常规测试集样本 |
| `seq_Redundancy_90` ~ `seq_Redundancy_30` | 序列同源性 OOD，步长 10% |
| `TM-score_0.9` ~ `TM-score_0.3` | 结构相似性 OOD，步长 0.1 |
| `OOD_Orphan` | 在 Train+Val 中无显著同源 hit |
| `OOD_ExtremeShort` / `OOD_ExtremeLong` | 极端长度 OOD（单链 < 60 / > 1000） |
| `OOD_IDR` | 高内在无序区域比例 |
| `OOD_FoldHoldout` | 折叠类型留一法（CATH 拓扑） |
| `OOD_SuperfamilyHoldout` | 超家族留一法（CATH 超家族） |
| `OOD_NewEC_L4` / `OOD_NewEC_L3` / `OOD_NewEC_L1` | 新 EC 功能类别，后缀表示 EC 层级：`L4`=4 级完整 EC 编号、`L3`=3 级、`L1`=一级大类 |
| `OOD_LongTail` | 训练集中低频出现的**预测标注**类别（该任务自身要预测的标签，如 CATH 架构、GO/EC 标签、二级结构、界面/结合位点比例等） |
| `OOD_LongTail_EC_L4` / `OOD_LongTail_EC_L3` | 训练集中低频出现的 **EC 功能**类别（用作非 EC 预测任务的 OOD 维度），后缀表示 EC 层级 |
| `OOD_Combinatorial` | 组合泛化 OOD（标签组合未在 Train+Val 出现） |
| `InD` | 纯 in-distribution 标识：所有 `OOD_*` + `seq_Redundancy_*` + `TM-score_*` 标记均为 False（NaN 视为 False）时为 True；旧名 `Pure_ID`，统一由 `.skills/ood-annotation-toolkit/scripts/finalize_ood_columns.py --mode ind` 在所有 OOD 列齐备后重算 |

> **OOD 标注 skills**：各 OOD 列的统一口径与通用脚本分别固化在四个 skill 中——`homology-ood-annotation`（seq_Redundancy_*/OOD_Orphan/TM-score_*）、`cath-ted-domain-annotation`（OOD_FoldHoldout/OOD_SuperfamilyHoldout）、`ec-function-ood-annotation`（OOD_NewEC_*/OOD_LongTail_EC_*）、`ood-annotation-toolkit`（OOD_IDR、值域 bin LongTail、Default=False→NaN 收尾、InD 重算、splits 体检）。新增或重算 OOD 列时必须调用对应 skill 的脚本，不得另起口径。

> **阈值变体**：个别任务对同一维度保留多个阈值列时以阈值后缀区分，例如 enzyme_optimal_ph 的 `OOD_LongTail_EC_L3_le5` / `OOD_LongTail_EC_L3_le10`（EC L3 频次 ≤5 / ≤10）与 `OOD_LongTail_pHbin_le50` / `OOD_LongTail_pHbin_le100`（pH 区间样本数 ≤50 / ≤100，属预测标签长尾）。

### 3.4 `Default` 列定义与处理逻辑

`Default` 列用于标识**常规测试集样本**（即未通过极端长度筛选提前移除的样本），是 OOD 评估的基准子集。

#### 定义规则（按任务类型区分）

| 任务类型 | `Default=True` 定义 | `Default=False` 定义 | 说明 |
|:---|:---|:---|:---|
| **LBA（ligand_binding_affinity）** | 所有测试集样本 | 无 | 训练集/验证集保留极端长度样本，测试集不做极端长度剔除，所有样本均为 Default=True |
| **PPI（ppi_prediction）** | 双链长度均 ≤1000 的常规样本 | 任一链 > 1000 的样本 | 划分前将含 >1000aa 单链的样本预挑出直入 test 作 `OOD_ExtremeLong`（`Default=False`）；train/val 单链 ≤1000 |
| **其他所有任务** | 长度在 [60, 1000] 范围内的正常长度样本 | 长度 < 60 或 > 1000 的极端长度样本 | 划分前或划分后已将极端长度样本从 train/val 移至 test 或单独作为 OOD 子集 |

#### 极端长度样本的处理要求（适用于所有含 `Default=False` 行的任务；LBA 无此类行）

当 `Default=False`（即样本为极端长度预挑出）时，必须满足以下约束：

1. **其他 OOD 标记列设为 `NaN`**：极端长度样本不参与序列同源性、结构相似性、IDR 等其他 OOD 维度的计算，除其极端长度原因列（`OOD_ExtremeShort` / `OOD_ExtremeLong`，即触发预挑出的那一列）外，对应列统一设为 `NaN`，与 `False` 区分。
2. **不参与 ID-Test 统计**：`Default=False` 的样本不计入 in-distribution 测试集。
3. **单独评估**：极端长度样本作为独立的 `OOD_ExtremeLength` / `OOD_ExtremeLong` / `OOD_ExtremeShort` 子集进行专门评估。

> **注意**：PPI 任务对含 >1000aa 单链的样本做"预挑出直入 test"处理（`Default=False`、`OOD_ExtremeLong=True`，1,325 行），其 `seq_Redundancy_*` / `TM-score_*` / `OOD_ExtremeShort` 等其余 OOD 列均为 NaN（2026-08-27 起与其他任务统一，此前曾对全行计算 seq/TM）；LBA 任务则全表 `Default=True`。

#### 超长蛋白硬性剔除（单链长度 > 2000）

除上述"标记/移至 test"的极端长度处理外，**单链氨基酸序列长度 > 2000 的样本一律直接从数据集中剔除**（train / val / test 全部移除，既不保留也不作为 OOD 子集）。理由：超长序列/结构会带来显存与处理成本问题，且属于分布极端尾部，不适合纳入常规评测。

- 判定单位：单条链的氨基酸序列长度（对 PPI 等 pair-level 任务，**任意一条链** > 2000 即剔除整个样本）。
- 适用范围：所有任务。

### 3.5 关键约束

1. **`label` 必须保存为字符串类型**：避免 `1.10`、`3.2.1.-` 等被解析为浮点数导致精度丢失。
2. **测试集 OOD 列必须完整**：每个测试样本的所有 OOD 标记都必须有明确值，不允许缺失。
3. **ID-Test 可识别**：未标记为任何 OOD 的测试样本即为 in-distribution 测试集（ID-Test）。
4. **Train / Val / Test 之间 `unique_id` 不得重复**。

## 4. 标准 OOD 测试集定义

| OOD 类型 | 定义 | 工具/方法 | 阈值/规则 | 科学意义 |
|:---|:---|:---|:---|:---|
| **OOD-LowHomology** | 测试样本与 train+val 的最大序列 identity 低于阈值 | mmseqs2 search | **< 90**、< 80、< 70、< 60、< 50、< 40、< 30 | 评估序列相似度降低时的泛化能力 |
| **OOD-Orphan** | 测试样本在 train+val 中无显著同源 hit | mmseqs2 search | max identity = 0（或 e-value > 1e-3） | 评估完全无同源训练样本的泛化能力 |
| **OOD-FamilyHoldout** | 按蛋白质家族留一法划分 | CATH / SCOPe / Pfam / ECOD | 留一法 | 评估对未见蛋白家族的泛化能力 |
| **OOD-SuperfamilyHoldout** | 按蛋白质超家族留一法划分 | 同上 | 留一法 | 评估远同源识别能力 |
| **OOD-NewFold** | 测试样本与 train+val 的最大结构相似度低于阈值 | Foldseek | **< 0.9**、< 0.8、< 0.7、< 0.6、< 0.5、< 0.4、< 0.3 | 评估结构折叠空间上的泛化能力 |
| **OOD-FoldHoldout** | 按折叠类型留一法划分 | CATH / SCOPe / ECOD | 留一法 | 评估对未见结构折叠的泛化能力 |
| **OOD-IDR** | 测试样本含有显著的内在无序区域 | metapredict v3 / flDPnn / IUPred | 连续 disorder propensity > 0.5 区域占序列长度 > 30% | 评估对富含无序区域蛋白的预测能力 |
| **OOD-ExtremeLength** | 长度处于分布两端的蛋白 | 长度统计 | length < 60 或 length > 1000 | 评估对极端长度样本的鲁棒性 |
| **OOD-NewFunction** | 按功能类别留一法划分 | EC / GO / InterPro | 留一法 | 评估对未见功能类别的泛化能力 |
| **OOD-LongTailFunction** | 测试样本带有 train 中低频的**功能类别**标签 | EC / GO / InterPro 标签频次统计 | EC 类别出现次数 **< 10**（可随任务调整） | 评估长尾功能类别预测能力 |
| **OOD-LongTail** | 测试样本带有 train 中低频的**预测标注**类别 | 标签频次统计 | 某标注类别样本数 **< 10** | 评估数据集预测标注的长尾泛化能力 |
| **OOD-CombinatorialFunction** | 测试样本标签组合未在 train+val 中出现 | 多标签组合统计 | 标签对 `(A, B)` 未在 train+val 中同时出现 | 评估组合泛化能力 |

## 5. DATA_PROCESS.md 写作标准

### 5.1 文档顶层结构

```markdown
# 标题（任务名称 + "数据集处理文档"）
> 简述文档目的与处理位置

---

## Stage N: 阶段名称（阶段目标简述）
> 阶段目标声明

### N.M 子步骤名称
#### N.M.K 更细粒度操作（如需要）

---

### X. 输出文件汇总
### X+1. 各维度统计汇总
### X+2. 注意事项
```

### 5.2 每个 Stage 必须包含

1. **项目/步骤概述**：核心任务和科学依据
2. **核心数据规模**：首次出现时必填，关键数字加粗
3. **数据来源**：输入文件表格（路径、用途、说明）
4. **处理步骤**：每个步骤包含目的、方法、脚本路径、结果
5. **验证与质控**：验证方法、指标分布、异常解释

### 5.3 表格使用规范

必须使用表格的场景：
- 数据规模统计
- 文件清单
- 筛选条件
- 对比结果
- OOD 划分结果

### 5.4 脚本与工具记录

每次调用脚本必须记录：
- 路径：`scripts/script_name.py`
- 用途：一句话说明
- 关键参数：并发数、阈值、模式
- 结果：成功/失败数量，耗时

使用外部工具时必须记录：
- 工具名称及版本
- 核心参数
- 比对方向
- 判定指标

### 5.5 数据质量检查

必须包含独立的**数据质量检查**章节：
1. **泄露检查**：Train ∩ Val、Train ∩ Test、Val ∩ Test 是否为空
2. **分布一致性**：长度、分辨率、类别分布等关键特征对比
3. **异常值处理**：列出被移除的异常样本及原因

### 5.6 输出文件记录

每个 Stage 末尾列出输出文件：
- 位置
- 文件名
- 行数/大小
- 格式
- 说明

## 6. 文本与数字风格

| 规范项 | 要求 |
|:---|:---|
| 语言 | 与项目主体语言保持一致 |
| 数字格式 | 千分位逗号分隔；百分比保留 1 位小数 |
| 加粗强调 | 关键结论、核心数字、推荐阈值、重要警告使用加粗 |
| 代码/路径 | 文件路径、脚本名、列名、字段值用反引号包裹 |
| 科学名称 | mmseqs2、Foldseek、PDB、UniProt 等保持原样 |

## 7. 快速检查清单

在 `DATA_PROCESS.md` 完成前，逐项确认：

- [ ] 文档头部包含处理位置说明
- [ ] 各 Stage 开头有目标声明
- [ ] 所有数据规模变化均有表格记录（输入 → 输出）
- [ ] 使用的脚本均已标注路径和参数
- [ ] 外部工具已记录版本和核心参数
- [ ] OOD 划分包含定义、工具、阈值、规模、科学意义
- [ ] 数据泄露检查已执行并记录结果
- [ ] 异常值/过滤样本已列出
- [ ] 输出文件清单完整（含路径、格式、大小）
- [ ] 关键字段/列已说明类型和含义
- [ ] 数据集划分文件符合 `{task_dir}/splits/` 格式规范
- [ ] `splits/` 目录下仅包含 train / val / test 三个最终划分文件
- [ ] 任务目录下 `data/` / `output/` / `logs/` / `pdbs/` / `splits/` / `scripts/` 组织符合规范
- [ ] 训练集/验证集不含 OOD 列，测试集包含完整 OOD 标记
- [ ] `Default` 列定义符合任务类型（LBA 全 True；PPI 含 >1000aa 链预挑出为 False；其他任务按极端长度区分）
- [ ] 极端长度样本（`Default=False`）的其他 OOD 列已设为 `NaN` 或 `False`
- [ ] `label` 列已强制保存为字符串类型
- [ ] 注意事项涵盖已知问题和特殊处理
- [ ] 设计意图对非显然策略做了说明
