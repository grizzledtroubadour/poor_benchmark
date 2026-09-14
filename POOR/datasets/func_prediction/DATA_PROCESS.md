# func_prediction 数据集处理文档（EC / GO-MF / GO-CC / GO-BP）

> 本文档记录 `datasets/func_prediction/` 从原始数据一次性复现的完整线性流程。
> 任务包含 4 个子任务：**EC**（`ec`）与 **GO** 三分支（`go_mf` / `go_cc` / `go_bp`），
> 最终产物为 `splits/` 下 12 个划分文件（每子任务 train / val / test 各一）。

## 任务定义与当前规模

从 PDB 实验结构出发，经 SIFTS 映射至 UniProt，为每个 UniProt 保留一条最优实验结构链，按子任务组织多标签分类数据集（EC 为 L3 级 `x.x.x.-` 标签；GO 为 term 级标签），按 deposition date 时间切分并追加多维 OOD 标注。

当前 `splits/` 规模（实测）：

| 子任务 | train | val | test | 其中 Default | test 列数 |
|:---|---:|---:|---:|---:|---:|
| **EC** | 15,389 | 1,687 | 5,047 | 4,243 | 31 |
| **GO-MF** | 18,179 | 2,084 | 5,609 | 4,609 | 34 |
| **GO-CC** | 13,878 | 1,730 | 5,613 | 4,256 | 34 |
| **GO-BP** | 15,232 | 1,742 | 5,142 | 4,121 | 34 |

处理流程分四个阶段（另有一个可选扩展阶段）：

- **阶段一**：源数据获取与解析（SIFTS → PDB 实验结构 → UniProt 映射 → 结构元数据与功能注释）
- **阶段二**：去冗余与数据精化（链级结构提取 → UniProt 级代表链筛选 → 长度过滤 → 序列聚类去冗余）
- **阶段三**：标签聚合与数据划分（按子任务划分 → 时间切分 → 测试集过滤 → 标签清洗）
- **阶段四**：OOD 标注（极端长度链并入 test → 多维度 OOD 标注 → InD/NaN 收尾）
- **阶段五（可选）**：AlphaFold DB 预测结构收集（仅用于训练集扩展）

---

## 阶段一：源数据获取与解析（SIFTS 结构-序列-功能映射）

> 目标：从 PDB 出发，通过 SIFTS 建立 **PDB chain 级别**与 UniProt 的精确映射，获取结构元数据及初步功能标注。这是数据集构建的**实验结构主入口**，可捕获 Swiss-Prot cross-reference 未覆盖的 PDB chains（包括 TrEMBL 来源及新提交结构）。

**脚本**：`scripts/01_parse_sifts.py`（SIFTS 解析、元数据批量查询、功能注释导出）。

### 1.1 SIFTS 核心映射下载（无脚本环节）

- **数据源**：EBI SIFTS (https://www.ebi.ac.uk/pdbe/docs/sifts/)
- **核心文件**（Flat File TSV，存放于项目级 `POOR/data/sifts/`）：
  - `sifts_chain_uniprot.tsv.gz` / `pdb_chain_uniprot.tsv.gz`：PDB chain ID ↔ UniProt Accession 的**残基级映射**
  - `sifts_chain_go.tsv.gz`：PDB chain ID ↔ GO ID 及证据代码（来自 UniProt-GOA）
  - `sifts_chain_ec.tsv.gz`：PDB chain ID ↔ EC 编号
  - `pdb_chain_cath.tsv.gz`：结构分类（供阶段四 4.6 OOD_FoldHoldout 使用）

### 1.2 PDB 链级元数据提取

- **数据源**：RCSB PDB GraphQL API (`https://data.rcsb.org/graphql`)
- **获取方式**：批量查询（batch size = 50 entries / request），断点续传（每 500 entries 保存 checkpoint）
- **GraphQL 查询字段**：`exptl.method`、`refine.ls_d_res_high`（X-ray 分辨率）、`refine.ls_R_factor_R_work` / `refine.ls_R_factor_R_free`、`em_3d_reconstruction.resolution`（EM 分辨率）
- **Resolution 取值逻辑**：优先取 `refine.ls_d_res_high`（X-ray），缺失时回退 `em_3d_reconstruction.resolution`（EM）
- **元数据覆盖度**（238,401 entries）：Resolution 94.9%，R-work 79.7%，R-free 78.6%

### 1.3 UniProt 参考序列获取

- **目的**：获取每条 PDB chain 映射到的 UniProt 全长序列，用于结构-序列比对验证与覆盖度计算
- **数据源**：UniProt REST API（按 100 accessions / request 批量查询 FASTA）
- **覆盖度**：~587,000 条序列，覆盖率 100%（产物在 `output/uniprot/` 与 `output/uniprot_reference_sequences.json`）

### 1.4 功能注释获取（通过 SIFTS 映射）

- **GO 标注**：直接采用 SIFTS 输出的 GO ID 及证据代码；SIFTS 的 GO 映射覆盖**整个 UniProtKB**（Swiss-Prot + TrEMBL）
- **EC 标注**：通过 SIFTS EC 映射获取
- **输出**：`data/sifts/sifts_go_annotations.json`（822M）、`data/sifts/sifts_ec_annotations.json`（43M）

### 1.5 初始筛选（长度 + 功能注释）

对原始约 95.8 万条 chains 执行初始筛选。本阶段**放宽结构质量限制**（暂不限制分辨率、方法、R-factor），仅保留两个核心条件：

| 筛选条件 | 阈值 | 说明 |
|:---|:---|:---|
| SIFTS 映射长度 | `60 <= mapped_length <= 1000` | 排除过短片段与超大蛋白 |
| 功能注释 | 必须有 EC **或** GO 标注 | 保证下游任务有监督信号 |

- **筛选结果**（以当前 `data/sifts/annotated_chains.csv` 重算）：**769,892** 条 chains、**56,059** 个 UniProt（历史执行为 769,730 / 56,057，差异为 SIFTS 再生漂移，见"注意事项"）
- **方法分布（历史执行）**：X-RAY DIFFRACTION 442,283；ELECTRON MICROSCOPY 319,654；SOLUTION NMR 6,670；SOLID-STATE NMR 650；ELECTRON CRYSTALLOGRAPHY 301；其他 172

> **设计意图**：过早限制分辨率/方法会大量过滤 EM 结构，而 EM 结构对大型复合物表征至关重要。质量指标在后续精筛阶段再应用。

### 1.6 初筛后结构文件下载

- **下载源**：RCSB PDB (`https://files.rcsb.org/download/`)，优先 `{pdb_id}.pdb.gz`，404 时回退 `{pdb_id}.cif.gz`
- **并发与容错**：10 线程；每文件最多重试 5 次，429 限流指数退避；checkpoint 每 100 个文件保存
- **存储位置**：`data/pdb/`（项目级共享目录，gzip 压缩，~255K 个文件）

### 1.7 输出文件

| 文件 | 说明 |
| :--- | :--- |
| `data/sifts/sifts_pdb_chains.csv` | 所有候选 PDB chains 及结构元数据（41M） |
| `data/sifts/sifts_uniprot_mapping.json` | PDB chain → UniProt 残基级映射（154M） |
| `data/sifts/sifts_go_annotations.json` | 链级 GO 标注（822M） |
| `data/sifts/sifts_ec_annotations.json` | 链级 EC 标注（43M） |
| `data/sifts/annotated_chains.csv` | 初筛候选链元数据（827,186 行；按 1.5 条件过滤即得候选池） |
| `output/uniprot/` | UniProt 参考序列与 Swiss-Prot 元数据 |

---

## 阶段二：去冗余与数据精化（PDB-UniProt 对齐）

> 目标：执行**链级结构提取**与 **UniProt 级去重**，为每个 UniProt 保留**一条最优 PDB 实验结构链**，并在序列层面聚类去冗余。阶段一保留"全量候选"（含同源多聚体全部拷贝），阶段二去粗取精。

### 2.1 UniProt 级代表链筛选

**脚本**：`scripts/02_select_representative_chains.py`（补缺口重实现；历史脚本与其产出 `representative_chains.csv` 未从旧工作目录迁移，本脚本按下列优先级表重实现）。

**策略**：按 UniProt ID 分组，每组仅保留**一条最优 PDB chain**。候选链按以下优先级逐层比较：

| 优先级 | 条件 | 排序方向 |
|:---|:---|:---|
| 1 | 提取模式 | A 档 > B 档 |
| 2 | 来源质量 | Swiss-Prot > TrEMBL（组内恒同，实际无区分力） |
| 3 | EC 优先 | 有 EC > 无 EC |
| 4 | 实验 GO 丰富度 | 降序 |
| 5 | GO 丰富度 | 降序 |
| 6 | 分辨率 | 升序 |
| 7 | R-free | 升序 |
| 8 | 映射长度 | 降序 |
| 9 | PDB ID + chain ID 字母序 | 升序（确定性 tiebreak） |

- **候选池构建**：`data/sifts/annotated_chains.csv`（1.5 初筛条件）+ `sifts_uniprot_mapping.json`（A/B 档判定：SIFTS 段存在 pdb_beg/pdb_end 均有效者 → A 档）+ `sifts_go_annotations.json`（实验 GO 数，ijson 流式）+ `output/uniprot_metadata.json`（reviewed）
- **输出**：`output/representative_chains.csv`（历史执行 **55,124** 条：A 档 31,770 / B 档 23,354；Swiss-Prot 23,374 / TrEMBL 31,750；有 EC 24,369；有实验 GO 14,100）
- **已知口径偏差**：历史执行在候选入池前对 B 档链做边界推断 identity ≥ 0.90 过滤（439,427 尝试 → 419,370 通过），其产物 `b_boundary_inference.jsonl` 未随仓库迁移；本脚本默认 B 档全量入池（可用 `--b-keep-list` 传入推断通过清单复现历史过滤）。2026-08-28 验证重跑（无 keep-list）：候选池 769,892 链 / 56,059 UniProt，选出 **56,057** 条代表链；与最终非冗余代表链比对，**38,835 / 47,852（81.2%）精确一致**，9,017 个 UniProt 选出了不同链（主因即 B 档过滤缺失与 A/B 档分类口径漂移，详见"注意事项"）。

**GO/EC 功能注释精筛（可选步骤，本次跳过）**：原计划对 TrEMBL 链要求"含 EC 或至少一个实验验证 GO 证据（EXP/IDA/IPI/IMP/IGI/IEP/TAS/IC/HTP/HDA）"，本次执行**跳过**，全部通过 1.5 初筛的 chains 参与排序。因跳过精筛，42.2%（23,258 个）代表链无实验验证功能注释；若执行精筛，总代表链约降至 39,000。

### 2.2 链级结构提取（三档策略）

**脚本**：`scripts/03_extract_representative_structures.py`
（从 `data/pdb/` 的原始 `.pdb.gz` / `.cif.gz` 中按链截取，输出到 `pdbs/`；脚本内路径为运行目录相对路径，需在含 `data/pdb/` 与 `pdbs/` 的工作根下执行）。

| 档位 | 条件 | 提取方式 | 优先级 |
|:---|:---|:---|:---|
| **A 档** | SIFTS segments 的 `pdb_beg` / `pdb_end` 有效 | 按官方残基范围精确截取 | 最高 |
| **B 档** | `pdb_beg` 或 `pdb_end` 为空，但边界推断成功 | 用 `sp_beg`/`sp_end` 子序列比对推断 PDB 边界后截取 | 中 |
| **C 档** | B 档推断失败或 identity 过低（<90%） | 回退提取整条链（或丢弃） | 最低 |

- **输入**：`output/representative_chains.csv`（2.1 产物）、`data/sifts/sifts_uniprot_mapping.json`、`output/b_boundary_inference.jsonl`（B 档边界推断中间产物，**未随仓库迁移**，重跑需先按下方"边界推断"方法重新生成）
- **输出路径**：`pdbs/{pdb_id}_{chain_id}.pdb`（单字符 chain ID）或 `pdbs/{pdb_id}_{chain_id}.cif`（多字符 chain ID，PDB 格式列宽无法容纳）
- **格式支持**：`.pdb.gz` 流式解析；`.cif.gz` 经 BioPython `MMCIFParser` 解析后保持 mmCIF 输出；按 `pdb_id` 分 batch，同一 cif.gz 只解析一次；多模型（NMR）条目截断于首个 `ENDMDL`（只保留 model 1，根因修复已内置）
- **提取结果**：成功提取 **55,124 个结构文件**（100% 覆盖代表链）：
  - A 档（完整边界精确截取）20,828 条（37.8%）+ A 档（部分边界推断）10,942 条（19.9%）
  - B 档（子序列比对推断）23,354 条（42.4%）；特殊链（插入码/递减范围）12 条手工修复；C 档 0 条
  - PDB 格式 47,621 个，mmCIF 格式 7,503 个

**B 档边界推断（序列比对恢复法）**：截取 `UniProt[sp_beg:sp_end]` 子序列与 PDB 链全长做局部比对（edlib HW 模式），映射回 PDB 残基号得到推断边界。以 identity ≥ 0.90 过滤后保留 **419,370 条** B 档候选链（剔除 19,906 条低质量 + 151 条无推断结果）。

> **历史说明**：当次执行使用全长 UniProt 序列比对（旧方法），identity/coverage 系统性偏高（identity ≥ 0.90 占 95.5%）；修正为子序列比对后分布会更严格。见"注意事项"。

### 2.3 长度过滤与序列聚类去冗余

#### 2.3.1 实际残基数统计与长度过滤（极端链单独保存）

**脚本**：`scripts/04_extract_seqs_and_count_residues.py`（从 `pdbs/` 提取序列并统计去重残基数，标准 20 种氨基酸 + MSE/SEC/PYL；产物写入 `output/`）

- **过滤规则**：正常长度区间为 **[60, 1000]** 残基；**长度 < 60 或 > 1000 的极端链不丢弃，而是单独保存**，供阶段四 4.1 并入测试集（由 `scripts/08_recover_extreme_length.py` 从初筛候选中收集，产出 `output/extreme_length_chains.csv`，**3,834 条**：ExtremeShort 3,146 + ExtremeLong 688，长度范围 3–1,997）；长度 **> 2000 的链全局剔除**（显存与成本考虑）
- **结果**：正常长度链 **53,902** 条（55,124 − 1,222 条出界，占 2.2%；A 档短链率 1.2%，B 档 3.6%）

#### 2.3.2 95% 序列同一性聚类（MMseqs2，无脚本环节）

- **工具**：MMseqs2 `mmseqs cluster`，参数 `--min-seq-id 0.95 -c 0.8 --cov-mode 0`
- **输入**：`output/filtered_sequences.fa`（53,902 条正常长度序列）
- **结果**：**47,852** 个 clusters（singleton 44,052 个，92.1%；最大 cluster 28 条序列）

#### 2.3.3 Cluster 内代表链筛选（无独立脚本环节）

对每个 cluster 应用与 2.1 相同的优先级排序保留最优一条，得 **47,852 条非冗余代表链**：

- A 档 28,266（59.1%）/ B 档 19,586（40.9%）；Swiss-Prot 22,515（47.1%）/ TrEMBL 25,337（52.9%）
- 方法分布：X-RAY 33,977（71.0%）、EM 11,868（24.8%）、NMR 1,978（4.2%）、其他 29（0.1%）
- 平均分辨率 2.42 Å；平均映射长度 300.5 残基；聚类缩减率 **11.2%**

#### 2.3.4 输出文件

| 文件 | 说明 |
| :--- | :--- |
| `output/filtered_sequences.fa` | 长度过滤后 53,902 条序列（~11 MB） |
| `output/chain_residue_counts.json` | 每链残基数统计（04 产物，重跑时再生） |
| `output/final_nonredundant_representatives.csv` | **最终非冗余代表链 47,852 条**（含 `extract_mode`、`cluster_rep`、`cluster_size` 等列） |
| `output/extreme_length_chains.csv` | 单独保存的极端长度链 3,834 条（含序列、deposition date、各任务标签） |
| `pdbs/` | 链级结构文件（当前 50,123 个 = 47,852 代表链 + 2,271 个并入 test 的极端链；`.pdb` 43,567 + `.cif` 6,556） |

---

## 阶段三：标签聚合与数据划分

> 目标：把 47,852 条代表链按子任务（EC / GO-MF / GO-CC / GO-BP）组织成数据集，按 **PDB deposition date** 时间切分 train/val/test，再清洗标签。中间产物在 `output/datasets_split_raw/`。

### 3.1 任务划分

**脚本**：`scripts/06_split_datasets.py`

| 数据集 | 分配规则 | 链数 |
|:---|:---|---:|
| **EC** | `ec_count > 0` | 22,095 |
| **GO-MF** | 至少 1 个 Molecular Function GO 标注 | 39,646 |
| **GO-CC** | 至少 1 个 Cellular Component GO 标注 | 27,972 |
| **GO-BP** | 至少 1 个 Biological Process GO 标注 | 34,607 |

> 同一条链可同时属于多个数据集，各数据集独立划分。GO namespace 由 `go-basic.obo` 解析（需从 Gene Ontology 下载，见"注意事项"）。

### 3.2 时间切分

**脚本**：`scripts/05_parse_deposition_dates.py`（从 `data/pdb/` 原始文件 HEADER / `_pdbx_database_status.recvd_initial_deposition_date` 解析日期，输出 `output/deposition_dates.json`）+ `scripts/06_split_datasets.py`（按日期升序分位数切分）。

| 数据集 | 训练:验证:测试 | 训练集截止 | 验证集区间 | 测试集起始 |
|:---|:---|:---|:---|:---|
| **EC** | 72 : 8 : 20 | ≤ 2019-06-25 | 2019-06-25 ~ 2021-04-29 | ≥ 2021-04-29 |
| **GO-MF** | 63 : 7 : 30 | ≤ 2019-05-19 | 2019-05-20 ~ 2020-10-07 | ≥ 2020-10-07 |
| **GO-CC** | 63 : 7 : 30 | ≤ 2020-10-14 | 2020-10-14 ~ 2021-11-22 | ≥ 2021-11-22 |
| **GO-BP** | 63 : 7 : 30 | ≤ 2020-01-15 | 2020-01-15 ~ 2021-04-20 | ≥ 2021-04-21 |

> 各数据集独立切分以保证内部比例精确，故 cutoff 不一致；时间 gap 保证测试集不泄漏到训练/验证。

### 3.3 测试集过滤

| 数据集 | 过滤规则 | 过滤前 | 过滤后 | 保留率 |
|:---|:---|---:|---:|---:|
| **EC** | **不过滤**（EC 分类可靠性不依赖 GO 式证据代码） | 4,419 | **4,419** | 100% |
| **GO-MF** | reviewed **OR** 任意实验 GO 证据 | 11,894 | **6,236** | 52.4% |
| **GO-CC** | 同上 | 8,392 | **4,692** | 55.9% |
| **GO-BP** | 同上 | 10,383 | **5,653** | 54.4% |

实验 GO 证据代码（10 种）：EXP, IDA, IPI, IMP, IGI, IEP, TAS, IC, HTP, HDA。GO 的 IEA 电子注释（约占 51.7%）噪声较高，评估时需过滤；EC 无独立证据代码体系，即使 TrEMBL 来源也有较高可靠性。

### 3.4 标签清洗

**脚本**：`scripts/07_clean_datasets.py`（输入/输出均为 `output/datasets_split_raw/`，清洗结果覆盖写回划分文件，另存 `{task}_labels.csv`、`*_meta.csv` 与 `cleaning_stats.json`）。

#### 3.4.1 EC 标签清洗（L4 → L3 归并）

1. **级别过滤**：只保留 3 级（`x.x.x.-`）和 4 级（`x.x.x.y`）EC 编号；
2. **全局归并**：所有 4 级编号统一提升为其 3 级前缀（`1.1.1.1` → `1.1.1.-`），避免多标签评估中的父子误判；
3. **训练集存在性过滤**：删除训练集中从未出现的 EC 编号；链上无剩余标签则删链。

| 数据集 | 清洗前 | 清洗后 | 训练集标签数 |
|:---|---:|---:|---:|
| EC train | 15,908 | **15,389** | **260**（全为 3 级标签） |
| EC val | 1,768 | **1,687** | 175 |
| EC test | 4,419 | **4,243** | 211 |

> 94.1% 样本为单标签，5.9% 为多标签。

#### 3.4.2 GO 标签清洗（频数过滤）

按 namespace 分别处理，只保留训练集中出现次数在 **[5, 5000]** 的 GO term（下限 5 过滤噪声低频 term，上限 5000 过滤过宽泛 term）；链上无剩余 term 则删链。

| 任务 | 清洗后 train | 清洗后 val | 清洗后 test | 标签数 |
|:---|---:|---:|---:|---:|
| **GO-MF** | **18,179** | 2,084 | 4,609 | 447 |
| **GO-CC** | **13,878** | 1,730 | 4,256 | 264 |
| **GO-BP** | **15,232** | 1,742 | 4,121 | 884 |

> GO-BP 删除比例最高，因 BP term 极度稀疏（67.6% 的 term 训练集仅出现 1–10 次）。

### 3.5 输出文件

| 文件 | 说明 |
|:---|:---|
| `output/datasets_split_raw/{task}_{train,val,test}.csv` | 清洗后标签（`chain_key, labels, num_labels`） |
| `output/datasets_split_raw/{task}_{train,val,test}_meta.csv` | 对应样本完整元数据 |
| `output/datasets_split_raw/split_summary.csv` | 划分统计摘要（含日期 cutoff） |
| `output/datasets_split_raw/cleaning_stats.json` | 清洗统计报告 |

> 阶段三到阶段四之间有一步**无独立脚本的 schema 落盘**：清洗后的划分经全项目统一 schema（`unique_id, aa_seq, struct_file, label`）写入 `splits/`（早期文档中的 `datasets_split/` 即现在的 `splits/`）。

---

## 阶段四：OOD 标注

> 目标：在清洗后的 train/val/test 基础上，对**测试集**追加多维度 OOD 标记列。训练/验证集保持完整不动；每个 OOD 维度一个布尔列（极端长度行除长度维度外为 NaN）。
>
> **统一口径**：各 OOD 维度的计算方法固化在四个 skill 中，重算时必须调用对应脚本，不得另起口径：
> - `.skills/homology-ood-annotation`：`seq_Redundancy_*` / `OOD_Orphan` / `TM-score_*`（mmseqs2 / Foldseek，支持 `--m8-cache` 重放）
> - `.skills/cath-ted-domain-annotation`：`OOD_FoldHoldout` / `OOD_SuperfamilyHoldout`（CATH+TED 域标注，经 `scripts/13_ted_augment_cath_ood.py` 调用）
> - `.skills/ec-function-ood-annotation`：`OOD_NewEC_*` / `OOD_LongTail_EC_*`（all-not-in 统一口径）
> - `.skills/ood-annotation-toolkit`：`OOD_IDR`、Default=False→NaN 兜底校验（NaN 由各生成脚本在流水线位置直接写入）、InD 重算、splits 体检

### 4.1 极端长度链并入测试集

阶段二 2.3.1 单独保存的极端长度链在此并入各子任务 test：

1. **候选收集**（`scripts/08_recover_extreme_length.py`）：从初筛候选中取 UniProt 不在 47,852 代表链中的链（每 UniProt 保留最优：has_ec → resolution → mapped_length → ID 字母序），从 `data/pdb/` 提取全链序列并解析 deposition date；`<60` → ExtremeShort，`1000 < len ≤ 2000` → ExtremeLong，`>2000` 丢弃。产出 `output/extreme_length_chains.csv`（3,834 条）。
2. **并入 test**（`scripts/09_integrate_extreme_to_splits.py`）：按子任务清洗标签（EC L4→L3 归并 + 须在 train 出现；GO 频数 [5,5000] + 须在 train 出现；GO test 额外要求 reviewed 或实验 GO 证据，与 3.3 一致），追加到 `splits/{task}_test.csv`，标记 `Default=False`、`OOD_ExtremeShort` / `OOD_ExtremeLong=True`，其余 OOD 列留待统一收尾。写前自动备份到 `output/splits_pre_extreme_backup/`。
3. **结构提取**（`scripts/10_extract_extreme_pdbs.py` / `scripts/11_extract_extreme_cifs.py`）：为并入的极端链从原始文件提取链级结构存入 `pdbs/`（全链提取；`.pdb` 文本过滤或 `.cif` 按 `auth_asym_id` 过滤；多模型条目同样截断于首个 `ENDMDL`）。

| 任务 | 原有 test | 新增 ExtremeShort | 新增 ExtremeLong | 更新后 test |
|:---|---:|---:|---:|---:|
| **EC** | 4,243 | 478 | 326 | **5,047** |
| **GO-MF** | 4,609 | 753 | 247 | **5,609** |
| **GO-CC** | 4,256 | 1,090 | 267 | **5,613** |
| **GO-BP** | 4,121 | 757 | 264 | **5,142** |

> **`Default` 与 NaN 惯例**：`Default=True` 表示长度在 [60, 1000] 的常规测试样本；极端长度行 `Default=False`，其除 `OOD_ExtremeShort` / `OOD_ExtremeLong` 外的所有 `OOD_*` 列**不参与计算**，由 `.skills/ood-annotation-toolkit/scripts/finalize_ood_columns.py --mode nan` 统一置为 **NaN**（空单元格，与"已计算、非 OOD"的 False 区分）。当前 4 个 test csv 均已按此处理（备份 `output/backup_before_nondefault_nan/`）。

### 4.2 OOD-TimeCutoff

数据集本身已按 deposition date 时间切分（阶段三 3.2），全部测试样本天然属于 time-cutoff OOD，作为基线参照，不单独设列。

### 4.3 序列同源 OOD（seq_Redundancy_* + OOD_Orphan）

- **方法**：mmseqs2 `easy-search -s 7` 高灵敏度搜索，**query = test，target = train+val**；每个测试样本取最大 `pident`，`seq_Redundancy_<T>` = max pident < T（T = 90…30）；**无任何显著 hit（e-value ≤ 1e-3）→ `OOD_Orphan=True`**。
- **重算路径**：`.skills/homology-ood-annotation/scripts/seq_homology_ood.py`。**现值为 2026-08-30 按本口径重算**（m8 缓存 `analysis/output/seq_ood_unify_assessment/func_{ec,go_bp,go_cc,go_mf}/mmseqs_easysearch.m8`，可 `--m8-cache` 重放；备份 `output/backup_before_seq_easysearch_unify/`，评估报告 `analysis/output/seq_ood_unify_assessment/REPORT.md`）。历史缓存 `output/test_vs_trainval.tsv`（6/9，`search`+`convertalis` 分数制 fident 格式）仅供溯源：实测 `OOD_Orphan` 旧列与它逐行一致，但 `seq_Redundancy_*` 旧列与它大面积不符（ec 即 1,018 行），且比现行口径系统性更严（高估 OOD），来源未能坐实（疑似低灵敏度/截断参数的历史 `mmseqs search` 或 UniRef50 口径叠加）。
- **[superseded]**：归档脚本 `orphan_mark.py` 的旧口径以 UniRef50 为参考集且要求 fident ≥ 0.30，与现行口径不同，不可用现行 skill 重放（见"注意事项"）。

| 任务 | Default test | 90 | 80 | 70 | 60 | 50 | 40 | 30 | OOD_Orphan |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **EC** | 4,243 | 94.0% | 85.5% | 78.2% | 67.9% | 52.9% | 34.2% | 17.2% | 498 (11.7%) |
| **GO-MF** | 4,609 | 93.3% | 83.7% | 74.5% | 63.8% | 51.4% | 37.6% | 24.4% | 902 (19.6%) |
| **GO-CC** | 4,256 | 92.7% | 82.5% | 74.0% | 65.1% | 55.2% | 43.9% | 31.6% | 1,059 (24.9%) |
| **GO-BP** | 4,121 | 93.0% | 83.3% | 73.9% | 63.0% | 51.3% | 38.8% | 25.7% | 834 (20.2%) |

> 上表为 2026-08-30 统一后实测（旧值见 `output/backup_before_seq_easysearch_unify/`）。统一使 seq_Redundancy 各档 True 数下降（旧列高估 OOD）、Orphan 上升（旧口径低估孤儿）、InD 上升（见 4.11）；旧口径的 3 行 Orphan⊈seq30 违例随之消解。

### 4.4 结构同源 OOD（TM-score_*）

- **方法**：Foldseek `easy-search`，query = test 结构，target = train+val 结构库；取最大 `alntmscore`，`TM-score_<T>` = max TM-score < T（T = 0.9…0.3）。
- **重算路径**：`.skills/homology-ood-annotation/scripts/struct_homology_ood.py`（同样支持 `--m8-cache` 重放）。

| 任务 | Default test | 0.9 | 0.8 | 0.7 | 0.6 | 0.5 | 0.4 | 0.3 |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| **EC** | 4,243 | 32.1% | 13.5% | 7.1% | 3.6% | 2.5% | 1.3% | 0.7% |
| **GO-MF** | 4,609 | 39.9% | 19.7% | 11.5% | 6.8% | 4.4% | 3.1% | 2.0% |
| **GO-CC** | 4,256 | 46.0% | 27.2% | 16.8% | 10.3% | 7.0% | 4.7% | 3.4% |
| **GO-BP** | 4,121 | 41.6% | 21.8% | 12.4% | 7.7% | 5.3% | 3.6% | 2.4% |

> Foldseek 对 train+val 的 hit 覆盖率极高（97.6%–99.5%），TM-score < 0.5 的样本占 2–7%，构成严格的"新折叠"OOD 子集。

### 4.5 OOD-IDR（高内在无序区域）

- **工具**：metapredict v3.0.2；预测对象为全部代表链序列（41 条含非标准残基的序列先删除 `X/B/Z/J/U`）
- **判定**：**绝对阈值法**——`IDR_ratio`（无序残基占比）> 0.3 → True
- **重算路径**：`.skills/ood-annotation-toolkit/scripts/idr_ood.py --mode residue-ratio --threshold 0.3 --clean strip`（本任务历史口径；缓存 `output/idr_features.csv`）

> **为什么不用相对阈值**：PDB 数据 86.3% 的链 IDR_ratio = 0，test 与 train+val 分布几乎一致，P90/P95 相对阈值会得到 0 个 OOD 样本。

| 任务 | Default test | OOD_IDR |
|:---|---:|---:|
| **EC** | 4,243 | 27 (0.6%) |
| **GO-MF** | 4,609 | 168 (3.6%) |
| **GO-CC** | 4,256 | 175 (4.1%) |
| **GO-BP** | 4,121 | 157 (3.8%) |

### 4.6 OOD_FoldHoldout / OOD_SuperfamilyHoldout（CATH+TED 域标注）

- **流程**：`scripts/13_ted_augment_cath_ood.py`（wrapper，循环四个子任务）：
  1. 从 splits 生成链清单 `output/func_{sub}_ted_aug_chain_list.csv`（脚本内置，原为无脚本的历史一次性产物；2026-08-28 验证重生成与现有文件逐行一致）
  2. 调用 `.skills/cath-ted-domain-annotation/scripts/cath_ted_annotate.py` 生成链级标注
  3. 调用同 skill 的 `cath_ted_holdout.py` 做 holdout 判定并写回 splits
- **口径**（链模式，`unique_id = PDB_chain`）：标签 = CATH v4.4 实验标注（名称匹配，域区间经 SIFTS 映射回链区间）∪ TED putative 标注（仅保留与链 UniProt 区间重叠 ≥50% 域长的域；H 级提供 C.A.T.H，T 级仅 C.A.T）；参考集 = 各子任务 **train** 的 CATH+TED 并集（train/val/test 对称增强）；any 语义：链任一 Topology / Superfamily 未在参考集出现 → True；无注释 → False；Default=False 行按 4.1 惯例置 NaN。
- **标注覆盖率（链级）**：仅 CATH 43.4%–59.8%，CATH+TED 并集 **86.2%–92.3%**（ec 91.1%、go_bp 90.9%、go_cc 86.2%、go_mf 92.3%）。

| 任务 | test 行数 | Default=True | OOD_FoldHoldout | OOD_SuperfamilyHoldout |
|:---|---:|---:|---:|---:|
| **EC** | 5,047 | 4,243 | 31 | 103 |
| **GO-BP** | 5,142 | 4,121 | 50 | 151 |
| **GO-CC** | 5,613 | 4,256 | 57 | 185 |
| **GO-MF** | 5,609 | 4,609 | 31 | 186 |

**输出文件**（每子任务一套，`{sub}` ∈ {ec, go_bp, go_cc, go_mf}）：`output/func_{sub}_chain_cath_ted_labels.csv`（链级标签）、`output/func_{sub}_chain_domains.csv`（域级明细）、`output/func_{sub}_ted_aug_chain_list.csv`（链清单）、`output/backup_before_ted_aug/{sub}_test.csv`（备份）。

> TED 为 foldseek/foldclass 计算预测标签（putative），CATH 实验标注优先、TED 仅补缺；曾修复 splits 大写 PDB ID 与标签表小写键不匹配、TED 逗号分隔双标签解析两个问题（详见 skill 文档）。

### 4.7 OOD_NewEC（新 EC 功能类别）

- **口径（统一 all-not-in）**：参考集 = train+val 出现的 EC 集合；样本有 EC 注释且其**所有**对应层级 EC 都不在参考集 → True；无注释 → False；Default=False 行 → NaN。EC 源为 **SIFTS 链级注释**（`data/sifts/sifts_chain_ec.tsv.gz`）。
- **重算路径**：`.skills/ec-function-ood-annotation/scripts/add_unified_newec.py`。
- **func_ec 仅保留 L4 一列**（`OOD_NewEC_L4`）：L3 是本任务的预测目标层级（label 为 `x.x.x.-` 形式），不设 `OOD_NewEC_L3`。
- **结果**：`OOD_NewEC_L4` True = **381**（Default 4,243 的 9.0%）。
- **go_\* 三任务**的 `OOD_NewEC_L3` / `OOD_NewEC_L4` 见 4.9（与 LongTailEC 同批落地）。

> **[superseded]** 旧口径（any 语义、json 源、要求四段全非 `-`）曾得 441，已被统一口径取代；旧脚本 `newfunction_mark.py` 已归档（见"归档清单"）。备份 `output/backup_before_newec_unify/`。

### 4.8 OOD_LongTail（预测标签长尾）

- **口径**：取各任务**训练集**标签频次 **≤ 10** 的标签为长尾集合；测试样本**任一标签**落入长尾集合 → True（any 语义，预测标签层级：EC 为 L3、GO 为 term 级）。
- **脚本**：`scripts/12_longtail_mark.py`（现行版本直接读写 `splits/`、仅标记 Default=True 行；历史执行读写的是中间文件 `output/datasets_split_raw/{task}_test_ood.csv`，见"注意事项"的路径偏差记录。2026-08-28 只读重算与当前 splits 完全一致）。

| 任务 | 训练集标签数 | 长尾标签数（≤10） | 测试集 OOD 样本 | 占比 |
|:---|---:|---:|---:|---:|
| **EC** | 260 | 104 (40.0%) | 163 | 3.8% |
| **GO-MF** | 447 | 127 (28.4%) | 295 | 6.4% |
| **GO-CC** | 264 | 68 (25.8%) | 235 | 5.5% |
| **GO-BP** | 884 | 372 (42.1%) | 684 | 16.6% |

### 4.9 OOD_LongTail_EC（EC 功能长尾，SIFTS L4/L3 注释级）

与 4.8 的"预测标签长尾"不同，本维度衡量样本携带的 **EC 功能注释**（SIFTS 链级真实 L4/L3 编号，**不是** label 列）在 train 中的频次：

- **口径（统一 all-not-in）**：样本有 EC 注释且其**所有** EC（对应层级）在 **train** 中频次都 < 5 / < 10（**含 0**）→ True；无注释 → False；Default=False 行 → NaN。含 0 使 NewEC ⊆ LongTailEC，两场景语义嵌套。
- **func_ec 仅保留 L4 两列**（`OOD_LongTail_EC_L4_le5` / `OOD_LongTail_EC_L4_le10`）：L3 级与预测目标层级冲突且与 4.8 的 `OOD_LongTail` 重复，已移除 L3 两列。重算路径：`.skills/ec-function-ood-annotation/scripts/add_unified_longtail_ec.py`。

| 列 | True 数（Default n=4,243） | 占比 |
|:---|---:|---:|
| `OOD_LongTail_EC_L4_le5` | 1,155 | 27.2% |
| `OOD_LongTail_EC_L4_le10` | 1,752 | 41.3% |

> 占比远高于 label 级长尾：SIFTS L4 注释远比 260 类的 label 空间细粒度。备份 `output/backup_before_longtail_unify/`、`output/backup_before_longtail_allnotin/`。

- **go_\* 三任务的 EC 功能 OOD 六列**（`OOD_NewEC_L3/L4` + `OOD_LongTail_EC_{L3,L4}_le{5,10}`）由 `.skills/ec-function-ood-annotation/scripts/add_go_ec_ood.py` 一次写出，口径与上述完全一致：

| 任务 | Default test | 无EC(→False) | NewEC_L3 | NewEC_L4 | LT_L3_le5 | LT_L3_le10 | LT_L4_le5 | LT_L4_le10 |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| **GO-MF** | 4,609 | 2,993 | 2 | 183 | 34 | 106 | 483 | 736 |
| **GO-BP** | 4,121 | 2,928 | 8 | 118 | 55 | 83 | 370 | 509 |
| **GO-CC** | 4,256 | 3,220 | 9 | 136 | 44 | 97 | 431 | 550 |

> **GO 版功能 OOD 不可行（分析结论）**：GO 标签词表封闭（test 标签全部 ⊆ train，标签级 NewGO 恒为 0）；改用 SIFTS 链级全量 GO 注释按 all-not-in 计算也为 0（每条链都携带 train 中的常见 IEA 高频 term）。LongTailGO 同样不可用（标签级 le5 恒 0——词表构建时已按频次 ≥5 过滤）。故不建 GO 版 OOD，改为扩展 EC 功能 OOD 六列。

### 4.10 OOD_Combinatorial（新标签组合，无独立留存脚本）

- **口径**：对多标签样本生成所有无序标签对 `(A, B)`（`A < B`）；测试样本含有**至少一个标签对**未在 train+val 组合集合中出现 → True。
- **设计意图**：即使每个单独标签都在训练集出现过，未见过的共现组合仍可能预测失败。
- **说明**：本维度为历史一次性计算落盘，无独立留存脚本；复现时按上述口径对 `splits/{task}_test.csv` 的 label 列计算即可。

| 任务 | Default test | 多标签样本 | OOD_Combinatorial | 占 test 比例 |
|:---|---:|---:|---:|---:|
| **EC** | 4,243 | 313 (7.4%) | 58 | 1.4% |
| **GO-MF** | 4,609 | 2,547 (55.3%) | 449 | 9.7% |
| **GO-CC** | 4,256 | 2,998 (70.4%) | 472 | 11.1% |
| **GO-BP** | 4,121 | 2,008 (48.7%) | 808 | 19.6% |

### 4.11 InD（纯 in-distribution 标识）

- **口径**：所有 `OOD_*` + `seq_Redundancy_*` + `TM-score_*` 均为 False（NaN 视为 False）→ `InD=True`。
- **重算路径**：`.skills/ood-annotation-toolkit/scripts/finalize_ood_columns.py --mode ind`（所有 OOD 列齐备后统一重算；备份 `output/backup_before_ind/`）。

| 任务 | InD=True | test 总行数 | 占比 |
|:---|---:|---:|---:|
| **EC** | 151 | 5,047 | 3.0% |
| **GO-MF** | 194 | 5,609 | 3.5% |
| **GO-CC** | 192 | 5,613 | 3.4% |
| **GO-BP** | 172 | 5,142 | 3.3% |

> 2026-08-30 序列同源列统一 easy-search 口径后重算（原 121 / 101 / 77 / 82；备份 `output/backup_before_seq_easysearch_unify/`）。

### 4.12 输出文件格式

最终数据集保存于 `splits/`（12 个文件，`{task}` ∈ {ec, go_mf, go_cc, go_bp}）：

| 文件 | 列 |
|:---|:---|
| `{task}_train.csv` / `{task}_val.csv` | `unique_id, aa_seq, struct_file, label`（4 列，无 OOD 列） |
| `ec_test.csv`（31 列） | 核心 4 列 + `seq_Redundancy_90…30`、`TM-score_0.9…0.3`、`OOD_IDR`、`OOD_Orphan`、`OOD_NewEC_L4`、`OOD_LongTail`、`OOD_Combinatorial`、`Default`、`OOD_ExtremeShort`、`OOD_ExtremeLong`、`OOD_FoldHoldout`、`OOD_SuperfamilyHoldout`、`OOD_LongTail_EC_L4_le5`、`OOD_LongTail_EC_L4_le10`、`InD` |
| `go_*_test.csv`（34 列） | 核心 4 列 + `seq_Redundancy_*`、`TM-score_*`、`OOD_IDR`、`OOD_Orphan`、`OOD_LongTail`、`OOD_Combinatorial`、`Default`、`OOD_ExtremeShort/Long`、`OOD_FoldHoldout`、`OOD_SuperfamilyHoldout`、`OOD_NewEC_L3`、`OOD_NewEC_L4`、`OOD_LongTail_EC_{L3,L4}_le{5,10}`（四列）、`InD` |

- `label` 强制保存为字符串，多标签以 `;` 分隔（EC 为 L3 级 `x.x.x.-` 形式）
- `struct_file` 仅文件名（`{pdb_id}_{chain}.pdb` 或 `.cif`），文件位于 `pdbs/`

### 4.13 OOD 划分小结

| OOD 类型 | 工具/skill | 核心指标 | 阈值/规则 | 难度特征 |
|:---|:---|:---|:---|:---|
| TimeCutoff | —（阶段三 3.2 天然） | deposition date | 单一基线 | 全部 test 为 OOD |
| seq_Redundancy_* + Orphan | mmseqs2（homology skill） | 序列 identity / 显著 hit | 90…30；无 hit | 阈值越低越难 |
| TM-score_* | Foldseek（homology skill） | 结构 TM-score | 0.9…0.3 | 阈值越低越难 |
| OOD_IDR | metapredict（toolkit） | IDR 残基占比 | > 0.3 绝对阈值 | 高无序蛋白 |
| FoldHoldout / SuperfamilyHoldout | CATH+TED（cath-ted skill） | Topology / Superfamily | 留一法，train 参考 | 未见折叠/超家族 |
| OOD_NewEC_L4 | SIFTS EC（ec skill） | EC L4 存在性 | all-not-in，train+val 参考 | 新酶功能（ec 9.0%） |
| OOD_LongTail | 标签频次（`scripts/12_longtail_mark.py`） | train 标签频次 ≤ 10 | any 语义 | 低频预测标签 |
| OOD_LongTail_EC_* | SIFTS EC（ec skill） | train EC 频次 < 5/10（含 0） | all-not-in | 所有功能均长尾 |
| OOD_Combinatorial | 标签共现（无留存脚本） | train+val 未见标签对 | any 语义 | 新标签组合 |
| ExtremeShort / ExtremeLong | 长度统计 | 残基数 < 60 / > 1000（≤2000） | 绝对阈值，Default=False | 极端长度 |

---

## 阶段五（可选扩展）：AlphaFold DB 预测结构收集

> **仅用于训练/验证集扩展**，不进入严格测试集。与阶段一/二（PDB 实验结构）互补，覆盖无 PDB 结构的 Swiss-Prot 条目。
> 相关脚本（`stage2.py` / `stage2_uniprot.py` / `stage2_finish.py`）路径已陈旧、与当前 splits 无耦合，已归档至 `output/scripts_archive/`，仅作历史参考；其产物（`output/alphafold_*`、`output/representative_proteins.csv` 等）仍保留。

流程要点：

1. **数据下载**：AlphaFold DB Swiss-Prot 全量预测结构 + 全局 pLDDT 摘要（仅 reviewed 条目）；
2. **质量解析**：残基级 pLDDT（B-factor 列）→ `output/alphafold_plddt_profiles.json`；全局 pLDDT < 70 仅用于训练；
3. **功能注释反推**：UniProt Swiss-Prot 的序列 / GO（含证据代码）/ EC / 元数据 → `output/alphafold_uniprot_annotations.json`；
4. **结构-序列一致性验证**：AlphaFold ATOM 序列 vs UniProt FASTA 全局比对（identity ≥ 0.95 且 coverage ≥ 0.90 通过）；
5. **与阶段二交叉去重**：已有 PDB 实验结构的 UniProt 不纳入（实验结构优先）；差集作为 AlphaFold 独有补充 → `output/source_overlap_map.json`。

| 场景 | 是否使用 AlphaFold | 说明 |
|:---|:---|:---|
| 严格基准测试 | 否 | 测试集仅 PDB 实验结构 |
| 训练集扩充 / 大规模预训练 | 是 | 无 PDB 结构的 Swiss-Prot 条目 |
| 预测结构专项评估 | 需单独报告 | 独立构建 AlphaFold 测试集 |

---

## 数据质量检查

以 `.skills/ood-annotation-toolkit/scripts/check_splits.py`（只读）对 4 个子任务的 splits 实测（2026-08-27）：

1. **泄露检查**：4 个子任务 train∩val、train∩test、val∩test 的 `unique_id` 交集**均为空**。
2. **列结构**：train/val 均为 4 列核心列、不含 OOD 列；test = 核心列 + 完整 OOD 块（ec 31 列、go_* 34 列）；无废弃列（`OOD-*` 连字符、`labels` 等）。
3. **NaN 惯例**：4 个 test csv 的 Default=False 行（极端长度样本）非极端 `OOD_*` 列**全部为空**，合规。
4. **嵌套违例已修复（2026-08-30）**：ec / go_mf / go_cc 曾各有 **1 行** `OOD_Orphan=True` 但 `seq_Redundancy_30=False`（`7QS4_A`、`7WFF_c`、`8YWA_A`）；实测三条序列对 train+val 在 e≤1e-3 下零 hit（Orphan=True 正确），其 seq30=False 来自宽松 e-value（0.01–1，pident 31–40）的历史检索轮次。序列同源列已于 2026-08-30 整体按 easy-search 口径重算（§4.3），违例消解，现行数据 Orphan ⊆ seq_Redundancy_30 成立。
5. **异常值处理**：长度 > 2000 的链全局剔除（阶段二 2.3.1 / 四 4.1）；B 档 identity < 0.90 的推断链剔除（19,906 条，阶段二 2.2）。

---

## 关键数据文件说明

| 文件/目录 | 说明 | 所属阶段 |
|:---|:---|:---|
| `data/sifts/sifts_*.json` / `*.tsv.gz` | SIFTS 链级映射与 GO/EC/CATH 注释（项目级共享） | 阶段一 |
| `data/sifts/annotated_chains.csv` | 初筛候选链元数据（827,186 行） | 阶段一 |
| `data/pdb/` | 原始 PDB/mmCIF 结构文件（项目级共享，~255K 个 .gz） | 阶段一 |
| `pdbs/` | 链级结构文件（50,123 个；`.pdb` 43,567 + `.cif` 6,556） | 阶段二 / 四 |
| `output/representative_chains.csv` | UniProt 级代表链清单（02 产物，重跑时再生） | 阶段二 |
| `output/filtered_sequences.fa` | 长度过滤后 53,902 条序列 | 阶段二 |
| `output/final_nonredundant_representatives.csv` | 最终非冗余代表链 47,852 条 | 阶段二 |
| `output/extreme_length_chains.csv` | 单独保存的极端长度链 3,834 条 | 阶段二 / 四 |
| `output/deposition_dates.json` | 各链 deposition date | 阶段三 |
| `output/datasets_split_raw/` | 划分与清洗中间产物（标签、元数据、统计） | 阶段三 |
| `output/test_vs_trainval.tsv`、`output/trainval_db*` | mmseqs2 比对缓存（可 `--m8-cache` 重放） | 阶段四 |
| `output/idr_features.csv` | metapredict IDR 特征缓存 | 阶段四 |
| `output/func_{sub}_chain_cath_ted_labels.csv` 等 | CATH+TED 标注产物（每子任务一套） | 阶段四 |
| `output/backup_before_*/` | 各次 splits 写回前的备份 | 阶段四 |
| `output/scripts_archive/` | 归档的旧口径/旧线脚本（见"归档清单"） | — |
| `splits/` | **最终划分文件（12 个）** | 阶段四 |
| `output/alphafold_*`、`output/representative_proteins.csv` 等 | AlphaFold 扩展产物（可选） | 阶段五 |

---

## 归档清单（`output/scripts_archive/`，不再使用）

| 脚本 | 归档原因 |
|------|----------|
| `fix_multimodel_pdbs.py` | 多模型（NMR）文件一次性修复已落入 `pdbs/` 全库（1,940 个文件截断修复，备份 `output/backup_before_multimodel_fix/`）；根因（首个 `ENDMDL` 即停）已修入 `scripts/03_extract_representative_structures.py` 与 `scripts/10_extract_extreme_pdbs.py` |
| `stage2.py` / `stage2_uniprot.py` / `stage2_finish.py` | 阶段五 AlphaFold 可选扩展脚本；路径陈旧（指向旧工作目录）、与当前 splits 无耦合，保留"可选扩展"说明（见阶段五） |
| `newfunction_mark.py` | 旧 NewFunction any 口径（json 源、四段全非 `-`），已被 ec skill 统一 all-not-in 口径取代（441 → 381） |
| `orphan_mark.py` | 旧 Orphan 口径（参考集 UniRef50 + fident ≥ 0.30），与现行 train+val 参考口径不同 |
| `fetch_uniprot_metadata.py` | 旧时间切分辅助脚本（硬编码已失效路径） |

## 注意事项

1. **func_ec 仅保留 L4 级 EC 功能 OOD 列**：L3（`x.x.x.-`）是本任务的预测目标层级，故不设 `OOD_NewEC_L3` / `OOD_LongTail_EC_L3_*`；L3 级长尾已由预测标签级 `OOD_LongTail`（4.8）覆盖。go_* 三任务的 label 是 GO term，与 EC 层级不冲突，故保留完整六列。
2. **`pdbs/` 双格式规则**：单字符 chain ID → `.pdb`，多字符 chain ID → `.cif`（PDB 格式列宽限制）；极端长度链中 cif 占 **20.4%**（784 / 3,834）。下游读取结构时须按 `struct_file` 扩展名分派解析器。
3. **代表链筛选的重建验证偏差（2026-08-28，新增）**：`scripts/02_select_representative_chains.py` 为缺口重实现。验证重跑（B 档全量入池）选出 56,057 条代表链，与最终非冗余代表链比对 **81.2%（38,835 / 47,852）精确一致**；9,017 个 UniProt 选出不同链。偏差来源：(a) 历史 B 档边界推断 identity ≥ 0.90 过滤的产物 `b_boundary_inference.jsonl` 未迁移，重跑时 B 档全量入池（554,852 vs 历史 419,370）；(b) A/B 档分类口径偏差（重跑 A 档 214,878 vs 历史 330,303，历史"A 档部分边界推断"的判定细节未文档化）；(c) SIFTS 再生漂移（候选池 769,892 vs 历史 769,730 条）。不一致属历史口径偏差，未改动 `pdbs/` 任何文件；如需精确复现，先重做 B 档边界推断并以 `--b-keep-list` 传入通过清单。
4. **B 档 identity 历史偏高**：当次边界推断使用全长 UniProt 序列比对（旧方法），identity/coverage 系统性偏高；修正为 `UniProt[sp_beg:sp_end]` 子序列比对后分布会更严格。重跑阶段二 2.2 时应注意此差异。
5. **Orphan ⊆ seq_Redundancy_30 的 3 行违例（已修复，2026-08-30）**：`7QS4_A`（ec）、`7WFF_c`（go_mf）、`8YWA_A`（go_cc）曾被标 Orphan=True 但 seq_Redundancy_30=False。机制：三条序列对 train+val 在 e≤1e-3（-s 7 / -s 7.5）下实测零 hit（Orphan=True 正确），而宽松的 e-value 上限（0.01–1）下存在 pident 31–40 的弱 hit——其 seq30=False 即来自该宽松口径的历史轮次；Orphan 列与 6/9 缓存 `output/test_vs_trainval.tsv` 逐行一致，seq 列则与该缓存大面积不符、来源未坐实。2026-08-30 四子任务序列同源列已按 skill easy-search（`-s 7`）整体重算写回（§4.3），嵌套关系恢复。
6. **无脚本环节**：MMseqs2 95% 聚类（`mmseqs cluster --min-seq-id 0.95 -c 0.8 --cov-mode 0`）、cluster 内代表链筛选（按 2.1 同优先级）、`sifts_chain_uniprot.tsv.gz` 等 SIFTS flat file 下载、go-basic.obo 下载、阶段三到四的 schema 落盘、OOD_Combinatorial 计算（4.10）、多模型 pdb 文件全库扫描（产物 `output/dup_multimodel_files.json`）。
7. **`go-basic.obo` 需自行下载**：`scripts/06_split_datasets.py` 依赖 GO OBO 文件解析 namespace（https://purl.obolibrary.org/obo/go/go-basic.obo，放在 `output/` 下），未随仓库分发；`output/go_id_to_namespace.json` 是既成的解析结果缓存（供 `scripts/07_clean_datasets.py` 使用）。
8. **`scripts/12_longtail_mark.py` 的 I/O 路径偏差**：历史执行在极端长度链并入之前运行，读写中间文件 `output/datasets_split_raw/{task}_train.csv` 与 `{task}_test_ood.csv`（未保留）；现行版本直接读写 `splits/` 且仅标记 Default=True 行，只读重算与当前 splits 数值一致（163 / 295 / 235 / 684）。
9. **旧路径映射**：早期文档中的 `datasets_split/`（最终划分）即现在的 `splits/`；阶段三中间产物统一在 `output/datasets_split_raw/`；项目级共享数据（SIFTS、原始 PDB、UniRef50、CATH、TED 缓存）在 `POOR/data/` 下；早期 `data/intermediate/` 中间目录未迁移，其产物路径已改为 `output/`（representative_chains.csv、b_boundary_inference.jsonl、filtered_sequences.fa 等）。
10. **多模型（NMR）修复已落数据**：全库曾发现 1,940 / 43,567 个 pdb 文件将全部模型无 `MODEL` 分隔地串联写入；已按"首个非连续重复残基处截断"修复（备份 `output/backup_before_multimodel_fix/`），重扫 0 残留；根因已修入提取脚本（首个 `ENDMDL` 即停），一次性修复脚本 `fix_multimodel_pdbs.py` 已归档。`.cif` 文件（6,556 个）扫描确认全部单模型。
11. **Default=False 行 seq/TM 列 NaN 统一（2026-08-28）**：四个子任务 test 的 Default=False（极端长度）行历史遗留 `seq_Redundancy_*` / `TM-score_*` 为 False（共 58,548 单元格），已统一为 NaN（"未计算"语义，与全项目惯例对齐；备份 `output/backup_nondefault_seqtm_nan/`，`InD` 不受影响已验证）。此后各生成脚本与 skill 脚本在流水线位置直接写入 NaN，`finalize_ood_columns.py --mode nan` 仅作兜底校验（当前数据上为 no-op）。
