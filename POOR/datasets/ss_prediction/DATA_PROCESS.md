# Secondary Structure Prediction (SSP) 数据集处理文档

> 记录从 PDB 实验结构到 SSP 数据集的完整处理流程，可按本文从原始数据一次性复现。
> 处理位置：`datasets/ss_prediction/`（二级结构预测任务主目录）
>
> **任务定义**：高质量 PDB 单链的二级结构预测，`label` 为逐残基 8-state DSSP
> 标注（整数映射 `H0 G1 I2 E3 B4 T5 S6 P7 C8`，`-1` = ignore index），
> `unique_id` 为 `PDB_chain`（如 `10GW_A`）。

**当前 splits 规模**（以 `splits/` 实际文件统计为准）：

| 数据集 | 数量 | 占比 | 沉积日期范围 | 长度范围 | 说明 |
|--------|------|------|-------------|---------|------|
| `splits/ssp_train.csv` | **29,622** | 71.7% | 1994-03 ~ 2019-06 | 31 – 998 | In-Distribution 训练集（含 834 条 31–59aa 短链） |
| `splits/ssp_val.csv` | **3,317** | 8.0% | 2019-06 ~ 2021-02 | 31 – 931 | In-Distribution 验证集（含 118 条短链） |
| `splits/ssp_test.csv` | **8,372** | 20.3% | 1999-07 ~ 2026-05 | 31 – 1,836 | Test 集（7,997 条时间切分 + 293 条 ExtremeShort 短链 + 82 条极长链），34 列（含完整 OOD 标记 + `InD` + `struct_label`） |

> 上表为当前最新状态（短链口径 31–59aa，见阶段一 1.3 与阶段四 4.1）。
> 全集 41,311 条 = ≥60 主池 39,984 + 31–59 短链 1,245 + 极长链 82。

## 设计原则

1. **直接从 PDB 出发**：不依赖 SIFTS→UniProt 映射，直接从 RCSB PDB 获取元数据并从原始结构文件提取链级信息。
2. **质量优先**：以 X-ray 晶体学为主力，严格限制分辨率（≤2.5 Å）和 R-free（≤0.25），确保 DSSP 标签的可靠性。
3. **可追溯**：每个阶段保留原始数据、中间结果和统计日志，支持断点续传和结果复现。

---

## 阶段一：源数据获取与解析

> 目标：从本地 PDB 镜像获取 entry 级元数据并做质量筛选，提取蛋白链序列。

### 1.1 PDB Entry 级元数据获取 —— `scripts/01_fetch_pdb_metadata_api.py`

**数据来源**：RCSB PDB GraphQL API (`https://data.rcsb.org/graphql`)

**查询字段**：

| GraphQL 字段 | 对应属性 | 说明 |
|-------------|---------|------|
| `exptl.method` | 实验方法 | X-RAY DIFFRACTION / SOLUTION NMR / ELECTRON MICROSCOPY 等 |
| `refine.ls_d_res_high` | X-ray 分辨率 | Å，X-ray 结构核心质量指标 |
| `refine.ls_R_factor_R_work` | R-work | 工作集 R-factor |
| `refine.ls_R_factor_R_free` | R-free | 自由集 R-factor，独立验证指标 |
| `em_3d_reconstruction.resolution` | EM 分辨率 | Cryo-EM 分辨率（X-ray 缺失时回退） |
| `rcsb_accession_info.deposit_date` | 沉积日期 | 用于后续时间切分 |

**处理逻辑**：Batch 查询 50 entries / request；每 500 entries 保存 checkpoint
（`output/api_checkpoint.json`）；单 batch 失败最多重试 5 次（指数退避）；
请求间隔 0.3s，遇 429 限流指数退避。

**输入**：`data/pdb/` 文件清单（255,227 entries，含 `.pdb.gz` + `.cif.gz`）
**输出**：

| 文件 | 说明 | 字段 |
|------|------|------|
| `output/pdb_metadata_api.csv` | 全量 PDB 元数据（254,975 行） | `pdb_id`, `method`, `resolution`, `r_work`, `r_free`, `deposition_date` |
| `output/api_checkpoint.json` | 断点检查点 | `processed` |

**实际执行规模**：5,105 个 batch，耗时 ~2.5 小时，API 成功率
254,975 / 255,227（**99.9%**）。

### 1.2 Entry 级质量筛选 —— `scripts/02_filter_entries.py`

**筛选条件**：

| 条件 | 阈值 | 理由 |
|------|------|------|
| **实验方法** | `X-RAY DIFFRACTION` | X-ray 坐标精度最高，DSSP 可靠性最佳 |
| **分辨率** | ≤ **2.5 Å** | 2.5 Å 是 DSSP 氢键判定的经验边界 |
| **R-free** | ≤ **0.25** | 自由集验证指标，≤0.25 属于高质量结构 |

**输入**：`output/pdb_metadata_api.csv`
**输出**：

| 文件 | 说明 |
|------|------|
| `output/stage1_filtered_entries.csv` | 筛选后的 entry 列表（121,160 行） |
| `output/stage1_filtering_stats.json` | 筛选统计 |

**实际筛选规模**：

| 步骤 | 实际数量 | 淘汰率 | 说明 |
|------|---------|--------|------|
| 全量 entries | **254,975** | — | API 成功返回的 entries |
| X-ray 筛选后 | **204,781** | -19.7% | X-ray 占比 80.3% |
| Resolution ≤2.5 Å | **160,940** | -21.4% | 分辨率均值 1.80 Å |
| **R-free ≤0.25** | **121,160** | **-24.7%** | **R-free 均值 0.211** |

> 筛选后数据质量极高：平均分辨率 **1.80 Å**，平均 R-free **0.211**。部分早期
> PDB entry 未记录 R-free（如 1UBQ），在筛选时被排除。

### 1.3 Chain 序列提取 —— `scripts/03_extract_chain_sequences.py`

**目标**：从 1.2 通过的 PDB 文件中提取所有蛋白链的氨基酸序列，仅用于后续
去重，不保存坐标文件。

**处理策略**：
- 对每个通过的 PDB entry（121,160 个），读取原始 `.pdb.gz` / `.cif.gz`
- 提取所有 ATOM/HETATM 记录中的蛋白链（仅 CA 原子去重）
- 将 3-letter 残基代码映射为 1-letter 序列（MSE→MET, SEC→CYS）
- NMR 结构仅取 model 1
- **长度下限 `MIN_RESIDUES=31`**：31–59aa 短链一并入池，≤30aa 在源头排除

**输入**：`output/stage1_filtered_entries.csv`、`data/pdb/`
**输出**：

| 文件 | 行数 | 说明 |
|------|------|------|
| `output/stage2_chain_sequences.fa` | 255,471 | 所有提取链的 FASTA 序列 |
| `output/stage2_chain_metadata.csv` | 255,471 | 链级元数据（pdb_id, chain_id, num_residues, has_mse, file_type, seq） |

> 现存磁盘上的这两个文件是 v1（`MIN_RESIDUES=60`）产物，为 249,018 行
> （不含短链）；按现行 03 重跑产出 255,471 行 = 249,018 + 31–59 短链
> 6,453。短链部分的实际提取记录为 `output/short_chains.fa` /
> `output/short_chains_meta.csv`（1–59aa 共 22,231 条，其中 31–59 段
> 6,453 条，见历史沿革）。

**提取结果统计**（≥60 主池实测；短链 6,453 条长度 31–59、中位约 32）：

| 指标 | 数值 |
|------|------|
| 处理的 PDB entries | 121,160 |
| 成功提取链的 entries | 118,182 (97.5%) |
| 提取的总链数 | **255,471**（≥60 的 249,018 + 31–59 的 6,453） |
| 单链 entries | 56,045 |
| 多链 entries | ~62,137 |
| 单 entry 最多链数 | 88（大分子复合物）|
| 平均链长度 | 271.3 残基（≥60 主池） |
| 长度范围 | 31 – 2,082 |

---

## 阶段二：去冗余与数据精化（全局序列聚类）

### 2.1 MMseqs2 全局去冗余 —— `scripts/04_run_mmseqs2_clustering.py`

**目标**：在 255,471 条链中聚类去冗余，得到 41,337 条代表链（40,067 主池
+ 1,270 短链）。

**工具与参数**：MMseqs2 `easy-linclust`，`--min-seq-id 0.95 -c 0.8 --cov-mode 0
--threads 32`，对 stage2 全部链（len≥31）**单阶段一次聚类**（直接喂
FASTA——mmseqs 18.8cc5c 的 easy-linclust 直接接受 FASTA 输入；预建 DB
在该版本会报 createdb 错误）；代表选择策略为簇内最长序列优先、再按
chain_key 字母序（04 脚本 `select_representatives`）。

> 复现说明：easy-linclust 对输入顺序敏感，边界簇的成员归属可能有 ±少量
> 差异。当前发布的数据由历史两阶段路线产出（主池 ≥60 单独聚类 + 短链
> <60 单独聚类并交叉去重 + 事后删除 ≤30 代表，见文末历史沿革）；单阶段
> 全量聚类在"簇内最长优先"规则下与两阶段路线等价（≥95%/0.8 双向覆盖下
> 短链只会并入 ≥60 簇并落选代表，等效于交叉去重丢弃），允许存在少量
> 偏差。短链聚类的历史实际值（在 1–59aa 池 22,231 条上运行）：9,229 簇
> → 交叉去重丢 150 → 幸存 9,079，其中 31–59 段 1,270 条。

**输入**：`output/stage2_chain_sequences.fa`、`output/stage2_chain_metadata.csv`
**输出**：

| 文件 | 行数 | 说明 |
|------|------|------|
| `output/stage4_nonredundant_sequences.fa` | 41,337 | 非冗余代表链 FASTA |
| `output/stage4_cluster_representatives.csv` | 41,337 | 代表链元数据（含 `chain_key`、`file_type`、`seq`；主池 40,067 + 短链 1,270） |
| `output/stage4_mmseqs_clusters.tsv` | — | 聚类结果 |
| `output/stage4_clustering_stats.json` | — | 聚类统计 |

> 现存磁盘上的 `stage4_*` 三个文件是 v1 产物（仅主池 40,067 行）；短链代表
> 的历史清单为 `output/short_chain_representatives.csv`（9,079 行，含 ≤30
> 段）；历史两阶段产物 `output/short_mmseqs_clusters.tsv`、
> `output/short_cross_dedup.m8` 保留备查。按现行 04 单阶段重跑产出合并的
> 41,337 行。

**聚类统计**（主池 ≥60 子集的历史实际值）：

| 指标 | 数值 |
|------|------|
| 输入总链数 | 249,018 |
| 总 clusters | **40,067** |
| Singleton clusters | 11,729 |
| Multi-member clusters | 28,338 |
| 去冗余率 | **83.91%** |

---

## 阶段三：标签聚合与计算（结构坐标提取 + DSSP 标注）

> 目标：为 41,337 条代表链（40,067 主池 + 1,270 短链）提取单链 PDB 坐标
> 文件，计算 DSSP 二级结构标签，并聚合为全量数据表。

### 3.1 单链 PDB 提取 —— `scripts/05_extract_representative_chains.py`

**策略**：从原始 `.pdb.gz` / `.cif.gz` 中提取代表链，保存为单链 PDB 文件
（多进程并行）。NMR 多模型条目仅保存 model 1（`io.set_structure(structure[0])`）。

**关键问题与解决**：
- DSSP v4 对 PDB 格式校验严格，缺少 `CRYST1` / `HEADER` / `EXPDTA` 或存在
  `ANISOU` 位置错误都会导致解析失败
- 使用 BioPython `PDBParser`/`MMCIFParser` 解析 + `PDBIO` 保存，可自动过滤掉
  `ANISOU` 等干扰记录
- 后处理添加最小化 PDB 头部（`HEADER`、`EXPDTA`、`CRYST1`、`END`）确保 DSSP v4 兼容
- 约 **367 条链**（0.9%）的链 ID 超过 1 字符（如 `AAA`、`1D`、`B2`），通过
  截断链 ID + 修改 BioPython chain 对象后成功提取。其中 55 条来自大型病毒
  结构（4Y4O、4YBB 等），因截断后链 ID 冲突导致残基类型不一致，后通过为
  每条链独立构建 `Structure` 对象修复

**输入**：`output/stage4_cluster_representatives.csv`、`data/pdb/`
**输出**：`pdbs/`（**41,327** 个 `.pdb` 文件 = 主池 40,067 + 短链 1,260；
短链 1,270 条代表中 23 条因 mmCIF 源多字符 chain ID 超出 PDB 格式限制
提取失败丢弃——与主池历史失败同类，另有 1 条短链在 DSSP 环节失败，见 3.2）

### 3.2 DSSP 标注 —— `scripts/06_run_dssp_parallel.py`

**工具**：`mkdssp` v4.6.1；24 进程并行，每个 worker 处理约 1,669 个文件。

**3-state 映射**：

| DSSP 原始状态 | 3-state 标签 | 结构类型 |
|--------------|-------------|---------|
| H (α-helix) | **H** | Helix |
| G (3₁₀-helix) | **H** | Helix |
| I (π-helix) | **H** | Helix |
| P (Polyproline II) | **H** | Helix |
| E (β-sheet) | **E** | Sheet |
| B (β-bridge) | **E** | Sheet |
| T (turn) | **C** | Coil |
| S (bend) | **C** | Coil |
| ` ` (coil) | **C** | Coil |

**输入**：`pdbs/`
**输出**：

| 文件 | 说明 |
|------|------|
| `output/dssp/` | 40,067 个 `.dssp` 原始输出文件（主池） |
| `output/dssp_short/` | 1,260 个短链 `.dssp` 原始输出文件 |
| `output/dssp_labels.csv` | 汇总表（40,067 行）：`chain_key`, `seq`, `ss8`, `ss3` |

**DSSP 标注统计**：主池成功标注 **40,067 / 40,067**（100%）；短链成功
**1,246 / 1,260**（1 条 `parse_fail` 丢弃，余 1,246 条含后被 QC 剔除的
7UYL_K，见 4.1）。主池总残基数 **10,407,197**；序列长度范围 55 – 2,082
残基（主池）；短链 31–59 残基。

**3-state 分布**：H 3,890,548（**37.4%**）/ E 2,510,675（**24.1%**）/
C 4,005,974（**38.5%**）。

**8-state 分布**：

| 标签 | 结构类型 | 残基数 | 占比 |
|------|---------|--------|------|
| H | α-helix | 3,199,406 | 30.7% |
| G | 3₁₀-helix | 419,045 | 4.0% |
| I | π-helix | 56,498 | 0.5% |
| P | Polyproline II (DSSP v4) | 215,599 | 2.1% |
| E | β-sheet | 2,386,784 | 22.9% |
| B | β-bridge | 123,891 | 1.2% |
| T | Turn | 1,181,359 | 11.4% |
| S | Bend | 863,036 | 8.3% |
| C | Coil | 1,961,579 | 18.8% |

### 3.3 全量数据聚合 —— `scripts/07_aggregate_ssp_full_data.py`

- **输入**：`output/stage4_cluster_representatives.csv`（行集、`aa_seq`、
  `file_type` 来源）、`output/dssp_labels.csv`（标签源 `ss8`）
- **输出**：`output/ssp_full_data.csv`（列
  `unique_id, aa_seq, file_type, labels`；`labels` 为 `[8 8 3 ...]` 风格
  空格分隔字符串，长度 == `len(aa_seq)`）
- **说明**：ss8 → 整数映射 `H0 G1 I2 E3 B4 T5 S6 P7 C8`。`aa_seq` 为链全长
  序列，比 DSSP 序列长（DSSP 丢弃缺骨架原子的残基）：序列一致时直接映射；
  不一致时（主池 5,737 条）用 `difflib.SequenceMatcher`（默认 autojunk=True）
  对齐——equal 块直接映射、replace 块按位置映射、未覆盖残基置 `-1`
  （ignore index）。该口径与现存文件逐行一致（已 diff 验证）。

> 现存磁盘上的 `ssp_full_data.csv` 是 v1 产物（40,067 行，仅主池，含后被
> 剔除的 4O9X_A）；短链标签按同一口径生成（difflib 对齐 235 条）。按现行
> 流程重跑 07 产出 **41,313 行** = 40,067 + 短链 1,246（含后被 QC 剔除的
> 7UYL_K）。

---

## 阶段四：数据划分与 OOD 标注

### 4.1 时间划分 —— `scripts/08_create_ssp_datasets.py`

- **输入**：`output/ssp_full_data.csv` + `output/stage4_cluster_representatives.csv`
  （链长）+ `output/pdb_metadata_api.csv`（沉积日期）
- **输出**：`splits/ssp_{train,val,test}.csv`；train/val 为
  `unique_id,aa_seq,struct_file,label`，test 追加 `Default,OOD_ExtremeLong`
- **划分策略**：
  - **超长链硬性剔除**：单链长度 > 2000 的样本直接剔除（全库仅 `4O9X_A`，
    2,082 残基），不进入任何集合
  - **QC 剔除**：`QC_EXCLUDE = {"7UYL_K"}`——该链 120 个结构残基中 71 个
    为 UNK（59% 未鉴定），序列提取按 `AA_MAP` 跳过 UNK（aa_seq 仅 49 残基）
    而结构提取保留 UNK，导致 label（49）与 struct_label（120）契约破坏、
    标注错帧（详见 4.4 质控记录）
  - **极长链**（`1,000 < seq_len ≤ 2000`）：**82 条**，不参与时间切分，直接
    并入 Test 集末尾（按 `unique_id` 排序），标记 `Default=False` /
    `OOD_ExtremeLong=True`
  - **剩余链**（`31 ≤ seq_len ≤ 1000`）：**41,229 条**，按 `deposition_date`
    升序稳定排序后按 **72:8:20 分位切分**（`n_train=int(n*0.72)`、
    `n_val=round(n*0.08)`，其余进 test），cutoff 为运行时对当前池动态计算。
    当前发布数据的实际 cutoff 为 2019-06-06 / 2021-02-09（系历史路线下对
    ≥60 主池分位计算、短链后按该边界补入，见历史沿革）；后续复现的 cutoff
    随池组成变化，不要求与本次一致

**划分结果**：train **29,622** / val **3,317** / test **8,372**
（= 7,997 条主池时间切分 + 293 条 31–59aa 短链（`OOD_ExtremeShort=True`）
+ 82 条极长链），合计 **41,311**。

**关键性质**：
- 按日期排序切分天然保证了**同一 PDB entry 的所有链落在同一集合**（无
  cross-set 泄漏）
- `Train/Val/Test` 之间 `unique_id` 无重叠（两两交集均为 0）
- 当前 splits 由历史路径产出（见历史沿革）；从头重跑与之的差异：
  (a) linclust 边界 tie-break；(b) **cutoff 随池组成移动**——重跑对含短链
  的全池重新分位，cutoff 日期与本次发布值（2019-06-06 / 2021-02-09）不同，
  边界附近的链归属会有少量差异；(c) 行序（历史路径将短链附于 train/val
  末尾，重跑按日期排序插入对应位置）。这些偏差不影响数据质量，属允许范围

### 4.2 OOD 标记总览

基于 Train+Val（32,939 条，含 952 条 31–59aa 短链）作为参考库，对 Test
（8,372 条）进行 OOD 标记（布尔列嵌入 `ssp_test.csv`）。

> **极端长度链的特殊处理（NaN 惯例）**：82 条极长链（`Default=False`）不参与
> 序列/结构比对与其他 OOD 维度计算，除 `OOD_ExtremeLong=True` 外其余所有
> OOD 列统一置 **NaN**（空单元格），`InD=False`。293 条 ExtremeShort 短链
> 正常参与全部计算（`Default=True`、`OOD_ExtremeShort=True`）。NaN 收尾口径由
> `.skills/ood-annotation-toolkit/scripts/finalize_ood_columns.py --mode nan`
> 统一执行。以下各维度统计均针对 **Default=True 的 8,290 条**。

| 列名 | 说明 |
|------|------|
| `Default` | 时间切分进入 Test 的 8,290 条为 `True`（含 293 条短链），极长链 82 条为 `False` |
| `OOD_ExtremeShort` | 长度 31–59 的 293 条短链为 `True`（正常参与 OOD 计算） |
| `OOD_ExtremeLong` | 长度 > 1000 的 82 条链为 `True` |
| `InD` | 纯 ID 样本：所有 `OOD_*` / `seq_Redundancy_*` / `TM-score_*` 标记均为 False（NaN 视为 False）。当前 **1,123 / 8,290（13.5%）**（2026-09-02 短链入池重算后；历史值见 4.4） |

#### 序列同源 OOD：`seq_Redundancy_90`~`30` + `OOD_Orphan`（skill）

方法学见 `.skills/homology-ood-annotation`（脚本
`.skills/homology-ood-annotation/scripts/seq_homology_ood.py`）：MMseqs2
`easy-search`（query=Test，target=Train+Val），取每个 test 样本对参考库的
最大序列同一性（`pident`），`seq_Redundancy_<T> = (max pident < T)`，无 hit
按 0 计；`OOD_Orphan` = 无任何 hit（保证 `OOD_Orphan ⊆ seq_Redundancy_30`）。
**现值为 2026-09-02 短链入池后重算**（easy-search `-s 7.5`、e≤1e-3、无 cov
过滤，m8 缓存 `output/ood_search/mmseqs_easysearch_test_vs_trainval_v3.m8`，
可 `--m8-cache` 重放）。更早的统一重算（2026-08-30，缓存
`analysis/output/seq_ood_unify_assessment/ssp/mmseqs_easysearch.m8`）与 v2 缓存
（`output/ood_search/mmseqs_easysearch_test_vs_trainval_v2.m8`）已被取代；
`output/ood_search/` 的历史中间结果与现值不符（旧实现系统性漏检 hit、高估
OOD），仅供溯源。

| 阈值 | OOD 数量 | OOD 占比 |
|------|---------|---------|
| < 90% | 6,894 | 83.2% |
| < 80% | 6,139 | 74.1% |
| < 70% | 5,548 | 66.9% |
| < 60% | 5,035 | 60.7% |
| < 50% | 4,327 | 52.2% |
| < 40% | 3,283 | 39.6% |
| < 30% | 2,016 | 24.3% |
| **OOD_Orphan** | **1,566** | **18.9%** |

> 上表为 2026-09-02 短链入池重算后 Default=True 子集实测（占比分母 8,290）。

#### 结构同源 OOD：`TM-score_0.9`~`0.3`（skill）

方法学见 `.skills/homology-ood-annotation`（脚本
`.skills/homology-ood-annotation/scripts/struct_homology_ood.py`）：Foldseek
**easy-search**（query=Test 结构，target=Train+Val 结构），取最大 `alntmscore`，
`TM-score_<T> = (max alntmscore < T)`，无 hit 按 0 计。target 含 31–59aa
短链结构。现值为 2026-09-02 短链入池后重算；m8 缓存：
`output/ood_search/foldseek_easysearch_test_vs_trainval_v3.m8`（可用
`--m8-cache` 重放判定）。参考集含大量短链后，部分 query 的 max alntmscore
被短链 target 抬高，可出现 >1 值（v2 重算时 305 个 query，最大 1.072；
为短 target 归一化所致，属预期现象）。

| 阈值 | OOD 数量 | OOD 占比 |
|------|---------|---------|
| < 0.9 | 2,863 | 34.5% |
| < 0.8 | 1,356 | 16.4% |
| < 0.7 | 696 | 8.4% |
| < 0.6 | 347 | 4.2% |
| < 0.5 | 188 | 2.3% |
| < 0.4 | 127 | 1.5% |
| < 0.3 | 102 | 1.2% |

> **口径修正记录（2026-08-28）**：此前四档列（0.9~0.6）实际源自旧式
> `foldseek search -e inf --max-seqs 1000`（target 误用整个 pdbs 目录含 test
> 自命中），且 m8 中 778 条 alnTMscore >1（最大 1.5，非法值），四档计数被
> 严重低估（如 0.6 档 19 → 252，0.9 档 2,194 → 2,713），列来源无法从现存
> 脚本复现。已按 skill easy-search 口径重算并补齐 0.5/0.4/0.3 三档（旧文档
> "TM-score < 0.5 样本 ≤5" 的结论随之失效）。重跑产物与对比见
> `output/foldseek_easysearch_compare/`，旧 test 备份于
> `output/backup_before_foldseek_easysearch_unify/`。旧搜索中间文件
> （`output/ood_search/foldseek_result*.m8` 等）仅作历史留存。

#### `OOD_IDR`（内在无序区域）—— `scripts/09_run_idr_prediction.py`

使用 **metapredict v3** 对全部代表链序列进行 IDR 预测。口径与
`.skills/ood-annotation-toolkit/scripts/idr_ood.py` 的
`--mode residue-ratio --threshold 0.3 --clean strip` 一致。

**处理流程**：
1. 序列预处理：删除非标准氨基酸（`X/B/Z/J/U`）后预测（本数据集中无此类残基）
2. IDR 定义：`disorder score > 0.5` 的残基
3. 核心特征：`IDR_ratio` = disorder > 0.5 残基总数 / 序列全长
4. 判定规则：`IDR_ratio > 0.3` → `OOD_IDR=True`

**输出**：`output/idr_predictions.csv`（40,067 行，主池）+ 写回
`splits/ssp_test.csv`。短链按同口径经 `idr_ood.py` 计算（缓存
`output/idr_cache_short.csv`，判定 `output/short_idr_flags.csv`）。

| 范围 | OOD-IDR 数量 | 占比 |
|------|-------------|------|
| 主池 40,067 条 | 474 | **1.18%** |
| Test（Default=True 8,290 条） | 159 | **1.9%** |

> 主池平均 `IDR_ratio` 仅 **2.05%**，远低于阈值 0.3。高 IDR 蛋白在实验结构
> 数据库中天然稀缺（晶体学难以解析高度无序区域）。31–59aa 短链的 IDR 比例
> 显著更高（短链全段易被判为无序），但入 test 的 293 条短链中仅 40 条
> True——大部分短链按日期 cutoff 分入了 train/val。

#### `OOD_LongTail`（8-state 长尾标签富集）—— `scripts/10_run_longtail_prediction.py`

基于 8-state DSSP 标签的全局分布，**I (π-helix, 0.54%)、B (β-bridge, 1.19%)、
P (Polyproline II, 2.07%)** 是占比最低的三类。对每条蛋白质计算
`LT_ratio = (I + B + P 残基数) / 有效残基总数`，取 **Top 5%** 标记为
`OOD_LongTail`。参考集为含短链的全体 41,311 条，当前阈值
**0.082914**（95 分位数，`output/longtail_threshold_v3.json`；历史阈值
0.081181（40,067 参考）→ 0.093333（49,007 参考）→ 0.082914）。

**输出**：`output/longtail_predictions.csv`（40,067 行，主池历史产物）+
写回 `splits/ssp_test.csv`。

| 范围 | OOD-LongTail 数量 | 占比 | LT_ratio 阈值 |
|------|------------------|------|--------------|
| 全量 41,311 条 | 2,066 | **5.0%** | ≥ 0.082914 |
| Test（Default=True 8,290 条） | 436 | **5.3%** | — |

> **注意**：Top 5% 阈值下 LT_ratio 最低为 8.29%，约为全局均值（3.65%）的
> 2.3 倍。Top 蛋白质几乎全部为 **P-rich**（PPII helix），I 和 B 在单蛋白中
> 极少高占比。

#### `OOD_NewEC_L3/L4`（新酶功能类别，skill）

**脚本**：`.skills/ec-function-ood-annotation/scripts/add_unified_newec.py`

- EC 数据源：SIFTS 链级注释（`data/sifts/sifts_chain_ec.tsv.gz`）
- 参考集 = train+val 出现的 EC 集合（对应层级）
- **all-not-in** 语义：样本有 EC 注释且其**所有** EC（对应层级）都不在
  train+val 中 → True；无 EC 注释 → False
- 层级解析：段数 ≥3 → L3=前三段；段数 ==4 → L4=完整编号（允许部分 "-"）
- Default=False 行两列置 NaN；备份：`output/backup_before_newec_unify/`

| 列名 | True 数（Default n=8,290，有 EC 注释 2,722） |
|------|------:|
| `OOD_NewEC_L4` | **245** |
| `OOD_NewEC_L3` | **1** |

> EC 注释对 Default 测试样本的覆盖率为 **32.8%**（2,722 / 8,290）。L3 新类别
> 规模过小，以 **L4 为主**。初版本地脚本 `run_ec_newcategory.py`（any 语义、
> 按日期截止推导划分）已被该 skill 取代并归档（superseded）。
> 现值为 2026-09-02 短链入池（train+val EC 参考集扩大）后 ss-only 重算。

#### `OOD_LongTail_EC_*`（EC 长尾功能，skill）

**脚本**：`.skills/ec-function-ood-annotation/scripts/add_unified_longtail_ec.py`

- 参考集：**train only**（逐样本计数：携带某 EC 类别的 train 样本数）
- 阈值：train 频次 < 5 / < 10（**含 0**）
- **all-not-in** 语义：样本所有 EC（对应层级）都满足频次条件 → True；无 EC
  注释 → False
- Default=False 行四列置 NaN；备份：`output/backup_before_longtail_unify/`

四列（追加到 test 末尾）：`OOD_LongTail_EC_L3_le5` / `OOD_LongTail_EC_L3_le10`
/ `OOD_LongTail_EC_L4_le5` / `OOD_LongTail_EC_L4_le10`。

| 列 | True 数（Default n=8,290，有 EC 2,722） |
|------|------:|
| `OOD_LongTail_EC_L3_le5` | 37 |
| `OOD_LongTail_EC_L3_le10` | 77 |
| `OOD_LongTail_EC_L4_le5` | 738 |
| `OOD_LongTail_EC_L4_le10` | 1,118 |

> L4 占比偏高（13.5%）是因为本任务按时间切分（train 为 2019-06 前条目），test
> 新条目携带的 EC 在旧 train 中天然低频。现值为 2026-09-02 短链入池后
> ss-only 重算。

#### `OOD_FoldHoldout` / `OOD_SuperfamilyHoldout`（CATH+TED 域标注）—— `scripts/11_ted_augment_cath_ood.py`

本任务驱动脚本依次完成三步（标注与 holdout 计算走
`.skills/cath-ted-domain-annotation`）：

0. **链清单生成**：从 `splits/ssp_{train,val,test}.csv` 的 `unique_id` 拆出
   小写 `pdb` + `chain` 并附 `splits` 列 → `output/ss_ted_aug_chain_list.csv`
   （**41,311** 条链 = 全库减去划分时剔除的 4O9X_A 与 7UYL_K）
1. **标注**：`.skills/cath-ted-domain-annotation/scripts/cath_ted_annotate.py`
   （`--offline`，TED 缓存基本命中项目级 `data/ted/ted_api_cache.json`，
   309 个未缓存 accession 记 none）→
   `output/ss_chain_cath_ted_labels.csv`、`output/ss_chain_domains.csv`
   （94,998 行域级明细）
2. **holdout 重算**：`.skills/cath-ted-domain-annotation/scripts/cath_ted_holdout.py`
   （`--split-col splits`；参考 = 新 train 29,622 行的 1,189 topos /
   4,214 sfs）→ 写回 `splits/ssp_test.csv`，Default=False 行两列置 NaN

**口径**（链模式，`unique_id=PDB_chain`）：

- 标签 = CATH v4.4 实验标注（名称匹配 `cath-domain-list.txt`，域区间经 SIFTS
  映射回链区间）∪ TED putative 标注（REST API 缓存；仅保留与链 UniProt 区间
  重叠 ≥50% 域长的域；H 级提供 C.A.T.H，T 级仅 C.A.T）
- 参考集 = **Train** 的 CATH+TED 并集（train/val/test 对称增强）
- any 语义：链任一 Topology/Superfamily 未在参考集出现 → True；无注释 →
  False；`Default=False` 行置 NaN

ss 数据源自全量 PDB 过滤去冗余，CATH 实验标注覆盖率仅 58.5%，TED 补缺收益
显著。标注来源分布（41,311 条链）：both 19,994 / ted-only 10,366 /
cath-only 3,852 / none 7,099。

**标注覆盖率（全部 41,311 条链）**：仅 CATH 58.5% → CATH+TED **82.8%**。

**Holdout 结果（Default=True n=8,290）**：

| 列 | True 数 |
|------|--------:|
| `OOD_FoldHoldout` | 20 |
| `OOD_SuperfamilyHoldout` | 96 |
| 无注释样本（→ False） | 3,511 |

#### finalize 校验与体检（skill）

1. **NaN 惯例与兜底校验**：Default=False（极长链）行除 `OOD_ExtremeLong` 外
   全部 OOD 标记列（含 `seq_Redundancy_*`）为 NaN、`InD=False`。2026-08-28
   起 NaN 由各生成脚本（`09`/`10`/`11` 及 skill 脚本）在流水线位置直接写入；
   `finalize_ood_columns.py --mode nan` 仅作兜底校验（当前数据上为 no-op）。
   同日已将 82 条极长链历史遗留的 `seq_Redundancy_*` 计算值（True）统一为
   NaN（574 单元格，备份 `output/backup_nondefault_seqtm_nan/`）
2. **InD 重算**：`finalize_ood_columns.py --mode ind`——所有 OOD 列齐备后按
   全并集重算 `InD`（当前 **1,123 / 8,290，13.5%**，2026-09-02 短链入池后；
   历史值见 4.4）
3. **体检**：`.skills/ood-annotation-toolkit/scripts/check_splits.py`

#### `ssp_test.csv` 列顺序（当前实际表头，34 列）

```
unique_id, aa_seq, struct_file, label, Default,
seq_Redundancy_90, seq_Redundancy_80, ..., seq_Redundancy_30,
TM-score_0.9, TM-score_0.8, TM-score_0.7, TM-score_0.6, TM-score_0.5, TM-score_0.4, TM-score_0.3,
OOD_Orphan, OOD_ExtremeShort, OOD_ExtremeLong, OOD_IDR, OOD_LongTail,
OOD_NewEC_L4, OOD_NewEC_L3,
OOD_FoldHoldout, OOD_SuperfamilyHoldout,
OOD_LongTail_EC_L3_le5, OOD_LongTail_EC_L3_le10,
OOD_LongTail_EC_L4_le5, OOD_LongTail_EC_L4_le10,
InD, struct_label
```

### 4.3 struct_label（结构残基级标注，skill 落盘）

train/val/test 三个文件均含 `struct_label` 列（末列）：对 `struct_file` 中的
**每个结构残基**（ATOM 记录按文件顺序、以 `(chain, resseq, icode)` 首次出现
枚举，含无 CA 原子的残基）给出对应的 DSSP 标注。结构模型只需解析 pdb + 读
`struct_label` 即可，二者等长、顺序一致，无需自行比对。

**落盘脚本**：`.skills/struct-label-alignment/scripts/add_struct_label.py`。
标注来源为已归档 struct_native 视图的 DSSP 直达标注（权威来源；视图文件现
归档于 `output/splits_struct_native_deprecated/`）——等长行（残基枚举一致，
约 98.7%）逐值照搬；含无 CA 残基的行（train 384 / val 12 / test 44）将全
ATOM 残基序列比对到 native 的 CA 序列继承标注，无 CA 残基置 **`-1`**。
`-1` 语义 = ignore index，训练/评估时跳过。备份：
`output/backup_before_struct_label_native/`。

- 验证：其它列逐单元格零变化；等长行与 native label **100% 一致**；
  `len(struct_label) == pdb 残基数` 抽样断言通过
- 含 `-1` 残基总数：train 871，val 62，test 176
- 结构模型统一使用主 `splits/` + `struct_label`（struct_native 视图因此冗余
  并已归档，见注意事项 11）

> 历史注记：`struct_label` 曾由本地脚本 `patch_struct_label_from_native.py`
> 移植写入（修复已落入 splits，脚本已归档）；现行重算路径统一走 skill
> `add_struct_label.py`。

### 4.4 数据质量检查

| 检查项 | 结果 |
|--------|------|
| Train ∩ Val / Train ∩ Test / Val ∩ Test（`unique_id`） | 均为 0（已核验） |
| `len(label) == len(aa_seq)` | 划分脚本按 `ssp_full_data.csv` 逐链生成（07 聚合时逐行断言） |
| Test OOD 列完整性 | 8,290 条 Default 样本（含 293 条 ExtremeShort 短链）所有 OOD 列均有 True/False 值；82 条极长链除 `OOD_ExtremeLong` 外统一 NaN（惯例） |
| `InD` 列重算 | 当前 = **1,123**（2026-09-02，与文件实值完全一致）；历史：1,135（2026-08-30 序列侧统一）→ 1,179（短链按 cutoff 划入）→ 1,123（≤30aa 短链删除后） |
| OOD 列嵌套性 | `OOD_Orphan ⊆ seq_Redundancy_30 ⊆ ... ⊆ seq_Redundancy_90`（计数单调，已核验） |

**质控记录（数据修正事件）**：

- **7UYL_K 剔除（2026-09-02）**：该短链 120 个结构残基中 71 个为 UNK
  （59% 未鉴定），序列提取按 `AA_MAP` 跳过 UNK（aa_seq 仅 49 残基）而结构
  提取保留 UNK，导致 label（49）与 struct_label（120）契约破坏、struct_label
  71 个位置错帧（对 DSSP 直接标注不一致）。已从 test 删除该行并移除
  `pdbs/7UYL_K.pdb`（备份 `output/backup_20260902/ssp_test.before_7uylk_removal.csv`）；
  现行 08 以 `QC_EXCLUDE` 常量钉死该剔除。短链中另有 7 条结构残基数 ≠
  aa_seq 长度（差 1–3 个无 CA 残基），经核对 struct_label 与 DSSP 完全一致，
  不受影响。
- **struct_label 恢复事件（2026-09-01/02）**：短链恢复流程中
  `add_struct_label.py --overwrite` 全表重算曾使 148 行（均为
  `aa_seq` 全长序列 ≠ 结构序列的行）的 struct_label 被改——这些行的新值
  是"结构序列比对到 aa_seq 再继承主 label"的串联比对结果，在低复杂度区
  发生错帧。以 DSSP 直接标注为裁判核对（可核对 125 行）：旧值一致率
  99.996%（46,206 残基仅 2 个差异），新值仅 71%。**该 148 行已从
  `output/backup_20260901/ssp_test.csv` 恢复为 struct_native 移植的旧值**
  （恢复前状态备份为 `output/backup_20260901/ssp_test.before_struct_label_restore.csv`）；
  短链新行不受影响（其 `aa_seq` 即 CA 提取序列，与结构一致，走直通路径，
  抽查 500 条与 DSSP 一致率 99.97%）。**教训：`aa_seq ≠ 结构序列` 的行不可用
  `add_struct_label.py --overwrite` 重算 struct_label。**

---

## 诊断与质控工具

- `.skills/ood-annotation-toolkit/scripts/check_splits.py`：splits 完整体检。
- `output/verify/`：`08_create_ssp_datasets.py --out-dir output/verify` 试跑
  产物（v1 口径 72:8:20 分位切分时点与 splits 逐行一致性的历史验证；
  现行 08 恢复 72:8:20 动态分位口径，该目录对应历史时点，仅供溯源）。
- 诊断脚本 `ss_struct_alignment_probe.py`（结构残基与 label 非 -1 位置对齐
  探查）已归档，见归档清单。

## 归档清单（`output/scripts_archive/`）

| 脚本 | 归档原因 |
|------|----------|
| `extract_subset.py` | 与 `05_extract_representative_chains.py` 冗余（断点续传期子集提取），且带**未修复的多模型 bug**（`io.set_structure(structure)` 会写出全部模型）；请勿复用，多模型正确口径见 05 |
| `run_ec_newcategory.py` | OOD-NewEC 初版（any 语义、按日期截止推导划分）；已被 `.skills/ec-function-ood-annotation/scripts/add_unified_newec.py` 统一口径取代（superseded） |
| `build_ss_struct_native.py` | struct_native 视图构建脚本；视图已弃用归档（`output/splits_struct_native_deprecated/`），功能由主 splits `struct_label` 列取代 |
| `patch_struct_label_from_native.py` | struct_label 一次性移植补丁；修复已落入 splits，重算改走 `.skills/struct-label-alignment/scripts/add_struct_label.py` |
| `ss_struct_alignment_probe.py` | 纯诊断探查脚本（结构-标签对齐核验），不进主流程 |
| `fix_multimodel_pdbs.py` | 一次性修复：26 个多模型（NMR）pdb 仅保留 model 1；修复已落数据，根因已修入 `05_extract_representative_chains.py`（`structure[0]`） |
| `add_longtail_ec_ood.py` | 初版单列 `OOD_LongTail_EC_L4` 补丁；已被 ec skill 统一四列取代 |
| `extract_pdb_metadata.py`、`extract_representative_structures.py`、`run_dssp_batch.py` | 早期管线遗留版本，已被 01/05/06 取代 |
| `12_extract_short_chains.py`、`13_cluster_short_chains.py`、`14_extract_short_structures.py`、`15_build_short_test_rows.py` | 短链恢复历史路径（2026-09-01）：<60 短链提取 / 聚类交叉去重 / 结构提取 / DSSP+追加 test；功能已分别线性化融入 03（MIN_RESIDUES=31）、04（单阶段全量聚类）、05/06（对代表清单通用），见历史沿革 |
| `16_redistribute_short_chains.py`、`17_apply_ood_updates.py` | 短链按原 cutoff 正常划分 + OOD 重算（2026-09-02）；划分逻辑已融入 08（72:8:20 分位切分 + QC_EXCLUDE） |
| `18_drop_short_chains_le30.py`、`19_apply_ood_updates_v3.py` | ≤30aa 短链删除 + 联动 OOD 重算（2026-09-02）；删除口径已钉入 03 的 MIN_RESIDUES=31 |

## 注意事项

1. **API 限流**：RCSB GraphQL API 有速率限制，脚本已内置退避逻辑。若频繁
   429，请增大 `REQUEST_DELAY`。
2. **断点续传**：`01_fetch_pdb_metadata_api.py` 支持断点续传，kill 后重新
   运行会自动从 checkpoint 继续。
3. **R-free 缺失**：部分早期 PDB entry 未记录 R-free（如 1UBQ），这些会在
   筛选时被排除。
4. **文件格式**：优先 `.pdb.gz`，缺失时回退 `.cif.gz`。mmCIF 中 chain ID 可能
   是多字符，提取时截断为单字符。
5. **DSSP v4 兼容性**：需确保 PDB 文件包含 `HEADER`、`EXPDTA`、`CRYST1`、`END`
   记录，且无 `ANISOU` 位置错误。BioPython `PDBIO` 提取可有效避免这些问题。
6. **极长链 83 → 82**：阶段一/二中长度 >1000 的代表链共 83 条，但 `4O9X_A`
   （2,082 残基）按"单链 >2000 一律剔除"规则在划分时被移除，最终 test 极长链
   为 **82 条**、test 总量 **8,372**。旧文档中"83 条 / 8,080 行"为剔除规则
   生效前的数字；短链历史路径中 test 总量曾为 17,020 / 10,456 / 10,455
   （见历史沿革）。
7. **`is_time_cutoff` 列已移除**：早期版本的辅助列，现行表头不含此列，等价
   信息即 `Default`。
8. **OOD-NewEC 口径演进**：初版（`run_ec_newcategory.py`，any 语义、按日期
   截止推导划分）结果为 L4=287 / L3=5；统一为 all-not-in 口径并以实际 splits
   为准后重算为 **L4=246 / L3=1**（`output/backup_before_newec_unify/`）。
9. **OOD-LongTailEC 口径演进**：早期单列 `OOD_LongTail_EC_L4`（初版脚本已
   归档）已被统一四列 `OOD_LongTail_EC_{L3,L4}_le{5,10}` 取代
   （`output/backup_before_longtail_unify/`）。
10. **TM-score 列方向修正**：早期 TM-score 列曾存在判定方向错误，已修正
    （备份 `output/ssp_test.before_tm_inversion.csv`）。此后 2026-08-28 又发现
    四档列实为旧式 `search -e inf` 口径且含非法 >1 值，已按 skill easy-search
    口径整体重算并补齐 0.5/0.4/0.3 档（见阶段四"结构同源 OOD"小节，
    备份 `output/backup_before_foldseek_easysearch_unify/`）。
11. **struct_native 视图已弃用归档**（2026-08-28）：`output/splits_struct_native_deprecated/`
    目录保留（analysis 脚本仍引用）。其 OOD 列止步于 `OOD_NewEC_L4/L3, InD`，
    缺少主 splits 后续新增的 FoldHoldout/SuperfamilyHoldout/LongTailEC 列，
    其 `InD`（871）也未按最新 OOD 列重算（主 splits 为 799）——仅作历史
    记录，勿用于 OOD 评估。功能由主 splits `struct_label` 列（见 4.3）取代。
12. **多模型（NMR）修复（归档注记）**：2026-08-28 全库扫描发现 26 个 pdb
    文件含多模型内容（23 个带 `MODEL`/`ENDMDL`，3 个模型串联无记录）。根因
    已修入 05（仅保存 model 1）；存量文件由一次性脚本修复（备份
    `output/backup_before_multimodel_fix/`），修复脚本已归档。修复后
    `struct_label` 已用 `--overwrite` 重算（残基枚举按 key 去重，标注值不受影响）。
13. **历史偏差记录（本次整理发现）**：现存 `output/ssp_full_data.csv` 的
    `file_type` 列有 **787 行误记为 `pdb`**（这些 entry 在 `data/pdb/` 仅有
    `.cif.gz` 源，`stage2_chain_metadata.csv` 记录为 `cif`）。该列下游
    （08 划分脚本）不使用，不影响 splits；重建脚本
    `07_aggregate_ssp_full_data.py` 按代表链元数据产出正确值，与现存文件
    diff 仅此 787 行差异（`unique_id/aa_seq/labels` 40,067 行逐行一致）。
14. **无脚本环节**：MMseqs2 去冗余（04）之前无额外手工步骤；`ssp_full_data.csv`
    与 `ss_ted_aug_chain_list.csv` 原为无脚本中间产物，本次已分别补建
    07 与并入 11，无遗留断链环节。

## 关键文件索引

| 文件路径 | 阶段 | 说明 |
|---------|------|------|
| `scripts/01_fetch_pdb_metadata_api.py` | 阶段一 | GraphQL 元数据获取 |
| `scripts/02_filter_entries.py` | 阶段一 | 质量筛选 |
| `scripts/03_extract_chain_sequences.py` | 阶段一 | 链序列提取 |
| `scripts/04_run_mmseqs2_clustering.py` | 阶段二 | MMseqs2 去冗余 |
| `scripts/05_extract_representative_chains.py` | 阶段三 | 单链 PDB 提取 |
| `scripts/06_run_dssp_parallel.py` | 阶段三 | DSSP 批量标注 |
| `scripts/07_aggregate_ssp_full_data.py` | 阶段三 | 全量数据聚合（ssp_full_data.csv） |
| `scripts/08_create_ssp_datasets.py` | 阶段四 | 数据集划分 |
| `scripts/09_run_idr_prediction.py` | 阶段四 | metapredict IDR 预测与 OOD-IDR 标记 |
| `scripts/10_run_longtail_prediction.py` | 阶段四 | 8-state LongTail 计算与 OOD-LongTail 标记 |
| `scripts/11_ted_augment_cath_ood.py` | 阶段四 | 链清单生成 + CATH+TED 域标注与 Fold/SuperfamilyHoldout 驱动 |
| `.skills/homology-ood-annotation/scripts/seq_homology_ood.py` | 阶段四 | seq_Redundancy_* / OOD_Orphan 统一口径 |
| `.skills/homology-ood-annotation/scripts/struct_homology_ood.py` | 阶段四 | TM-score_* 统一口径 |
| `.skills/ec-function-ood-annotation/scripts/add_unified_newec.py` | 阶段四 | OOD_NewEC_L3/L4 统一口径 |
| `.skills/ec-function-ood-annotation/scripts/add_unified_longtail_ec.py` | 阶段四 | OOD_LongTail_EC 四列统一口径 |
| `.skills/cath-ted-domain-annotation/scripts/cath_ted_annotate.py` | 阶段四 | CATH+TED 链级标注 |
| `.skills/cath-ted-domain-annotation/scripts/cath_ted_holdout.py` | 阶段四 | Fold/SuperfamilyHoldout 重算 |
| `.skills/ood-annotation-toolkit/scripts/idr_ood.py` | 阶段四 | OOD_IDR 统一工具（residue-ratio 口径，与 09 一致） |
| `.skills/ood-annotation-toolkit/scripts/finalize_ood_columns.py` | 阶段四 | Default=False 行 NaN 收尾 + InD 重算 |
| `.skills/struct-label-alignment/scripts/add_struct_label.py` | 阶段四 | struct_label 落盘（重算路径） |
| `output/pdb_metadata_api.csv` | 阶段一 | 全量元数据（254,975 行） |
| `output/stage1_filtered_entries.csv` | 阶段一 | 筛选后 entries（121,160 行） |
| `output/stage2_chain_sequences.fa` | 阶段一 | 全量链 FASTA（磁盘现存为 v1 的 249,018 条 ≥60 链；按现行 03 重跑为 255,471 条，含 31–59 短链 6,453） |
| `output/stage4_cluster_representatives.csv` | 阶段二 | 代表链列表（磁盘现存为 v1 的 40,067 行主池；按现行 04 重跑为 41,337 行，含短链 1,270） |
| `pdbs/` | 阶段三 | 单链 PDB 文件（41,327 个 = 主池 40,067 + 31–59 短链 1,260） |
| `output/dssp_labels.csv` | 阶段三 | DSSP 标签汇总（40,067 行主池；短链 dssp 在 `output/dssp_short/`，1,260 个） |
| `output/ssp_full_data.csv` | 阶段三 | 全量数据（磁盘现存为 v1 的 40,067 行，含已剔除的 4O9X_A；重跑为 41,313 行） |
| `output/idr_predictions.csv` | 阶段四 | IDR 预测结果（40,067 行主池；短链缓存 `output/idr_cache_short.csv`） |
| `output/longtail_predictions.csv` | 阶段四 | LongTail 计算结果（40,067 行主池历史产物；当前阈值见 `output/longtail_threshold_v3.json`） |
| `output/ss_chain_cath_ted_labels.csv` | 阶段四 | 链级 CATH+TED 标签（41,311 行） |
| `output/ss_chain_domains.csv` | 阶段四 | 域级明细（94,998 行） |
| `output/ss_ted_aug_chain_list.csv` | 阶段四 | 链清单（41,311 行，11 生成） |
| `splits/ssp_train.csv` | 阶段四 | 训练集（29,622 行 = 28,788 主池 + 31–59aa 短链 834） |
| `splits/ssp_val.csv` | 阶段四 | 验证集（3,317 行 = 3,199 主池 + 短链 118） |
| `splits/ssp_test.csv` | 阶段四 | 测试集（8,372 行 = 7,997 主池 + 293 短链 + 82 极长链，含完整 OOD 标记 + InD + struct_label） |
| `output/short_chains.fa` / `output/short_chains_meta.csv` | 历史 | 短链恢复路径的短链序列与元数据（1–59aa 22,231 条；其中 31–59 段 6,453 条入池） |
| `output/short_chain_representatives.csv` | 历史 | 短链幸存代表链（9,079 行，含 ≤30 段；31–59 段 1,270 条并入代表清单） |
| `output/short_clustering_stats.json` | 历史 | 短链聚类与交叉去重统计（1–59 池：9,229 簇 → 丢 150 → 9,079） |
| `output/short_dssp_failures.csv` | 历史 | 短链 DSSP/结构失败清单（138 行） |
| `output/backup_20260901/ssp_test.csv` | 历史 | 短链恢复路径中 ExtremeShort 追加前的 test 备份（8,079 行） |
| `output/backup_20260902/` | 历史 | 短链重新划分前的 splits 三件套备份（test 17,020 行） |
| `output/backup_20260902_min30/` | 历史 | ≤30aa 删除前的 splits 三件套备份（test 10,455 行） |
| `output/ood_search/mmseqs_easysearch_test_vs_trainval_v2.m8` | 历史 | 短链划分后 seq 同源搜索缓存（现行缓存为 v3） |
| `output/ood_search/foldseek_easysearch_test_vs_trainval_v2.m8` | 历史 | 短链划分后结构同源搜索缓存（现行缓存为 v3） |
| `output/ood_search/mmseqs_easysearch_test_vs_trainval_v3.m8` | 阶段四 | 现行 seq 同源搜索缓存（含短链参考集，min30 后） |
| `output/ood_search/foldseek_easysearch_test_vs_trainval_v3.m8` | 阶段四 | 现行结构同源搜索缓存（含短链参考集，min30 后） |
| `output/idr_cache_short.csv` / `output/short_idr_flags.csv` | 历史 | 短链 IDR 计算缓存与判定（8,941 条，含已删除的 ≤30 段） |
| `output/longtail_threshold_v2.json` / `output/longtail_threshold_v3.json` | 历史/阶段四 | LongTail 阈值记录（0.081181 → 0.093333 → 0.082914，v3 为现行） |
| `output/splits_struct_native_deprecated/` | 历史 | 结构原生 split（已弃用归档；功能由主 splits `struct_label` 列取代） |

## 历史沿革

短链（ExtremeShort）处理的实际演进路径（现行文档已将其线性化融入
阶段一/二/四，下列脚本均归档于 `output/scripts_archive/`）：

1. **2026-09-01 短链恢复**（12–15）：03 曾以 `MIN_RESIDUES=60` 过滤短链；
   恢复流程从原始 PDB 重提取 <60 链 22,231 条（1–59aa）→ easy-linclust
   9,229 簇 → 代表 9,229 → 对现有 splits 交叉去重丢 150 → 幸存 9,079 →
   结构提取失败 123 + DSSP 失败 15 → 8,941 条以 `Default=False` /
   `OOD_ExtremeShort=True` 全部追加进 test（17,020 行中间态）。
2. **2026-09-02 按原 cutoff 正常划分**（16/17）：重建 08 的日期 cutoff
   （train ≤ 2019-06-06、val ≤ 2021-02-09），短链分流 train 5,672 /
   val 892 / test 2,377 并按扩大后的参考集重算全部 OOD 标注；同日 QC
   剔除 7UYL_K（见 4.4 质控记录）。
3. **2026-09-02 ≤30aa 删除**（18/19）：短链口径收紧为 31–59aa，删除
   ≤30 共 7,695 条（train 4,838 / val 774 / test 2,083）及孤立结构，
   OOD 再次联动重算，达当前状态（41,311 = 29,622 / 3,317 / 8,372）。
4. **2026-09-02 文档线性化**：短链处理融入 03（`MIN_RESIDUES=31`）、
   04（单阶段全量聚类）、08（72:8:20 分位切分 + `QC_EXCLUDE`），
   12–19 归档。

**等价性说明**：当前数据与从头按本文档重跑的差异有三：(a) easy-linclust
边界 tie-break（linclust 对输入顺序敏感，边界簇成员归属可能有 ±少量
差异）；(b) **cutoff 随池组成移动**（重跑对含短链的全池重新分位，cutoff
日期与本次发布值 2019-06-06 / 2021-02-09 不同，边界附近的链归属有少量
差异，见 4.1）；(c) 行序（历史路径将短链附于 train/val 末尾，重跑按
日期排序插入对应位置）。这些偏差不影响数据质量，属允许范围。
此外磁盘上部分中间产物（stage2/stage4/ssp_full_data 等）仍为 v1 版本
（行数以各节注记为准，重跑即更新）。
