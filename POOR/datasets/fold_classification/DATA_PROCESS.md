# Fold Classification 数据集处理文档

> 记录从 CATH S95 列表与本地 PDB 结构数据到最终 OOD 划分数据集的完整处理流程，
> 按"从原始数据一次性复现"的线性顺序组织。
> 处理位置：`datasets/fold_classification/`（折叠分类任务主目录）
>
> **任务定义**：CATH v4.4 S95 非冗余域的折叠分类，`label` = C.A.T（Topology）三级
> 分类号（如 `3.40.50`），共 **1,163** 类；`unique_id` 为 CATH 域 ID（如 `2da1A01`）。

**当前 splits 规模**（以 `splits/` 实际文件统计为准）：

| 集合 | 行数 | 说明 |
|------|-----:|------|
| `splits/cath_train.csv` | **43,473** | `unique_id,aa_seq,struct_file,label`，含全部 1,163 个 Topology |
| `splits/cath_val.csv` | **6,132** | 同上 |
| `splits/cath_test.csv` | **12,911** | 核心 4 列 + 27 个标注列（31 列，含完整 OOD 标记与 `InD`） |

---

## 阶段一：源数据获取与解析（CATH S95 结构域提取）

> 目标：从 CATH v4.4 S95 列表和本地 PDB 实验结构中提取全部 62,915 个结构域的
> ATOM 坐标，生成逐域 PDB 文件与元数据汇总表。

### 1.1 数据来源

| 数据 | 位置 | 用途 |
|------|------|------|
| `cath-classification-data/cath-domain-list-S95.txt` | `data/CATHv44/` | **核心入口**：S95 非冗余域列表（62,915 条） |
| `cath-classification-data/cath-domain-boundaries.txt` | `data/CATHv44/` | 结构域边界（PDB 编号，Domall CDF 2.0） |
| `cath-classification-data/cath-domain-boundaries-seqreschopping.txt` | `data/CATHv44/` | 备选边界（SEQRES 连续编号） |
| `sequence-data/cath-domain-seqs-S95.fa` | `data/CATHv44/` | S95 域参考序列（FASTA），划分池 `aa_seq` 的来源 |
| PDB 实验结构（`.pdb.gz` 为主，`.cif.gz` 补全） | `data/pdb/` | 结构源，下载自 RCSB（`https://files.rcsb.org/download/`） |

S95 列表每行 12 列：`domain_id`（7 字符：`PDBID`+`Chain`+`DomainNum`）、`C.A.T.H`
四级分类号、`S35/S60/S95/S100` 序列聚类编号、`count`、`length`、`resolution`
（999=NMR，1000=obsolete）。注意 **S95 编号不是全局唯一的**，仅在 `(S35, S60)`
上下文中唯一，唯一标识一个 S95 聚类需完整键 `(C, A, T, H, S35, S60, S95)`。

### 1.2 缺失 PDB 补全 —— `scripts/01_download_missing_for_s95.py`

- **输入**：S95 列表 + `data/pdb/` 本地索引（脚本内自动对比，生成
  `output/missing_for_s95.txt`，无手工步骤）
- **输出**：补齐的 `data/pdb/*.pdb.gz|.cif.gz`；日志 `output/download_missing_s95.log`；
  失败列表 `output/download_failed_s95.txt`（如有）
- **说明**：10 线程并发，优先 `.pdb.gz`，404 回退 `.cif.gz`。
  实际执行：S95 列表 41,403 个唯一 PDB ID，缺失 246 个，**246 / 246 全部
  成功下载**，覆盖率 100%。

### 1.3 全量结构域提取 + 汇总表 —— `scripts/02_extract_cath_s95_domains.py`

- **输入**：S95 列表、域边界文件、`data/pdb/`
- **输出**：`output/cath_s95_domains/{domain_id}.pdb`（后整体迁移为 `pdbs/`）；
  `output/cath_s95_summary.csv`（**62,901 行**，字段 `domain_id, pdb_id, chain,
  dom_num, type, C, A, T, H, S35, S60, S95, S100, length_expected, n_residues,
  n_atoms, resolution`）；断点 `output/extract_checkpoint.pkl`
- **说明**：
  - `ProcessPoolExecutor` 多进程（`cpu_count()-1` workers），按域类型选择策略：
    whole-chain（`dom_num=00`）直接提取整条链；chopped 域按
    `cath-domain-boundaries.txt` 的 segment 范围过滤 ATOM 记录
  - `.cif.gz` 提取后统一转为 PDB 格式输出，原子序号重排；NMR 多模型条目在首个
    `ENDMDL` 即停（仅保留 model 1）
  - 每 500 个域保存 checkpoint，支持断点续传
  - **汇总表在提取结束后由脚本从产物文件全量重建**（逐文件统计
    `n_residues`/`n_atoms` 并与 S95 元数据合并），一次性产出完整
    `cath_s95_summary.csv`。早期版本在并行循环内按 checkpoint 追加写 CSV 有缺陷
    （仅保存了部分批次），曾由一次性补丁 `fix_summary.py` 重建，该逻辑已融入
    本脚本，补丁已归档（见归档清单）
  - 结果：success **62,901** / empty（无原子坐标）14；耗时约 15–20 分钟

### 1.4 结构库落位 `pdbs/`

提取产物 62,901 个域 PDB 文件整体迁移为任务最终结构库 `pdbs/`
（`struct_file = {domain_id}.pdb` 即对应此目录）。后续 foldseek 结构搜索与
下游模型训练均直接引用 `pdbs/`。

### 1.5 CATH 层级分布统计（基于 62,901 个提取域）

| 层级 | 类别数 | 说明 |
|------|------|------|
| Class (C) | 5 | 二级结构组成大类 |
| Architecture (C.A) | 43 | 二级结构大尺度排列 |
| **Topology (C.A.T)** | **1,472** | 二级结构连接方式（**本任务分类粒度**） |
| Superfamily (C.A.T.H) | 6,631 | 进化同源超家族（过细，3,207 个仅 1 域，未采用） |

Class 分布：1 Mainly Alpha 13,417（21.3%）、2 Mainly Beta 16,227（25.8%）、
3 Alpha Beta 31,613（**50.3%**）、4 Few Secondary Structures 675（1.1%）、
6 Special 969（1.5%）。最大超家族为 `2.60.40.10`（4,542 域），分布高度长尾。

---

## 阶段二：标签聚合与计算（样本池构建）

> **脚本**：`scripts/03_build_domain_pool.py`
>
> **分类标签**：`label = "C.A.T"` 字符串（如 `3.40.50`）。Superfamily（6,631 类，
> 81% 类别样本过少）过细、Architecture（43 类）过粗，均不采用。

- **输入**：`output/cath_s95_summary.csv`（域清单与 C.A.T.H 分类号）；
  `data/CATHv44/sequence-data/cath-domain-seqs-S95.fa`（`aa_seq` 来源）
- **输出**：`output/domain_pool.csv`（**62,896 行**，列
  `unique_id, aa_seq, struct_file`）
- **说明**：
  - `aa_seq` 取 CATH 域参考序列（域级权威序列；PDB ATOM 记录常缺失 loop
    残基，不作序列源），`struct_file = {domain_id}.pdb` 并断言存在于 `pdbs/`
  - 剔除 `length_expected > 1000` 的超长域 **5 个**（1ej6B00、1u6gC00、
    3a6pA00、3egwA02、3gjxD00；按 CATH 列表长度判定）
  - **2026-09-01 补剔**：3w3uA01（`length_expected=999` 未触发上述剔除，但实际
    参考序列 1,047 残基）按实际序列长度口径从 Train 补删（Train 43,474 → 43,473；
    其 Topology 1.25.10 在 Train 仍有 91 条代表），其 `pdbs/3w3uA01.pdb` 一并移除
  - 该池与 v1（C.A 粒度）划分备份 `output/splits_ca_backup/` 三文件并集
    **行级等价**（已 diff 验证一致）；v1 备份保留作历史对照与回退

---

## 阶段三：去冗余与数据精化

- **序列去冗余**：由数据源本身完成——S95 列表即 CATH 官方按 95% 序列同一性
  聚类的非冗余代表域集，本任务不再额外聚类。
- **划分前过滤**（在 `04_build_topology_split.py` 内执行，属数据精化）：

| 操作 | 数量 | 说明 |
|------|------|------|
| 删除单例 Topology | **302** 个 Topology（302 域） | 仅 1 个样本，无法同时进 Train 和 Test |
| EC7 先验预留 | **413** 域 | 携带 EC7（Translocases）L1 大类的域**全部移出 Train/Val 池**，直入 Test（多 EC 域按"携带即中"） |
| 删除纯 EC7 零池 Topology | 6 个 Topology（**77** 域） | 全部成员均为 EC7 先验，永不出现于 Train，删除 |
| 保留的 EC7 先验 | **336** 域 | 413 − 77，全部进入 Test（占 Test **2.6%**） |

被删除的 379 个域清单存于 `output/dropped_domains.csv`（含 topology、
superfamily、is_newec 标记）。

**选取 EC7 作为 NewEC 先验的理由**（各 EC 大类 holdout 对 Test 分布的冲击对比）：

| EC holdout | 域数 | 涉及 Topology | 占 Test 比例 |
|------------|-----:|--------------:|:----------:|
| **EC7（采用）** | **413** | **40** | **3.2%** |
| EC6 | 1,584 | 93 | 11.5% |
| EC5 | 1,818 | 125 | 13.0% |
| EC1/2/3 | 5,364–9,595 | 215–389 | 31.9–47.5% |

EC1/2/3 是酶的主体，holdout 会把近半数据集推进 Test；EC7 冲击最小，且为
2018 年新增大类，与经典代谢酶机制差异明显，适合"新功能泛化"评估。

> **口径说明**：本任务的 `OOD_NewEC` 是**划分前预留的 EC7（L1）先验**，
> 与 `.skills/ec-function-ood-annotation` 的 NewEC 口径（all-not-in、train+val
> 参考、L3/L4 层级）**不同，不适用该 skill 重算**（skill 注意事项亦明确
> fold 为例外）。

**多模型（NMR）文件修复（归档注记）**：2026-08-28 全库扫描发现 4,076 / 62,901
个 pdb 文件将全部模型无 `MODEL` 分隔地串联写入。根因已修入
`02_extract_cath_s95_domains.py`（首个 `ENDMDL` 即停）；存量文件由一次性脚本
截断修复（保留 model 1，原件备份于 `output/backup_before_multimodel_fix/`），
修复脚本已归档（见归档清单）。修复后全库重扫 0 残留重复。

---

## 阶段四：数据划分与 OOD 标注

### 4.1 Topology 分层 7:1:2 划分 —— `scripts/04_build_topology_split.py`（seed=42）

- **输入**：`output/domain_pool.csv`（默认，03 产物；缺失时回退
  `output/splits_ca_backup/` v1 池，二者行级等价）；`output/cath_s95_summary.csv`
  （C.A.T.H 分类号来源）；`data/sifts/sifts_ec_annotations.json`（SIFTS 链级
  EC 注释，NewEC 先验依据）
- **输出**：`splits/cath_{train,val,test}.csv`（test 携带 `OOD_NewEC` 列，其余
  OOD 列由 4.2 追加）；`output/dropped_domains.csv`
- **说明**：样本过滤按阶段三表格执行；非 EC7 池内样本按 **Topology 分层、
  7:1:2** 随机划分（seed=42）；池内仅 1 个样本的 Topology 强制分入 Train
  （其 Test 覆盖由先验样本保证）；EC7 先验样本全部并入 Test。

| 集合 | Domain 数 | 说明 |
|------|----------:|------|
| **Train** | **43,473** | 含全部 **1,163** 个 Topology（train 类别频次 min/median/max = 1/6/5,584） |
| **Val** | **6,132** | — |
| **Test** | **12,911** | 含 336 个 EC7 先验样本（`OOD_NewEC=True`） |
| 合计 | **62,516** | = 62,896（池）− 379（过滤）− 1（3w3uA01 补剔） |

**脚本内置断言校验**：三方 `unique_id` 两两无交集；全部 1,163 个 Topology
同时出现在 Train 和 Test；`label` 全部为 `C.A.T` 字符串格式；Train/Val 中
EC7 样本数为 0。

> `label` 列强制保存为字符串类型；加载时须显式 `dtype={'label': str}`，
> 避免 `1.10` 被解析为浮点导致类别合并。

### 4.2 Test 集 OOD 标注 —— `scripts/05_annotate_test_ood.py`（一站式）

> **总口径**：**全部 12,911 行完整参与所有依赖型 OOD 计算**（含 336 个
> NewEC 先验行与 767 个 ExtremeShort 行）；`Default` 全表为 True（本任务无
> "划分前单独预留"的长度样本，短序列与所有样本一样走正常分层划分）。
> 因此 Default=False 行的 NaN 惯例在本任务无适用对象。

各轴分述如下（统计均以当前 `splits/cath_test.csv` 实值为准）：

#### 序列同源 OOD：`seq_Redundancy_90`~`30` + `OOD_Orphan`

- **工具**：mmseqs2 `easy-search -s 7.5 --format-output query,target,pident`
  （easy-search 默认 e-value ≤ 1e-3 显著性门槛）
- **方向**：query = Test 全表，target = Train+Val（49,606 条）
- **判定**：每个 Test 样本取 max pident（无 hit 按 0.0 计）；
  `seq_Redundancy_<T>` = max pident < T（T ∈ 90,80,…,30）；
  `OOD_Orphan` = 无任何显著 hit
- **缓存**：`output/ood_mmseqs_test_vs_trainval.m8`（重跑时自动复用）
- **方法学**：与 `.skills/homology-ood-annotation`（`seq_homology_ood.py`）
  口径一致，该 skill 以本脚本为基准；重放重算可用其 `--m8-cache` 直接复现

| 阈值 | Domain 数 | 占 Test |
|------|----------:|--------:|
| < 90% | 10,257 | 79.4% |
| < 80% | 9,118 | 70.6% |
| < 70% | 8,114 | 62.8% |
| < 60% | 7,055 | 54.6% |
| < 50% | 5,586 | 43.3% |
| < 40% | 3,876 | 30.0% |
| **< 30%** | **2,353** | **18.2%** |

**OOD_Orphan**：**1,978** domains（15.3%）。子集关系 `OOD_Orphan ⊆
seq_Redundancy_30 ⊆ … ⊆ seq_Redundancy_90` 成立（check_splits 体检确认）。

#### 结构同源 OOD：`TM-score_0.9`~`0.3`

- **工具**：foldseek `easy-search --format-output query,target,alntmscore`；
  query/target 目录用符号链接组建（避免拷贝 ~10 GB 结构），query id 取文件名
  去路径与 `.pdb` 后缀
- **判定**：每个 Test 结构取 max alntmscore（无 hit 按 0.0 计）；
  `TM-score_<T>` = max alntmscore < T（T ∈ 0.9,…,0.3）
- **缓存**：`output/ood_foldseek_test_vs_trainval.m8`
- **方法学**：同 `.skills/homology-ood-annotation`（`struct_homology_ood.py`）

| 阈值 | Domain 数 | 占 Test |
|------|----------:|--------:|
| < 0.9 | 4,022 | 31.2% |
| < 0.8 | 1,755 | 13.6% |
| < 0.7 | 839 | 6.5% |
| < 0.6 | 495 | 3.8% |
| < 0.5 | 348 | 2.7% |
| < 0.4 | 273 | 2.1% |
| **< 0.3** | **237** | **1.8%** |

12,681 / 12,911 个 Test 结构有 foldseek hit。

#### `OOD_SuperfamilyHoldout`（post-split 超家族未见标注）

**定义**：正常 7:1:2 划分**之后**，Test 中 superfamily（C.A.T.H，取自
`output/cath_s95_summary.csv`）未在 Train/Val 中出现的样本。不预留任何样本，
划分分布完全不受影响；这些样本是普通划分样本，完整参与所有依赖型 OOD 轴。

**规模**：**761** domains（5.9%），涉及 **697** 个 superfamily、**223** 个
Topology。

> **语义注意**：自然落入 Test 的"未见超家族"几乎全部是极小超家族——大超家族
> 所有成员恰好全部落入 20% Test 的概率随大小指数衰减。本轴测的是"训练时从未
> 见过的稀有超家族"（超家族级别的长尾/孤儿）。
>
> **fold 任务使用 CATH 原生标签**（域 ID 自带 C.A.T.H），域级分类 OOD 直接由
> CATH 分类号推导，**不适用** `.skills/cath-ted-domain-annotation`（TED 域
> 标注流程）；`OOD_FoldHoldout` 因与"Topology 即分类标签、每个 Topology 必须
> 进 Train"的设计在定义上不兼容，本任务不设置该列。

#### `OOD_ExtremeShort` 与 `Default`

- **定义**：`len(aa_seq) < 60` → `OOD_ExtremeShort=True`，**对全部行如实标注**
  （含 24 条 <60 的 NewEC 先验行），长度属性与 `Default` 语义独立
- **规模**：**767** domains（5.9%）
- **`Default`**：全表 True。短序列在划分前未被单独预留，与所有样本一样走
  7:1:2 分层划分，故 ES 样本同为 Default=True（与 PPI 等任务约定一致）；
  `Default=False` 仅用于"划分前单独预留"的场景，本任务不存在
- ES 行完整参与所有依赖型 OOD 计算，`idr_ratio` 保留；全 Test 无 >1,000
  残基样本（03 建池时已剔除 5 个超长域），无 `OOD_ExtremeLong` 列

#### `OOD_LongTail`（预测标签长尾）

**定义**：**Train 中出现频次 ≤ 10** 的 Topology（`label`）类别，其 Test 样本
标记为长尾 OOD（对全部行计算）。

**规模**：**753** 个长尾 Topology，**1,081** 个 Test 样本（8.4%）。

#### `OOD_IDR`（内在无序区域）

- **工具**：metapredict v3，残基 disorder propensity > 0.5 判无序；预测前过滤
  非 20 种标准氨基酸字符（strip 口径）
- **判定**：`idr_ratio` = 无序残基占比（residue-ratio），> **0.3** →
  `OOD_IDR=True`
- **缓存**：`output/idr_ratio_full.csv`（按序列增量缓存；另从
  `output/splits_ca_backup/cath_test.csv` 继承旧算结果，避免重算）
- **方法学**：与 `.skills/ood-annotation-toolkit` 的 `idr_ood.py
  --mode residue-ratio --threshold 0.3 --clean strip` 口径一致

**规模**：**401** domains（3.1%）；`idr_ratio` 全 12,911 行有值（范围
0.0–1.0，无缺失）。

### 4.3 skill 链收尾（EC 长尾四列与 InD 重算）

> 以下步骤脚本均在 `.skills/` 下，按序执行：

1. **EC 长尾四列**：`.skills/ec-function-ood-annotation/scripts/add_unified_longtail_ec.py`
   （全项目批量驱动，内置 fold 的 EC 源适配，直接写回 `splits/`）。
   口径（all-not-in）：样本有 EC 注释且其**所有** EC（对应层级）在 **Train**
   中频次都 < 5 / < 10（**含 0**）→ True；无 EC 注释 → False。EC 数据源与
   NewEC 先验一致（`data/sifts/sifts_ec_annotations.json`，取完整 EC 编号）。

   | 列 | True 数 | 占 Test |
   |------|--------:|--------:|
   | `OOD_LongTail_EC_L3_le5` | 368 | 2.9% |
   | `OOD_LongTail_EC_L3_le10` | 415 | 3.2% |
   | `OOD_LongTail_EC_L4_le5` | 1,459 | 11.3% |
   | `OOD_LongTail_EC_L4_le10` | 2,199 | 17.0% |

   > 全部 336 行 NewEC 先验必然命中 LongTailEC——EC7 在 Train 频次为 0（含 0
   > 口径），语义嵌套符合预期。

2. **InD 重算**：`.skills/ood-annotation-toolkit/scripts/finalize_ood_columns.py
   --mode ind`。`InD` = 所有 `OOD_*` + `seq_Redundancy_*` + `TM-score_*` 标记
   均为 False（NaN 视为 False）→ True；并集包含 `OOD_LongTail_EC_*` 四列，
   须在所有 OOD 列齐备后统一重算（旧名 `Pure_ID` 原位更名重算；
   `--mode nan` 对本任务无对象——全表 Default=True）。
   **规模**：**1,805** domains（**14.0%**）。

3. **体检**：`.skills/ood-annotation-toolkit/scripts/check_splits.py`
   （全部通过，2026-08-27）。

### 4.4 最终产物

位置：`datasets/fold_classification/splits/`

| 文件 | 行数 | 列 |
|------|------|------|
| `cath_train.csv` | 43,473 | `unique_id, aa_seq, struct_file, label` |
| `cath_val.csv` | 6,132 | 同上 |
| `cath_test.csv` | 12,911 | 核心 4 列 + 27 个标注列（31 列） |

`cath_test.csv` 完整列结构（31 列）：

| 列名 | 类型 | 说明 |
|------|------|------|
| `Default` | bool | 全表 True（本任务无划分前长度预留样本） |
| `InD` | bool | 纯 ID-Test（1,805 个，14.0%） |
| `seq_Redundancy_90` ~ `seq_Redundancy_30` | bool | 序列同源 OOD（7 个阈值） |
| `TM-score_0.9` ~ `TM-score_0.3` | bool | 结构同源 OOD（7 个阈值） |
| `OOD_Orphan` | bool | Train+Val 中无显著同源 hit |
| `OOD_SuperfamilyHoldout` | bool | superfamily 未在 Train/Val 出现（post-split） |
| `OOD_NewEC` | bool | EC7（转位酶）划分前先验预留 |
| `OOD_ExtremeShort` | bool | < 60 残基（含 24 条 NewEC 先验行） |
| `OOD_LongTail` | bool | Train 频次 ≤ 10 的 Topology |
| `idr_ratio` | float | metapredict 无序残基占比 |
| `OOD_IDR` | bool | idr_ratio > 0.3 |
| `OOD_LongTail_EC_L3_le5/le10`、`OOD_LongTail_EC_L4_le5/le10` | bool | Train 频次 < 5 / < 10（含 0）的 EC L3 / L4 |

### 4.5 各 OOD 维度统计汇总（Test = 12,911，全行口径）

| OOD 维度 | 阈值/说明 | Domain 数 | 占 Test |
|----------|----------|----------:|--------:|
| **OOD-Prior（NewEC）** | EC7 holdout（划分前预留） | **336** | **2.6%** |
| seq_Redundancy_90 ~ 30 | max pident < 90%…30% | 10,257 / 9,118 / 8,114 / 7,055 / 5,586 / 3,876 / 2,353 | 79.4%…18.2% |
| TM-score_0.9 ~ 0.3 | max alntmscore < 0.9…0.3 | 4,022 / 1,755 / 839 / 495 / 348 / 273 / 237 | 31.2%…1.8% |
| OOD_Orphan | mmseqs2 无显著 hit | 1,978 | 15.3% |
| OOD_SuperfamilyHoldout | post-split 未见超家族 | 761 | 5.9% |
| OOD_ExtremeShort | < 60 aa（含 24 条先验行） | 767 | 5.9% |
| OOD_LongTail | Train 频次 ≤ 10 的 Topology | 1,081 | 8.4% |
| OOD_LongTail_EC_L3_le5 / le10 | Train 频次 < 5 / < 10（含 0） | 368 / 415 | 2.9% / 3.2% |
| OOD_LongTail_EC_L4_le5 / le10 | 同上（L4 层级） | 1,459 / 2,199 | 11.3% / 17.0% |
| OOD_IDR | idr_ratio > 0.3 | 401 | 3.1% |
| **InD（纯 ID-Test）** | 无任何 OOD 标记 | **1,805** | **14.0%** |

**子集统计（备查）**

**ExtremeShort 子集（767 行）**：seq_Redundancy_90…30 = 534/491/440/397/323/
270/258；TM-score_0.9…0.3 = 410/284/224/189/166/151/140；Orphan 258；
SF-Holdout 151；LongTail 110；IDR 122；LongTailEC L3le5/L3le10/L4le5/L4le10
= 24/25/54/72；NewEC 24（即 24 条 <60 的 EC7 先验行）。

**NewEC 先验子集（336 行）**：seq_Redundancy_90…30 = 334/330/322/313/292/249/
163；TM-score_0.9…0.3 = 220/115/56/33/24/21/21；Orphan 132；SF-Holdout 36；
LongTail 43；IDR 18。

### 4.6 数据质量检查（全部通过）

| 检查项 | 结果 |
|--------|------|
| Train / Val / Test `unique_id` 两两无交集 | ✅ |
| 全部 1,163 个 Topology 同时出现在 Train 和 Test | ✅ |
| `label` 全部为 `C.A.T` 字符串格式 | ✅ |
| Train/Val 中 EC7（NewEC）样本数 = 0（SIFTS 反查） | ✅ |
| `OOD_SuperfamilyHoldout=True` 行的 superfamily 均不在 Train/Val | ✅ |
| `OOD_Orphan ⊆ seq_Redundancy_30` 阶梯单调 | ✅ |
| ExtremeShort 行 `Default=True`、依赖型 OOD 列有真实值、`idr_ratio` 保留 | ✅ |
| 全部行所有 OOD 列无缺失值 | ✅ |
| `InD` 与"无任何 OOD 标记"双向一致 | ✅ |
| 全部 `struct_file` 对应的 PDB 文件存在于 `pdbs/` | ✅ |
| `.skills/ood-annotation-toolkit/scripts/check_splits.py` 体检 | ✅ 全部通过（2026-08-27） |

---

## 诊断与质控工具

- `.skills/ood-annotation-toolkit/scripts/check_splits.py`：splits 完整体检
  （列完整性、嵌套单调性、InD 一致性等），见 4.6。
- `04_build_topology_split.py` / `05_annotate_test_ood.py` 内置断言与汇总打印
  （见 4.1/4.2）。

## 归档清单（`output/scripts_archive/`）

| 脚本 | 归档原因 |
|------|----------|
| `fix_summary.py` | 一次性补丁：修复 02 早期版本 CSV 追加写缺陷、重建完整汇总表；逻辑已融入 `02_extract_cath_s95_domains.py`（提取结束全量重建） |
| `fix_multimodel_pdbs.py` | 一次性修复：4,076 个多模型（NMR）pdb 截断保留 model 1；修复已落数据，根因已修入 02（首个 `ENDMDL` 即停） |
| `extract_cath_domains_test.py` | 开发期验证脚本（50 域抽样比对，mean identity 0.969），边界解析逻辑已确认，不进主流程 |
| `fill_es_ood.py`、`fill_prior_ood.py` | 2026-08-27 补丁期一次性补齐脚本（ES 行/先验行补齐依赖型 OOD）；逻辑已并入 `05_annotate_test_ood.py` 现行口径 |
| `add_longtail_ec_ood.py` | 初版单列 `OOD_LongTail_EC_L4` 补丁；已被 ec-function-ood-annotation skill 的统一四列取代 |

## 注意事项

1. **Missing Residues**：PDB ATOM 记录常缺失 loop 区域，`n_residues` <
   `length_expected` 属正常现象
2. **表达标签**：少数参考序列含 His-tag 等载体序列，PDB 结构中不存在
3. **Chain '0'**：CATH 中 chain='0' 表示 PDB 无 chain field，SIFTS 查询时按
   `['0','A',' ']` 回退匹配
4. **mmCIF 转换**：从 `.cif.gz` 提取的域统一输出为 PDB 格式，原子序号重排
5. **S95 编号非全局唯一**：需完整键 `(C,A,T,H,S35,S60,S95)` 标识聚类
6. **OOD 重叠**：同一 domain 可同时属于多个 OOD 维度，为 OOD 评估标准做法
7. **EC 注释覆盖**：仅部分链有 SIFTS EC 注释；无 EC 注释的样本在
   `OOD_LongTail_EC_*` 上记 False（非 NaN）
8. **fold 的 OOD_NewEC 为 EC7 先验口径**：划分前预留、单布尔列，与
   ec-function-ood-annotation skill 的 NewEC 口径（all-not-in、L3/L4 两列）
   不同，**不适用** skill 重算；LongTailEC 四列则走 skill 统一口径（见 4.3）
9. **不适用 cath-ted skill**：fold 域自带 CATH 原生 C.A.T.H 标签，无需 TED
   补缺；不设置 `OOD_FoldHoldout`（与 Topology 分层设计定义不兼容，见 4.2）
10. **无脚本环节**：>1,000 残基剔除（5 个域）在 `03_build_domain_pool.py` 内
    按 `length_expected` 口径完成；v1（C.A 粒度）→ v2（C.A.T）的池迁移已由
    03 补建并 diff 验证等价，无遗留手工步骤
11. **历史偏差记录**：本次整理对 03 新建池与 v1 备份池做行级 diff，**完全一致**，
    未发现偏差
12. **历史演进（简要）**：v1（2026-06，C.A 粒度，43 类区分度过低）划分备份于
    `output/splits_ca_backup/`；v2（2026-08）标签改为 Topology 分层 7:1:2，
    NewEC 先验从 EC5+EC6 收窄为仅 EC7，SuperfamilyHoldout 改为 post-split
    标注，FoldHoldout 移除；2026-08-27 起 ES 行与先验行完整参与所有依赖型
    OOD，`Pure_ID` 更名 `InD` 并按全 OOD 并集重算。各次写回前的 splits 备份
    见 `output/backup_*/`
