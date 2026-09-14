# 蛋白质配体结合位点预测数据集处理文档

> 记录从 PDB 结构数据出发，构建蛋白质配体结合位点（Ligand Binding Site, LBS）残基级预测任务数据集的完整处理流程，可按本文从原始数据一次性复现。
> 处理位置：`datasets/ligand_binding_site/`（配体结合位点预测任务主目录）
> 数据源：RCSB PDB（`https://www.rcsb.org/`），本地镜像位于 `data/pdb/`（项目根）

## 任务定义与当前规模

本任务定义蛋白质配体结合位点为：**配体重原子 4.5 Å 范围内的蛋白质残基**。配体范围包括小分子辅因子、药物、金属离子复合物等，但排除水分子、核酸、标准氨基酸、常见结晶添加剂以及未解析的 UNK/UNL（过滤规则见阶段二）。对每个蛋白链-配体对，生成与 `aa_seq` 逐残基对齐的二进制标签数组（`1` = 结合残基，`0` = 非结合残基）。

当前 splits 规模（以 `splits/` 实际文件统计为准）：

| 指标 | Train | Val | Test |
|------|------:|------:|------:|
| **样本数** | **71,093** | **7,938** | **17,225** |
| 唯一蛋白链数 | 65,462 | 6,718 | 14,328 |
| 唯一 PDB 数 | 60,136 | 5,174 | 11,379 |
| 唯一配体种类数 | 31,940 | 2,826 | 5,125 |

Test 内部组成（三者互斥且覆盖全部 Test）：

| 子集 | 定义 | 样本数 | 占比 |
|------|------|------:|-----:|
| `Default` | 60 ≤ length ≤ 1000 | **13,723** | 79.7% |
| `OOD_ExtremeShort` | length < 60 | **1,974** | 11.5% |
| `OOD_ExtremeLong` | length > 1000 | **1,528** | 8.9% |

结构文件：`pdbs/` 共 96,979 个 `{unique_id}.pdb`（含 723 个对应已丢弃无 SMILES 样本的冗余文件，见注意事项）。

---

## 阶段一：源数据获取与解析（PDB → chain-ligand 标签）

> 目标：遍历本地 PDB/mmCIF 结构文件，识别蛋白-配体复合物并计算结合残基标签，聚合为统一的 chain-ligand 级数据集。

### 核心数据规模

| 指标 | 数值 |
|------|------|
| **PDB 候选结构总数** | 255,227 |
| 成功处理结构数 | 128,725 |
| 处理超时结构数 | 1,099 |
| 解析失败结构数 | 77 |
| 聚合前总对数（含重复） | 763,605 |
| **按 `unique_id` 去重后 chain-ligand 对数** | **645,774** |
| 唯一配体种类数 | 46,869 |

### Step 1: 批量提取蛋白-配体复合物

**脚本**：`scripts/01_extract_binding_sites.py`（核心解析逻辑在库文件 `scripts/extract_binding_sites_core.py`，被 01 等 import，不编号）

- 输入：`data/pdb/`（项目根）下所有 `.pdb.gz` 和 `.cif.gz`
- 输出：`output/binding_sites_batch_{NNNNN}.csv`
- 关键参数：重原子距离阈值 `4.5 Å`（`scipy.spatial.cKDTree`）；配体最小重原子数 `5`；最小结合残基数 `> 3`（代码实现为至少 3 个）
- 结构解析使用 BioPython `PDBParser` / `MMCIFParser`

### Step 2: 聚合批次文件

**脚本**：`scripts/02_aggregate_batches.py`

- 输入：`output/binding_sites_batch_*.csv`（共 221 个 batch）
- 按 `unique_id` 去重（保留首次出现记录）：763,605 → **645,774** 对
- 输出：`output/ligand_binding_sites_raw.csv`、`output/ligand_binding_sites_std.csv`

> 归档注记：主批量进程曾因进程池崩溃中断，一次性诊断/恢复脚本 `03_diagnose_remaining_files.py`（单文件诊断分类 OK/Timeout/Error/Crash，超时 10 s）与 `05_reprocess_timeouts.py`（60 s 超时重试，前 40 个最小文件仍全部超时，零产出后终止）已完成使命，均归档至 `output/scripts_archive/`；1,099 个超时大文件（约占全部 PDB 的 0.4%）接受跳过。断点续跑脚本 `02_resume_extract_binding_sites.py` 同属一次性归档。

### 关键分布统计

成功处理的 128,725 个结构中，Top 10 配体（按 chain-ligand 对数统计）：

| 配体 | 对数 |
|------|------|
| CLA | 61,981 |
| HEM | 16,022 |
| ADP | 15,397 |
| ATP | 12,849 |
| BCR | 12,013 |
| CDL | 10,082 |
| SF4 | 9,406 |
| BCL | 9,100 |
| FAD | 8,118 |
| NAD | 7,482 |

### 异常值/过滤样本

| 类型 | 数量 | 原因 |
|------|------|------|
| 处理超时 | 1,099 | 单文件在 10 s（诊断）/ 60 s（重试）内无法完成解析 |
| 解析错误 | 77 | mmCIF 缺少必要列（`_atom_site.id`、`_atom_site.pdbx_PDB_ins_code`）或 B factor 异常 |

---

## 阶段二：标签聚合与计算（标签语义与配体过滤规则）

> 标签在阶段一提取时即按以下口径计算（实现于 `scripts/extract_binding_sites_core.py`），本阶段无额外脚本，仅固化语义。

### 标签定义

- 结合残基：任一重原子与配体任一重原子距离 ≤ 4.5 Å 的蛋白残基
- `label`：与 `aa_seq` 逐残基对齐的二进制数组，序列化为方括号包裹、单空格分隔的整数字符串（如 `[0 0 1 0]`），**强制保存为字符串**且 `len(label) == len(aa_seq)`

### 配体过滤规则

| 过滤项 | 规则 |
|--------|------|
| 水分子 | 排除 `HOH`, `DOD`, `WAT` 等 |
| 核酸 | 排除 `A/T/G/C/U/DA/DT/...` |
| 标准氨基酸 | 排除 20 种标准 AA |
| 修饰氨基酸 | 34 种映射回标准 AA（如 `MSE→M`），并视为蛋白链一部分 |
| 结晶添加剂 | 排除 127 种常见添加剂（清单见 `extract_binding_sites_core.py`） |
| 最小配体大小 | 保留重原子数 ≥ 5 的配体 |
| 最小结合位点 | 保留结合残基数 > 3 的 chain-ligand 对 |

### 验证与质控

- 所有 `label` 数组长度与 `aa_seq` 长度一致
- 对每个 chain-ligand 对记录 `num_binding_residues` 与 `num_ligand_heavy_atoms`
- 按 `unique_id`（`{pdb_id}_{chain_id}_{ligand_name}_{ligand_chain}_{ligand_resnum}`）去重，避免同一蛋白链-配体对重复

---

## 阶段三：去冗余与数据精化

> 目标：配体-wise 序列去冗余、提取独立蛋白链结构文件并校验序列-结构一致性、注释配体 SMILES/ECFP4 并丢弃无法解析的样本。

### 3.1 配体-wise 95% 序列去冗余

为避免模型在同质化样本上过拟合并保证测试集评估公平性，对 `ligand_binding_sites_std.csv` 做序列层面去冗余。去冗余以**配体名称**为分组单位：相同配体下的蛋白链按 95% 序列相似度聚类，每簇仅保留一条代表序列；不同配体之间的序列即使高度相似也不互相删除，以保留跨配体的结合多样性。

| 指标 | 数值 |
|------|------|
| 去冗余前 chain-ligand 对数 | 645,774 |
| **去冗余后 chain-ligand 对数** | **109,954** |
| 去冗余前 PDB 数 | 128,725 |
| 去冗余后 PDB 数 | 88,208 |
| 唯一配体种类数 | 46,869（保持不变） |

**Step 1: mmseqs2 迭代聚类去冗余** — `scripts/03_dedup_by_ligand_95.py`

- 输入：`output/ligand_binding_sites_std.csv`；输出：`output/ligand_binding_sites_std_dedup95.csv`
- 为每个 `ligand_name` 生成 150 aa 的配体特异性前缀标签，确保 mmseqs2 不会把不同配体的序列聚到同一簇
- `mmseqs easy-cluster --min-seq-id 0.95 -c 0.95 --cov-mode 1 --cluster-mode 1`，迭代执行直至本轮不再移除样本
- 每簇代表按质量优先规则保留：最多结合残基 → 最高分辨率 → 最长蛋白链 → 最多配体重原子

**Step 2: 精确后处理去冗余** — `scripts/04_exact_cleanup_dedup95.py`

- 对样本数 ≤ 50 的配体，用 **edlib** 精确全局比对（identity ≥ 0.95）构建相似图，每连通分量保留一条序列
- 对样本数 > 50 的配体，单独运行 `mmseqs easy-cluster` 精确聚类
- 末尾经 subprocess 调用 `03_dedup_by_ligand_95.py` 做最终收敛检查，确保无残留 ≥95% 相似对

**验证结果**：数量最多的前 5 种配体分别重跑 `mmseqs easy-cluster`（相同参数），所有序列各自成簇（`max_cluster_size = 1`，`redundant = 0`）；再次迭代去冗余立即收敛（第二轮移除数为 0）。

Top 10 配体去冗余前后对比：

| 配体 | 去冗余前 | 去冗余后 |
|------|---------|---------|
| CLA | 61,981 | 795 |
| HEM | 16,022 | 981 |
| ADP | 15,397 | 1,727 |
| ATP | 12,849 | 1,363 |
| BCR | 12,013 | 433 |
| CDL | 10,082 | 589 |
| SF4 | 9,406 | 652 |
| BCL | 9,100 | 103 |
| FAD | 8,118 | 948 |
| NAD | 7,482 | 927 |

### 3.2 蛋白质链结构提取

> 注：结构提取针对划分后的样本清单执行（即阶段四时间划分 + 配体过滤产出的 splits），此处按主题归入数据精化；线性复现时它在阶段四 Step 2 之后运行。

**脚本**：`scripts/07_extract_protein_chains.py`

- 按 PDB ID 分组，每个源结构文件只解析一次；对 `.pdb.gz` 使用手动快速 ATOM 解析，对 `.cif.gz` 解析 `_atom_site` 循环
- 仅保留**蛋白质残基**（标准氨基酸 + 常见修饰氨基酸），跳过核酸、配体、结晶添加剂；修饰氨基酸（MSE、SEP、TPO 等）即使以 `HETATM` 记录出现也会保留
- 输出链 ID 截断为单字符（符合标准 PDB 格式），原子序号按链重新编号，避免大结构原子序号超 5 列导致坐标错位
- 多模型（NMR）结构只保留首个模型：`parse_pdb_file` 首个 `ENDMDL` 即停、`parse_cif_file` 按 `pdbx_PDB_model_num` 仅保留 model 1；mmCIF 带引号原子名（如 `"C2'"`）已做去引号处理（根因修复，历史上的 549 个多模型文件已由一次性脚本重提取，见归档清单）
- 源文件名同时支持大写/小写两种形式；多进程并发（16 进程），约 6 分钟完成
- 输出：`pdbs/{unique_id}.pdb`（当时 96,979 个样本全部提取成功，缺失源结构 0 个）

### 3.3 序列-结构一致性校验（QC）

**脚本**：`scripts/13_verify_seq_structure_consistency.py`（诊断/质控工具，作为 QC 收尾步骤保留）

- 批量提取后约 0.80% 样本的 `aa_seq` 与 `pdbs/` 中提取的序列不一致（修饰残基命名、重复残基编号、非蛋白原子等原因）
- 历史差异详情见 `output/seq_structure_mismatches.csv`；当前 splits 全部样本 `aa_seq` 与结构文件序列一致率 **100.00%**

> 归档注记：不一致样本的一次性修复脚本 `15_patch_inconsistent_from_structure.py`（用 core 的蛋白识别逻辑重算 `aa_seq`/`label`）与 `16_reextract_inconsistent_biopython.py`（Biopython 重提取并对重复残基编号重编号）已将 96,201/96,979（99.20%）修复至 100%，修复已落数据、根因在 core/07，均归档。

### 3.4 配体 SMILES / ECFP4 注释与缺失样本丢弃

**脚本**：`scripts/08_add_ligand_smiles_ecfp4.py`、`scripts/09_drop_missing_smiles_and_rebuild.py`

- 从 `unique_id` 提取 `ligand_name`，从 PDB CCD（`data/pdb/components.cif.gz`）的 `_pdbx_chem_comp_descriptor` 段读取 SMILES / InChI
- 用 RDKit 计算 canonical SMILES 与 ECFP4（Morgan, radius=2, nBits=2048）；对无法直接 kekulize 的芳香/带电 SMILES 使用 `sanitize=False` + 跳过 kekulization 的 fallback
- ECFP4 以稀疏位点列表存储，如 `[12,45,200]`
- 09 丢弃 `ligand_smiles` 为空的样本：Train 548 条、Test 20 条（Val 0 条），并就地刷新 NewFunction/LongTail 增量列（默认快速模式；`--recompute_full_ood` 会调用 11 全量重算）

| 指标 | Train | Val | Test |
|------|------:|------:|------:|
| 写入 SMILES / ECFP4 前 | 71,641 | 7,938 | 17,400 |
| 丢弃无 SMILES 后 | **71,093** | **7,938** | **17,380** |

验证：`ligand_smiles` 缺失（Train / Val / Test）= 0 / 0 / 0；`ligand_ecfp4` 空列表 = 0；标签长度与序列长度不一致 = 0。

> 注：09 号脚本写出的 `has_ec_annotation` / `OOD_NewFunction` 为统一前的历史列名，现已被 `OOD_NewEC_L4` 等取代（见阶段四）；脚本内 2026-08 已适配统一列名（`label`/`struct_file`/`OOD_Extreme*`）。

---

## 阶段四：数据划分与 OOD 标注

> 目标：按时间顺序划分 Train/Val/Test，配体相似度过滤与极端长度重组，并在测试集上构建完整的 OOD 标记列。
> 各维度 OOD 的统一口径由对应 skill 固化（`.skills/` 下），重算时必须调用 skill 脚本，不得另起口径。

### Step 1: 时间顺序划分（Train / Val / Test）

**脚本**：`scripts/05_create_time_based_splits.py`

- 输入：`output/ligand_binding_sites_std_dedup95.csv`
- 按 PDB 结构沉积日期（deposition date）对所有唯一 `pdb_id` 排序，按 70% / 10% / 20% 时间切分
- 同一 PDB 的所有 chain-ligand 对只出现在一个集合中，避免数据泄露
- 无法解析沉积日期的 37 个 PDB（共 49 对）被分配到 **Test** 集

| 集合 | PDB 数 | chain-ligand 对数 | 配体种类数 | 时间范围 |
|------|--------|-------------------|-----------|---------|
| **Train** | 61,731 | 73,608 | 32,630 | 1974-06-06 ~ 2021-04-28 |
| **Val** | 8,818 | 12,487 | 6,480 | 2021-04-28 ~ 2022-12-22 |
| **Test** | 17,659 | 23,859 | 12,316 | 2022-12-23 ~ 2026-05-14（含 20 个无日期 PDB） |

泄露检查（划分时点）：Train ∩ Val、Train ∩ Test、Val ∩ Test（按 `pdb_id`）均为 0。

> 历史注记：早期按 `pdb_id` 随机打乱（种子 42）的 80/10/10 方案（`06_create_splits.py`）已被时间划分完全取代，脚本已归档。

### Step 2: 验证集/测试集配体相似度过滤

**脚本**：`scripts/06_filter_test_by_ligand_similarity.py`

- 从 CCD 解析配体 SMILES，RDKit 生成 ECFP4 指纹（Morgan, radius=2, nBits=2048）
- 对 Val/Test 中每个配体，计算其与训练集所有配体的最大 Tanimoto 相似度，仅保留最大相似度 **> 0.6** 的样本（训练集不变化）

| 划分 | 过滤前 | 过滤后 | 缺失 SMILES | 指纹失败 | 被阈值筛掉 |
|------|--------|--------|-------------|----------|------------|
| Val  | 12,487 | **8,555** | 31 | 234 | 3,667 |
| Test | 23,859 | **14,816** | 35 | 475 | 8,533 |

### Step 3: 极端长度子集重组与超长剔除

**脚本**：`scripts/10_rebuild_splits_with_extreme_subset.py`

- 所有 `protein_length < 60` → `OOD_ExtremeShort`；`protein_length > 1000` → `OOD_ExtremeLong`
- train/val/test 中的极端长度样本统一移至 **test**；长度正常的 test 样本标记 `Default = True`
- Default 子集的 OOD 列经 subprocess 调用 `scripts/11_build_ood_splits.py` 重算
- 重组后规模（结合阶段三 3.4 的 SMILES 丢弃）：Train 71,093 / Val 7,938 / Test 17,380（Default 13,723 + Short 1,974 + Long 1,683）
- **超长剔除（无独立脚本）**：按全项目统一规则，单链长度 > 2000 的样本一律从数据集剔除。本任务受影响样本全部在 Test 的 `OOD_ExtremeLong` 中，共 **155** 条（长度 2,003–4,584），Test 17,380 → **17,225**；剔除前备份见 `output/len2000_removal_backup/`

> 注：10 号脚本输出的 `ExtremeShort`/`ExtremeLong`/`is_time_cutoff` 等为历史列名，后续统一为 `OOD_ExtremeShort`/`OOD_ExtremeLong` 等；脚本内 2026-08 已适配统一核心列名（`label`/`struct_file`），并修复了原脚本未定义 `ROOT` 的路径缺陷。

### Step 4: 序列同源性 OOD（seq_Redundancy_* / OOD_Orphan，仅 Default 子集）

**正式口径**：`.skills/homology-ood-annotation/scripts/seq_homology_ood.py`——query = Default test，target = train+val；**mmseqs2 `easy-search -s 7`**（easy-search 默认 e-value ≤ 1e-3 显著性门槛），`--format-output query,target,pident,evalue`；取 max pident 按 90/80/…/30 阶梯判定 `seq_Redundancy_*`；`OOD_Orphan` = 无显著 hit（实现为 best e-value > 1e-3，与"无 hit"等价）。query 按 chain_key（`{pdb}_{chain}`）去重。

- 本地脚本 `scripts/11_build_ood_splits.py` 内嵌 mmseqs 段已标注 **[superseded]**（重算以 skill 为准）；历史旧式 `mmseqs search -e 100 -c 0.0` 口径曾导致 Orphan 与 seq_Redundancy_30 计数倒挂，由一次性修复脚本（`24_fix_seq_homology_easysearch.py`，已归档）按 easy-search -s 7 重算修复
- 缓存：`output/ood/mmseqs_max_identity_easysearch.csv`（chain_key 键）；注意重放前需确认键空间为 chain_key 而非 unique_id
- Default=False（极端长度）行的 `seq_Redundancy_*` / `OOD_Orphan` 列一律为 NaN（不参与计算）
- 单调性自检：`OOD_Orphan` ⊆ `seq_Redundancy_30` ⊆ … ⊆ `seq_Redundancy_90` 成立

### Step 5: 结构相似性 OOD（TM-score_*，仅 Default 子集）

**脚本**：`scripts/11_build_ood_splits.py`（内嵌 foldseek 调用）；统一口径见 `.skills/homology-ood-annotation/scripts/struct_homology_ood.py`

- 当前缓存为 **easy-search 统一口径**（2026-08-28 重算）：`foldseek easy-search query target out tmp --format-output query,target,alntmscore`（默认 e-value ≤ 10、sensitivity 9.5；query = test 全部唯一 chain_key 14,328 条，target = train+val 唯一 chain_key 72,180 条，经符号链接组目录避免拷贝），取 max alntmscore 按 0.9/0.8/0.7/0.6/0.5/0.4/0.3 阶梯判定 `TM-score_*`
- 旧口径为 `foldseek search -a -e inf --max-seqs 1000` + convertalis：对比实验发现其缓存中 2,433/11,406 条（21%）max_tmscore > 1（非法值，疑似当年 format-output 取错列），且合法值整体系统性偏高（937 条中位偏高 0.074），导致各档 OOD 计数被严重低估（如 0.5 档 4 → 241）。对比脚本与报告存档于 `output/foldseek_easysearch_compare/`（run_compare.py / compare_report.txt / easysearch_max_tmscore.csv），旧缓存与旧 test 备份于 `output/backup_before_foldseek_easysearch_unify/`
- 缓存：`output/ood/foldseek_max_tmscore.csv`（chain_key 键，14,328 条 = test 全部唯一链）；Default=False 行各列置 NaN

### Step 6: IDR 无序 OOD（idr_ratio / OOD_IDR，仅 Default 子集）

**统一工具**：`.skills/ood-annotation-toolkit/scripts/idr_ood.py`（metapredict v3，实测 3.0.2）；本地实现于 `scripts/11_build_ood_splits.py`

- 序列清洗：**map 模式**（B→N, Z→Q, J→L, U→C, O→K，其余非标准字符→A）
- `idr_ratio` = disorder score > 0.5 的**无序残基占序列全长比例**（residue-ratio）
- 判定：**`OOD_IDR` = `idr_ratio` > 0.3**
- 重算命令对应 `idr_ood.py --mode residue-ratio --threshold 0.3 --clean map`（与 ss/func/fold 一致）
- 缓存：脚本默认缓存路径 `output/ood/idr_predictions.csv`（按清洗后序列缓存；当前未保留该文件，仅有中间产物 `output/ood/idr_fastas/`）；Default=False 行 `idr_ratio` 留空、`OOD_IDR` 置 NaN

### Step 7: EC 功能新颖性 OOD（OOD_NewEC_L4 / OOD_NewEC_L3，仅 Default 子集）

**统一口径与脚本**：`.skills/ec-function-ood-annotation/scripts/`（`ec_annotate.py` + `ec_ood.py`，批量驱动 `add_unified_newec.py`）

- EC 数据源：SIFTS 链级 `data/sifts/sifts_chain_ec.tsv.gz`（`unique_id` 前两段 → `PDB_CHAIN`）；train+val 唯一 EC 编号 3,380 个
- 参考集 = train+val 出现的 EC 集合；**all-not-in** 语义：样本有 EC 注释且其所有 EC（对应层级）都不在 train+val 中 → True；无 EC 注释 → False；层级解析允许部分 `-`（如 `3.4.22.-` → L3=`3.4.22`）
- `OOD_NewEC_L4` 由 ec skill 统一驱动脚本就地重算；`OOD_NewEC_L3` 为统一时新增列
- Default 子集有 EC 注释的样本：**6,307 / 13,723（46.0%）**；Default=False 行两列置 NaN
- 本地增量脚本 `22_add_newfunction_longtail_ood.py` 的 EC 部分已被 skill 取代（脚本已归档）

### Step 8: EC 长尾功能 OOD（OOD_LongTail_EC_*，仅 Default 子集）

**统一口径与脚本**：`.skills/ec-function-ood-annotation/scripts/add_unified_longtail_ec.py`（批量驱动；初版单列脚本 `26_add_longtail_ec_ood.py` 已归档）

- 四列：`OOD_LongTail_EC_L3_le5` / `OOD_LongTail_EC_L3_le10` / `OOD_LongTail_EC_L4_le5` / `OOD_LongTail_EC_L4_le10`
- 判定：样本有 EC 注释且其**所有** EC（对应层级）在 **train** 中频次都 < 5 / < 10（**含 0**）→ True（all-not-in）；无 EC 注释 → False
- EC 数据源与键推导同 Step 7；Default=False 行四列置 NaN

### Step 9: 结合位点占比长尾 OOD（OOD_LongTail，仅 Default 子集）

**脚本**：`scripts/11_build_ood_splits.py`

- 解析 `label`（与 `aa_seq` 逐残基对齐的二进制数组），计算 `binding_site_ratio = 结合残基数 / 蛋白链长度`
- `binding_site_ratio < 0.01` → `OOD_LongTail = True`；Default=False 行 `binding_site_ratio` 留空、`OOD_LongTail` 置 NaN

### Step 10: CATH+TED 域标注 OOD（OOD_FoldHoldout / OOD_SuperfamilyHoldout，仅 Default 子集）

**脚本**：`scripts/12_ted_augment_cath_ood.py`（wrapper，标注与 holdout 判定分别调用 `.skills/cath-ted-domain-annotation/scripts/cath_ted_annotate.py` 与 `.skills/cath-ted-domain-annotation/scripts/cath_ted_holdout.py`）

- 样本粒度：链级（`unique_id` 前两段为 PDB/链，如 `10BL_A_ADP_A_503` → `10bl_A`；同一链的不同配体样本共享同一标注）
- 标签 = CATH v4.4 实验标注（名称匹配 `cath-domain-list.txt`）∪ TED putative 标注（REST API，缓存于 `data/ted/ted_api_cache.json`；仅保留与链 UniProt 区间重叠 ≥50% 域长的域；H 级提供 C.A.T.H，T 级仅 C.A.T）
- 参考集 = **Train** 的 CATH+TED 并集（train/val/test 对称增强）；any 语义：链任一 Topology/Superfamily 未在参考集出现 → True；无注释 → False；Default=False 行置 NaN

数据规模与覆盖率：

- 链清单：86,508 条链（`output/lbs_ted_aug_chain_list.csv`，历史一次性生成），25,336 个唯一 UniProt accession；TED 缓存全部命中（离线运行）；3,902 条链无 SIFTS UniProt 映射
- 标注来源分布：both 41,953 / ted-only 29,238 / cath-only 4,674 / none 10,643
- Topology 覆盖率（全部链）：仅 CATH 53.9% → CATH+TED **87.7%**
- 输出：`output/lbs_chain_cath_ted_labels.csv`（链级标签）、`output/lbs_chain_domains.csv`（域级明细，228,136 行）、`output/lbs_chain_holdout_flags.csv`（holdout 判定结果，重跑 12 号 wrapper 时生成）

> 历史注记：CATH 旧标注路线（`20_annotate_cath_domains.py`、`21_annotate_cath_domains_mmseqs2.py`）已被 CATH+TED 并集取代，脚本已归档。

### Step 11: struct_label 结构残基级标注

**脚本**：`.skills/struct-label-alignment/add_struct_label.py`（先 `--dry-run` 再写入）

- train/val/test 三个文件均含 `struct_label` 列（末列）：对 `struct_file` 中的**每个结构残基**（ATOM 记录按文件顺序、以 `(chain, resseq, icode)` 首次出现枚举，含无 CA 原子的残基）给出对应的结合位点标注
- 映射方式：pdb 残基序列与 `aa_seq` 全局比对，结构残基继承比对位置的 `label`；无法映射回序列的结构残基记 **`-1`**（ignore index，训练/评估时跳过）
- 验证：其它列逐单元格零变化；`len(struct_label) == pdb 残基数` 全量断言通过；pdb 残基数 == `len(aa_seq)` 的行 `struct_label` 与 `label` 完全一致；含 `-1` 的行数 / `-1` 残基总数：train 0 / 0，val 0 / 0，test 0 / 0
- 备份：`output/backup_before_struct_label/`

### Step 12: NaN 校验与 InD 重算

**工具**：`.skills/ood-annotation-toolkit/scripts/finalize_ood_columns.py --mode both`

- **NaN 惯例**：Default=False（极端长度）行除 `Default` / `OOD_ExtremeShort` / `OOD_ExtremeLong` 外，所有 OOD 维度列（`seq_Redundancy_*`、`TM-score_*`、`OOD_Orphan`、`OOD_IDR`、`OOD_NewEC_*`、`OOD_LongTail*`、`OOD_FoldHoldout`、`OOD_SuperfamilyHoldout`）一律为 NaN（空单元格，表"未计算"）；`idr_ratio` / `binding_site_ratio` 同样仅 Default 行有值。2026-08-28 起该 NaN 由各生成脚本（`10`/`11`/`12` 及 skill 脚本）在流水线位置直接写入，`finalize_ood_columns.py --mode nan` 仅作兜底校验（当前数据上为 no-op）
- **InD**：所有 `OOD_*` + `seq_Redundancy_*` + `TM-score_*` 均为 False（NaN 视为 False）时 `InD = True`，由 `finalize_ood_columns.py --mode ind` 统一重算
- 体检：`.skills/ood-annotation-toolkit/scripts/check_splits.py`

### OOD 统计（以当前 splits 实测为准；占比基于 Default n=13,723）

| OOD 列 | 定义 | True 数 | 占比 |
|:-------|:-----|--------:|-----:|
| `seq_Redundancy_90` | max pident < 90% | 7,908 | 57.6% |
| `seq_Redundancy_80` | max pident < 80% | 7,115 | 51.8% |
| `seq_Redundancy_70` | max pident < 70% | 6,418 | 46.8% |
| `seq_Redundancy_60` | max pident < 60% | 5,547 | 40.4% |
| `seq_Redundancy_50` | max pident < 50% | 4,562 | 33.2% |
| `seq_Redundancy_40` | max pident < 40% | 3,191 | 23.3% |
| `seq_Redundancy_30` | max pident < 30% | 1,859 | 13.5% |
| `TM-score_0.9` | max TM-score < 0.9 | 3,017 | 22.0% |
| `TM-score_0.8` | max TM-score < 0.8 | 1,310 | 9.5% |
| `TM-score_0.7` | max TM-score < 0.7 | 691 | 5.0% |
| `TM-score_0.6` | max TM-score < 0.6 | 376 | 2.7% |
| `TM-score_0.5` | max TM-score < 0.5 | 241 | 1.8% |
| `TM-score_0.4` | max TM-score < 0.4 | 189 | 1.4% |
| `TM-score_0.3` | max TM-score < 0.3 | 136 | 1.0% |
| `OOD_Orphan` | 无显著同源 hit（e-value > 1e-3） | 1,310 | 9.5% |
| `OOD_IDR` | 无序残基占比 > 30% | 140 | 1.0% |
| `OOD_NewEC_L4` | 全部 L4 EC 未在 train+val 出现 | 271 | 2.0% |
| `OOD_NewEC_L3` | 全部 L3 EC 未在 train+val 出现 | 19 | 0.1% |
| `OOD_LongTail` | 结合位点占比 < 1% | 330 | 2.4% |
| `OOD_FoldHoldout` | 任一 Topology 未在 train 出现 | 15 | 0.1% |
| `OOD_SuperfamilyHoldout` | 任一 Superfamily 未在 train 出现 | 135 | 1.0% |
| `OOD_LongTail_EC_L3_le5` | 全部 L3 EC train 频次 < 5 | 39 | 0.3% |
| `OOD_LongTail_EC_L3_le10` | 全部 L3 EC train 频次 < 10 | 68 | 0.5% |
| `OOD_LongTail_EC_L4_le5` | 全部 L4 EC train 频次 < 5 | 803 | 5.9% |
| `OOD_LongTail_EC_L4_le10` | 全部 L4 EC train 频次 < 10 | 1,142 | 8.3% |

- `InD`（全部 OOD 标记均为 False）：**5,209 / 17,225（占 Test 30.2%）**
- CATH+TED 无注释样本（两列判 False）：4,234 条 Default 样本
- 除 `Default` / `OOD_ExtremeShort` / `OOD_ExtremeLong` 外，所有 OOD 列仅在 Default 子集计算，极端长度行（3,502 行）统一为 NaN

### 验证与质控

| 检查项 | 结果 |
|--------|------|
| Train / Val / Test `unique_id` 两两交集 | 0 |
| Default 子集与 Train / Val 的 `pdb_id` 交集 | 0 |
| Test 极端长度行与 Train 共享 `pdb_id`（见注意事项 5） | 516 行 / 284 个 PDB |
| Test 极端长度行与 Val 共享 `pdb_id` | 284 行 / 127 个 PDB |
| `label` 长度与 `aa_seq` 长度不一致 | 0 |
| `Default` + `OOD_ExtremeShort` + `OOD_ExtremeLong` == Test 总数 | 17,225 == 17,225 |
| 极端长度行其他 OOD 列全为 NaN | 是（3,502 行） |
| `OOD_Orphan` ⊆ `seq_Redundancy_30` ⊆ … ⊆ `seq_Redundancy_90` | 是 |
| `TM-score_0.3` ⊆ … ⊆ `TM-score_0.9`（Default 子集） | 是 |
| 缺失 `ligand_smiles`（Train / Val / Test） | 0 / 0 / 0 |
| splits 样本缺失 `pdbs/` 结构文件 | 0 |
| `aa_seq` 与结构文件序列一致率 | 100.00% |

---

## 输出文件汇总

### 目录结构

```
datasets/ligand_binding_site/
├── DATA_PROCESS.md              # 本文档
├── data/                        # 任务专属数据（当前为空；原始结构/SIFTS 在项目根 data/）
├── logs/                        # 运行日志
├── pdbs/                        # 提取的蛋白质链 PDB 结构文件（96,979 个）
├── splits/                      # 最终划分文件（仅 train/val/test 三个 CSV）
│   ├── ligand_binding_site_train.csv
│   ├── ligand_binding_site_val.csv
│   └── ligand_binding_site_test.csv
├── scripts/                     # 处理脚本（按流程编号）
│   ├── 01_extract_binding_sites.py            # 阶段一：批量提取
│   ├── extract_binding_sites_core.py          # 库文件：解析与结合位点检测逻辑（被 01 等 import）
│   ├── 02_aggregate_batches.py                # 阶段一：批次聚合
│   ├── 03_dedup_by_ligand_95.py               # 阶段三：配体-wise 95% 去冗余
│   ├── 04_exact_cleanup_dedup95.py            # 阶段三：精确后处理
│   ├── 05_create_time_based_splits.py         # 阶段四：时间划分
│   ├── 06_filter_test_by_ligand_similarity.py # 阶段四：配体相似度过滤
│   ├── 07_extract_protein_chains.py           # 阶段三：链结构提取（快速版）
│   ├── 08_add_ligand_smiles_ecfp4.py          # 阶段三：SMILES/ECFP4 注释
│   ├── 09_drop_missing_smiles_and_rebuild.py  # 阶段三：缺失丢弃与重建
│   ├── 10_rebuild_splits_with_extreme_subset.py # 阶段四：极端长度重组
│   ├── 11_build_ood_splits.py                 # 阶段四：OOD 主体（foldseek/IDR/LongTail；mmseqs 段 superseded）
│   ├── 12_ted_augment_cath_ood.py             # 阶段四：CATH+TED 域标注 wrapper
│   └── 13_verify_seq_structure_consistency.py # 诊断/QC：序列-结构一致性校验
├── output/                      # 中间产物、缓存、备份
│   ├── ligand_binding_sites_raw.csv / _std.csv / _std_dedup95.csv
│   ├── binding_sites_batch_*.csv              # 阶段一批次输出
│   ├── ood/                                   # mmseqs/foldseek/IDR 缓存
│   ├── intermediate_splits/                   # 历史中间划分文件
│   ├── lbs_chain_cath_ted_labels.csv / lbs_chain_domains.csv / lbs_chain_holdout_flags.csv
│   ├── backup_*/                              # 各次 splits 写回前备份
│   ├── len2000_removal_backup/                # >2000 剔除前备份
│   └── scripts_archive/                       # 归档的一次性/被取代脚本
```

### 关键输出文件

**最终划分文件**（位置：`datasets/ligand_binding_site/splits/`）

| 文件 | 行数 | 格式 | 说明 |
|------|------|------|------|
| `ligand_binding_site_train.csv` | 71,093 | CSV | 最终训练集（无 OOD 列，7 列） |
| `ligand_binding_site_val.csv` | 7,938 | CSV | 最终验证集（无 OOD 列，7 列） |
| `ligand_binding_site_test.csv` | 17,225 | CSV | 最终测试集（6 个核心列 + 29 个 OOD/标注列 + `struct_label`，共 36 列） |

**结构文件**（位置：`datasets/ligand_binding_site/pdbs/`）

- 共 96,979 个 `{unique_id}.pdb` 文件，每个文件仅包含对应蛋白链的 ATOM 记录
- 当前 splits 全部 96,256 个样本均有对应结构文件，且 `aa_seq` 与结构文件序列完全一致
- 其中 723 个 PDB 文件对应已被丢弃的无 SMILES 样本（不在 splits 中），保留备用

**中间产物**

- `output/ligand_binding_sites_raw.csv` / `ligand_binding_sites_std.csv` / `ligand_binding_sites_std_dedup95.csv`：阶段一/三主数据
- `output/intermediate_splits/`：历史中间划分（`*_dedup95*`、`*_ligsim60*`、时间划分映射表等）
- `output/ood/`：OOD 计算缓存（`mmseqs_max_identity_easysearch.csv`、`foldseek_max_tmscore.csv`、`idr_predictions.csv`，query 键均为 chain_key）

#### 字段说明（train / val，7 列）

| 字段 | 类型 | 说明 |
|------|------|------|
| `unique_id` | string | 样本唯一标识，格式 `{pdb_id}_{chain_id}_{ligand_name}_{ligand_chain}_{ligand_resnum}` |
| `aa_seq` | string | 蛋白链氨基酸序列 |
| `struct_file` | string | 结构文件名（仅文件名），即 `{unique_id}.pdb`，位于 `pdbs/` |
| `ligand_smiles` | string | 配体 canonical SMILES（无 SMILES 样本已丢弃） |
| `ligand_ecfp4` | string | 配体 ECFP4（Morgan, radius=2, nBits=2048）稀疏位点列表，如 `[12,45,200]` |
| `label` | string | **强制字符串**，与 `aa_seq` 逐残基对齐的结合位点二进制数组，方括号包裹、单空格分隔（如 `[0 0 1 0]`） |
| `struct_label` | string | 按 **pdb 结构残基顺序** 对齐的标注（见阶段四 Step 11）；无法映射回序列的结构残基记 `-1` |

#### 字段说明（test 追加列，按当前 `splits/ligand_binding_site_test.csv` 表头）

| 字段 | 类型 | 说明 |
|------|------|------|
| `Default` | bool | 长度正常（60 ≤ length ≤ 1000）的 Test 样本 |
| `OOD_ExtremeShort` / `OOD_ExtremeLong` | bool | length < 60 / > 1000（与 Default 互斥） |
| `seq_Redundancy_90` ~ `seq_Redundancy_30` | bool / NaN | max pident 低于对应阈值（仅 Default 计算，极端行 NaN） |
| `TM-score_0.9` ~ `TM-score_0.3` | bool / NaN | max TM-score 低于对应阈值（仅 Default 计算，极端行 NaN） |
| `OOD_Orphan` | bool / NaN | 在 Train+Val 中无显著 mmseqs2 hit |
| `idr_ratio` | float / 空 | 无序残基占序列全长比例（仅 Default 有值） |
| `OOD_IDR` | bool / NaN | `idr_ratio` > 0.3 |
| `OOD_NewEC_L4` / `OOD_NewEC_L3` | bool / NaN | 全部 EC（L4/L3 层级）未在 train+val 出现（all-not-in） |
| `binding_site_ratio` | float / 空 | 结合残基数 / 蛋白链长度（仅 Default 有值） |
| `OOD_LongTail` | bool / NaN | `binding_site_ratio` < 0.01 |
| `OOD_FoldHoldout` / `OOD_SuperfamilyHoldout` | bool / NaN | 任一 CATH+TED Topology / Superfamily 未在 train 出现 |
| `OOD_LongTail_EC_L3_le5` / `_L3_le10` / `_L4_le5` / `_L4_le10` | bool / NaN | 全部 EC（对应层级）train 频次 < 5 / < 10（含 0） |
| `InD` | bool | 纯 ID 样本：所有 `OOD_*` / `seq_Redundancy_*` / `TM-score_*` 均为 False（NaN 视为 False），由 `finalize_ood_columns.py --mode ind` 统一重算 |

---

## 诊断与质控工具

| 脚本 | 用途 |
|------|------|
| `scripts/13_verify_seq_structure_consistency.py` | 校验 splits 中 `aa_seq` 与 `pdbs/*.pdb` 提取序列的一致性（QC 收尾步骤，当前一致率 100.00%） |

---

## 归档清单（`output/scripts_archive/`）

| 脚本 | 归档原因 |
|------|----------|
| `02_resume_extract_binding_sites.py` | 一次性断点续跑 |
| `03_diagnose_remaining_files.py` | 一次性崩溃恢复诊断（单文件分类 OK/Timeout/Error/Crash） |
| `05_reprocess_timeouts.py` | 一次性超时重试（前 40 个最小文件仍全部超时，零产出后终止） |
| `06_create_splits.py` | 随机 80/10/10 划分，已被时间划分（05）取代 |
| `11_extract_protein_chains.py` | 慢版结构提取，已被 07（快速版）取代 |
| `13_move_intermediate_files.py` | 一次性中间文件迁移 |
| `15_patch_inconsistent_from_structure.py` | 序列-结构不一致一次性修复，修复已落数据，根因在 core/07 |
| `16_reextract_inconsistent_biopython.py` | 同上（Biopython 重提取） |
| `20_annotate_cath_domains.py` / `21_annotate_cath_domains_mmseqs2.py` | CATH 旧标注路线，已被 CATH+TED 并集（12 + cath-ted skill）取代 |
| `22_add_newfunction_longtail_ood.py` | EC 列已被 ec skill 统一口径取代、LongTail 与 11 重复 |
| `23_restore_idr_ood.py` | 一次性 IDR 灾难恢复 |
| `24_fix_seq_homology_easysearch.py` | 一次性 seq 同源口径修复（easy-search -s 7），重算用 homology skill |
| `26_add_longtail_ec_ood.py` | LongTailEC 初版单列，已统一为四列 all-not-in 口径（ec skill） |
| `26_fix_multimodel_pdbs.py` | 549 个多模型（NMR）pdb 一次性重提取，修复已落数据、根因已修入 07；注意其内部 `import_module('11_extract_protein_chains_fast')` 为改名前模块名（现 `07_extract_protein_chains.py`） |
| `fix_lbs_aaseq.py` | 文档未提及的一次性 aa_seq 修复，列名失效（`pdb_file`/`labels`），修复已落数据 |

---

## 注意事项

1. **PDB 数据版本**：使用本地 RCSB PDB 镜像（项目根 `data/pdb/`），建议记录下载时间以保障可复现性。
2. **结合位点定义**：配体重原子 4.5 Å 范围内的蛋白残基；配体范围排除水、核酸、标准氨基酸、127 种结晶添加剂（阶段二）。
3. **IDR 口径细节**：判定列为 residue-ratio > 0.3（clean=map）；`11_build_ood_splits.py` 内部还计算了"最长连续无序区占比"但未用于判定——早期文档中"连续无序区域 > 30%"的描述不准确，以阶段四 Step 6 为准。
4. **TM-score 档位与口径（2026-08-28 已统一）**：0.9~0.3 七档齐全，与 skill 标准阶梯一致；foldseek 缓存已按 easy-search 统一口径重算（旧 `search -a -e inf` 缓存含 21% 非法 >1 值，曾严重低估各档计数，对比细节见阶段四 Step 5 与 `output/foldseek_easysearch_compare/compare_report.txt`）。
5. **pdb_id 泄漏口径**：时间划分时点 train/val/test 按 `pdb_id` 零交集；极端长度重组（阶段四 Step 3）把 train/val 的极端长度链移入 test，使 516 / 284 个 test 极端长度行与 train / val 共享 `pdb_id`（不同链或同 PDB 的不同配体样本）。**Default 子集与 train/val 仍保持 `pdb_id` 零交集**，常规 ID/OOD 评估不受影响。
6. **`NA` 链 ID**：部分蛋白链 ID 为字符串 `"NA"`，读取 CSV 时必须 `keep_default_na=False`（各写回脚本统一 `pd.read_csv(dtype=str, keep_default_na=False)` + `to_csv(index=False, lineterminator="\n")`）。
7. **配体 ECFP4 存储格式**：`ligand_ecfp4` 为稀疏列表字符串（如 `[12,45,200]`），可用 `ast.literal_eval` 或 `json.loads` 转为列表再转 numpy 数组。
8. **OOD 比对方向**：序列/结构 OOD 均按 **Default Test vs Train+Val** 计算；验证集参与参考但不作为 OOD 评估对象；极端长度样本不参与这些 OOD 计算（对应列 NaN）。
9. **超时文件**：1,099 个大文件在 10 s/60 s 超时内无法完成解析已跳过（约占全部 PDB 的 0.4%），其潜在样本不在数据集中。
10. **pdbs/ 冗余文件**：723 个 PDB 文件对应已丢弃的无 SMILES 样本，不在 splits 中；如需严格对齐可按 `struct_file` 列过滤。
11. **>2000 剔除（无脚本环节）**：155 个超长样本（2,003–4,584 aa，全部原属 `OOD_ExtremeLong`）按全项目统一规则从 test 剔除，备份在 `output/len2000_removal_backup/`。
12. **pdbs 多模型（NMR）修复（2026-08-28，已落数据）**：全库扫描曾发现 549 个 pdb 文件将约 20–25 个模型无分隔串联写入同一文件（残基重复、原子序号超 99,999 列错位），其中 7 个文件因列错位产生幻影残基，曾导致 `struct_label` 出现 `-1`。根因已修入 `07_extract_protein_chains.py`（首个 `ENDMDL` 即停 / mmCIF 仅保留 model 1 / 带引号原子名去引号），受影响文件由归档脚本 `26_fix_multimodel_pdbs.py` 重提取（原文件备份于 `output/backup_before_multimodel_fix/`），`struct_label` 已重算、`-1` 归零。
13. **本次整理（2026-08）脚本适配记录**：09/10/11 三个脚本内部旧列名 `pdb_file`/`labels`/`file_type` 已适配为统一 schema（`struct_file`/`label`）；同时修复了 10 号脚本未定义 `ROOT` 的路径缺陷与 09 号脚本备份目录的失效路径（`ROOT/"ligand_binding_site"` → 任务目录 `output/`），09 号脚本补充了 CSV 读回布尔列的字符串归一化。这些脚本属历史一次性步骤，重放时应以阶段四 skill 链为准。
