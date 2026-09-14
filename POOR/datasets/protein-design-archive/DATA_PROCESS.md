# Protein Design Archive (PDA) 额外测试集构建文档（v2）

> 从 PDA 数据库提取人工设计蛋白链，为 POOR 任务构建**额外测试集**。
> 当前覆盖 **fold prediction**（`splits/cath_design.csv`）和 **func_prediction 的 EC 任务**（`splits/ec_design.csv`）。
>
> 数据源：`POOR/data/protein-design-archive/`（PDA 网站仓库快照，Nature Biotechnology 2025，核心数据 `backend/scripts/data.json`，1,879 条目 / 2,389 链）
>
> **v2 设计（2026-08 重建，替换 v1）**：
> 1. 样本池 = PDA 全部条目的 **chain_type D+U** 链（2,207 条），不再仅限 synthetic-construct 来源；**N（天然来源）链完全不使用**。
> 2. 设计集内部**仅精确去重**（完全相同序列合并），不做 95% 序列同一性聚类——微小序列差异也可能对应不同功能。
> 3. train/val 重叠剔除仍为 **PDB id 匹配 OR mmseqs max pident ≥ 95**。
> 4. EC 标签**只用该链自身的标注**（PDB 沉积者 ∪ SIFTS），不做聚类内转移。
>
> v1（synthetic-only + 95% 去重）产物归档于 `output/backup_v1/`。

**当前 splits 规模**（实测）：

| 文件 | 行数 | 说明 |
|------|------|------|
| `splits/cath_design.csv` | **122** | fold 额外测试集（29 个 Topology 标签），列结构与 `fold_classification/splits/cath_test.csv` 一致 |
| `splits/ec_design.csv` | **31** | EC 额外测试集（15 个 L3 标签，1 条多标签），列与 `func_prediction/splits/ec_test.csv` 完全一致 |

---

## Stage 1: 源数据获取与解析

### 1.1 设计蛋白链提取 — `scripts/01_extract_chains.py`

- 输入：`POOR/data/protein-design-archive/backend/scripts/data.json`
- 提取 `chain_type ∈ {D, U}` 的全部链记录 → **2,207 条**（D 1,223 / U 984），涉及 1,840 个 PDB 条目；N 链（182 条）跳过
- 逗号合并 `chain_id`（同一 entity 的相同拷贝，1,017 条记录）取第一条为代表
- `unique_id = {pdb}_{chain}`；序列取 `chain_seq_nat`
- 条目级 PDA CATH 标注一并带出：`cath_codes` / `cath_cat_set` / `label_candidate`（仅条目级 distinct C.A.T 唯一时无歧义，250 条）

**输出**：`output/design_chains_all.csv`

### 1.2 单链结构提取 — `scripts/02_extract_structures.py`

- 源文件：`POOR/data/pdb/{PDB}.pdb.gz`（优先）或 `{PDB}.cif.gz`（文件名大写回退）
- `.pdb` 源：按链列过滤 ATOM/HETATM → `.pdb`；`.cif` 源：纯文本行过滤 `_atom_site.auth_asym_id` → `.cif`（避免 BioPython 解析大 CIF 卡死）
- 多模型（NMR）条目在首个 `ENDMDL` 即停（只保留 model 1）
- 结果：**2,207 / 2,207 全部成功**

**输出**：`output/chains_all/`（暂存，去重前）

---

## Stage 2: 标签聚合与计算

### 2.1 Fold（CATH）标签

直接使用 §1.1 带出的**条目级** PDA CATH 标注。

> **粒度限制（重要）**：PDA 的 `cath_full` 来自其对设计结构做 Foldseek→CATH 搜索后的**条目级码集合**（上游 `merizo_to_json.py`），**不保留链级/结构域级对应关系**。因此 fold 测试集只使用**条目级 distinct C.A.T == 1** 的链（`label_candidate` 非空），多 C.A.T 条目不做消歧。

### 2.2 EC 标签（严格、链级、无转移）

由 `scripts/06_build_ec_testset.py` 前半段完成。每条链的 EC **只来自其自身条目**的标注：

- `.pdb`：COMPND 段 `MOLECULE:/CHAIN:/EC:` 子字段按链归属
- `.cif`：`_entity.pdbx_ec`（entity 级）经 `_atom_site.(auth_asym_id, label_entity_id)` 映射到链
- ∪ SIFTS `data/sifts/sifts_ec_annotations.json`（`{PDB}_{chain}` 大写）

**完全相同的序列**跨条目重复出现时，要求所有标注成员的 L3 标签集合一致，不一致则整组剔除（本次 0 组冲突，冲突清单会写 `output/ec_dropped_conflict.csv`）。

D+U 链中命中 **69 条**（`output/design_ec_annotations.csv`）→ 精确去重 **43 条**（0 冲突）。

---

## Stage 3: 去冗余与数据精化

### 3.1 精确去重 — `scripts/03_dedup_exact.py`

- **仅合并完全相同的序列**（`aa_seq` 精确匹配）：2,207 → **1,726 条**
- 组内代表：优先结构提取成功者，其次首条
- 代表链结构**复制**（非软链）到 `pdbs/`（1,726 个文件，供后续任务复用）

**输出**：`output/design_chains_dedup_exact.csv`、`pdbs/`

---

## Stage 4: 数据划分与 OOD 标注

> 本任务产出的是**额外测试集**，无 train/val 划分；本阶段 = 针对 fold / EC 两个任务分别做标签筛选、与对应 train/val 的重叠剔除、OOD 标注与 skill 统一收尾。

### 4.1 Fold 测试集筛选 — `scripts/04_build_fold_testset.py`

| 过滤 | 规则 | 数量 |
|------|------|------|
| 标签无歧义 | 条目级 distinct C.A.T 唯一（`label_candidate`） | 213 |
| 标签可用 | C.A.T ∈ fold train 的 1,163 个 Topology 标签 | 211 |
| PDB id 重叠剔除 | pdb_id ∈ fold train/val | −32 |
| 序列重叠剔除 | mmseqs2 `easy-search -s 7.5` vs train+val，max pident ≥ 95 | −88（并集 −89） |
| **最终** | | **122**（29 个 Topology 标签） |

**输出**：`output/cath_design_base.csv`、`output/removed_overlap.csv`、`output/design_mmseqs_vs_trainval.m8`（Step 4.2 复用）

### 4.2 Fold OOD 标注 — `scripts/05_annotate_ood.py`

列结构与 `fold_classification/splits/cath_test.csv` 一致。口径与 skill 链对齐：

- 序列同源 `seq_Redundancy_*` / `OOD_Orphan`（mmseqs2 easy-search vs fold train+val）与结构相似 `TM-score_*`（foldseek easy-search）：与 `.skills/homology-ood-annotation/`（`seq_homology_ood.py` / `struct_homology_ood.py`）口径一致
- `OOD_IDR`（metapredict v3、residue-ratio @0.3）：与 `.skills/ood-annotation-toolkit/idr_ood.py` fold 口径一致
- `OOD_SuperfamilyHoldout`：链候选 superfamily（PDA `cath_codes` 中 C.A.T == label 的 C.A.T.H）均未在 train/val 出现
- `OOD_NewEC`：SIFTS EC7 先验口径（fold 任务惯例；设计蛋白多无 EC7 注释，基本不适用）
- `OOD_ExtremeShort`（<60 aa）；`OOD_LongTail`（label train 频次 ≤ 10）

> **`Default` 全 True**：本表是额外测试集，但构建时不存在"划分前预挑出"的样本，按全项目 Default 本意整表标记 True，与 fold 任务"ES 行也是 Default=True"的惯例一致（2026-08-27 由全 False 改判，备份 `output/backup_before_default_true/`）。仅 ES 行（`OOD_ExtremeShort=True`，45 条）的其他依赖型 OOD 列未计算（留空）——**已知口径滞后点**：fold 现行口径已补齐 ES 行的依赖型 OOD，此处未同步补齐。

| OOD 轴 | True 数 | 说明 |
|--------|--------:|------|
| OOD_ExtremeShort | 45 | 37%（设计蛋白富含迷你蛋白/短肽） |
| OOD_Orphan | 37 | 占 77 条依赖行 48% |
| seq_Redundancy_30 | 39 | 51% |
| TM-score_0.9 / 0.5 | 44 / 1 | 新折叠极少 |
| OOD_SuperfamilyHoldout / OOD_NewEC / OOD_LongTail | 0 | — |
| OOD_IDR | 2 | — |
| **InD** | **8** | 6.6% |

**校验**：PDB id 零重叠、max pident 94.5 < 95、无完全相同序列、label ⊆ train、struct_file 100% 存在、InD 双向一致、Default 全 True、unique_id 唯一。

**输出**：`splits/cath_design.csv`、`output/design_foldseek_vs_trainval.m8`

### 4.3 EC 测试集筛选 — `scripts/06_build_ec_testset.py`（后半段）

- L4→L3 归一化；多标签 `;` 连接；要求所有 L3 ∈ EC train 260 标签：43 → **43 条**
- 重叠剔除 vs EC train/val（PDB id 匹配 OR mmseqs2 `easy-search -s 7` max pident ≥ 95）：pdb_id −2、pident≥95 −12，并集 −12 → **最终 31 条**（15 个 L3 标签）

**输出**：`output/ec_design_base.csv`、`output/ec_removed_overlap.csv`、`output/design_ec_mmseqs_vs_trainval.m8`

### 4.4 EC OOD 标注 — `scripts/07_annotate_ec_ood.py`

列与 `func_prediction/splits/ec_test.csv` **完全一致**（无 idr_ratio 列）：

- `seq_Redundancy_*` / `OOD_Orphan`（mmseqs2 `-s 7` vs EC train+val，复用 §4.3 的 .m8）、`TM-score_*`（foldseek vs EC train+val 结构）：与 `.skills/homology-ood-annotation/` 口径一致
- `OOD_IDR`（metapredict v3 residue-ratio @0.3）：与 `.skills/ood-annotation-toolkit/idr_ood.py` func 口径一致
- `OOD_NewEC_L4`（任一 4 级 EC ∉ train+val L4 集合，无 L4 者 False）；`OOD_LongTail`（train 频次 ≤10）；`OOD_Combinatorial`（标签对未在 train+val 共现）
- `Default` 全 True（同 §4.2）；本批无极端长度样本，全部行的 OOD 列均已计算

| OOD 轴 | True 数 | 占 31 |
|--------|--------:|------:|
| seq_Redundancy_90 / 30 | 26 / 2 | 84% / 6% |
| OOD_Orphan | 2 | 6% |
| TM-score_0.9 / 0.5 | 7 / 4 | 23% / 13% |
| OOD_IDR | 1 | 3% |
| NewEC_L4 / LongTail / Combinatorial / ES / EL | 0 | — |
| **InD** | **4** | 12.9% |

**校验（全部通过）**：与 EC train/val PDB id 零重叠、无完全相同序列、max pident 94.6 < 95、label ⊆ train 标签集（15 个 L3）、struct_file 100% 存在、列与 ec_test.csv 完全一致、布尔列无空值、unique_id 唯一、Orphan 行 seq_Redundancy 全 True、输入池仅 D+U。

**输出**：`splits/ec_design.csv`、`output/design_ec_foldseek_vs_trainval.m8`

### 4.5 skill 统一收尾（脚本均在 `.skills/` 下）

1. **LongTailEC 四列**（`ec-function-ood-annotation/add_unified_longtail_ec.py`，批量驱动）：作用于 `cath_design.csv`——设计链有 EC 注释且其**所有** EC（L3/L4）在 fold `cath_train.csv` 中频次都 < 5 / < 10（含 0，all-not-in）→ True；ES 行置 NaN；非 ES 行计算（含 NewEC prior 行，同 OOD_LongTail 惯例）。设计蛋白极少有 SIFTS EC 注释：77 条非 ES 行中仅 3 条有 EC，**仅 `OOD_LongTail_EC_L4_le10` True = 1**，其余三列全 0（统计 `output/longtail_ec_unified_summary.json`）。train 参考频次与 fold 任务一致（精确键，3,053 个 L4 类别）。
2. **InD 重算**（`ood-annotation-toolkit/finalize_ood_columns.py --mode ind`）：两个 splits 文件的 `InD` 均由其统一重算（纯 ID：所有 `OOD_*` / `seq_Redundancy_*` / `TM-score_*` 均为 False）。`--mode nan` 对本任务无影响（只作用于 Default=False 行，本任务全 True）。
3. cath_design 的 `InD`（8）并集含 `OOD_LongTail_EC_*` 列，故较旧 Pure_ID=9 少 1。

---

## 目录结构

```
protein-design-archive/
├── DATA_PROCESS.md                      # 本文档（v2）
├── scripts/
│   ├── 01_extract_chains.py             # D+U 全量链清单
│   ├── 02_extract_structures.py         # 单链结构提取
│   ├── 03_dedup_exact.py                # 精确去重（仅完全相同序列）
│   ├── 04_build_fold_testset.py         # fold 筛选
│   ├── 05_annotate_ood.py               # fold OOD 标注
│   ├── 06_build_ec_testset.py           # EC 标注提取 + 筛选
│   └── 07_annotate_ec_ood.py            # EC OOD 标注
├── output/
│   ├── design_chains_all.csv            # D+U 链清单（2,207）
│   ├── chains_all/                      # 去重前全部单链结构（2,207）
│   ├── design_chains_dedup_exact.csv    # 精确去重总表（1,726，后续任务复用）
│   ├── cath_design_base.csv / removed_overlap.csv
│   ├── design_mmseqs_vs_trainval.m8 / design_foldseek_vs_trainval.m8
│   ├── design_ec_annotations.csv        # 链级 EC 命中（69）
│   ├── ec_design_base.csv / ec_removed_overlap.csv
│   ├── design_ec_mmseqs_vs_trainval.m8 / design_ec_foldseek_vs_trainval.m8
│   ├── cath_design_longtail_ec_summary.json   # 初版单列 LongTailEC 统计
│   ├── longtail_ec_unified_summary.json       # 统一口径 LongTailEC 四列统计
│   ├── dup_multimodel_files.json        # 多模型修复时的受影响文件清单
│   ├── s04_tmp/ ... s07_tmp/            # Step 4–7 中间缓存（foldseek 软链等）
│   ├── scripts_archive/                 # 归档脚本（见归档清单）
│   ├── backup_before_longtail_aug/      # ↓ LongTailEC 统一口径多轮迭代与
│   ├── backup_before_longtail_unify/    #   InD 重算前的 splits 备份
│   ├── backup_before_longtail_allnotin/
│   ├── backup_before_ind/
│   ├── backup_before_multimodel_fix/    # 多模型修复前原始 pdb 备份
│   └── backup_v1/                       # v1（synthetic-only + 95% 去重）归档
├── pdbs/                                # 1,726 个精确去重后设计链结构
└── splits/
    ├── cath_design.csv                  # fold 额外测试集（122）
    └── ec_design.csv                    # EC 额外测试集（31）
```

## 归档清单（`output/scripts_archive/`，含原因）

| 脚本 | 归档原因 |
|------|---------|
| `fix_multimodel_pdbs.py` | 一次性修复：203 个 NMR 多模型 pdb 将全部模型无 `MODEL` 分隔串联写入，在首个非连续重复残基处截断（备份 `output/backup_before_multimodel_fix/`）；修复已落数据，根因已修入 `02_extract_structures.py`（首个 `ENDMDL` 即停）。验证：修复后全库重扫 0 残留重复、0 格式错误行，203 个文件均可被 Bio.PDB 解析且为单模型，52 个 `.cif` 全部单模型 |
| `08_add_longtail_ec_ood.py` | 初版单列 LongTailEC（any 语义、仅 `OOD_LongTail_EC_L4` 一列），已被 ec skill 统一 all-not-in 四列口径取代（§4.5） |

## 注意事项

1. **条目级 CATH 粒度**：PDA 标注无链级对应，多 C.A.T 条目（多数）未纳入 fold 测试集。
2. **标签是计算预测值/沉积者标注**：fold 的 CATH 标签来自 PDA 的 Foldseek→CATH 搜索；EC 标签来自 PDB 沉积者 ∪ SIFTS，未必经独立实验验证。
3. **U 链语义**：unknown 来源的链按"PDA 收录即设计蛋白"处理，可能混有少量未标注清楚的天然链。
4. **重叠剔除语义**：PDB id 重叠 OR 与 train/val max pident ≥ 95% 剔除；与 fold/EC **test** 集的样本不去重（都是测试集）。
5. **设计集内部仅精确去重**：≥95% 同源的设计变体（如同一系列的多代设计）会同时保留在测试集中——评估时需注意这些样本间高度相关，不应视为独立样本。
6. **OOD_NewEC 对设计蛋白基本不适用**（fold 任务；设计蛋白多无 UniProt/EC7 注释）。
7. **cath_design ES 行口径滞后（已知偏差）**：45 条 ES 行的其他依赖型 OOD 列未计算（留空）；fold 现行口径已补齐 ES 行依赖型 OOD，此处未同步补齐；`finalize_ood_columns.py --mode nan` 对其无影响（只作用于 Default=False 行）。
8. **`Default` 全 True**：两个额外测试集构建时不存在"划分前预挑出"的样本，按 Default 本意整表标记 True（2026-08-27 由全 False 改判，备份 `output/backup_before_default_true/`）。
9. **历史演进**：初版 LongTailEC 单列脚本（any 语义）已归档；2026-08-27 起全项目统一为 all-not-in 四列口径（§4.5），InD 由 `finalize_ood_columns.py --mode ind` 统一重算。
