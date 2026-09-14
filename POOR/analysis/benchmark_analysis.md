# POOR 论文数据分析与实验规划

## 第一部分：无需训练模型的数据分析（已完成 & 待完成）

### 一、已完成分析

#### 1. 数据集规模统计表

**状态**: ✅ 已完成（2026-09-01 重算）

**文件**: `analysis/output/dataset_scale_stats.json`

**脚本**: `analysis/scripts/dataset_scale_stats.py`

**核心结果**:

### Classification Tasks

| Task | Train | Val | Test | Total | Labels | Seq Length (median) | Extreme in Train/Val? |
|------|-------|-----|------|-------|--------|---------------------|----------------------|
| fold_classification | 43,473 | 6,132 | 12,911 | 62,516 | 1,163 classes (C.A.T) | 9-960 (med=133) | ⚠️ YES (2,940, 5.9%，全部 <60) |
| func_prediction (EC) | 15,389 | 1,687 | 5,047 | 22,123 | 260 classes | 3-1997 (med=298) | ✅ No |
| func_prediction (GO-MF) | 18,179 | 2,084 | 5,609 | 25,872 | 447 classes | 3-1997 (med=236) | ✅ No |
| func_prediction (GO-CC) | 13,878 | 1,730 | 5,613 | 21,221 | 264 classes | 3-1997 (med=222) | ✅ No |
| func_prediction (GO-BP) | 15,232 | 1,742 | 5,142 | 22,116 | 884 classes | 3-1997 (med=228) | ✅ No |
| ppi_prediction | 58,968 | 8,394 | 16,264 | 83,626 | train pos=29,484 / neg=29,484 | 3-1997 (med=280) | ⚠️ YES (1,803, 2.7%) |

### Regression Tasks

| Task | Train | Val | Test | Total | Label Range | Seq Length (median) | Extreme in Train/Val? |
|------|-------|-----|------|-------|-------------|---------------------|----------------------|
| enzyme_kinetics_prediction | 15,479 | 2,150 | 5,002 | 22,631 | [-6.00, 6.00] (μ=0.88) | 16-1928 (med=379) | ✅ No |
| enzyme_optimal_ph | 7,794 | 866 | 2,955 | 11,615 | [1.50, 12.50] (μ=7.13) | 23-1476 (med=361) | ✅ No |
| ligand_binding_affinity | 8,366 | 930 | 2,367 | 11,663 | [0.40, 15.22] (μ=6.26) | 1-1504 (med=231) | ⚠️ YES (1,212, 13.0%) |

### Per-Residue Tasks

| Task | Train | Val | Test | Total | Label Type | Seq Length (median) | Extreme in Train/Val? |
|------|-------|-----|------|-------|------------|---------------------|----------------------|
| ligand_binding_site | 71,093 | 7,938 | 17,225 | 96,256 | binary per-residue | 4-1991 (med=291) | ✅ No |
| ppis_prediction | 6,043 | 672 | 1,980 | 8,695 | binary per-residue | 7-1484 (med=220) | ✅ No |
| ss_prediction | 29,622 | 3,317 | 8,372 | 41,311 | 8-state + 3-state DSSP | 31-1836 (med=222) | ⚠️ YES (952, 2.9%，全部 <60) |

**关键发现**:
- 总样本量：**429,645** 条（Train 303,516 / Val 37,642 / Test 88,487）
- 序列长度范围：Fold 序列最短（med≈133aa）；回归与结合位点类任务偏长；剔除 >2000 超长样本后，全体单链长度上限为 2000aa
- 8 个任务的 train/val 无极端长度样本（<60 或 >1000）；4 个任务（fold_classification / ppi_prediction / ligand_binding_affinity / ss_prediction）的 train/val 仍含极端长度单链，其中 **LBA 占比最高**（train/val 的 13.0% 样本含极端长度单链，几乎全为多链复合物中的 >1000aa 长链）；fold 的 2,940 条与 SSP 的 952 条全部为 <60 短样本（只标注不剔除 / 正常划分的设计使然）；PPI 重划分后极端长度占比已降至 2.7%
- func_prediction 各任务 test 的最短序列（3aa）来自 2026-08 找回的 ExtremeShort 测试样本（仅划入 test，train/val 不受影响）；ss_prediction 的 <60 样本（8,941 条，最短 1aa）为 2026-09-01 恢复、09-02 按原时间 cutoff 正常划分到 train/val/test（详见下文 2026-09-01/02 注记）
- PPI 唯一的空链样本（O75746 同源二聚体）已于 2026-08-28 从 UniProt 补全序列并修复（见 `datasets/ppi_prediction/DATA_PROCESS.md` Step 5 注记），最小长度 3 为真实短链；LBA 长度按 `|` 拆链统计，最小值 1 为复合物中沉积的短肽/占位链片段（原始 PDB 沉积如此，见 `datasets/ligand_binding_affinity/` 核查记录）
- 所有任务均包含结构文件引用（pdb/cif）
- **2026-09-01/02 极端长度口径修正**：① SSP 恢复 <60 链（22,231 条 → 95% 去重 + 与现有 splits 交叉去重 → 8,941 条），并于 09-02 按原时间 cutoff（train≤2019-06-06、val≤2021-02-09）正常划分到 train/val/test（+5,672 / +892 / 2,377），test 短链行 `Default=True` 且 seq/TM/Orphan/IDR/LongTail/NewEC/LongTailEC/Fold-SuperfamilyHoldout 全部按新参考集重算；另剔除 1 条 59% UNK 残基的低质量短链 7UYL_K（序列提取跳过 UNK 而结构保留导致标注错帧）；随后按用户要求仅保留 31–59aa 短链，删除 ≤30aa 共 7,695 条（train 4,838 / val 774 / test 2,083），OOD 各列随之再次按新参考集重算；② pH train 的 3 条 aa_seq<60 样本（P05959/B1KRG2/P0C7A9）移入 test（ES 24→27）；③ fold train 删除唯一 >1000 样本 3w3uA01（1,047aa，建池时按 CATH `length_expected=999` 口径未触发剔除的边界残留）

---

#### 2. OOD场景覆盖度矩阵

**状态**: ✅ 已完成（2026-09-01 重算）

**文件**: `analysis/output/ood_coverage_matrix.csv` / `.md`

**脚本**: `analysis/scripts/ood_coverage_matrix.py`

**核心结果**（各 OOD 场景在 Test 集中的计数与占比；占比基于全测试集 N_test，`-` 表示该任务无此场景。）:

| OOD Scenario | kcat | pH | CATH | EC | GO-BP | GO-CC | GO-MF | LBA | LBS | PPI | PPIS | SSP |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Test set size | 5,002 | 2,955 | 12,911 | 5,047 | 5,142 | 5,613 | 5,609 | 2,367 | 17,225 | 16,264 | 1,980 | 8,372 |
| Default (ID) subset | 4412 (88.2%) | 2,717 (91.9%) | 12911 (100.0%) | 4243 (84.1%) | 4121 (80.1%) | 4256 (75.8%) | 4609 (82.2%) | 2367 (100.0%) | 13723 (79.7%) | 14939 (91.9%) | 1685 (85.1%) | 8,290 (99.0%) |
| InD subset (no OOD flag) | 366 (7.3%) | 50 (1.7%) | 1805 (14.0%) | 151 (3.0%) | 172 (3.3%) | 192 (3.4%) | 194 (3.5%) | 1091 (46.1%) | 5209 (30.2%) | 1944 (12.0%) | 134 (6.8%) | 1,123 (13.4%) |
| Extreme short (<60 aa) | 59 (1.2%) | 27 (0.9%) | 767 (5.9%) | 478 (9.5%) | 757 (14.7%) | 1090 (19.4%) | 753 (13.4%) | 15 (0.6%) | 1974 (11.5%) | 715 (4.4%) | 261 (13.2%) | 293 (3.5%) |
| Extreme long (>1000 aa) | 531 (10.6%) | 211 (7.1%) | - | 326 (6.5%) | 264 (5.1%) | 267 (4.8%) | 247 (4.4%) | 303 (12.8%) | 1528 (8.9%) | 1325 (8.1%) | 34 (1.7%) | 82 (1.0%) |
| Seq identity < 90% | 3936 (78.7%) | 2,540 (86.0%) | 10257 (79.4%) | 3989 (79.0%) | 3834 (74.6%) | 3945 (70.3%) | 4300 (76.7%) | 575 (24.3%) | 7908 (45.9%) | 10655 (65.5%) | 1457 (73.6%) | 6,894 (82.3%) |
| Seq identity < 80% | 3415 (68.3%) | 2,396 (81.1%) | 9118 (70.6%) | 3627 (71.9%) | 3434 (66.8%) | 3511 (62.6%) | 3856 (68.7%) | 498 (21.0%) | 7115 (41.3%) | 9856 (60.6%) | 1358 (68.6%) | 6,139 (73.3%) |
| Seq identity < 70% | 3063 (61.2%) | 2,227 (75.4%) | 8114 (62.8%) | 3318 (65.7%) | 3044 (59.2%) | 3150 (56.1%) | 3433 (61.2%) | 454 (19.2%) | 6418 (37.3%) | 9027 (55.5%) | 1298 (65.6%) | 5,548 (66.3%) |
| Seq identity < 60% | 2645 (52.9%) | 2,025 (68.5%) | 7055 (54.6%) | 2880 (57.1%) | 2595 (50.5%) | 2769 (49.3%) | 2939 (52.4%) | 409 (17.3%) | 5547 (32.2%) | 8098 (49.8%) | 1202 (60.7%) | 5,035 (60.1%) |
| Seq identity < 50% | 2145 (42.9%) | 1,674 (56.6%) | 5586 (43.3%) | 2246 (44.5%) | 2113 (41.1%) | 2348 (41.8%) | 2371 (42.3%) | 368 (15.5%) | 4562 (26.5%) | 6982 (42.9%) | 1066 (53.8%) | 4,327 (51.7%) |
| Seq identity < 40% | 1218 (24.4%) | 1,197 (40.5%) | 3876 (30.0%) | 1449 (28.7%) | 1599 (31.1%) | 1869 (33.3%) | 1735 (30.9%) | 294 (12.4%) | 3191 (18.5%) | 5496 (33.8%) | 879 (44.4%) | 3,283 (39.2%) |
| Seq identity < 30% | 386 (7.7%) | 644 (21.8%) | 2353 (18.2%) | 729 (14.4%) | 1059 (20.6%) | 1346 (24.0%) | 1124 (20.0%) | 241 (10.2%) | 1859 (10.8%) | 3653 (22.5%) | 602 (30.4%) | 2,016 (24.1%) |
| Orphan (no significant hit) | 173 (3.5%) | 422 (14.3%) | 1978 (15.3%) | 498 (9.9%) | 834 (16.2%) | 1059 (18.9%) | 902 (16.1%) | 181 (7.6%) | 1310 (7.6%) | 3017 (18.6%) | 458 (23.1%) | 1,566 (18.7%) |
| TM-score < 0.9 | 956 (19.1%) | 926 (31.3%) | 4022 (31.2%) | 1361 (27.0%) | 1714 (33.3%) | 1958 (34.9%) | 1839 (32.8%) | 767 (32.4%) | 3017 (17.5%) | 8172 (50.2%) | 770 (38.9%) | 2,863 (34.2%) |
| TM-score < 0.8 | 339 (6.8%) | 411 (13.9%) | 1755 (13.6%) | 574 (11.4%) | 898 (17.5%) | 1156 (20.6%) | 909 (16.2%) | 662 (28.0%) | 1310 (7.6%) | 5345 (32.9%) | 426 (21.5%) | 1,356 (16.2%) |
| TM-score < 0.7 | 106 (2.1%) | 171 (5.8%) | 839 (6.5%) | 300 (5.9%) | 510 (9.9%) | 717 (12.8%) | 529 (9.4%) | 447 (18.9%) | 691 (4.0%) | 3771 (23.2%) | 224 (11.3%) | 696 (8.3%) |
| TM-score < 0.6 | 13 (0.3%) | 40 (1.4%) | 495 (3.8%) | 151 (3.0%) | 319 (6.2%) | 438 (7.8%) | 313 (5.6%) | 306 (12.9%) | 376 (2.2%) | 2843 (17.5%) | 99 (5.0%) | 347 (4.1%) |
| TM-score < 0.5 | 6 (0.1%) | 8 (0.3%) | 348 (2.7%) | 105 (2.1%) | 217 (4.2%) | 297 (5.3%) | 204 (3.6%) | 201 (8.5%) | 241 (1.4%) | 2186 (13.4%) | 47 (2.4%) | 188 (2.2%) |
| TM-score < 0.4 | 6 (0.1%) | 4 (0.1%) | 273 (2.1%) | 57 (1.1%) | 147 (2.9%) | 201 (3.6%) | 142 (2.5%) | 135 (5.7%) | 189 (1.1%) | 1510 (9.3%) | 28 (1.4%) | 127 (1.5%) |
| TM-score < 0.3 | 1 (0.0%) | 3 (0.1%) | 237 (1.8%) | 31 (0.6%) | 99 (1.9%) | 145 (2.6%) | 90 (1.6%) | 79 (3.3%) | 136 (0.8%) | 943 (5.8%) | 20 (1.0%) | 102 (1.2%) |
| IDR-containing | 48 (1.0%) | 39 (1.3%) | 401 (3.1%) | 27 (0.5%) | 157 (3.1%) | 175 (3.1%) | 168 (3.0%) | 93 (3.9%) | 140 (0.8%) | 3743 (23.0%) | 28 (1.4%) | 159 (1.9%) |
| New EC L4 | 853 (17.1%) | 647 (21.9%) | 336 (2.6%) | 381 (7.5%) | 118 (2.3%) | 136 (2.4%) | 183 (3.3%) | 46 (1.9%) | 271 (1.6%) | 203 (1.2%) | 144 (7.3%) | 245 (2.9%) |
| New EC L3 | 3 (0.1%) | 20 (0.7%) | - | - | 8 (0.2%) | 9 (0.2%) | 2 (0.0%) | 8 (0.3%) | 19 (0.1%) | 10 (0.1%) | 9 (0.5%) | 1 (0.0%) |
| Long-tail (label-frequency) | - | - | 1081 (8.4%) | 163 (3.2%) | 684 (13.3%) | 235 (4.2%) | 295 (5.3%) | - | 330 (1.9%) | - | 48 (2.4%) | 436 (5.2%) |
| Long-tail EC L3 ≤5 | 28 (0.6%) | 118 (4.0%) | 368 (2.9%) | - | 55 (1.1%) | 44 (0.8%) | 34 (0.6%) | 24 (1.0%) | 39 (0.2%) | 17 (0.1%) | 62 (3.1%) | 37 (0.4%) |
| Long-tail EC L3 ≤10 | 75 (1.5%) | 221 (7.5%) | 415 (3.2%) | - | 83 (1.6%) | 97 (1.7%) | 106 (1.9%) | 37 (1.6%) | 68 (0.4%) | 72 (0.4%) | 137 (6.9%) | 77 (0.9%) |
| Long-tail EC L4 ≤5 | 1524 (30.5%) | 1,763 (59.7%) | 1459 (11.3%) | 1155 (22.9%) | 370 (7.2%) | 431 (7.7%) | 483 (8.6%) | 189 (8.0%) | 803 (4.7%) | 444 (2.7%) | 515 (26.0%) | 738 (8.8%) |
| Long-tail EC L4 ≤10 | 2163 (43.2%) | 2,208 (74.7%) | 2199 (17.0%) | 1752 (34.7%) | 509 (9.9%) | 550 (9.8%) | 736 (13.1%) | 344 (14.5%) | 1142 (6.6%) | 698 (4.3%) | 659 (33.3%) | 1,118 (13.4%) |
| Long-tail pH bin ≤50 | - | 25 (0.8%) | - | - | - | - | - | - | - | - | - | - |
| Long-tail pH bin ≤100 | - | 58 (2.0%) | - | - | - | - | - | - | - | - | - | - |
| Long-tail Aff bin ≤50 | - | - | - | - | - | - | - | 19 (0.8%) | - | - | - | - |
| Long-tail Aff bin ≤100 | - | - | - | - | - | - | - | 65 (2.7%) | - | - | - | - |
| Long-tail kcat bin ≤50 | 7 (0.1%) | - | - | - | - | - | - | - | - | - | - | - |
| Long-tail kcat bin ≤100 | 133 (2.7%) | - | - | - | - | - | - | - | - | - | - | - |
| Combinatorial (new label pair) | - | - | - | 58 (1.1%) | 808 (15.7%) | 472 (8.4%) | 449 (8.0%) | - | - | - | - | - |
| Fold holdout | 43 (0.9%) | 91 (3.1%) | - | 31 (0.6%) | 50 (1.0%) | 57 (1.0%) | 31 (0.6%) | 38 (1.6%) | 15 (0.1%) | 44 (0.3%) | 65 (3.3%) | 20 (0.2%) |
| Superfamily holdout | 160 (3.2%) | 311 (10.5%) | 761 (5.9%) | 103 (2.0%) | 151 (2.9%) | 185 (3.3%) | 186 (3.3%) | 73 (3.1%) | 135 (0.8%) | 265 (1.6%) | 235 (11.9%) | 96 (1.1%) |

> 注：① 列名已按统一 schema 重命名（`OOD_NewFunction` → `OOD_NewEC_L4/L3`、`OOD_LongTailFunction` → `OOD_LongTail_EC_L*` 等；CATH 任务的 `OOD_NewEC` 归入 New EC L4 行）。② 序列同源性 / 结构相似性 / Orphan 等默认在 `Default=True` 子集内判定，非 Default（极端长度先验划分）样本对应列为 NaN，计数按 False 处理、占比分母仍为全测试集。③ `InD subset` 为不属于任何 OOD 场景的样本（显式 `InD` 列）。④ NewEC / LongTailEC 为 all-not-in 口径（以 train 或 train+val 为参考集，EC L3/L4 层级，频次阈值 ≤5/≤10 且含 0）。⑤ `PPI` 为 pair-level（16,264 条），`PPIS` 为残基级 interface（1,980 条），二者 IDR 等比例差异源于 pair-level 判定规则（任一链触发即 True）。⑥ Fold/Superfamily holdout 在非 CATH 任务上基于 CATH+TED 域标注判定（见 `.skills/cath-ted-domain-annotation/`）。

---

#### 3. OOD 正交性 / OOD 场景标记重叠度

**状态**: ✅ 已完成（2026-09-01 重算）

**文件**: `analysis/output/ood_overlap_matrix.csv`

**脚本**: `analysis/scripts/ood_coverage_matrix.py`（同一脚本输出）

**方法**: 对每对 OOD 轴在同一任务的 test 集上计算 Jaccard 相似度
`J(A,B) = |A∩B| / |A∪B|`。`J=0` 表示两个 OOD 集合完全不相交；`J` 接近 1
表示二者高度重合。重叠度过高的轴不宜同时作为独立评估维度。

**核心结果**:

- 共计算 301 对 OOD 轴组合；整体 Jaccard 均值 **0.071**。
- **完全正交（J = 0.0）** 的有 **121 对**，例如各任务的 `OOD_ExtremeShort` 与 `OOD_ExtremeLong`、极端长度与 EC 功能新颖性、CATH holdout 与组合泛化等。
- **高重叠（J > 0.5）** 共 13 对，其中 11 对出现在 **Orphan（无同源命中）** 与 **seq_Redundancy_30（序列同一性 <30%）** 之间——这是预期的，因为 Orphan 是 seq_Redundancy_30 的极限子集：

| 任务 | 重叠轴对 | Jaccard | Orphan 且 <30% | 仅 Orphan | 仅 <30% |
|------|----------|---------|----------------:|----------:|--------:|
| CATH | OOD_Orphan ∩ seq_Redundancy_30 | 0.841 | 1,978 | 0 | 375 |
| PPI  | OOD_Orphan ∩ seq_Redundancy_30 | 0.826 | 3,017 | 0 | 636 |
| GO-MF | OOD_Orphan ∩ seq_Redundancy_30 | 0.802 | 902 | 0 | 222 |
| GO-BP | OOD_Orphan ∩ seq_Redundancy_30 | 0.788 | 834 | 0 | 225 |
| GO-CC | OOD_Orphan ∩ seq_Redundancy_30 | 0.787 | 1,059 | 0 | 287 |
| PPIS | OOD_Orphan ∩ seq_Redundancy_30 | 0.761 | 458 | 0 | 144 |
| SSP  | OOD_Orphan ∩ seq_Redundancy_30 | 0.777 | 1,566 | 0 | 450 |
| LBA  | OOD_Orphan ∩ seq_Redundancy_30 | 0.751 | 181 | 0 | 60 |
| LBS  | OOD_Orphan ∩ seq_Redundancy_30 | 0.705 | 1,310 | 0 | 549 |
| EC   | OOD_Orphan ∩ seq_Redundancy_30 | 0.683 | 498 | 0 | 231 |
| pH   | OOD_Orphan ∩ seq_Redundancy_30 | 0.655 | 422 | 0 | 222 |

另外 2 对高重叠出现在 LBA（复合物 multimer TM 口径下，结构新颖性与序列低同源度的耦合增强）：

| 任务 | 重叠轴对 | Jaccard | 交集 | 仅前者 | 仅后者 |
|------|----------|---------|-----:|-------:|-------:|
| LBA  | seq_Redundancy_30 ∩ TM-score_0.5 | 0.524 | 152 | 89 | 49 |
| LBA  | OOD_Orphan ∩ TM-score_0.5 | 0.522 | 131 | 50 | 70 |

- 除上述 13 对外，其余 288 对 Jaccard **均 ≤ 0.448**（最高为 kcat 的
  Orphan–seq30，同为子集关系但未过 0.5 线），绝大多数
  （75% 分位数）≤ 0.055。

**按任务平均 Jaccard**（数值越低表示该任务内部 OOD 轴越相互独立）：

| 任务 | 平均 Jaccard | 最高 Jaccard | 说明 |
|------|-------------:|-------------:|------|
| kcat | 0.042 | 0.448 | 序列侧 2026-08-30 统一后 Orphan–seq30 降至 0.448（<0.5），任务内最高重叠已低于高重叠线 |
| LBS  | 0.047 | 0.705 | 最高为 Orphan–seq30（子集关系）；TM-score 改 easy-search 口径后各档间重叠上升，平均值由 0.033 升至 0.048（现 0.047） |
| PPIS | 0.055 | 0.761 | 受 Orphan–seq30 驱动；TM-score 改 easy-search 口径后平均值由 0.077 降至 0.055；2026-08-29 序列侧（seq_Redundancy/Orphan）统一 easy-search 后均值/极值不变（0.055 / 0.761） |
| GO-CC | 0.060 | 0.787 | 功能预测；2026-08-30 序列侧统一后 Orphan–seq30 进入高重叠（子集关系） |
| SSP  | 0.070 | 0.777 | 受 Orphan–seq30 驱动；≤30aa 短链删除后，此前的短链耦合轴组（ES/IDR/TM0.5/Orphan）消失，仅剩 Orphan–seq30 子集关系 |
| GO-MF | 0.060 | 0.802 | 功能预测；同上 |
| EC   | 0.064 | 0.683 | 功能预测；同上 |
| pH   | 0.067 | 0.655 | 受 Orphan–seq30 驱动 |
| GO-BP | 0.069 | 0.788 | 功能预测；同上 |
| CATH | 0.105 | 0.841 | 受 Orphan–seq30 驱动 |
| PPI  | 0.106 | 0.826 | 受 Orphan–seq30 驱动 |
| LBA  | 0.128 | 0.751 | 受 Orphan–seq30 驱动；复合物 multimer TM 口径后平均值由 0.116 升至 0.128 |

**关键发现**: 去掉 Orphan–低同源度之间的包含关系后，各 OOD 维度基本相互独立，意味着它们捕获的是**互补的分布偏移类型**，适合作为多轴 OOD 评估基准。

**论文用途**: Figure 2（OOD Overlap Matrix）或 Supplementary Table

---

#### 4. 标签长尾分布

**状态**: ✅ 已完成（2026-09-01 重算（数值不变））

**文件**: `analysis/output/label_distribution.json`

**脚本**: `analysis/scripts/label_distribution.py`（train 集标签组合频率；SSP 为逐残基 token）

**核心结果**:

| 任务 | 标签数 | Head 20%覆盖 | Tail 50%覆盖 | Singletons比例 | Zipf系数 |
|------|--------|-------------|-------------|---------------|----------|
| EC | 558 | 88.0% | 2.6% | 31.7% | -1.652 |
| GO-MF | 3,436 | 81.0% | 9.5% | 63.6% | -0.995 |
| GO-CC | 2,506 | 84.8% | 9.0% | 75.6% | -0.814 |
| GO-BP | 3,191 | 81.1% | 10.5% | 70.1% | -0.924 |
| SSP | 9 | 53.4% | 7.7% | 0.0% | -1.762 |
| Fold | 1,163 | 86.6% | 3.8% | 14.0% | -1.445 |

**类别不平衡（二分类任务）**:

| 任务 | Train Pos% | Test Pos% |
|------|-----------|-----------|
| PPI | 50.0% | 50.0% |
| PPIS（残基级） | 22.8% | 22.9% |

**回归任务标签分布**:

| 任务 | Train Mean±Std | Test Mean±Std | Range |
|------|---------------|---------------|-------|
| LBA | 6.26 ± 1.97 | 6.11 ± 1.80 | [0.40, 15.22] |
| kcat | 0.88 ± 1.67 | 0.93 ± 1.62 | [-6.00, 6.00] |
| pH | 7.13 ± 1.28 | 7.39 ± 1.08 | [1.50, 12.50] |

**关键发现**:
- GO任务极度长尾：>60%标签仅出现1次（singletons）
- GO-CC最极端：75.6% singletons，Zipf=-0.814（最平缓）
- EC相对均衡：31.7% singletons，Zipf=-1.652（最陡峭）
- Fold 重划分为 C.A.T 粒度后类别数从 41 增至 1,163，singletons 14.0%，长尾程度介于 EC 与 GO 之间
- SSP 标签为 9 个有效 DSSP 状态（0–8，旧版"16"为未剥离 `[`/`]` 的 token 计数瑕疵，已修正）
- PPIS 残基级阳性率 22.8%（旧版 99.2% 为统计口径错误——误按"含任意阳性残基的样本比例"计算，已修正）， PPI 样本级严格 1:1 平衡

**论文用途**: Figure 3（Label Long-Tail Distribution / Zipf Curves）

---

#### 5. 蛋白质实验结构数据占比

**状态**: ✅ 已完成（2026-09-01 重算）

**文件**: `analysis/output/structure_source_stats.json`

**脚本**: `analysis/scripts/compute_structure_source_stats.py`

**方法**: 依据最终划分文件中的 `struct_file`（或 `struct_file1`/`struct_file2`）命名规则判定结构来源：
- `AF-*` 或 `af__*` 视为 **AlphaFold 预测结构**；
- 其余 PDB/CATH/SIFTS/PDBbind 来源视为 **实验结构**；
- `enzyme_optimal_ph` 使用 `metadata.csv` 中的 `structure_source`（`AF`/`PDB`）字段；
- `ppi_prediction` 为 pair-level，若两条链中**任一链**为 AlphaFold 则该 pair 计入 Predicted。

**核心结果**:

| Task | Train (N / Exp% / Pred%) | Val (N / Exp% / Pred%) | Test (N / Exp% / Pred%) | Source |
|:-----|:-------------------------|:-----------------------|:------------------------|:-------|
| kcat | 15,479 / 0.0% / 100.0% | 2,150 / 0.0% / 100.0% | 5,002 / 0.0% / 100.0% | AlphaFold v6 (UniProt) |
| pH | 7,794 / 0.0% / 100.0% | 866 / 0.0% / 100.0% | 2,955 / **93.0%** / **7.0%** | AlphaFold v6 + RCSB PDB (metadata) |
| Fold (CATH) | 43,473 / 100.0% / 0.0% | 6,132 / 100.0% / 0.0% | 12,911 / 100.0% / 0.0% | RCSB PDB (CATH S95) |
| EC | 15,389 / 100.0% / 0.0% | 1,687 / 100.0% / 0.0% | 5,047 / 100.0% / 0.0% | RCSB PDB (SIFTS) |
| GO-BP | 15,232 / 100.0% / 0.0% | 1,742 / 100.0% / 0.0% | 5,142 / 100.0% / 0.0% | RCSB PDB (SIFTS) |
| GO-CC | 13,878 / 100.0% / 0.0% | 1,730 / 100.0% / 0.0% | 5,613 / 100.0% / 0.0% | RCSB PDB (SIFTS) |
| GO-MF | 18,179 / 100.0% / 0.0% | 2,084 / 100.0% / 0.0% | 5,609 / 100.0% / 0.0% | RCSB PDB (SIFTS) |
| LBA | 8,366 / 100.0% / 0.0% | 930 / 100.0% / 0.0% | 2,367 / 100.0% / 0.0% | PDBbind v2020/v2024 |
| LBS | 71,093 / 100.0% / 0.0% | 7,938 / 100.0% / 0.0% | 17,225 / 100.0% / 0.0% | RCSB PDB local mirror |
| PPI | 58,968 / 10.7% / 89.3% | 8,394 / 4.8% / 95.2% | 16,264 / **3.5%** / **96.5%** | PINDER (pair-level; 任一链为 AF 即计 Predicted) |
| PPIS | 6,043 / 100.0% / 0.0% | 672 / 100.0% / 0.0% | 1,980 / 100.0% / 0.0% | DIPS-PLUS (RCSB PDB) |
| SSP | 29,622 / 100.0% / 0.0% | 3,317 / 100.0% / 0.0% | 8,372 / 100.0% / 0.0% | RCSB PDB (X-ray/NMR/EM filtered) |

**按样本汇总**（含 pair-level 任务的样本计数）:

| 类别 | 样本数 | 占比 |
|:-----|------:|-----:|
| 实验结构 (Experimental) | 321,792 | 74.9% |
| 预测结构 (Predicted) | 107,853 | 25.1% |
| **合计** | **429,645** | **100.0%** |

> 注：PPI 按 pair 计数；若把 PPI 的两条链拆开为链级统计，则实验链数会进一步上升，但 pair-level 评估以 pair 为单位，因此按 pair 汇报。PPI 负样本中的 XTAL 晶体 pair 全部为实验结构，因此重划分后实验结构占比（尤其 train）较旧版有所上升。

**关键发现**:
- **100% 实验结构**：Fold、EC、GO-BP/GO-CC/GO-MF、LBA、LBS、PPIS、SSP。
- **100% AlphaFold 预测**：kcat（train/val/test 全部为 AF v6）。
- **混合来源**：optimal_ph 的 train/val 全为 AF，test 中 93.0% 为 PDB 实验结构（来自 Enzyme Optimum pH 数据库），7.0% 为 AF；PPI 的 test 中 96.5% 的 pair 至少含一条 AF 链，实验 pair 仅占 3.5%。
- 整体上，**74.9%** 的样本使用实验结构，**25.1%** 使用 AlphaFold 预测结构；预测结构主要集中在 PPI 和 kcat 任务。

**论文用途**: Table 1 或 Supplementary Table（Dataset / Structure Source Overview）

---

#### 7. 数据质量报告

**状态**: 🔄 完成

> **快照说明**：§7.1–§7.4 的数字为 2026-07 初快照（结构/序列级重核对未随 08 月数据集更新重跑，fold/PPI 等已重划分任务的绝对计数会有出入，但结论性质不变）；§7.5–§7.6 为 2026-08 新增内容。

**7.1 aa_seq 与结构提取序列一致性**

对全部任务的 train/val/test 划分，逐样本解析结构文件（`pdbs/` 或 `pdb_file` 列引用）的 ATOM 记录、逐链提取氨基酸序列（含常见修饰残基映射），与 `aa_seq` 比对。PPI 为 pair-level，按 A/B 两条链分别核对（下表 TOTAL 为链数）。

- **脚本**: ~~`analysis/scripts/verify_aa_seq_vs_structure.py`~~（2026-09 analysis/scripts 清理时移除，结果快照保留）
- **结果**: `analysis/output/aa_seq_structure_consistency.json`
- **规模**: 518,470 条待核对序列，318,673 个去重结构文件

分类定义：`EXACT`=逐字节一致；`结构多出`=aa_seq 是结构序列子串（结构解析出更多残基）；`aa_seq更长`=结构序列是 aa_seq 子串（结构有未解析残基/缺失密度）；`MISMATCH`=既不相等也非子串（多为内部单残基 indel 或前导 `X` 占位导致的错位，抽样比对显示仍为同一蛋白，中位一致度 ~0.98）。

| 任务 | 核对数 | EXACT | 完全一致率 | 结构多出 | aa_seq更长 | MISMATCH | 空结构 |
|------|------:|------:|------:|------:|------:|------:|------:|
| kcat | 22,631 | 22,631 | **100.0%** | 0 | 0 | 0 | 0 |
| optimal_ph | 11,615 | 11,615 | **100.0%** | 0 | 0 | 0 | 0 |
| LBS | 96,256 | 96,045 | **99.8%** | 0 | 0 | 211 | 0 |
| ppis | 8,695 | 8,682 | **99.9%** | 0 | 0 | 13 | 0 |
| ss | 40,066 | 34,384 | **85.8%** | 280 | 501 | 4,901 | 0 |
| func | 87,150 | 73,697 | **84.6%** | 394 | 10,029 | 3,030 | 0 |
| PPI（按链） | 177,344 | 144,623 | **81.5%** | 2,385 | 9,953 | 20,383 | 0 |
| LBA | 11,663 | 5,264 | **45.1%** | 0 | 6,395 | 4 | 0 |
| fold | 62,896 | 24,764 | **39.4%** | 10 | 16,618 | 21,502 | 2 |
| **合计** | **518,470** | **421,704** | **81.3%** | 3,070 | 43,496 | 50,038 | 2 |

**结论与解读**：
- **完全一致**：`kcat`、`optimal_ph`（AlphaFold 全长模型，序列即结构）；`LBS`（99.8%）、`ppis`（99.9%）实质一致，非一致仅为个别点差异/单残基 indel。
- **系统性差异（非错配，同一蛋白）**：
  - `fold`：`aa_seq` 为**整链序列**（常含前导 `X` 占位），结构为 **CATH 域**——长度/起点不同；抽样比对中位一致度 0.98。
  - `LBA`：绝大多数为 `aa_seq更长`——沉积 PDB 有大量未解析残基（晶体缺失密度），序列高度一致但结构不完整。
  - `func`：非一致主要是 `aa_seq更长`（SIFTS 子序列 vs 结构缺失残基）。
  - `PPI`：81.5% 完全一致，其余多为结构缺失残基或 AlphaFold 全长 vs 截断链，同一蛋白。
  - `ss`：非一致的 ~14% 抽样中位一致度 0.989（100% ≥0.95），为单残基 indel；**对逐残基任务意味着这些样本"标签↔结构"存在 1 残基级错位，需评估**。
- **确凿的坏样本**：`fold` 有 **2 个空结构文件**（`4bpe000`、`4v19000`，0 残基），应单独修复或剔除。
- **总体**：未发现"序列-结构张冠李戴"式错配；全部非一致均为同一蛋白的细节差异（缺失残基 / indel / 前导 X / 域 vs 整链）。

**7.2 残基级任务三方一致性（LBS / ppis / ss）**

对三个 per-residue 任务，核对 `aa_seq` ↔ 结构提取序列 ↔ 逐残基标注 三者一致性。

- **脚本**: ~~`analysis/scripts/residue_level_consistency.py`~~（2026-09 analysis/scripts 清理时移除，结果快照保留）
- **结果**: `analysis/output/residue_level_consistency.json`

**关键结论：`len(labels) == len(aa_seq)` 在三任务全部 145,017 个样本上 100% 成立**——标注与 `aa_seq` 严格逐残基对齐，无一例外。因此"标注 ↔ 结构"的一致性等价于"`aa_seq` ↔ 结构"的一致性。

| 任务 | 标注长度==aa_seq | aa_seq==结构 | 主要不一致原因 | 非一致中位一致度 |
|------|:---:|------:|------|:---:|
| LBS | 100% (96,256) | 100% | 211 例已由结构序列替换 `aa_seq`（单残基末端替换），标注与结构完全对齐 | 1.000 |
| ppis | 100% (8,695) | 99.9% | 154 例低质量/占位符结构已删除，61 例已按全局对齐裁剪为代表链；剩余 13 例结构与 aa_seq 等长（6 例纯替换 + 7 例内部 indel 抵消） | 0.989 |
| ss | 100% (40,066) | 85.8% | 4,737 例结构内部缺失残基 + 501 端部缺失（晶体未解析密度），已用 `-1` 掩码 | 0.990 |

**标注编码与掩码机制**：
- **LBS**：标注 `{0,1}`（结合位点二分类），无掩码；`aa_seq==结构` 达 100%，标注逐残基对齐结构。
- **ppis**：标注 `{0,1,-1}`，`-1`（140 个残基）= 该残基在 PDB 中存在但原始界面标注缺失（文档已定义）。结构-序列-标注长度已统一。
- **ss**：标注 `{0..8, -1}`，其中 **`-1`（30,051 个残基）= aa_seq 中未被结构解析的残基（缺失密度）的掩码**。验证：`#(标注≥0) == 结构残基数` 占 97.8%、`#(-1) == (aa_len − 结构长度)` 占 98.9%——即 ss 的 `aa_seq` 是完整序列、结构仅含已解析残基、缺失部分在标注中以 `-1` 掩码。**故 ss 的 aa_seq 与结构长度差异并非错误，而是被 `-1` 完整记账**。

**逐任务原因诊断**：
- **LBS（211 例）**：已将 `aa_seq` 替换为结构提取序列，消除末端单残基替换；内部残基与标注完全对齐，影响可忽略。
- **ppis（13 例）**：已删除 154 例无法恢复的低质量/占位符结构，并将 61 例不一致结构按全局对齐裁剪为代表链（跳过结构插入、保留与 `aa_seq` 一一对应的残基）。剩余 13 例结构长度与 `aa_seq` / 标注完全一致（6 例纯点替换 + 7 例内部插入-缺失长度抵消），标注逐残基对齐。
- **ss（5,238 例）**：绝大多数为结构未解析残基（内部/端部缺失密度），`aa_seq` 完整而结构不完整，标注已用 `-1` 掩码处理，属设计性一致；非一致中位一致度 0.99。

**ss 结构输入的专用视图**：`splits/` 的全序列 + `-1` 掩码适合序列输入模型；但对结构输入模型，实测约 2.6%（207/8,079 test）样本的 `-1` 掩码未如实记录所有内部缺口，"丢弃 `-1`" 会静默错帧。为此新增**结构原生视图** `datasets/ss_prediction/splits/splits_struct_native/`（脚本已归档：`datasets/ss_prediction/output/scripts_archive/build_ss_struct_native.py`）：`aa_seq` = `pdbs/` 的 CA 序列（模型实际输入），`labels` = 以 DSSP `ss8` 为源经全局比对映射到每个结构残基，`len(labels)==len(aa_seq)==M(pdb)` 在全部 40,066 条上成立、结构-标签严格 1:1（DSSP 精确 98.8% / 比对 1.2% / 回退 6，identity 均值 0.99992，无隔离样本）。该视图后续已弃用（见下注）。

> **2026-08-28 弃用**：主 `splits/` 新增 `struct_label` 列（与 struct_native 的 label 同源同值、残基枚举为全 ATOM 口径，且 OOD 列保持同步）后，该视图已冗余，归档至 `datasets/ss_prediction/output/splits_struct_native_deprecated/`。结构模型请使用主 `splits/` + `struct_label`。

**7.3 序列合法性检查（非法氨基酸字符）**

对全部 12 个任务视图共 **552,093 条序列 / 约 1.68 亿残基** 逐字符扫描（脚本 ~~`analysis/scripts/sequence_legality.py`~~，2026-09 清理时移除；结果 `analysis/output/sequence_legality.json`）。字符分四类：标准 20 种 `ACDEFGHIKLMNPQRSTVWY`、IUPAC 歧义/稀有码 `X B Z J U O`、链分隔符 `|`、其余一律判为非法。

**核心结论：0 条非法字符、0 条空/NaN 序列**——无小写、数字、间隔符（`-./*`）、空白或乱码。所有非标准字符均为生物学合法编码或格式约定，无需清洗。

| 任务 | #序列 | 空/NaN | 含歧义码 | 含分隔符 | 非法 | 非标准字符明细 |
|------|------:|:---:|------:|------:|:---:|------|
| kcat | 22,631 | 0 | 0 | 0 | 0 | — |
| optimal_ph | 11,615 | 0 | 0 | 0 | 0 | — |
| fold(CATH) | 62,517 | 0 | 10,922 | 0 | 0 | `X`×47,977（前导占位/未解析残基） |
| func_EC | 22,123 | 0 | 15 | 0 | 0 | `X`×611 |
| func_GO-BP | 22,116 | 0 | 13 | 0 | 0 | `X`×478 |
| func_GO-CC | 21,221 | 0 | 21 | 0 | 0 | `X`×618 |
| func_GO-MF | 25,872 | 0 | 22 | 0 | 0 | `X`×735 |
| LBA | 11,663 | 0 | 0 | 6,393 | 0 | `\|`×12,657（多链复合物分隔符） |
| LBS | 96,256 | 0 | 0 | 0 | 0 | — |
| PPI | 167,252 | 0 | 272 | 0 | 0 | `X`×409, `B`×101, `U`×58（硒代半胱氨酸）, `Z`×41, `O`×2（吡咯赖氨酸） |
| PPIS | 8,695 | 0 | 0 | 0 | 0 | — |
| SSP | 40,066 | 0 | 0 | 0 | 0 | — |
| SSP(struct-native) | 40,066 | 0 | 16 | 0 | 0 | `X`×161（CA 提取的非标准残基） |

> 注：`#序列` 中 PPI 按两条链分别计数（`aa_seq1` + `aa_seq2`），LBA 的 `|` 用于分隔同一复合物的多条链（如同源四聚体=4×468aa）。下游 dataloader 需按 `|` 拆链、并对 `X/B/Z/J/U/O` 做 tokenizer 兜底（多数蛋白语言模型已内置 `X`/`<unk>`）。

**7.4 格式一致性审查（列名 schema + 标签格式）**

在完成全库列名统一（`label`/`aa_seq(1,2)`/`struct_file(1,2)`、`OOD_` 前缀、`OOD_NewEC_L*`/`OOD_LongTail(_EC_L*)` 等）后，对全部 13 个任务视图、39 个 split 文件做结构化审查（脚本 `analysis/scripts/format_consistency.py`，结果 `analysis/output/format_consistency.json`）。

审查维度：核心列齐全、无弃用列残留（`labels`/`file_type`/`pdb_file*`/`is_time_cutoff`/`OOD-*`/`OOD_NewFunction*`/`OOD_LongTailFunction*`/`OOD_CATH_*`/`protein_*_seq`）、`label` 非空、train/val 无 OOD 列而 test 含 `Default`+OOD、per-residue 任务 `len(label)==len(aa_seq)`、train/val/test 间 `unique_id` 无泄漏、train/val 列集一致。

**核心结论：194/194 项结构化检查全部通过。**

| 任务 | Train/Val/Test | label 类型 | 核心列齐 | 弃用列 | 跨集泄漏 | per-res 长度 |
|------|---------------:|:---:|:---:|:---:|:---:|:---:|
| kcat | 15479/2150/5002 | float（回归） | ✅ | 0 | 无 | — |
| optimal_ph | 7797/866/2952 | float（回归） | ✅ | 0 | 无 | — |
| fold(CATH) | 34963/4994/22939 | str（C.A 码，需 dtype=str） | ✅ | 0 | 无 | — |
| func_EC | 15389/1687/4243 | str | ✅ | 0 | 无 | — |
| func_GO-BP | 15232/1742/4121 | str | ✅ | 0 | 无 | — |
| func_GO-CC | 13878/1730/4256 | str | ✅ | 0 | 无 | — |
| func_GO-MF | 18179/2084/4609 | str | ✅ | 0 | 无 | — |
| LBA | 8366/930/2367 | float（回归） | ✅ | 0 | 无 | — |
| LBS | 71093/7938/17225 | str | ✅ | 0 | 无 | ✅ |
| PPI | 65020/7195/16457 | int（0/1） | ✅ | 0 | 无 | — |
| PPIS | 6043/672/1980 | str | ✅ | 0 | 无 | ✅ |
| SSP | 28788/3199/8079 | str | ✅ | 0 | 无 | ✅ |
| SSP(struct-native) | 28788/3199/8079 | str | ✅ | 0 | 无 | ✅ |

其余 12 个任务的 `label` 类型均符合语义：回归任务（kcat/optimal_ph/LBA）为 float、PPI 交互为 int、多标签（EC/GO）与 per-residue（LBS/PPIS/SSP）为 str。

**论文用途**: Supplementary Material（数据质量）

**7.5 struct_label 列（残基级任务的结构对齐标注，2026-08 新增）**

三个残基级任务（LBS / PPIS / SSP）的全部 split 新增末列 `struct_label`：按 pdb 文件全 ATOM 残基（含无 CA 残基，按 `(chain, resseq, icode)` 首次出现顺序）逐残基对齐的标注，`-1` 为 ignore index。结构输入模型可直接以 `pdbs/` 结构 + `struct_label` 训练评估，无需再做序列-结构比对。

- **SSP**：直接移植原 `splits_struct_native` 的 DSSP 直达标注（权威源）；主 `splits/` 全序列 label 的缺口放置存在约 8% train 行局部错帧瑕疵，不用于结构模型。原 `splits_struct_native/` 视图已退役归档至 `datasets/ss_prediction/output/splits_struct_native_deprecated/`。
- **LBS / PPIS**：由结构接触提取流程原生对齐生成，等长行与全序列 label 100% 一致、`-1` 归零。
- **工具**：`.skills/struct-label-alignment/`（`add_struct_label.py`，支持 `--dry-run`/`--overwrite`）。

**7.6 结构库（pdbs/）修复记录（2026-08）**

对全部 10 个任务的结构库做全量扫描与修复，修复后全库重扫 0 残留问题、所有文件 Bio.PDB 可解析为单模型、无"结构长于序列"样本：

| 任务 | 问题 | 规模 | 修复 |
|------|------|------|------|
| PPIS | pdb 列偏移（resname 后多一空格致 chain/resseq 右移，全库 Bio.PDB 不可读）+ 5 个链拼接截断 | 8,833 文件 | `ppis_prediction/output/scripts_archive/18_fix_pdb_column_layout.py`（坐标逐原子验证不变；已归档），提取脚本（现 `scripts/08_extract_representative_structures.py`）格式串已修补 |
| LBS | NMR 多模型无 ENDMDL 分隔串联（同残基重复出现） | 549 文件 | `output/scripts_archive/26_fix_multimodel_pdbs.py` 重提取 model 1（已归档）；`scripts/07_extract_protein_chains.py` 打补丁（ENDMDL 即停 + cif 首模型 + 带引号原子名去引号） |
| SSP | 多模型串联（23 个带 MODEL 记录 + 3 个无记录） | 26 文件 | `output/scripts_archive/fix_multimodel_pdbs.py`（已归档）；`scripts/05_extract_representative_chains.py` 改 `structure[0]` |
| Fold / func / PDA | 同款 NMR 串联 | 4,076 / 1,940 / 203 文件 | 各任务 `output/scripts_archive/fix_multimodel_pdbs.py`（已归档），提取脚本同步补丁 |
| PPI / LBA / kcat / pH | 扫描确认无问题 | — | kcat 7,057 个文件含 MODEL 记录但为单模型，无需处理 |

各任务修复前备份保存在对应 `output/backup_*/` 目录。

---

#### 8. protein-design-archive 额外测试集（2026-08 新增）

**状态**: ✅ 已完成

**文件**: `datasets/protein-design-archive/splits/cath_design.csv`、`ec_design.csv`

**文档**: `datasets/protein-design-archive/DATA_PROCESS.md`

**内容**: 从 Protein Design Archive（1,879 个条目）构建的人工设计蛋白额外测试集，用于评估模型在"自然界不存在"的蛋白上的泛化。使用全部条目（不限纯人工设计），仅去除完全相同序列（保留微小序列差异），不做 95% 同一性聚类；pdbs/ 保存去重后的完整单链结构。

| 测试集 | 样本数 | 对应任务 | 说明 |
|--------|------:|----------|------|
| `cath_design.csv` | 122 | fold_classification | 保留 CATH 可标注且落在现有 label 空间的设计蛋白，剔除与 fold train/val 重叠（id 或 ≥95% 序列同一性）的样本；OOD 列与主测试集同 schema |
| `ec_design.csv` | 31 | func_prediction (EC) | EC 标注不做聚类内转移，保证每条样本标注准确可信 |

- 全部行 `Default=True`（设计蛋白不参与主数据集的正常划分，不适用 Default=False 的先验划分语义）；
- OOD 标注（seq/TM/Orphan/NewEC/LongTailEC 等）以主数据集对应任务的 train/val 为参考集计算；
- 后续可复用 pdbs/ 与流程构建其他任务的设计蛋白额外测试集。

**论文用途**: Supplementary / Figure（De novo 设计蛋白泛化评估）

---

### 二、待完成的无需模型分析

#### 6. 与现有Benchmark的对比表

**状态**: ⏳ 待完成

**内容**: 与CASP、CAFA、PDBbind、CATH、EnzyBase等现有基准的定量对比

**维度**: 任务数、OOD场景、样本量、结构可用性、评估维度

**论文用途**: Table 2（Related Work / Motivation）

**工作量**: 需手动整理文献数据，无需代码

---

## 第二部分：需要训练模型的实验

### 一、基线模型评估（核心实验）

#### 1. ID（In-Distribution）性能基准

**目标**: 在train/val上训练，在test的ID子集上评估

**方法**: 筛选test集中所有OOD标记均为False的样本作为ID test

**模型**: ESM-2、ProstT5、AlphaFold2-derived等预训练模型+fine-tune

**评估指标**: 任务相关（ACC、F1、MCC、MAE、Spearman等）

**论文用途**: Table 4（ID性能基准）

**计算需求**: 高（需GPU训练）

---

#### 2. OOD性能衰减分析

**目标**: 对比ID vs各OOD场景的性能下降

**方法**: 在同一模型上，分别评估test集的ID子集和每个OOD子集

**预期结果**: OOD场景性能显著低于ID，衰减幅度与OOD难度相关

**论文用途**: Figure 4（OOD性能衰减图）, Table 5（各OOD场景性能对比）

**计算需求**: 高（基于实验1的模型，额外推理）

---

#### 3. OOD阈值敏感性分析

**目标**: 验证OOD阈值选择的合理性

**方法**: 变化SeqRedundancy（40%, 30%, 20%）和TM-score（0.5, 0.3, 0.2）阈值，观察性能衰减曲线

**论文用途**: Figure 6（阈值敏感性曲线）

**计算需求**: 高（需多次重新划分OOD子集并评估）

**注意**: 此分析需要模型在不同OOD子集上的性能输出

---

#### 4. OOD场景组合效应（难度分析）

**目标**: 分析多个OOD标记同时出现的样本的模型性能

**方法**: 基于模型输出，对比单一OOD vs组合OOD的性能衰减

**论文用途**: Figure 7（组合OOD效应分析）

**计算需求**: 高（基于实验2的数据，额外分析）

**注意**: 此分析需要模型在组合OOD子集上的性能数据，属于模型依赖分析

---

#### 5. 模型规模与OOD性能关系

**目标**: 验证更大模型是否对OOD更鲁棒

**方法**: 对比ESM-2 (8M/35M/150M/650M/3B)在不同OOD场景的表现

**论文用途**: Figure 8（模型规模 vs OOD性能）

**计算需求**: 极高（需训练/推理多个模型规模）

---

### 二、OOD检测实验（基于模型输出）

#### 6. 无监督OOD检测

**目标**: 不依赖OOD标签，检测test样本是否为OOD

**方法**:
- 基于预训练模型的embedding密度（如MSP、Mahalanobis距离）
- 基于重建误差（如VAE、AE）
- 基于梯度/不确定性（如MC-Dropout）

**评估**: AUROC、AUPR、FPR@95%TPR

**论文用途**: Figure 5（OOD检测ROC曲线）, Table 6（检测性能对比）

**计算需求**: 中高（需训练检测器或提取embedding）

**注意**: 此实验需要模型输出的embedding或概率分布，属于模型依赖分析

---

### 三、OOD缓解策略实验（基于模型训练）

#### 7. 数据增强与重采样

**目标**: 缓解OOD性能衰减

**方法**:
- 对长尾标签进行过采样
- 对ExtremeLength/IDR样本进行数据增强
- 基于序列相似度的训练集扩展

**评估**: ID性能保持 + OOD性能提升

**论文用途**: Table 7（缓解策略效果对比）

**计算需求**: 高

---

#### 8. 域自适应/迁移学习

**目标**: 提升模型对OOD域的泛化能力

**方法**:
- DANN（Domain-Adversarial Neural Network）
- CORAL（特征对齐）
- Mixup/CutMix等数据混合策略

**评估**: 同实验5

**论文用途**: 与实验5对比

**计算需求**: 高

---

#### 9. 集成方法与不确定性量化

**目标**: 提升模型对OOD样本的识别能力

**方法**:
- Deep Ensemble（多模型集成）
- SWA（Stochastic Weight Averaging）
- Evidential Deep Learning

**评估**: ID性能 + OOD检测AUROC + 校准性（ECE）

**论文用途**: Table 8（集成方法对比）

**计算需求**: 极高（需训练多个模型）

---

## 论文图表规划

### 无需模型的图表（可直接生成）

| 图表 | 内容 | 来源 |
|------|------|------|
| Table 1 | 数据集规模统计 | 结构序列统计 |
| Figure 1 | OOD覆盖热力图 | OOD覆盖度矩阵 |
| Figure 2 | OOD重叠矩阵 | OOD重叠分析 |
| Figure 3 | 标签长尾/Zipf曲线 | 标签分布分析 |
| Table 2 | Benchmark对比 | 待完成（文献整理） |

### 需要模型的图表（实验后生成）

| 图表 | 内容 | 来源 |
|------|------|------|
| Table 3 | ID性能基准 | 实验1 |
| Figure 4 | OOD性能衰减 | 实验2 |
| Figure 5 | OOD检测ROC | 实验3 |
| Table 4 | OOD检测对比 | 实验3-4 |
| Table 5 | 缓解策略对比 | 实验5-7 |
| Figure 6 | 阈值敏感性 | 实验8 |
| Figure 7 | 组合OOD效应 | 实验9 |
| Figure 8 | 模型规模vsOOD | 实验10 |

---

### 最小可行论文（MVP）

- **必须完成**: 阶段一全部 + 实验1-4
- **可选增强**: 实验5（缓解策略）+ 实验6（敏感性）+ 实验7-9（组合效应）
- **理想完整**: 全部实验

---

## 文件位置

分析脚本与中间结果组织如下（相对仓库根 `POOR/`）：

```
analysis/
├── scripts/                                  # 分析脚本（仅保留数据统计/论文作图工具）
│   ├── dataset_scale_stats.py                # §1 数据集规模统计（2026-08 新建）
│   ├── ood_coverage_matrix.py                # §2 OOD 覆盖度矩阵 + §3 重叠度（2026-08 更新场景）
│   ├── label_distribution.py                 # §4 标签长尾分布（2026-08 新建）
│   ├── compute_structure_source_stats.py     # §5 结构来源占比（2026-08 修正 ROOT）
│   ├── plot_figure2_ood_orthogonality.py     # Figure 2 绘图
│   └── plot_ood_jaccard_heatmaps.py          # OOD Jaccard 热力图绘图
└── output/                                   # 分析中间结果
    ├── dataset_scale_stats.json              # §1（2026-08-28）
    ├── ood_coverage_matrix.csv / .md         # §2（2026-08-28）
    ├── ood_overlap_matrix.csv                # §3（2026-08-28）
    ├── label_distribution.json               # §4（2026-08-28）
    ├── structure_source_stats.json           # §5（2026-08-28）
    ├── aa_seq_structure_consistency.json     # §7.1（2026-07 快照）
    ├── residue_level_consistency.json        # §7.2（2026-07 快照）
    ├── sequence_legality.json                # §7.3（2026-08-28 重算）
    ├── format_consistency.json               # §7.4（2026-07 快照）
    └── ood_jaccard_heatmaps/                 # Figure 2 图片
```

> 注（2026-08-28）：任务特定的一次性脚本已下沉至各任务 `scripts/` 后，各任务 scripts 已于同日完成全面整理：全部流程脚本按执行顺序重编号（`NN_*.py`），一次性补丁脚本（`fix_multimodel_pdbs.py`、`18_fix_pdb_column_layout.py`、`fix_lbs_aaseq.py`、`build_ss_struct_native.py`、`ss_struct_alignment_probe.py`、`patch_struct_label_from_native.py` 等）已归档至各任务 `output/scripts_archive/`（修复均已落入数据、根因已修入正常流程脚本），`crop_ppis_structures.py` 现为 `ppis_prediction/scripts/12_crop_ppis_structures.py`；各任务处理流程以对应 `DATA_PROCESS.md`（四阶段结构）为准。

可复用的数据处理流程已沉淀为项目 skill（`.skills/`）：

```
.skills/
├── cath-ted-domain-annotation/    # CATH+TED 域标注（topology/superfamily），供 FoldHoldout/SuperfamilyHoldout OOD 判定
├── ec-function-ood-annotation/    # EC 标注 + NewEC / LongTailEC OOD 统一口径（all-not-in，L3/L4，阈值 5/10）
├── homology-ood-annotation/       # 序列同源性 / Orphan / TM-score OOD 标注
├── ood-annotation-toolkit/        # OOD 标注通用工具（Default/InD/NaN 语义）
├── pood-benchmark-process/        # 数据集构建总流程
└── struct-label-alignment/        # 残基级任务 struct_label 生成与校验
```

> 注（2026-09-01）：`analysis/scripts/` 只保留数据统计/论文作图工具（上表 6 个）；QA 审计脚本（`comprehensive_audit.py`、`struct_audit.py`、`audit_foldseek_easysearch.py`、`verify_aa_seq_vs_structure.py`、`residue_level_consistency.py`、`sequence_legality.py`）与一次性调研脚本（`combinatorial_ec_feasibility.py`）已删除，其产出（`output/comprehensive_audit_20260830/` 等）保留备查。`ppis_prediction/scripts/12_crop_ppis_structures.py` 对 `verify_aa_seq_vs_structure.py` 的 import 已改为内联 Biopython 映射。

> **注意**：仓库根 `.gitignore` 忽略 `pdbs/`、`splits/`、`output/`、`check/`，因此 `analysis/output/` 及各任务的 `output/` **不纳入版本控制**；`analysis/scripts/` 与 `analysis/benchmark_analysis.md` 已入库。
