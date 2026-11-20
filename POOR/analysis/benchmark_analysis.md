# POOR: Data Analysis and Experiment Planning

## Part I: Analyses without model training (completed and planned)

### I. Completed analyses

#### 1. Dataset scale statistics

**Status:** ✅ Complete (recomputed 2026-09-01)

**File**: `analysis/output/dataset_scale_stats.json`

**Script**: `analysis/scripts/dataset_scale_stats.py`

**Results**:

### Classification Tasks

| Task | Train | Val | Test | Total | Labels | Seq Length (median) | Extreme in Train/Val? |
|------|-------|-----|------|-------|--------|---------------------|----------------------|
| fold_classification | 43,473 | 6,132 | 12,911 | 62,516 | 1,163 classes (C.A.T) | 9-960 (med=133) | ⚠️ YES (2,940, 5.9%, all <60) |
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
| ss_prediction | 29,622 | 3,317 | 8,372 | 41,311 | 8-state + 3-state DSSP | 31-1836 (med=222) | ⚠️ YES (952, 2.9%, all <60) |

**Key findings:**
- Total: **429,645** samples (train 303,516 /val 37,642 /test 88,487).
- Fold has the shortest sequences (median ≈133 aa); regression and binding-site tasks tend to be longer. After removing >2000-aa chains, the single-chain upper bound is 2000 aa.
- Eight task views have no extreme lengths (<60 or >1000) in train/val. Four retain them: fold_classification, ppi_prediction, ligand_binding_affinity, ss_prediction. **LBA has the highest fraction**, 13.0%, almost entirely >1000-aa chains in multichain complexes. All 2,940 fold and 952 SSP extremes are short chains retained/partitioned by design; PPI repartitioning reduces extremes to 2.7%.
- Function-task minimum test length 3 aa comes from ExtremeShort samples recovered in 2026-08, added only to test. SSP <60-aa chains (8,941, minimum 1 aa) were restored on 2026-09-01 and split by the original time cutoffs on 09-02; see the notes below.
- The only empty PPI chain, an O75746 homodimer, was completed from UniProt on 2026-08-28 (ppi_prediction/DATA_PROCESS.md, Step 5). Minimum length 3 is a real short chain. LBA lengths are counted after splitting on `|`; minimum 1 reflects deposited short peptides/placeholders, as documented in its audit.
- Every task references structure files (PDB/CIF).
- **Extreme-length corrections, 2026-09-01/02:** (1) Restore SSP <60-aa chains: 22,231 →95% deduplication and cross-deduplication against existing splits →8,941; partition on 09-02 using original cutoffs train≤2019-06-06, val≤2021-02-09 (+5,672 /+892 /2,377). Short test chains have Default=True; recompute seq/TM/Orphan/IDR/LongTail/NewEC/LongTailEC/fold-superfamily holdouts against new references. Remove low-quality 7UYL_K (59% UNK), whose sequence extraction skipped UNK retained in the structure and misaligned labels. Then retain only 31–59-aa short chains as requested, removing 7,695 ≤30-aa rows (train 4,838 /val 774 /test 2,083) and recomputing OOD again. (2) Move three short pH training samples P05959/B1KRG2/P0C7A9 to test, ES 24→27. (3) Remove the only >1000-aa fold training sample 3w3uA01 (1,047 aa), missed because its CATH length_expected=999.

---

#### 2. OOD scenario coverage matrix

**Status**: ✅ Complete (recomputed 2026-09-01)

**File**: `analysis/output/ood_coverage_matrix.csv` / `.md`

**Script**: `analysis/scripts/ood_coverage_matrix.py`

**Results:** counts and percentages per test OOD scenario; percentages use the full N_test; `-` means not applicable.

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

> Notes: (1) Unified column names replace OOD_NewFunction with OOD_NewEC_L4/L3 and OOD_LongTailFunction with OOD_LongTail_EC_L*. CATH OOD_NewEC is placed in the New EC L4 row. (2) Sequence/structure/Orphan flags are normally evaluated on Default=True only; prior-held-out extreme rows receive NaN, counted as False, while percentages still use the full test denominator. (3) InD means no OOD flag (explicit InD column). (4) NewEC/LongTailEC use all-not-in against train or train+val, L3/L4, thresholds ≤5/≤10 including zero. (5) PPI is pair-level (16,264), PPIS residue-interface level (1,980); IDR differences reflect PPI any-chain semantics. (6) Non-CATH fold/superfamily holdouts use CATH+TED domain annotations; see .skills/cath-ted-domain-annotation/.

---

#### 3. OOD orthogonality / scenario-flag overlap

**Status**: ✅ Complete (recomputed 2026-09-01)

**File**: `analysis/output/ood_overlap_matrix.csv`

**Script**: `analysis/scripts/ood_coverage_matrix.py`(same script output)

**Method:** compute pairwise test-set Jaccard similarity within each task:
`J(A,B) = |A∩B| / |A∪B|`. J=0 means disjoint sets; J near 1 means high overlap.
Highly overlapping axes should not be treated as independent evaluation dimensions.

**Results**:

- Evaluate 301 OOD-axis pairs; mean Jaccard **0.071**.
- **121 pairs have J=0**, including ExtremeShort versus ExtremeLong, extreme length versus EC novelty, and CATH holdout versus combinatorial generalization.
- **13 pairs have J>0.5**, including 11 Orphan versus seq_Redundancy_30 pairs, expected because Orphan is a limiting subset of low-identity proteins:

| Task | Axis pair | Jaccard | Orphan and <30% | Orphan only | <30% only |
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

The remaining two high-overlap pairs occur in LBA, where complex-level multimer TM strengthens sequence–structure novelty coupling:

| Task | Axis pair | Jaccard | Intersection | First only | Second only |
|------|----------|---------|-----:|-------:|-------:|
| LBA  | seq_Redundancy_30 ∩ TM-score_0.5 | 0.524 | 152 | 89 | 49 |
| LBA  | OOD_Orphan ∩ TM-score_0.5 | 0.522 | 131 | 50 | 70 |

- The remaining 288 pairs all have Jaccard **≤0.448**; the maximum is kcat Orphan–seq30, also nested but below 0.5.
  The vast majority (75th percentile)
  are ≤0.055.

**Mean Jaccard by task:** lower values indicate less overlapping task-specific OOD axes.

| Task | Mean Jaccard | Max Jaccard | Notes |
|------|-------------:|-------------:|------|
| kcat | 0.042 | 0.448 | Sequence unification on 2026-08-30 reduces Orphan–seq30 to 0.448, below the high-overlap threshold. |
| LBS | 0.047 | 0.705 | Orphan–seq30 is largest (nested). TM easy-search increased mean overlap from 0.033 to 0.048; currently 0.047. |
| PPIS | 0.055 | 0.761 | Driven by Orphan–seq30. TM easy-search reduces the mean from 0.077 to 0.055; sequence easy-search unification on 2026-08-29 leaves mean/max unchanged at 0.055 /0.761. |
| GO-CC | 0.060 | 0.787 | Function prediction; sequence unification on 2026-08-30 makes nested Orphan–seq30 a high-overlap pair. |
| SSP | 0.070 | 0.777 | Driven by Orphan–seq30. Removing ≤30-aa chains removes the former ES/IDR/TM0.5/Orphan coupling, leaving the nested pair. |
| GO-MF | 0.060 | 0.802 | Function prediction; same pattern as above. |
| EC | 0.064 | 0.683 | Function prediction; same pattern as above. |
| pH | 0.067 | 0.655 | Driven by Orphan–seq30. |
| GO-BP | 0.069 | 0.788 | Function prediction; same pattern as above. |
| CATH | 0.105 | 0.841 | Driven by Orphan–seq30. |
| PPI | 0.106 | 0.826 | Driven by Orphan–seq30. |
| LBA | 0.128 | 0.751 | Driven by Orphan–seq30; multimer TM increases mean overlap from 0.116 to 0.128. |

**Finding:** after excluding the Orphan–low-identity nesting, OOD dimensions are largely independent in this overlap analysis, capturing **complementary shift types** for multi-axis evaluation.

**Paper use**: Figure 2(OOD Overlap Matrix) or Supplementary Table

---

#### 4. Label long-tail distribution

**Status:** ✅ Complete; corrected 2026-09-14 to per-label frequencies matching paper Table 1. Old label-combination counts remain as combo_* references.

**File**: `analysis/output/label_distribution.json`

**Script:** `analysis/scripts/label_distribution.py`: training per-label frequencies, splitting multilabel strings on `;`; SSP counts residue tokens. Also outputs combo_* statistics.

**Definition:** each EC number /GO term is a class, matching Table 1. Old counts such as EC 558 and GO-MF 3,436 counted whole label strings as combinations; retained in JSON combo_* fields.

**Results**(per-label definition):

| Task | Labels | Top-20% coverage | Bottom-50% coverage | Singleton fraction | Zipf exponent |
|------|--------|-------------|-------------|---------------|----------|
| EC | 260 | 75.4% | 4.9% | 7.7% | -1.596 |
| GO-MF | 447 | 77.9% | 7.0% | 0.0% | -1.277 |
| GO-CC | 264 | 82.6% | 5.7% | 0.0% | -1.308 |
| GO-BP | 884 | 68.8% | 11.9% | 0.0% | -0.995 |
| SSP | 9 | 53.6% | 7.7% | 0.0% | -1.763 |
| Fold | 1,163 | 86.6% | 3.8% | 14.0% | -1.445 |

**Class imbalance (binary tasks)**:

| Task | Train Pos% | Test Pos% |
|------|-----------|-----------|
| PPI | 50.0% | 50.0% |
| PPIS(residue-level) | 22.8% | 22.9% |

**Regression target distributions**:

| Task | Train Mean±Std | Test Mean±Std | Range |
|------|---------------|---------------|-------|
| LBA | 6.26 ± 1.97 | 6.11 ± 1.80 | [0.40, 15.22] |
| kcat | 0.88 ± 1.67 | 0.93 ± 1.62 | [-6.00, 6.00] |
| pH | 7.13 ± 1.28 | 7.39 ± 1.08 | [1.50, 12.50] |

**Key findings**:
- GO has no singleton labels but remains strongly skewed: top 20% covers 68.8%–82.6%, with shallow Zipf slopes -1.00 to -1.31.
- GO-BP has the largest vocabulary (884), shallowest Zipf slope (-0.995), and lowest head coverage (68.8%).
- EC has 260 L3 classes, 7.7% singletons, Zipf -1.596.
- C.A.T-level fold repartitioning increases classes from 41 to 1,163, with 14.0% singletons, the heaviest singleton tail among classification tasks.
- Reference combo_* counts: EC 558 combinations (31.7% singletons), GO-MF 3,436 (63.6%), GO-CC 2,506 (75.6%), GO-BP 3,191 (70.1%), showing highly unique multilabel combinations.
- SSP has nine valid DSSP states 0–8. The old count 16 resulted from unstripped `[`/`]` tokens and is corrected.
- PPIS residue-positive rate is 22.8%; the old 99.2% incorrectly counted samples containing any positive residue. PPI remains exactly 1:1 at sample level.

**Paper use**: Figure 3(Label Long-Tail Distribution / Zipf Curves)

---

#### 5. Experimental structure coverage

**Status**: ✅ Complete (recomputed 2026-09-01)

**File**: `analysis/output/structure_source_stats.json`

**Script**: `analysis/scripts/compute_structure_source_stats.py`

**Method:** infer structure source from final struct_file or struct_file1/struct_file2 naming:
- `AF-*` or `af__*`: **AlphaFold predictions**.
- Other PDB/CATH/SIFTS/PDBbind sources: **experimental structures**.
- enzyme_optimal_ph uses metadata.csv structure_source (`AF`/`PDB`).
- PPI is pair-level: count a pair as predicted if **either chain** is AlphaFold.

**Results**:

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
| PPI | 58,968 / 10.7% / 89.3% | 8,394 / 4.8% / 95.2% | 16,264 / **3.5%** / **96.5%** | PINDER (pair-level; predicted if either chain is AF) |
| PPIS | 6,043 / 100.0% / 0.0% | 672 / 100.0% / 0.0% | 1,980 / 100.0% / 0.0% | DIPS-PLUS (RCSB PDB) |
| SSP | 29,622 / 100.0% / 0.0% | 3,317 / 100.0% / 0.0% | 8,372 / 100.0% / 0.0% | RCSB PDB (X-ray/NMR/EM filtered) |

**Sample-level totals**, including pair-level task samples:

| Category | Samples | Fraction |
|:-----|------:|-----:|
| Experimental structures | 321,792 | 74.9% |
| Predicted structures | 107,853 | 25.1% |
| **Total** | **429,645** | **100.0%** |

> PPI counts pairs, not individual chains. Chain-level experimental coverage would be higher, but pair-level reporting matches evaluation. All XTAL crystal-packing negatives use experimental structures; repartitioning increases experimental fractions, especially in training.

**Key findings**:
- **100% experimental:** Fold, EC, GO-BP/CC/MF, LBA, LBS, PPIS, SSP.
- **100% AlphaFold:** kcat, all splits use AF v6.
- **Mixed:** pH train/val are all AF; test is 93.0% experimental PDB from the Enzyme Optimum pH database and 7.0% AF. PPI test pairs are 96.5% at least one AF chain, 3.5% fully experimental.
- Overall: **74.9%** experimental, **25.1%** AlphaFold, concentrated in PPI and kcat.

**Paper use**: Table 1 or Supplementary Table (Dataset / Structure Source Overview)

---

#### 7. Data quality report

**Status**: 🔄 Complete

> **Snapshot:** §7.1–7.4 use early-2026-07 counts; structure/sequence audits were not rerun for August dataset updates. Absolute counts for repartitioned fold/PPI may differ, while qualitative findings remain. §7.5–7.6 were added in 2026-08.

**7.1 Agreement between aa_seq and structure-derived sequences**

Parse ATOM records from each train/val/test structure (pdbs/ or pdb_file references), extract chains with modified-residue mapping, and compare to aa_seq. Check PPI chains A/B separately; TOTAL below counts chains.

- **Script:** ~~`analysis/scripts/verify_aa_seq_vs_structure.py`~~, removed during 2026-09 cleanup; result snapshots retained.
- **Results**: `analysis/output/aa_seq_structure_consistency.json`
- **Scale:** 518,470 sequences, 318,673 deduplicated structure files.

Categories: EXACT =byte-identical; structure-longer =aa_seq is a substring of extracted structure; aa_seq-longer =structure is a substring of aa_seq (unresolved residues/density); MISMATCH =neither equal nor substring, usually internal single-residue indels or leading X placeholders. Sampled alignments confirm the same protein, median identity ≈0.98.

| Task | Checked | EXACT | Exact fraction | Structure longer | aa_seq longer | MISMATCH | Empty structures |
|------|------:|------:|------:|------:|------:|------:|------:|
| kcat | 22,631 | 22,631 | **100.0%** | 0 | 0 | 0 | 0 |
| optimal_ph | 11,615 | 11,615 | **100.0%** | 0 | 0 | 0 | 0 |
| LBS | 96,256 | 96,045 | **99.8%** | 0 | 0 | 211 | 0 |
| ppis | 8,695 | 8,682 | **99.9%** | 0 | 0 | 13 | 0 |
| ss | 40,066 | 34,384 | **85.8%** | 280 | 501 | 4,901 | 0 |
| func | 87,150 | 73,697 | **84.6%** | 394 | 10,029 | 3,030 | 0 |
| PPI(per chain) | 177,344 | 144,623 | **81.5%** | 2,385 | 9,953 | 20,383 | 0 |
| LBA | 11,663 | 5,264 | **45.1%** | 0 | 6,395 | 4 | 0 |
| fold | 62,896 | 24,764 | **39.4%** | 10 | 16,618 | 21,502 | 2 |
| **Total** | **518,470** | **421,704** | **81.3%** | 3,070 | 43,496 | 50,038 | 2 |

**Interpretation:**
- **Exact:** kcat/pH, full-length AlphaFold; LBS 99.8% and PPIS 99.9% are essentially consistent, with isolated substitutions/indels.
- **Systematic differences within the same protein:**
  - Fold aa_seq is a whole chain, often with leading X placeholders, while structures are CATH domains; length/start differs, sampled median identity 0.98.
  - LBA mostly has longer aa_seq because deposited structures omit unresolved residues; sequence agreement is high but structures are incomplete.
  - Function mismatches are mainly longer aa_seq: SIFTS subsequences versus unresolved structural residues.
  - PPI exact agreement is 81.5%; remaining cases mostly unresolved residues or full-length AF versus truncated chains of the same protein.
  - SSP: sampled median identity 0.989 among ≈14% non-exact cases (all ≥0.95), reflecting single-residue indels. **Residue tasks can have one-residue label–structure offsets requiring assessment.**
- **Definite bad samples:** two empty fold structures, 4bpe000 and 4v19000, with zero residues; repair or remove separately.
- **Overall:** no evidence of unrelated sequence–structure pairing; differences concern the same protein: missing residues, indels, leading X, domains versus whole chains.

**7.2 Three-way residue-task consistency (LBS / ppis / ss)**

For three residue-level tasks, compare aa_seq, structure-derived sequences, and residue annotations.

- **Script:** ~~`analysis/scripts/residue_level_consistency.py`~~, removed during 2026-09 cleanup; snapshots retained.
- **Results**: `analysis/output/residue_level_consistency.json`

**Finding: len(labels)==len(aa_seq) for all 145,017 samples**, so annotations align to aa_seq residue by residue. Label–structure consistency therefore reduces to aa_seq–structure consistency.

| Task | Label length==aa_seq | aa_seq==structure | Main discrepancy | Median non-exact identity |
|------|:---:|------:|------|:---:|
| LBS | 100% (96,256) | 100% | Replace aa_seq with structural sequences in 211 terminal-substitution cases; labels align fully. | 1.000 |
| ppis | 100% (8,695) | 99.9% | Remove 154 low-quality/placeholders and align/crop 61 structures; 13 remaining equal-length cases: six substitutions and seven cancelling internal indels. | 0.989 |
| ss | 100% (40,066) | 85.8% | 4,737 internal and 501 terminal unresolved-residue cases, masked with -1. | 0.990 |

**Label encoding and masking**:
- **LBS:** labels {0,1}, binding-site binary classification without masks; aa_seq matches structures 100%, with aligned residue annotations.
- **PPIS:** labels {0,1,-1}; -1 at 140 residues denotes structural residues lacking source interface labels. Structure/sequence/label lengths are unified.
- **SSP:** labels {0..8,-1}. **-1 masks 30,051 unresolved aa_seq residues.** Nonnegative-label count equals structural residues in 97.8%; -1 count equals aa_len minus structure length in 98.9%. Full aa_seq and resolved-only structures differ by design, with missing positions accounted for by -1.

**Task-specific diagnosis**:
- **LBS, 211 cases:** replace aa_seq with extracted structural sequences, removing terminal substitutions; internal labels align and effects are negligible.
- **PPIS, 13 cases:** remove 154 unrecoverable low-quality/placeholders; align/crop 61 structures, skipping insertions and retaining residues matched one-to-one to aa_seq. The remaining 13 are equal-length (six substitutions, seven cancelling indels), with aligned labels.
- **SSP, 5,238 cases:** mostly internal/terminal unresolved density, with complete sequences and incomplete structures. -1 masks preserve consistency by design; median non-exact identity 0.99.

**Structure-input SSP view:** full-sequence splits with -1 masks suit sequence models, but ≈2.6% (207/8,079 test) masks miss internal gaps, so simply dropping -1 silently misaligns labels. The structure-native view at datasets/ss_prediction/splits/splits_struct_native/ (archived script: datasets/ss_prediction/output/scripts_archive/build_ss_struct_native.py) uses structural CA sequences as aa_seq and globally maps DSSP ss8 labels to every structural residue. All 40,066 rows satisfy len(labels)==len(aa_seq)==M(pdb), strict 1:1 alignment: 98.8% exact DSSP, 1.2% alignment, six fallbacks, mean identity 0.99992, no quarantines. Later deprecated below.

> **Deprecated 2026-08-28:** main splits now include struct_label from the same source/values as struct_native, enumerating all ATOM residues and retaining synchronized OOD columns. Archive: datasets/ss_prediction/output/splits_struct_native_deprecated/. Structure models should use main splits plus struct_label.

**7.3 Sequence legality (invalid amino-acid characters)**

Scan **552,093 sequences /≈168 million residues** across 12 task views character by character. Script ~~analysis/scripts/sequence_legality.py~~ was removed in 2026-09; results: analysis/output/sequence_legality.json. Classes: standard 20 ACDEFGHIKLMNPQRSTVWY; IUPAC ambiguity/rare X B Z J U O; chain separator |; all others invalid.

**Zero invalid characters and zero empty/NaN sequences:** no lowercase, digits, punctuation (-./*), whitespace, or corrupt characters. Nonstandard characters are biologically valid or formatting conventions; no cleaning needed.

| Task | #Sequences | Empty/NaN | Ambiguous codes | Separators | Invalid | Nonstandard character details |
|------|------:|:---:|------:|------:|:---:|------|
| kcat | 22,631 | 0 | 0 | 0 | 0 | — |
| optimal_ph | 11,615 | 0 | 0 | 0 | 0 | — |
| fold(CATH) | 62,517 | 0 | 10,922 | 0 | 0 | `X`×47,977(leading placeholders/unresolved residues) |
| func_EC | 22,123 | 0 | 15 | 0 | 0 | `X`×611 |
| func_GO-BP | 22,116 | 0 | 13 | 0 | 0 | `X`×478 |
| func_GO-CC | 21,221 | 0 | 21 | 0 | 0 | `X`×618 |
| func_GO-MF | 25,872 | 0 | 22 | 0 | 0 | `X`×735 |
| LBA | 11,663 | 0 | 0 | 6,393 | 0 | `\|`×12,657(multichain complex separators) |
| LBS | 96,256 | 0 | 0 | 0 | 0 | — |
| PPI | 167,252 | 0 | 272 | 0 | 0 | `X`×409, `B`×101, `U`×58(selenocysteine), `Z`×41, `O`×2(pyrrolysine) |
| PPIS | 8,695 | 0 | 0 | 0 | 0 | — |
| SSP | 40,066 | 0 | 0 | 0 | 0 | — |
| SSP(struct-native) | 40,066 | 0 | 16 | 0 | 0 | `X`×161(nonstandard residues from CA extraction) |

> PPI sequence counts include aa_seq1 and aa_seq2 separately. LBA | separates chains within a complex, e.g. a homotetramer is 4×468 aa. Dataloaders must split chains and handle X/B/Z/J/U/O through tokenizer fallbacks (most pLMs already provide X/<unk>).

**7.4 Schema and label-format consistency**

After project-wide schema unification (label/aa_seq(1,2)/struct_file(1,2), OOD_ prefixes, OOD_NewEC_L*/OOD_LongTail(_EC_L*)), audit 13 task views and 39 split files. Script: analysis/scripts/format_consistency.py; results: analysis/output/format_consistency.json.

Check complete core columns; no obsolete labels/file_type/pdb_file*/is_time_cutoff/OOD-*/OOD_NewFunction*/OOD_LongTailFunction*/OOD_CATH_*/protein_*_seq; nonempty labels; OOD only in test with Default; residue label/sequence length agreement; no cross-split unique_id leakage; matching train/val column sets.

**All 194/194 structured checks pass.**

| Task | Train/Val/Test | Label type | Complete core columns | Obsolete columns | Cross-split leakage | Residue lengths |
|------|---------------:|:---:|:---:|:---:|:---:|:---:|
| kcat | 15479/2150/5002 | float(regression) | ✅ | 0 | None | — |
| optimal_ph | 7797/866/2952 | float(regression) | ✅ | 0 | None | — |
| fold(CATH) | 34963/4994/22939 | str(C.A codes; require dtype=str) | ✅ | 0 | None | — |
| func_EC | 15389/1687/4243 | str | ✅ | 0 | None | — |
| func_GO-BP | 15232/1742/4121 | str | ✅ | 0 | None | — |
| func_GO-CC | 13878/1730/4256 | str | ✅ | 0 | None | — |
| func_GO-MF | 18179/2084/4609 | str | ✅ | 0 | None | — |
| LBA | 8366/930/2367 | float(regression) | ✅ | 0 | None | — |
| LBS | 71093/7938/17225 | str | ✅ | 0 | None | ✅ |
| PPI | 65020/7195/16457 | int(0/1) | ✅ | 0 | None | — |
| PPIS | 6043/672/1980 | str | ✅ | 0 | None | ✅ |
| SSP | 28788/3199/8079 | str | ✅ | 0 | None | ✅ |
| SSP(struct-native) | 28788/3199/8079 | str | ✅ | 0 | None | ✅ |

The other 12 task views use semantically appropriate labels: float for kcat/pH/LBA regression, int for PPI, strings for EC/GO multilabel and LBS/PPIS/SSP residue targets.

**Paper use**: Supplementary Material(data quality)

**7.5 struct_label: structure-aligned residue labels (added 2026-08)**

All splits of LBS/PPIS/SSP add a final struct_label column, aligned to all PDB ATOM residues (including residues without CA), ordered by first (chain, resseq, icode) occurrence; -1 is ignored. Structure models can train/evaluate directly from pdbs/ plus struct_label without realigning sequences.

- **SSP:** import authoritative DSSP-native labels from the former splits_struct_native. Full-sequence main labels have local gap-placement offsets in ≈8% of training rows and are not used by structure models. Archive: datasets/ss_prediction/output/splits_struct_native_deprecated/.
- **LBS/PPIS:** structure-contact extraction generates native alignment; equal-length rows match full-sequence labels 100%, with -1 mapped to zero.
- **Tool:** .skills/struct-label-alignment/add_struct_label.py, supports --dry-run/--overwrite.

**7.6 Structure library (pdbs/) repairs, 2026-08**

Scan and repair structures across all ten tasks. Post-repair rescans find zero remaining issues; all files parse as single models in Bio.PDB; no structure-longer-than-sequence samples:

| Task | Issue | Scale | Repair |
|------|------|------|------|
| PPIS | PDB column shift: extra space after resname shifts chain/resseq, making the library unreadable in Bio.PDB; five truncated chain concatenations | 8,833 files | Archived ppis_prediction/output/scripts_archive/18_fix_pdb_column_layout.py; atomwise coordinates unchanged. Patch the format string in scripts/08_extract_representative_structures.py. |
| LBS | NMR models concatenated without ENDMDL; duplicate residues | 549 files | Archived output/scripts_archive/26_fix_multimodel_pdbs.py re-extracts model 1. Patch scripts/07_extract_protein_chains.py: stop at ENDMDL, first CIF model, unquote atom names. |
| SSP | Concatenated models: 23 with MODEL records, three without | 26 files | Archived output/scripts_archive/fix_multimodel_pdbs.py; extraction uses structure[0]. |
| Fold /func /PDA | Same NMR concatenation | 4,076 /1,940 /203 files | Per-task archived output/scripts_archive/fix_multimodel_pdbs.py plus extraction-script patches. |
| PPI /LBA /kcat /pH | No issues found | — | kcat has 7,057 files with MODEL records but only one model; no changes needed. |

Pre-repair backups are in each task output/backup_*/ directory.

---

#### 8. protein-design-archive additional test sets (added 2026-08)

**Status**: ✅ Complete

**File**: `datasets/protein-design-archive/splits/cath_design.csv`, `ec_design.csv`

**Documentation**: `datasets/protein-design-archive/DATA_PROCESS.md`

Additional designed-protein test sets from 1,879 PDA entries evaluate generalization to proteins absent from nature. Use all entries, not only purely artificial designs; merge only exact sequences, retaining small differences without 95% identity clustering. pdbs/ stores complete deduplicated single-chain structures.

| Test set | Samples | Task | Notes |
|--------|------:|----------|------|
| cath_design.csv | 122 | fold_classification | CATH-annotatable designs in the current label space, excluding fold train/val ID or ≥95% identity overlap; same OOD schema as main test. |
| ec_design.csv | 31 | func_prediction (EC) | No EC transfer within clusters, preserving per-sample label provenance. |

- All Default=True: designs do not undergo main-dataset splitting, so prior-reserved Default=False semantics do not apply.
- Compute seq/TM/Orphan/NewEC/LongTailEC against the main task train/val references.
- Reuse pdbs/ and the pipeline for additional design-test tasks later.

**Paper use**: Supplementary / Figure(generalization to de novo designs)

---

### II. Planned analyses without models

#### 6. Comparison with existing benchmarks

**Status**: ⏳ Pending

**Content**: Quantitative comparison with CASP, CAFA, PDBbind, CATH, EnzyBase, and other benchmarks

**Dimensions**: Task counts, OOD scenarios, sample sizes, structure availability, and evaluation dimensions

**Paper use**: Table 2(Related Work / Motivation)

**Effort**: Manual literature compilation; no code required

---

## Part II: Experiments requiring model training

### I. Baseline evaluation (core experiments)

#### 1. In-distribution performance baseline

**Goal**: Train using train/val and evaluate on the ID test subset

**Method**: Use test samples whose OOD flags are all False as the ID test subset

**Models**: Fine-tune pretrained models such as ESM-2, ProstT5, and AlphaFold2-derived models

**Metrics**: Task-specific: ACC, F1, MCC, MAE, Spearman, etc.

**Paper use**: Table 4(ID baseline)

**Compute**: High (GPU training)

---

#### 2. OOD performance degradation

**Goal**: Compare ID performance with each OOD scenario

**Method**: Evaluate ID and each OOD test subset using the same model

**Expected result**: Lower performance under OOD, with degradation related to shift difficulty

**Paper use**: Figure 4(OOD degradation plot), Table 5(per-scenario performance)

**Compute**: High (additional inference using Experiment 1 models)

---

#### 3. OOD threshold sensitivity

**Goal**: Assess the choice of OOD thresholds

**Method**: Vary SeqRedundancy (40%, 30%, 20%) and TM-score (0.5, 0.3, 0.2), then examine degradation curves

**Paper use**: Figure 6(threshold-sensitivity curves)

**Compute**: High (repeated OOD subset selection and evaluation)

**Note**: Requires model performance on different OOD subsets

---

#### 4. Combined OOD scenarios (difficulty analysis)

**Goal**: Analyze model performance on samples carrying multiple OOD flags

**Method**: Compare single versus combined OOD degradation using model outputs

**Paper use**: Figure 7(combined OOD analysis)

**Compute**: High (additional analysis of Experiment 2 outputs)

**Note**: Requires performance on combined OOD subsets; model-dependent analysis

---

#### 5. Model scale and OOD performance

**Goal**: Test whether larger models are more robust to OOD

**Method**: Compare ESM-2 scales (8M/35M/150M/650M/3B) across OOD scenarios

**Paper use**: Figure 8(model scale versus OOD performance)

**Compute**: Very high (training/inference at multiple model scales)

---

### II. OOD detection using model outputs

#### 6. Unsupervised OOD detection

**Goal**: Detect OOD test samples without OOD-label supervision

**Method**:
- Use pretrained embeddings/scores, e.g. MSP or Mahalanobis distance
- Reconstruction errors, e.g. VAE or AE
- Gradient/uncertainty methods, e.g. MC-Dropout

**Evaluation**: AUROC, AUPR, FPR@95%TPR

**Paper use**: Figure 5(OOD detection ROC), Table 6(detection comparison)

**Compute**: Medium-high (train detectors or extract embeddings)

**Note**: Requires model embeddings or probability distributions; model-dependent

---

### III. OOD mitigation through model training

#### 7. Data augmentation and resampling

**Goal**: Mitigate OOD performance degradation

**Method**:
- Oversample long-tail labels
- Augment ExtremeLength/IDR samples
- Expand training data using sequence similarity

**Evaluation**: Preserve ID performance and improve OOD performance

**Paper use**: Table 7(mitigation strategy comparison)

**Compute**: High

---

#### 8. Domain adaptation / transfer learning

**Goal**: Improve generalization to OOD domains

**Method**:
- DANN (Domain-Adversarial Neural Network)
- CORAL(feature alignment)
- Data mixing such as Mixup/CutMix

**Evaluation**: Same as Experiment 5

**Paper use**: Compare with Experiment 5

**Compute**: High

---

#### 9. Ensembles and uncertainty quantification

**Goal**: Improve recognition of OOD samples

**Method**:
- Deep Ensemble(multiple-model ensemble)
- SWA (Stochastic Weight Averaging)
- Evidential Deep Learning

**Evaluation**: ID performance + OOD detection AUROC + calibration (ECE)

**Paper use**: Table 8(ensemble comparison)

**Compute**: Very high (multiple model training runs)

---

## Planned paper figures and tables

### Model-free figures/tables (ready to generate)

| Figure/table | Content | Source |
|------|------|------|
| Table 1 | Dataset scale statistics | Structure/sequence statistics |
| Figure 1 | OOD coverage heatmap | OOD coverage matrix |
| Figure 2 | OOD overlap matrix | OOD overlap analysis |
| Figure 3 | Label long-tail/Zipf curves | Label distribution analysis |
| Table 2 | Benchmark comparison | Pending(literature compilation) |

### Model-dependent figures/tables (after experiments)

| Figure/table | Content | Source |
|------|------|------|
| Table 3 | ID baseline | Experiment 1 |
| Figure 4 | OOD degradation | Experiment 2 |
| Figure 5 | OOD detection ROC | Experiment 3 |
| Table 4 | OOD detection comparison | Experiment 3-4 |
| Table 5 | Mitigation comparison | Experiment 5-7 |
| Figure 6 | Threshold sensitivity | Experiment 8 |
| Figure 7 | Combined OOD effects | Experiment 9 |
| Figure 8 | Model scale versus OOD | Experiment 10 |

---

### Minimum viable paper (MVP)

- **Required**: All Stage 1 work +Experiments 1–4
- **Optional extensions**: Experiment 5 (mitigation) +6 (sensitivity) +7–9 (combined effects)
- **Full version**: All experiments

---

## File locations

Analysis scripts and intermediate outputs relative to the POOR/ root:

```
analysis/
├── scripts/                                  # Analysis scripts: statistics and paper figures only
│   ├── dataset_scale_stats.py                # §1 Dataset scale statistics (added 2026-08)
│   ├── ood_coverage_matrix.py                # §2 OOD coverage matrix +§3 overlap (scenarios updated 2026-08)
│   ├── label_distribution.py                 # §4 Label long-tail distribution (added 2026-08)
│   ├── compute_structure_source_stats.py     # §5 Structure-source coverage (ROOT fixed 2026-08)
│   ├── plot_figure2_ood_orthogonality.py     # Figure 2 Plotting
│   └── plot_ood_jaccard_heatmaps.py          # OOD Jaccard heatmaps
└── output/                                   # Intermediate analysis outputs
    ├── dataset_scale_stats.json              # §1 (2026-08-28)
    ├── ood_coverage_matrix.csv / .md         # §2 (2026-08-28)
    ├── ood_overlap_matrix.csv                # §3 (2026-08-28)
    ├── label_distribution.json               # §4 (2026-08-28)
    ├── structure_source_stats.json           # §5 (2026-08-28)
    ├── aa_seq_structure_consistency.json     # §7.1(2026-07 snapshot)
    ├── residue_level_consistency.json        # §7.2(2026-07 snapshot)
    ├── sequence_legality.json                # §7.3(recomputed 2026-08-28)
    ├── format_consistency.json               # §7.4(2026-07 snapshot)
    └── ood_jaccard_heatmaps/                 # Figure 2 Images
```

> Note (2026-08-28): after moving task-specific one-off scripts into task directories, all workflows were renumbered NN_*.py in execution order. One-off fixes (fix_multimodel_pdbs.py, 18_fix_pdb_column_layout.py, fix_lbs_aaseq.py, build_ss_struct_native.py, ss_struct_alignment_probe.py, patch_struct_label_from_native.py, etc.) are archived in per-task output/scripts_archive/. Repairs are applied to data and root causes patched in normal extraction scripts. crop_ppis_structures.py is now ppis_prediction/scripts/12_crop_ppis_structures.py. Each four-stage DATA_PROCESS.md is authoritative.

Reusable processing workflows are provided as project skills (.skills/):

```
.skills/
├── cath-ted-domain-annotation/    # CATH+TED domain annotations for topology/superfamily holdouts
├── ec-function-ood-annotation/    # EC annotation and canonical NewEC/LongTailEC (all-not-in, L3/L4, thresholds 5/10)
├── homology-ood-annotation/       # Sequence homology / Orphan / TM-score OOD annotation
├── ood-annotation-toolkit/        # Generic OOD tools with Default/InD/NaN semantics
├── pood-benchmark-process/        # Overall dataset construction workflow
└── struct-label-alignment/        # Residue-task struct_label generation and validation
```

> Note (2026-09-01): analysis/scripts/ retains only the six statistics/plotting tools listed above. QA scripts (comprehensive_audit.py, struct_audit.py, audit_foldseek_easysearch.py, verify_aa_seq_vs_structure.py, residue_level_consistency.py, sequence_legality.py) and one-off combinatorial_ec_feasibility.py were removed; outputs such as output/comprehensive_audit_20260830/ remain for inspection. ppis_prediction/scripts/12_crop_ppis_structures.py now inlines Biopython mappings instead of importing verify_aa_seq_vs_structure.py.

> **Note:** root .gitignore excludes pdbs/, splits/, output/, and check/. Thus analysis/output/ and per-task output/ are **not version controlled**; analysis/scripts/ and analysis/benchmark_analysis.md are tracked.
