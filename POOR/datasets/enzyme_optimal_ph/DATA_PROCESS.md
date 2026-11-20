# Enzyme Optimal pH Prediction (optimal_ph_prediction): Dataset Processing

> Task: optimal_ph_prediction (regression: enzyme pH optimum).
> Source: Zenodo [18405148](https://zenodo.org/records/18405148), *Dataset EnzyBase12k: A Curated Dataset for Predicting Enzyme pH Optima Using pHoptNN*.
> Location: `datasets/enzyme_optimal_ph/`.
>
> Complete four-stage linear workflow to reproduce current `splits/` and `pdbs/` from raw data. Statistics and column names were recomputed and checked read-only against current `splits/*.csv` headers.

## Task and current scale

Predict optimal pH from enzyme sequence/structure. `label` stores the pH optimum as a string.

| Split | File | Samples | Description |
|:---|:---|---:|:---|
| Train | `splits/optimal_ph_prediction_train.csv` | 7,794 | AF structures; four columns |
| Val | `splits/optimal_ph_prediction_val.csv` | 866 | AF structures; four columns |
| Test | `splits/optimal_ph_prediction_test.csv` | 2,955 | PDB plus extreme-length samples; 34 columns including full OOD annotations |
| **Total** | — | **11,615** | 11,615 `.cif` files in `pdbs/` |

Test includes 2,717 `Default=True` original PDB samples, 211 `OOD_ExtremeLong` (`seq_length >1000`), and 27 `OOD_ExtremeShort` (24 `seq_length <60` plus three moved from train on 2026-09-01; §4.2.1). Both extreme groups have Default=False.

Scripts in execution order:

| No. | Script | Stage |
|:---|:---|:---|
| 01 | `01_verify_structures.py` | ① Extract and verify structures |
| 02 | `02_split_dataset.py` | ④ Initial split: extreme-length holdout, PDB test, AF 9:1 |
| 03 | `03_prepare_splits_and_pdbs.py` | ④ Unified schema and structure export |
| 04 | `04_compute_ood_markers.py` | ④ Sequence/structure homology, orphan, and IDR OOD |
| 05 | `05_add_phbin_longtail_ood.py` | ④ pH-bin long-tail OOD |
| 06 | `06_ted_augment_cath_ood.py` | ④ CATH+TED fold/superfamily holdouts |
| 07 | `07_move_short_train_rows_to_test.py` | ④ Move aa_seq<60 training rows to test (§4.2.1) |
| Diagnostics | `diag_analyze_metadata.py`, `diag_validate_sequence_structure.py`, `diag_cluster_sequences_mmseqs.py` | Diagnostics and quality control (§7) |

---

## Stage 1: Obtain and parse sources

### 1.1 Zenodo files and downloads

| File | Size reported by API | Description |
|:---|---:|:---|
| `metadata.csv` | ~6.4 MB | 11,615 rows ×13 columns |
| `DATASET.zip` | ~1.87 GB | PDB/PQR structures |
| `Dataset_mmcifs.zip` | ~559 MB | CIF structures |

```bash
# Metadata
curl -L -o data/metadata.csv "https://zenodo.org/api/records/18405148/files/metadata.csv/content"

# Large structure archives: resumable downloads recommended.
wget -c -t 0 --timeout=300 -O data/DATASET.zip \
  "https://zenodo.org/api/records/18405148/files/DATASET.zip/content"

wget -c -t 0 --timeout=300 -O data/Dataset_mmcifs.zip \
  "https://zenodo.org/api/records/18405148/files/Dataset_mmcifs.zip/content"
```

Store downloads in `data/`; download logs are in `logs/`.

### 1.2 Extract and verify structures

**Script:** `scripts/01_verify_structures.py`<br>
**Input:** `data/DATASET.zip`, `data/Dataset_mmcifs.zip`, `data/metadata.csv`<br>
**Output:** `data/dataset/` (`DATASET/EnzyBase12k_pbs|pqrs`), `data/dataset_mmcifs/` (`EnzyBase12k_mmcifs`), `output/structure_verification.json`

- Extract both archives unless already extracted, count files, and match names to uniprot_id /pdb_id_final.
- Verification: **11,615** files for each PDB/PQR/CIF format, matching metadata rows; all filenames match uniprot_id.

### 1.3 Metadata fields

| Field | Type | Description |
|:---|:---|:---|
| `uniprot_id` | str | UniProt accession, later unique_id |
| `ph_optimum` | float | Enzyme optimal pH target |
| `structure_source` | str | PDB /AF |
| `structure_mode` | str | AF2_downloaded /AF3 /AF3_template |
| `pdb_id_final` | str | Final PDB ID, PDB source only |
| `resolution` | float | Å, PDB source only |
| `pdb_chain` | str | PDB chain ID, PDB source only |
| `ec_id` | str | EC numbers; **782 UniProts contain multiple ECs separated by `; `** (§9) |
| `cut_to_uniprot_domain` | int | Whether truncated to a UniProt domain |
| `systematic_name` | str | Systematic name |
| `organism` | str | Source organism |
| `seq_length` | int | Sequence length |
| `uniprot_seq_cut` | str | Truncated amino-acid sequence |

Metadata statistics from diag_analyze_metadata.py: 11,615 samples, no repeated UniProts; mean pH 7.20, median 7.50, range 1.5–12.5; AF 76.3%, PDB 23.7%; median PDB resolution 1.95 Å (n=2,735); median length 386, range 32–1500; EC coverage 98.6%.

---

## Stage 2: Aggregate and compute labels

Use source ph_optimum from metadata.csv directly as label, without multi-source aggregation. ec_id also comes directly from metadata.csv and supplies Stage 4 EC OOD; retain this historical source (§9).

---

## Stage 3: Redundancy reduction and refinement

> **No samples were removed by quality/redundancy filters.** The following steps are diagnostics only.

### 3.1 Sequence–structure validation: diagnostic

**Script:** `scripts/diag_validate_sequence_structure.py`<br>
**Input:** metadata.csv plus extracted structures, preferring CIF<br>
**Output:** `output/sequence_structure_validation.csv` and `_summary.json`

Globally align first-chain structure sequences to uniprot_seq_cut using edlib. All 11,615 samples parse successfully. Mean/median identity 0.938/0.959; 77.6% have identity ≥0.90, none <0.50. Structure truncation places 22.4% between 0.70 and 0.90; no severe discrepancies. **No identity filtering applied.**

### 3.2 MMseqs2 95% clustering: diagnostic, not applied

**Script:** `scripts/diag_cluster_sequences_mmseqs.py`<br>
**Input:** metadata.csv uniprot_seq_cut<br>
**Output:** `output/sequences.fa`, `output/mmseqs_cluster_95/`, `output/mmseqs_cluster_95_summary.json`

Parameters: min-seq-id 0.95, coverage 0.8, cov-mode 0. Clustering 11,615 sequences gives 11,058 clusters (10,622 singletons, maximum size nine), 557 redundant rows (4.8%). **Diagnostic only**, demonstrating low redundancy; all 11,615 samples remain.

### 3.3 Candidate filters

| Criterion | Suggested threshold | Applied? | Notes |
|:---|:---|:---|:---|
| Length | 60–1000 | ✅ at splitting | Out-of-range samples become extreme-length test OOD (§4.1) |
| PDB resolution | ≤3.0 Å | ❌ | Candidate only |
| pH outliers | 1.0–13.0 | ❌ | Actual 1.5–12.5, no filtering needed |
| Sequence–structure agreement | Identity ≥0.90 | ❌ | Candidate only |
| Redundancy | MMseqs2 95% identity | ❌ | Diagnostic only (§3.2) |

---

## Stage 4: Split and annotate OOD

### 4.1 Initial split

**Script:** `scripts/02_split_dataset.py`<br>
**Input:** `data/metadata.csv`<br>
**Output:** intermediate `splits/{train,val,test,ood_extremelength}.csv`, consumed and removed next; `output/split_summary.json`

1. **Extreme-length holdout:** reserve 235 samples with seq_length <60 or >1000.
2. **Test:** all remaining PDB-source samples, 2,717.
3. **Train/val:** randomly split remaining AF samples 9:1, random_seed=42.

### 4.2 Unified schema and structure export

**Script:** `scripts/03_prepare_splits_and_pdbs.py`<br>
**Input:** initial split CSVs and extracted structures<br>
**Output:** `splits/optimal_ph_prediction_{train,val,test}.csv`, 11,615 `pdbs/*.cif`, `output/split_mapping.json`

- Merge ood_extremelength.csv into test; remove initial intermediate files.
- Original PDB test rows have Default=True; extremes have Default=False.
- Add OOD_ExtremeLong (seq_length >1000) and OOD_ExtremeShort (<60).
- Extract aa_seq from structures, preferring CIF; store pH label as string.
- Copy structures to pdbs/, named by unique_id; prefer CIF, fall back to PDB.
- **Current schema:** train/val: unique_id, aa_seq, struct_file, label; test additionally has Default and two length flags. Historically: unique_id, file_type, aa_seq, labels, migrated project-wide (§9).

Recomputed current splits:

| Split | Samples | Structures | Mean pH | Median aa_seq length |
|:---|---:|:---|---:|---:|
| Train | 7,794 | AF | 7.129 | 367 |
| Val | 866 | AF | 7.139 | 363 |
| Test | 2,955 | PDB +extremes | 7.390 | 347 |

> Means/medians use current label/aa_seq. Structure-derived aa_seq differs from metadata seq_length.

### 4.2.1 Extreme-length correction, 2026-09-01

**Script:** `scripts/07_move_short_train_rows_to_test.py`

Initial holds use metadata seq_length (UniProt-domain length), whereas aa_seq is structure-derived. Three metadata-normal training samples have structural sequences <60: P05959 55 aa, B1KRG2 59 aa, P0C7A9 46 aa, all label 8.0. Move them to test with Default=False, OOD_ExtremeShort=True, OOD_ExtremeLong=False, other OOD columns NaN, InD=False. Structures already reside in pdbs/. Train 7,797→7,794; test 2,952→2,955; ES 24→27. Backup: `output/backup_20260901/`.

> Pending discrepancy: approximately 90 OOD_ExtremeLong=True rows have aa_seq lengths 747–999. EL still follows metadata seq_length; unifying to len(aa_seq) requires recomputing EL flags and updating statistics/documentation.

### 4.3 General OOD conventions

Except Default, length flags, and InD, OOD columns including seq_Redundancy_* and TM-score_* are computed **only for Default=True**. Default=False extremes receive **NaN**, meaning uncomputed, finalized in §4.10.

Canonical methods reside in .skills/homology-ood-annotation, ood-annotation-toolkit, cath-ted-domain-annotation, and ec-function-ood-annotation. Each step below gives local and canonical tools.

### 4.4 Sequence homology: seq_Redundancy_* and OOD_Orphan

**Script:** `scripts/04_compute_ood_markers.py`<br>
**Input:** final train/val/test CSVs<br>
**Output:** update test.csv; `output/ood_markers_summary.json`

- MMseqs2 easy-search (`-s 7`), Default test versus train+val; use maximum query pident.
- seq_Redundancy_90…30: max pident **<** threshold, 10% steps. **Changed ≤ to < on 2026-08-29** for project consistency. At one-decimal precision, 35 boundary rows had extra True flags (3/1/3/4/9/8/7 across thresholds); recomputed and corrected. Backup `output/backup_before_threshold_lt_unify/`; audit `analysis/output/ph_threshold_audit/`; cache `output/ood_mmseqs_test_vs_trainval.m8`. P14618 at exactly 90.0 loses its sole OOD flag, increasing InD 49→50.
- OOD_Orphan: no significant train+val hit.
- Canonical tool: `.skills/homology-ood-annotation/scripts/seq_homology_ood.py`, supports easy-search and --m8-cache replay.

### 4.5 Structural similarity: TM-score_*

**Script:** same `04_compute_ood_markers.py` step.

- Foldseek easy-search, Default test against train+val CIF structures in pdbs/; maximum alntmscore.
- TM-score_0.9…0.3: maximum **<** threshold, step 0.1. Changed ≤ to < on 2026-08-29 with sequence flags; no exact-boundary cases, so values unchanged. Cache: `output/ood_foldseek_test_vs_trainval.m8`.
- Canonical: `.skills/homology-ood-annotation/scripts/struct_homology_ood.py`.

### 4.6 Disorder: OOD_IDR

**Script:** `scripts/04_compute_ood_markers.py`; correct IDR logic is integrated, replacing archived second-pass compute_ood_idr.py.

- metapredict v3 predict_disorder_fasta batches; residue disorder propensity >0.5 denotes disorder.
- **Residue ratio** (disordered residues /sequence length) **>0.1**. Default rows only; others receive NaN before saving (§4.10).
- Historical task threshold 0.1 differs from SSP/func/fold/PPI 0.3; see the definition table in `.skills/ood-annotation-toolkit/scripts/idr_ood.py`.
- Intermediate `output/test_sequences_for_idr.fasta`; **39** positives among 2,717 Default rows, reflecting mostly structured enzymes.

### 4.7 Fold/superfamily holdouts

**Wrapper:** `scripts/06_ted_augment_cath_ood.py`<br>
**Input:** test.csv and accession list `output/optimal_ph_ted_aug_acc_list.csv`<br>
**Output:** `output/optimal_ph_acc_cath_ted_labels.csv`; update two test columns

Use the sole canonical implementations `.skills/cath-ted-domain-annotation/scripts/cath_ted_annotate.py` and cath_ted_holdout.py.

**UniProt mode**, with unique_id equal to accession:

- Labels: experimental CATH v4.4 unioned through reverse SIFTS UniProt→PDB-chain mapping ∪ putative full-length TED without domain filtering. REST API cache shared at data/ted/ted_api_cache.json; H supplies C.A.T.H, T only C.A.T.
- Reference: **training** CATH+TED union, with symmetric train/val/test augmentation.
- Any semantics: any topology/superfamily absent from training →True; unannotated →False.

**Coverage:**

- 11,615 accession keys, all TED cache hits, offline, zero new queries.
- Sources: ted-only 8,167 /both 2,787 /none 595 /cath-only 66.

| Definition | Topology coverage across 11,615 accessions |
|------|----------------|
| CATH only | 24.6% |
| CATH+TED | **94.9%** |

**Holdouts, Default n=2,717:** fold 91, superfamily 311, unannotated →False 65. Backup `output/backup_before_ted_aug/`.

### 4.8 EC function OOD: NewEC and LongTailEC

**Canonical skill:** .skills/ec-function-ood-annotation; batch drivers scripts/add_unified_newec.py and add_unified_longtail_ec.py write outputs.

- EC source: metadata.csv ec_id; partial hierarchy allowed, e.g. 3.4.22.- →L3 3.4.22.
- **Multiple ECs:** split `; ` strings in 782 UniProts before evaluation.
- OOD_NewEC_L3/L4: **all-not-in**; annotated and all ECs at that level absent from **train+val** →True; unannotated →False.
- OOD_LongTail_EC_{L3,L4}_le{5,10}: **all-not-in**; all ECs have **training** frequency <5 /<10, **including zero** →True; unannotated →False.
- Summaries: output/newec_unified_summary.json and longtail_ec_unified_summary.json. Backups: output/backup_before_newec_unify/ and backup_before_longtail_unify/.

| Column | True, Default n=2,717 |
|------|------:|
| `OOD_NewEC_L3` | 20 |
| `OOD_NewEC_L4` | 647 |
| `OOD_LongTail_EC_L3_le5` | 118 |
| `OOD_LongTail_EC_L3_le10` | 221 |
| `OOD_LongTail_EC_L4_le5` | 1,763 |
| `OOD_LongTail_EC_L4_le10` | 2,208 |

> With only 7,794 training samples and UniProt-counted EC frequencies, L4 tails cover 64.9%–81.3%, limiting discrimination. **L3 is more informative.**

> Superseded: archived compute_additional_ood.py used single EC strings and ≤ thresholds in OOD_NewFunction /OOD_LongTailFunction_L3_le{5,10}. The canonical skill replaces it.

### 4.9 pH-bin long-tail OOD

**Script:** `scripts/05_add_phbin_longtail_ood.py`, extracted from archived compute_additional_ood.py with unified names<br>
**Input:** final train/test CSVs<br>
**Output:** update test.csv and output/ood_markers_summary.json

- pd.cut bins: **[1.5, 13.5), width 0.5**.
- Default test samples are True when their **training** bin frequency ≤threshold.
- Canonical replay: `.skills/ood-annotation-toolkit/scripts/value_bin_longtail_ood.py`; parameter table includes this task.

| Column | Threshold | True, Default n=2,717 |
|------|------|------:|
| `OOD_LongTail_pHbin_le50` | Training bin count ≤50 | 25 |
| `OOD_LongTail_pHbin_le100` | Training bin count ≤100 | 58 |

### 4.10 Validate non-Default NaN and finalize InD

**Tool:** `.skills/ood-annotation-toolkit/scripts/finalize_ood_columns.py`, --mode nan then --mode ind; validate with check_splits.py in the same skill.

- **NaN:** Default=False rows retain only length flags; all computed OOD_*, seq_Redundancy_*, TM-score_* flags are empty/NaN, meaning not computed rather than False. Since 2026-08-28, scripts 04/05/06 and skills write these at generation time. --mode nan is a **fallback check**, currently no-op. Historical backup: output/backup_before_nondefault_nan/.
- **InD:** all OOD_*, seq_Redundancy_*, TM-score_* False, treating NaN as False. **50** rows (1.7% of test). Backup: output/backup_before_ind/.
- I/O: `pd.read_csv(dtype=str, keep_default_na=False)` and `to_csv(index=False, lineterminator="\n")`. Back up first; verify non-target columns unchanged afterward.

### 4.11 Current OOD counts: Default n=2,717

| Column | Definition | True |
|:---|:---|---:|
| `seq_Redundancy_90`–`30` | Max pident <90/80/70/60/50/40/30 | 2,540 /2,396 /2,227 /2,025 /1,674 /1,197 /644 |
| `TM-score_0.9`–`0.3` | Max TM <0.9/0.8/0.7/0.6/0.5/0.4/0.3 | 926 /411 /171 /40 /8 /4 /3 |
| `OOD_Orphan` | No significant train+val hit | 422 |
| `OOD_IDR` | Disordered-residue ratio >0.1 | 39 |
| `OOD_FoldHoldout` | CATH+TED topology absent from training | 91 |
| `OOD_SuperfamilyHoldout` | CATH+TED superfamily absent from training | 311 |
| `OOD_NewEC_L3` /`L4` | All ECs absent from train+val | 20 /647 |
| `OOD_LongTail_EC_L3_le5` /`le10` | All L3 training frequencies <5 /<10, including zero | 118 /221 |
| `OOD_LongTail_EC_L4_le5` /`le10` | All L4 training frequencies <5 /<10, including zero | 1,763 /2,208 |
| `OOD_LongTail_pHbin_le50` /`le100` | Training bin count ≤50 /≤100 | 25 /58 |
| `InD` | All OOD flags False, NaN treated as False | 50 |

Full 34-column test header:

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

## 5. Quality checks

### 5.1 Leakage

| Check | Result |
|:---|---:|
| Train ∩Val | 0 |
| Train ∩Test | 0 |
| Val ∩Test | 0 |
| Duplicate unique_id within any split | 0 |

### 5.2 Completeness

| Check | Result |
|:---|---:|
| Final files in splits/ | 3: train/val/test |
| Structures in pdbs/ | 11,615, all CIF |
| Every split sample has a structure | ✅ |
| label stored as string | ✅ |
| No train/val OOD columns; full 34-column test | ✅ |
| Default=False computed OOD columns NaN | ✅ |
| OOD_Orphan ⊆seq_Redundancy_30, 422≤644 | ✅ |

### 5.3 Outliers

- OOD_ExtremeLong (211) and OOD_ExtremeShort (24) have Default=False; evaluate separately, without other OOD annotations (NaN).
- IDR positives: 39 (1.4% of Default), rare among mostly structured enzymes.
- NewEC_L3 positives: 20 (0.7% of Default), reflecting substantial train+val/test L3 overlap.
- LongTail_EC_L3_le10: 221 (8.1% of Default), showing L3 long tails. L4 covers 64.9%–81.3% and is less discriminative (§4.8).

---

## 6. Outputs and layout

```
enzyme_optimal_ph/
├── data/                       # Raw and intermediate data
│   ├── metadata.csv
│   ├── DATASET.zip / Dataset_mmcifs.zip
│   ├── dataset/                # Extracted PDB/PQR: DATASET/EnzyBase12k_pbs|pqrs
│   └── dataset_mmcifs/         # Extracted CIF: EnzyBase12k_mmcifs
├── output/                     # Reports, intermediates, backups
│   ├── metadata_analysis.json / METADATA_ANALYSIS_REPORT.md / figures/
│   ├── structure_verification.json
│   ├── sequence_structure_validation.csv (+ _summary.json)
│   ├── sequences.fa / mmseqs_cluster_95/ (+ _summary.json)
│   ├── split_summary.json / split_mapping.json
│   ├── ood_markers_summary.json
│   ├── newec_unified_summary.json / longtail_ec_unified_summary.json
│   ├── optimal_ph_ted_aug_acc_list.csv / optimal_ph_acc_cath_ted_labels.csv
│   ├── test_sequences_for_idr.fasta
│   ├── backup_before_*/        # Backups before split updates
│   └── scripts_archive/        # Archived one-off/obsolete scripts (§8)
├── pdbs/                       # 11,615 final CIF structures
├── splits/                     # Three final train/val/test CSVs
├── scripts/                    # Numbered scripts listed above
└── logs/                       # Download logs
```

| Output | Location | Description |
|:---|:---|:---|
| Final splits | `splits/optimal_ph_prediction_{train,val,test}.csv` | 7,794 /866 /2,955 rows; 34 test columns |
| Structures | `pdbs/*.cif` | 11,615 files named by unique_id |
| CATH+TED labels | `output/optimal_ph_acc_cath_ted_labels.csv` | Accession-level union, domains, and sources |
| Stage summaries | `output/*.json` | See layout |

---

## 7. Diagnostic tools

These scripts remain as diag_* outside the main workflow:

| Script | Purpose | Output |
|:---|:---|:---|
| `diag_analyze_metadata.py` | Metadata statistics and plots | output/metadata_analysis.json, METADATA_ANALYSIS_REPORT.md, figures/ |
| `diag_validate_sequence_structure.py` | edlib sequence–structure alignment (§3.1) | output/sequence_structure_validation.csv |
| `diag_cluster_sequences_mmseqs.py` | 95% clustering redundancy analysis (§3.2) | output/mmseqs_cluster_95*/ |

---

## 8. Archived scripts: output/scripts_archive/

| Script | Reason archived |
|:---|:---|
| `compute_ood_idr.py` | IDR@0.1 integrated into 04_compute_ood_markers.py, eliminating the second pass |
| `compute_additional_ood.py` | pH bins extracted into 05_add_phbin_longtail_ood.py; EC definitions replaced by canonical skill |
| `fix_labels_dtype.py` | One-off string-label conversion; current 03 writes strings directly |
| `restrict_ood_to_default.py` | Historical one-off False values for non-Default OOD replaced by project-wide NaN since 2026-08-27 (§4.10) |
| `update_ood_summary.py` | Companion summary updater archived with those repairs |

---

## 9. Historical notes and known deviations

1. **Schema:** original 03 output unique_id/file_type/aa_seq/labels migrated to unique_id/aa_seq/struct_file/label. The current script writes the unified schema directly. Audit 2026-08-28: all columns in a 320-row rerun match current splits exactly.
2. **Non-Default semantics:** restrict_ood_to_default.py historically set seq/TM/Orphan to False. Since 2026-08-27 use NaN (§4.10); archive that script. Since 2026-08-28, 04/05/06 and skills write NaN during generation; finalize --mode nan is fallback-only and currently no-op.
3. **EC columns:** obsolete OOD_NewFunction /OOD_LongTailFunction_L3_le5 /le10 no longer exist, replaced by OOD_NewEC_L3/L4 and four LongTail_EC columns (§4.8).
4. **IDR threshold:** **0.1**, versus 0.3 for most tasks. Formerly corrected by second-pass compute_ood_idr.py; now integrated into 04. See the toolkit idr_ood.py cross-task table.
5. **95% clustering not applied:** 557 redundant rows are diagnostic only; all 11,615 samples retained.
6. **Multiple EC strings:** split ec_id on `; ` for 782 UniProts before all EC statistics (§4.8).
7. **EC drift:** recomputation must use historical metadata.csv ec_id, not current UniProt annotations; see the EC skill source table.
8. **Putative TED:** prioritize experimental CATH, fill gaps with TED. Approximately 8% of accessions are not_found in TED and treated as unannotated; see the CATH+TED skill notes.
