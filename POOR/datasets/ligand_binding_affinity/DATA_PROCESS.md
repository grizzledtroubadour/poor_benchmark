# Ligand Binding Affinity Dataset Processing

> Task: ligand_binding_affinity (regression: pK = -log10(Kd/Ki))
> Source: free PDBbind+ v2020.R1 distribution (https://www.pdbbind-plus.org.cn/download)
> Location: `datasets/ligand_binding_affinity/`
>
> This document describes the complete four-stage, linear workflow for reproducing the current `splits/` and `pdbs/` from raw data.
> Statistics and column names were recomputed and checked against the current `splits/*.csv` headers in a read-only audit.

## Task and current scale

Predict binding affinity pK from protein sequence/structure and ligand SMILES/ECFP4. The target `label` is stored as a string with six decimal places.

| Split | File | Samples | Kd / Ki | Notes |
|:---|:---|---:|:---|:---|
| Train | `splits/ligand_binding_affinity_train.csv` | 8,366 | 4,638 / 3,728 | 6 columns |
| Validation | `splits/ligand_binding_affinity_val.csv` | 930 | 516 / 414 | 6 columns |
| Test | `splits/ligand_binding_affinity_test.csv` | 2,367 | 1,684 / 683 | 38 columns, including all OOD annotations |
| **Total** | — | **11,663** | — | 11,663 `{pdb_id}.pdb` files in `pdbs/` |

Overall ratio: approximately **7.2 : 0.8 : 2.0**; train:val ≈ **9 : 1**. Test `Default` is always True (the LBA/PPI convention; see §9).

Distribution: **PDBbind v2020.R1**, the free PDBbind+ Demo package. Structures were reprocessed with the PDBbind v2024 pipeline and have improved quality over the old PDBbind-CN v2020 distribution. The index is `INDEX_general_PL.2020R1.lst`, containing 19,037 protein–ligand complexes without precomputed pK; pK is calculated during parsing.

Scripts in execution order:

| No. | Script | Stage |
|:---|:---|:---|
| 01 | `01_download_pdbbind2020.sh` | ① Download and extract archives |
| 02 | `02_parse_index.py` | ① Parse index |
| 03 | `03_clean_metadata.py` | ② Clean labels and retain Kd/Ki |
| 04 | `04_parse_ligands.py` | ② Parse ligands: SMILES + ECFP4 |
| 05 | `05_extract_protein_sequences.py` | ③ Extract sequences by chain, separated by `\|` |
| 06 | `06_build_dataset_csv.py` | ③ Write sample-level dataset |
| 07 | `07_copy_whole_protein_pdbs.py` | ③ Write whole-protein structures |
| 08 | `08_cluster_time_split.py` | ④ Ligand clustering, within-cluster temporal splitting, and ratio balancing |
| 09 | `09_format_splits_to_skill.py` | ④ Unified schema and sequence-homology OOD |
| 10 | `10_compute_structure_idr_ood.py` | ④ Structural similarity and IDR OOD |
| 11 | `11_ted_augment_cath_ood.py` | ④ CATH+TED fold/superfamily holdouts |
| 12 | `12_add_value_bin_longtail_ood.py` | ④ Affinity-bin long-tail OOD |

---

## Stage 1: Obtain and parse source data

### 1.1 Download and extract archives

**Script:** `scripts/01_download_pdbbind2020.sh`<br>
**Input:** manually downloaded `data/index.tar.gz` and `data/P-L.tar.gz`<br>
**Output:** `data/index/INDEX_general_PL.2020R1.lst` and `data/P-L/{1981-2000,2001-2010,2011-2019}/{pdb_id}/`

```bash
cd ligand_binding_affinity
bash scripts/01_download_pdbbind2020.sh
```

Each extracted sample directory contains `{pdb_id}_protein.pdb`, `{pdb_id}_ligand.mol2`, `{pdb_id}_ligand.sdf`, and `{pdb_id}_pocket.pdb`.

### 1.2 Parse index

**Script:** `scripts/02_parse_index.py`<br>
**Input:** `data/index/INDEX_general_PL.2020R1.lst`<br>
**Output:** `output/pdbbind_v2020_metadata.csv`: 19,037 records, zero parsing failures; affinity types Kd/Ki/IC50. Compute pK = -log10(affinity_value) here.

---

## Stage 2: Aggregate and compute labels

### 2.1 Clean labels: retain Kd/Ki

**Script:** `scripts/03_clean_metadata.py`<br>
**Input:** `output/pdbbind_v2020_metadata.csv`<br>
**Output:** `output/pdbbind_v2020_metadata_cleaned.csv` and `_stats.json`

Exclude records with:

1. `affinity_type == IC50`.
2. Non-exact `affinity_operator` values such as `<`, `<=`, or `~`.
3. `[Incomplete ligand]`.
4. `[Covalent complex]`.
5. `[Different protein in assay]`.
6. `[Different ligand in assay]`.

> `[Uncommon element]` indicates uncommon ligand elements, such as boron or metals. These records are retained: 96 samples. Further handling may be needed for element-sensitive models.

| Item | Value |
|------|------|
| Before cleaning | 19,037 |
| After cleaning | **11,664** (Kd 6,839 / Ki 4,825) |
| Retention | 61.3% |

### 2.2 Parse ligands: SMILES + ECFP4

**Script:** `scripts/04_parse_ligands.py`<br>
**Input:** `output/pdbbind_v2020_metadata_cleaned.csv`<br>
**Output:** `output/pdbbind_v2020_metadata_ligand.csv`, `_stats.json`, and `ligand_parse_failures.csv`

- Prefer `{pdb_id}_ligand.sdf`; fall back to `{pdb_id}_ligand.mol2`.
- RDKit generates `ligand_smiles` (isomeric canonical SMILES), `ligand_ecfp4` (ECFP4/Morgan fingerprint, radius=2, nBits=2048, stored as ON-bit indices), and `ligand_num_heavy_atoms`.

| Item | Value |
|------|------|
| Input | 11,664 |
| Successfully parsed | **11,663** (`2pll` failed and was removed) |
| Parsing failures | 1 |

---

## Stage 3: Redundancy reduction and refinement

> Three sequence/ligand deduplication approaches were evaluated (95% sequence identity, 95% ligand similarity, and 95% SMILES similarity) but **not adopted**. All 11,663 cleaned samples remain. Scripts are archived (§8); `output/pdbbind_v2020_metadata_dedup95*.csv` files are historical artifacts only.
> Refinement comprises sequence–structure alignment and whole-protein structure export in the following three steps.

### 3.1 Extract protein sequences by chain

**Script:** `scripts/05_extract_protein_sequences.py`<br>
**Input:** `output/pdbbind_v2020_metadata_ligand.csv` and `data/P-L/**/{pdb_id}_protein.pdb`<br>
**Output:** `output/pdbbind_v2020_metadata_seq.csv`, `_stats.json`, and `sequence_extraction_failures.csv`

- Extract standard and common modified amino acids **by chain** from source protein PDBs.
- Join chains with `|` (no separator for single chains). Also generate `protein_length` (sum of chain lengths, excluding `|`), `protein_chains` (`;`-separated), and `num_protein_chains`.
- Deduplicate by (chain, res_seq, res_name), retaining protein residues only. Extraction from source files therefore agrees with extraction from `pdbs/` after step 07. Historically, archived `17_reextract_sequences_from_pdbs.py` corrected chain handling afterward; that logic is now incorporated here.
- Successful: **11,663**; failures: 0. Multichain: **6,393**; single-chain: 5,270. Mean length: 532.8 aa; range: 24–5,280 aa.

### 3.2 Write sample-level dataset

**Script:** `scripts/06_build_dataset_csv.py`<br>
**Input:** `output/pdbbind_v2020_metadata_seq.csv`<br>
**Output:** `output/pdbbind_v2020_dataset.csv` and `output/protein_chain_extraction_stats.json`

- Keep the 11,663 successfully extracted samples and add `protein_pdb_file` (`pdbs/{pdb_id}.pdb`).
- Replace the chain-level route in archived `08_extract_protein_chains.py`; individual-chain PDBs are no longer written. Retain `chain_extraction_success` / `chain_extraction_error` for downstream compatibility, matching sequence-extraction success.

Core fields in `pdbbind_v2020_dataset.csv`:

| Field | Description |
|------|------|
| `pdb_id` | Four-character PDB code; later used as `unique_id` |
| `resolution` / `release_year` | Resolution (Å) / release year for temporal splitting |
| `affinity_type` / `affinity_value` / `pK` | Kd or Ki / affinity in M / regression target |
| `ligand_smiles` / `ligand_ecfp4` / `ligand_num_heavy_atoms` | Ligand features |
| `protein_sequence` / `protein_length` | Protein sequence (chains separated by `\|`) / total amino acids |
| `protein_chains` / `num_protein_chains` | Chain IDs / number of chains |
| `protein_pdb_file` | Whole-protein PDB path: `pdbs/{pdb_id}.pdb` |
| `protein_file` / `ligand_sdf_file` / `pocket_file` | Relative source-file paths |

### 3.3 Write whole-protein structures

**Script:** `scripts/07_copy_whole_protein_pdbs.py`<br>
**Input:** `output/pdbbind_v2020_dataset.csv` and `data/P-L/**/{pdb_id}_protein.pdb`<br>
**Output:** 11,663 `pdbs/{pdb_id}.pdb` files and `output/whole_protein_pdb_stats.json`

- Filenames use **PDB ID only**. Filter non-protein `HETATM` residues (water, ligands, cofactors, etc.); keep standard and common modified amino acids.
- Stop at the first `ENDMDL` to avoid duplicates from NMR models.

---

## Stage 4: Split data and annotate OOD

### 4.1 Ligand clustering, within-cluster temporal splitting, and balancing

**Script:** `scripts/08_cluster_time_split.py`<br>
**Input:** `output/pdbbind_v2020_dataset.csv`<br>
**Output:** intermediate `splits/{train,val,test}.csv` (consumed and removed by 09) and `output/cluster_split_stats.json`

1. **Clique partition:** greedily cluster ECFP4 fingerprints, requiring pairwise Tanimoto ≥0.6 within each cluster. Isolated samples form singleton clusters: 5,394 clusters, 3,205 singletons, maximum size 271, mean size 2.16.
2. **Within-cluster temporal split:** sort each cluster by `release_year`, breaking ties by `pdb_id`, then split train/val/test 7:1:2.
3. **Similarity fallback:** iteratively move val/test samples with maximum training Tanimoto <0.6 into train: **2,866** samples moved.
4. **Ratio balancing:** incorporates former `12_balance_val_from_test.py` / `13_adjust_val_to_train_9to1.py`. After fallback, move the earliest third of test (1,183 samples) into val, then the earliest 442 val samples into train, giving train:val ≈9:1 (actual 8.996).

Final split: train 8,366 (71.7%), val 930 (8.0%), test 2,367 (20.3%).

### 4.2 Unified schema and sequence-homology OOD

**Script:** `scripts/09_format_splits_to_skill.py`<br>
**Input:** `output/pdbbind_v2020_dataset.csv` and `splits/{train,val,test}.csv`<br>
**Output:** `splits/ligand_binding_affinity_{train,val,test}.csv` (intermediate files removed), `output/{test,train_val}_for_mmseqs.fasta`, and `output/mmseqs_test_vs_trainval.tsv`

- Unified schema: `unique_id, aa_seq, struct_file, ligand_smiles, ligand_ecfp4, label`. Store `label` (pK) as a string with six decimal places and `struct_file` as a filename (`{pdb_id}.pdb`).
- Add test columns: `Default` (all True), `OOD_ExtremeShort` / `OOD_ExtremeLong` (total protein length <60 / >1,000; flag without removing), `seq_Redundancy_90..30`, and `OOD_Orphan`.
- Compute maximum test-to-train+val sequence identity with mmseqs2. Remove `|` in FASTA because mmseqs2 accepts only standard amino acids. `--mmseqs-cache` replays existing convertalis results.
- **Superseded:** this script contains the historical first-pass mmseqs2 implementation. The canonical implementation/replay tool is `.skills/homology-ood-annotation/scripts/seq_homology_ood.py` (easy-search / `--m8-cache`).

| Field | Description |
|------|------|
| `unique_id` | Unique sample ID, equal to `pdb_id` |
| `aa_seq` | Protein amino-acid sequence; multichain sequences separated by `\|` |
| `struct_file` | Structure filename only: `{pdb_id}.pdb`, located in `pdbs/` |
| `ligand_smiles` / `ligand_ecfp4` | Canonical ligand SMILES / ECFP4 ON-bit indices |
| `label` | Affinity pK, **stored as a string with six decimal places** |

### 4.3 Structural similarity and IDR OOD

**Script:** `scripts/10_compute_structure_idr_ood.py`<br>
**Input:** `splits/ligand_binding_affinity_{train,val,test}.csv` and `pdbs/`<br>
**Output:** update test.csv in place with `TM-score_0.9..0.3`, `idr_ratio`, and `OOD_IDR`; update `output/cluster_split_stats.json`

- **Foldseek-Multimer complex alignment (since 2026-08-29):** `easy-multimersearch`, using 2,367 test complexes as queries and 9,296 train+val structures as targets. Take the maximum symmetric mean of each complex pair's `complexqtmscore`/`complexttmscore`. Foldseek createdb splits multichain PDBs (`<file>_<chain>`); normalize hits back to the complex stem (bare name for single-chain files). Complex-level TM matches the complex–ligand affinity task. Previously: single-chain `search -a` with symmetric `qtmscore/ttmscore` means (different from other tasks' `alntmscore`). See `analysis/output/foldseek_easysearch_audit/lba/multimer/three_way_compare.txt` for old symmetric mean / single-chain alntm max / complex multimer comparisons. Old test backup: `output/backup_before_multimer_tm_unify/`; per-stem cache: `output/foldseek_multimer_max_tm_per_stem.csv`.
- metapredict estimates intrinsic disorder; `idr_ratio >0.1` sets `OOD_IDR`. Unified tool: `.skills/ood-annotation-toolkit/scripts/idr_ood.py`.

### 4.4 CATH fold/superfamily holdouts: final CATH+TED union

**Script:** `scripts/11_ted_augment_cath_ood.py`<br>
**Input:** `splits/ligand_binding_affinity_{train,val,test}.csv`, project-level `data/sifts/` and `data/CATHv44/`<br>
**Output:** `output/ligand_ted_aug_chain_list.csv`, `ligand_chain_cath_ted_labels.csv`, `ligand_chain_domains.csv`, `ligand_pdb_cath_ted_labels.csv`, and `ligand_pdb_cath_ted_holdout.csv`; update two test columns in place

Use the **CATH+TED union**: prioritize CATH v4.4 experimental annotations; fill gaps with putative TED annotations. Unified workflow: `.skills/cath-ted-domain-annotation`. Its `scripts/cath_ted_holdout.py` is the sole holdout implementation (training-only reference, any semantics, unannotated →False).

1. **Chain list:** 11,663 PDBs × chains present in either SIFTS file or `cath-domain-list.txt`. SIFTS-only seeding misses **944 CATH-annotated chains** across 342 PDBs, including 278 chains in 87 test PDBs.
2. Run `cath_ted_annotate.py` with `--offline`. TED requires 50% domain overlap; H level provides C.A.T.H, T level only C.A.T. API cache: project-level `data/ted/ted_api_cache.json`.
3. Union chain labels by PDB. **Symmetric augmentation:** use the CATH+TED union for training references as well, avoiding falsely unseen TED-only test labels.
4. Run `cath_ted_holdout.py`, updating only the two test holdout columns.

**Scale:** 27,352 chains, 11,663 PDBs. Sources: both 18,990 / cath-only 3,219 / ted-only 2,115 / none 3,028. TED cache hits: 3,811 accessions; zero new queries.

**Test results (n=2,367):** `OOD_FoldHoldout` **38**, `OOD_SuperfamilyHoldout` **73**, unannotated samples (silently False) 149. Backup: `output/backup_before_ted_aug/`.

> Superseded: archived CATH-only `16_cath_ec_ood.py` left 27.5% of test PDBs unannotated and yielded only 24/45 fold/superfamily holdouts. Its EC section is also superseded by §4.5.

### 4.5 EC function OOD: `OOD_NewEC_*` / `OOD_LongTail_EC_*`

> Canonical definitions, hierarchy parsing, and drivers: `.skills/ec-function-ood-annotation` (`scripts/add_unified_newec.py`, `scripts/add_unified_longtail_ec.py`). EC source: `data/sifts/sifts_chain_ec.tsv.gz`, unioned across all chains per PDB.

- `OOD_NewEC_L4` / `OOD_NewEC_L3`: reference = EC labels in **train+val**; **all-not-in** semantics. True if annotated and all ECs at the given level are absent from the reference; unannotated →False.
- `OOD_LongTail_EC_{L3,L4}_le{5,10}`: **training-only** per-sample frequencies, **including zero**. All-not-in: True if every EC at that level has training frequency <5 / <10; unannotated →False.
- Legacy `OOD_NewFunction` (any semantics, 54 samples) and `OOD_LongTailFunction` (181 samples) were **removed and replaced**. Backups: `output/backup_before_newec_unify/`, `output/backup_before_longtail_allnotin/`.

### 4.6 Affinity-bin long-tail OOD

**Script:** `scripts/12_add_value_bin_longtail_ood.py`<br>
**Input:** `splits/ligand_binding_affinity_{train,test}.csv`<br>
**Output:** update test.csv with `OOD_LongTail_AffBin_le50` / `le100`; write `output/ligand_longtail_affbin_summary.json`

- Bin pK labels with **width 0.5**, range [-0.5, 16.5); count **training** samples per bin.
- True if a test sample's bin has training count ≤50 / ≤100. All test samples have Default=True.
- Training pK (n=8,366) is approximately normal (mean=6.26, std=1.97), with symmetric long tails. Positive test labels span (0.66, 12.30), representing extreme affinities at both ends.
- Labels are not clipped (unlike kcat's upstream CatPred-DB clipping to [-6, 6]).
- Generic tool: `.skills/ood-annotation-toolkit/scripts/value_bin_longtail_ood.py`.

### 4.7 Finalize InD and validate

> Finalizer: `.skills/ood-annotation-toolkit/scripts/finalize_ood_columns.py`; validator: `scripts/check_splits.py` in the same skill.

- `InD`: True when all `OOD_*`, `seq_Redundancy_*`, and `TM-score_*` flags are False (NaN treated as False). Recompute with `finalize_ood_columns.py --mode ind` after all OOD columns are ready. Recorded count: **1,091 / 2,367 (46.1%)**; backup: `output/backup_before_ind/`.
- NaN convention (`--mode nan`): all test `Default` values are True (LBA/PPI convention), so there are no Default=False rows requiring NaN conversion.
- `check_splits.py` passes leakage and subset-relation checks (§5).

### 4.8 Test OOD statistics: final recorded splits

| OOD dimension | Flag | True count |
|----------|--------|----------|
| Default test set | `Default` | 2,367 |
| Extreme short / long | `OOD_ExtremeShort` / `OOD_ExtremeLong` | 15 / 303 |
| Sequence identity <90%–30% | `seq_Redundancy_90`–`seq_Redundancy_30` | 575 / 498 / 454 / 409 / 368 / 294 / 241 |
| No homologous hit | `OOD_Orphan` | 181 |
| Structural novelty (complex TM <0.9–0.3) | `TM-score_0.9`–`TM-score_0.3` | 767 / 662 / 447 / 306 / 201 / 135 / 79 |
| High disorder (idr >0.1) | `OOD_IDR` | 93 |
| CATH+TED superfamily / fold holdout | `OOD_SuperfamilyHoldout` / `OOD_FoldHoldout` | 73 / 38 |
| New EC L4 / L3 | `OOD_NewEC_L4` / `OOD_NewEC_L3` | 46 / 8 |
| Long-tail EC L3 (le5 / le10) | `OOD_LongTail_EC_L3_le5` / `le10` | 24 / 37 |
| Long-tail EC L4 (le5 / le10) | `OOD_LongTail_EC_L4_le5` / `le10` | 189 / 344 |
| Long-tail affinity bins (le50 / le100) | `OOD_LongTail_AffBin_le50` / `le100` | 19 / 65 |
| Pure InD samples | `InD` | 1,247 (52.7%) |

Current test header (38 columns):

```
unique_id, aa_seq, struct_file, ligand_smiles, ligand_ecfp4, label,
Default, OOD_ExtremeShort, OOD_ExtremeLong,
seq_Redundancy_90, ..., seq_Redundancy_30,
TM-score_0.9, ..., TM-score_0.3,
idr_ratio, OOD_IDR, OOD_Orphan,
OOD_SuperfamilyHoldout, OOD_FoldHoldout,
OOD_NewEC_L4, OOD_LongTail_AffBin_le50, OOD_LongTail_AffBin_le100,
OOD_LongTail_EC_L3_le5, OOD_LongTail_EC_L3_le10,
OOD_LongTail_EC_L4_le5, OOD_LongTail_EC_L4_le10,
OOD_NewEC_L3, InD
```

### 4.9 Ligand chemical similarity: Tanimoto

Maximum Tanimoto similarity of each val/test sample to all training samples:

| Split | Mean | Median | Min | Max |
|------|------|--------|--------|--------|
| Val | 0.858 | 0.873 | 0.600 | 1.000 |
| Test | 0.846 | 0.835 | 0.600 | 1.000 |

| Similarity interval | Val | Test |
|------------|-----|------|
| 0.6–0.7 | 179 | 495 |
| 0.7–0.8 | 198 | 551 |
| 0.8–0.9 | 103 | 290 |
| ≥0.9 | 450 | 1,031 |

> Every val/test sample has maximum training Tanimoto ≥0.6, the intended chemical-support floor.

---

## 5. Data quality checks

1. **Leakage:** pairwise Train/Val/Test `unique_id` intersections are empty. `check_splits.py` passes all checks, including `OOD_Orphan ⊆ seq_Redundancy_30`.
2. **Split design:** samples are split chronologically **within** each ligand-similarity cluster (pairwise Tanimoto ≥0.6). Clusters may span train/val/test and are not held out as groups. The val/test maximum-training Tanimoto floor of 0.6 deliberately preserves chemical support. This is not a single global temporal cutoff.
3. **Distribution consistency:** length and Kd/Ki proportions are consistent across subsets; test pK range (0.66, 12.30) lies within training support.
4. **Sequence–structure consistency:** `aa_seq` aligns sample by sample with `pdbs/`, using chain-wise extraction (§3.1). Historical re-extraction found zero mismatches.
5. **Outliers:** remove `2pll` after ligand parsing failure (§2.2). IC50, inequalities, and problematic annotations are excluded during cleaning (§2.1).

---

## 6. Output files

```
ligand_binding_affinity/
├── DATA_PROCESS.md              # This document
├── data/                        # Raw data
│   ├── index/
│   └── P-L/
├── pdbs/                        # 11,663 whole-protein PDBs: {pdb_id}.pdb
├── scripts/                     # Numbered processing scripts listed above
├── output/
│   ├── pdbbind_v2020_metadata.csv
│   ├── pdbbind_v2020_metadata_cleaned.csv
│   ├── pdbbind_v2020_metadata_ligand.csv
│   ├── pdbbind_v2020_metadata_seq.csv
│   ├── pdbbind_v2020_dataset.csv          # Sample-level dataset; Stage 4 input
│   ├── pdbbind_v2020_metadata_dedup95*.csv # Historical unused deduplication route
│   ├── ligand_parse_failures.csv
│   ├── whole_protein_pdb_stats.json
│   ├── cluster_split_stats.json
│   ├── {test,train_val}_for_mmseqs.fasta / mmseqs_test_vs_trainval.tsv
│   ├── ligand_ted_aug_chain_list.csv      # PDB × chain list with split membership
│   ├── ligand_chain_cath_ted_labels.csv   # Chain-level CATH+TED labels
│   ├── ligand_chain_domains.csv           # Domain details: 81,029 rows
│   ├── ligand_pdb_cath_ted_labels.csv / ligand_pdb_cath_ted_holdout.csv
│   ├── ligand_longtail_affbin_summary.json
│   ├── scripts_archive/                   # Archived scripts (§8)
│   └── backup_before_*/                   # Backups before split updates
└── splits/
    ├── ligand_binding_affinity_train.csv  # 8,366
    ├── ligand_binding_affinity_val.csv    # 930
    └── ligand_binding_affinity_test.csv   # 2,367, including all OOD flags
```

---

## 7. Diagnostics and quality-control tools

No standalone `diag_*` scripts are used. Quality control relies on:

- `.skills/ood-annotation-toolkit/scripts/check_splits.py`: leakage, column completeness, and subset relations.
- `output/cluster_split_stats.json`, generated by `08_cluster_time_split.py`: clustering, split statistics, and max_train_tanimoto distributions.
- Per-step `*_stats.json` and failure-list CSVs (§6).

---

## 8. Archived scripts: `output/scripts_archive/`

| Script | Reason archived |
|:---|:---|
| `05_dedup_sequences_95.py` | Evaluated but unused 95% sequence-identity deduplication |
| `06_dedup_by_ligand_95.py` | Evaluated but unused 95% ligand-similarity deduplication |
| `07_dedup_by_smiles_95.py` | Evaluated but unused 95% SMILES-similarity deduplication |
| `08_extract_protein_chains.py` | Chain-level structures replaced by whole proteins; `06_build_dataset_csv.py` now generates `pdbbind_v2020_dataset.csv` |
| `10_time_split_and_tanimoto.py` | Replaced by `08_cluster_time_split.py` |
| `12_balance_val_from_test.py` | Balancing incorporated into `08_cluster_time_split.py` |
| `13_adjust_val_to_train_9to1.py` | Train:val 9:1 adjustment incorporated into 08 |
| `16_cath_ec_ood.py` | CATH-only and legacy EC definitions replaced by `11_ted_augment_cath_ood.py` and the ec-function-ood-annotation skill |
| `17_reextract_sequences_from_pdbs.py` | Chain-wise `\|` sequence extraction incorporated into `05_extract_protein_sequences.py` |

---

## 9. Notes: history, non-script steps, and known deviations

1. **LBA/PPI convention: Default always True.** Train/val retain extreme lengths: training has 67 samples <60 aa and 1,008 with total length >1,000 aa. Test does not exclude extreme lengths; all 2,367 test samples have `Default=True`. Extreme-length flags annotate without removing; no Default=False rows require NaN conversion.
2. **Multichain complexes:** `aa_seq` separates chains with `|` (1,264 multichain test samples, 4,591 train). Remove separators for downstream tools as needed. Length is the sum across chains.
3. **Long single chains near the boundary:** single-chain lengths >1,000 cluster around 1,014–1,016 aa. Of 25 long single-chain training samples, 16 fall in this interval; all three validation samples are 1,015–1,016. Test has only four, at 1,029–1,195. Other long samples exceed the threshold only in combined multichain length; no single chain exceeds 2,000 aa.
4. **≤3-aa fragments are source-deposition features (audit 2026-08-28):** seven train, two val, and two test samples contain such fragments. Individual checks against PDBbind sources in `data/P-L/` and the RCSB mirror confirm they are not extraction artifacts: peptide inhibitor chains (1eb1 retains only PRO 2), UNK placeholders/short stubs (2j9n GLN+UNK; 3e4a/6byz/4hva/3t6y/3t64/3t70 poly-Ala/His), collapsed alternative conformations (5jjm Ser-Asp ×12), and modified HETATM residues (4aze SEP→S; note the overflow record `HETATM17267` without spaces). Chain-wise ATOM+HETATM extraction and modified-residue mapping faithfully retain these chains; no repair is needed.
5. **Floating-point reserialization drift (identified 2026-08-28):** repeated intermediate reads/writes introduce a 1-ulp difference in `3o56`'s `affinity_value`: current `7.939999999999993e-10`, regenerated by 05+06 `7.939999999999995e-10`, original `metadata.csv` `7.939999999999999e-10`. The pK string remains identical, with no label or split effect. Retain this historical deviation without correction.
6. **Reproduction audit (2026-08-28):** regenerate `pdbbind_v2020_dataset.csv` from `pdbbind_v2020_metadata_ligand.csv` using 05+06. All 11,663 rows ×36 columns agree except the single 3o56 floating-point deviation above. A full rerun of consolidated `08_cluster_time_split.py` (5,394 clusters; moves 2,866 +1,183 +442) reproduces current train/val/test membership **exactly**.
7. **Historical changes:** CATH-only fold/superfamily counts 24/45 become 38/73 after TED augmentation. NewEC/LongTailEC use unified all-not-in semantics, replacing `OOD_NewFunction`/`OOD_LongTailFunction`; affinity-bin flags were added later; the toolkit recomputes InD. Backups precede each update in `output/backup_before_*/`.
