# Enzyme Kinetics Prediction (kcat): Dataset Processing

> Complete linear reproduction workflow for `datasets/enzyme_kinetics_prediction/`, from raw CatPred-DB data to final splits. Statistics are based on recomputation from the current `splits/*.csv` files.

## Task and current scale

The **CatPred-DB kcat** regression task uses `label = log10(kcat)`, AlphaFold v6 predicted structures, and UniProt accession keys. An accession may have multiple rows for different substrates.

- Zenodo: `https://zenodo.org/records/14775076` (approximately 1.1 GB; processed datasets, PDB records, and raw BRENDA/SABIO CSVs).
- GitHub: `https://github.com/maranasgroup/CatPred-DB` (approximately 397 MB).

Current measured split sizes:

| File | Rows | UniProts | Columns | Notes |
|------|-----:|-----------:|-----:|------|
| `splits/kcat_train.csv` | 15,479 | 4,888 | 6 | Training; six core columns |
| `splits/kcat_val.csv` | 2,150 | 662 | 6 | Validation; same schema as train |
| `splits/kcat_test.csv` | 5,002 | 1,486 | **37** | Six core plus 31 OOD/annotation columns; 4,412 Default=True rows and 590 extreme-length rows |

Four processing stages:

- **Stage 1:** obtain and parse sources: CatPred-DB, AlphaFold structures, PDB filtering, and Q9Y233 truncation.
- **Stage 2:** aggregate and compute labels: reactant filtering and ECFP4 fingerprints.
- **Stage 3:** redundancy reduction and refinement: no sequence-clustering deduplication for this task; see the refinement measures below.
- **Stage 4:** split and annotate OOD: length-based splits → standard columns → all OOD dimensions → NaN/InD finalization.

## Directory layout

```
enzyme_kinetics_prediction/
├── data/                                   # Task sources, caches, and processing logs (*.log)
│   ├── catpred-db-github/                  # GitHub repository contents
│   ├── zenodo/                             # Extracted Zenodo files
│   ├── kcat_filtered_with_pdb.csv          # After PDB filtering: 22,736 rows
│   └── kcat_with_reactant_ecfp4.csv        # Reactant filtering + ECFP4: 22,690 rows
├── scripts/                                # Processing scripts 01–08, in execution order
├── splits/                                 # Three final split files
├── pdbs/                                   # 7,058 decompressed AlphaFold v6 .pdb structures
├── output/                                 # Intermediate outputs, summaries, OOD subsets, backups
│   ├── scripts_archive/                    # Archived historical scripts
│   └── backup_*/                           # Backups before split updates
└── DATA_PROCESS.md                         # This document
```

> Processing logs are in `data/*.log`, not a separate `logs/` directory. `pdbs/` and `splits/` contain the final training/evaluation artifacts.

---

# Stage 1: Obtain and parse source data

> Goal: a structure-supported CatPred-DB kcat dataset: **22,736** rows and 7,057 UniProts.

## 1.1 Download data: manual step

- Clone GitHub into `data/catpred-db-github/`; log: `data/catpred-db-github_clone.log`.
- Extract Zenodo into `data/zenodo/`.

## 1.2 Obtain AlphaFold structures

- **Script:** `scripts/01_download_alphafold_missing.py`.
- **Input:** `data/kcat_missing_alphafold_uniprots.txt`.
- **Output:** project-level `data/alphafold_structures/AF-*-F1-model_v6.pdb.gz`; log: `data/alphafold_download_progress.log`.
- **Details:** kcat contains **7,232** unique UniProts. The shared structure library initially covered 4,522. Resumable AlphaFold DB API downloads with eight threads increase coverage to **7,057** (97.6%).

## 1.3 Filter and prepare PDBs

- **Script:** `scripts/02_process_kcat_with_pdbs.py`.
- **Input:** `data/catpred-db-github/datasets/processed/kcat_max_wt_singleSeqs_wpdbs.csv` and project-level `data/alphafold_structures/`.
- **Output:** `data/kcat_filtered_with_pdb.csv` (22,736 rows), `pdbs/*.pdb.gz`, `data/kcat_uniprots_removed_no_pdb.txt`, `data/kcat_pdb_sequence_consistency.csv`, and `data/kcat_substrate_stats.json`.
- **Details:**
  - Retain **22,736** rows / **7,057** UniProts with local structures; exclude 461 UniProts without structures.
  - Exact sequence agreement: 7,045 / 7,057 (99.83%).
  - Valid `reactant_smiles` coverage: 99.8%; 12,114 unique `reaction_smiles`. Major EC classes: Hydrolases (3), Oxidoreductases (1), Transferases (2), Lyases (4). Zero exact enzyme + complete-reactant-SMILES duplicate rows.

- **Script:** `scripts/03_decompress_pdbs_and_truncate_q9y233.py`.
- **Input:** `pdbs/*.pdb.gz` and `data/kcat_filtered_with_pdb.csv`.
- **Output:** decompressed `pdbs/*.pdb`, truncated Q9Y233 structure, and `data/q9y233_truncation_mapping.json`.
- **Details:** decompress all `.pdb.gz` files. Q9Y233's AlphaFold model differs from the dataset sequence; locate the dataset sequence with a sliding window and retain its **779** residues. The truncated structure is `AF-Q9Y233-F1-model_v6.pdb`; retain the full-length backup as `AF-Q9Y233-F1-model_v6_full.pdb`. Total: 7,058 files =7,057 UniProts +one full-length backup.

---

# Stage 2: Aggregate and compute labels

> Goal: **22,690** clean rows with reactant SMILES and ECFP4 fingerprints.

## 2.1 Filter reactants and compute ECFP4

- **Script:** `scripts/04_compute_reactant_ecfp4.py`.
- **Input:** `data/kcat_filtered_with_pdb.csv`.
- **Output:** `data/kcat_with_reactant_ecfp4.csv` (22,690 rows, `ecfp4_bit_0`–`ecfp4_bit_2047` and the binary-string `ecfp4` column), `output/kcat_ecfp4_summary.json`; log: `data/compute_reactant_ecfp4.log`.
- Use only `reactant_smiles`; exclude 46 rows without reactants, leaving **22,690** rows /7,051 UniProts.
- Compute RDKit Morgan fingerprints (radius=2, 2048 bits) for each reactant and combine molecules with bitwise OR. Mean active bits per combined fingerprint: **57.6**.
- Historical Step 3 updated archived `kcat_random_*` splits; this dead code has been removed.

## 2.2 Targets and length distribution: 22,690 rows

- `label = log10_value` (log10(kcat)), clipped/censored upstream by CatPred-DB to **[-6, 6]**. Raw BRENDA range: 2.24e-11 to 9.3e8 s⁻¹.
- Mean **0.89**, standard deviation **1.65**, minimum -6.00, maximum 6.00. See §4.6 for censoring effects on long-tail bins.

| Length | Count | Fraction |
|----------|-----:|-----:|
| <60 | 60 | 0.3% |
| 60–1000 | 22,041 | 97.1% |
| >1000 | 589 | 2.6% |

---

# Stage 3: Redundancy reduction and refinement

> **No sequence-clustering deduplication:** multiple substrates per UniProt are part of the task. Stage 4 uses **UniProt-level splitting** to prevent sample leakage. Refinement measures are distributed across the workflow:

| Measure | Location | Count |
|----------|------|------|
| Exclude missing structures | §1.3 | 461 UniProts |
| Exclude missing reactant SMILES | §2.1 | 46 rows |
| Move entire normal-length val/test UniProts with reactant similarity <0.6 into train | §4.1 | 69 UniProts |
| Remove entire extreme-length UniProts with training reactant similarity <0.6 | §4.1 | 21 rows /5 UniProts |
| Remove single chains >2000 aa, a project-wide rule without a task-specific script | §4.7 | 38 rows /10 UniProts |
| Schema migration: `pdb_file`→`struct_file`, `labels`→`label`, remove `file_type`; project-wide step without a task script | After initial OOD pass | — |

---

# Stage 4: Split data and annotate OOD

> Goal: UniProt-disjoint 7:1:2 splits; extreme-length samples assigned to test as OOD candidates; complete all OOD dimensions using canonical definitions.
>
> **Canonical implementations:** use the corresponding skill scripts for recomputation; do not introduce parallel definitions.
> - `.skills/homology-ood-annotation`: `seq_Redundancy_*` / `OOD_Orphan` (mmseqs2 easy-search) / `TM-score_*` (Foldseek; supports `--m8-cache` replay).
> - `.skills/ood-annotation-toolkit`: `idr_ratio` / `OOD_IDR` (metapredict), value-bin long-tail flags, NaN/InD finalization, and split validation.
> - `.skills/cath-ted-domain-annotation`: `OOD_FoldHoldout` / `OOD_SuperfamilyHoldout`, invoked through `scripts/07_ted_augment_cath_ood.py`.
> - `.skills/ec-function-ood-annotation`: unified all-not-in `OOD_NewEC_*` / `OOD_LongTail_EC_*`.

## 4.1 Length-based OOD train/validation/test splits

- **Script:** `scripts/05_create_length_ood_split.py`.
- **Input:** `data/kcat_with_reactant_ecfp4.csv`.
- **Output:** `splits/kcat_{train,val,test}.csv`, `output/kcat_split_summary.json`; log: `data/create_length_ood_split.log`.
- **Strategy:**
  1. Separate 649 extreme-length rows (<60 or >1000 aa).
  2. Randomly split normal-length **UniProts** 7:1:2 (`seed=42`), keeping each protein within one split.
  3. Move entire normal-length val/test UniProt groups whose maximum reactant ECFP4 Tanimoto to train is <**0.6** into train: 26 val +43 test =**69** UniProts.
  4. Remove extreme-length UniProt groups with maximum similarity <0.6 to expanded train: **5** UniProts /**21** rows. Assign remaining extreme-length rows to test as OOD candidates.

**Script-time split results:**

| Split | Rows | UniProts | Length range | Extreme-length rows |
|-------|-----:|-----------:|----------|-----------:|
| Train | 15,479 | 4,888 | 60–1000 | 0 |
| Val | 2,150 | 662 | 60–1000 | 0 |
| Test | 5,040 | 1,496 | Mixed | 628 |
| Removed | 21 | 5 | Extreme | — |

> §4.7 later removes 38 rows >2000 aa, leaving **5,002** test rows /**1,486** UniProts, including 590 extreme-length rows.

## 4.2 Standard POOD and length-OOD columns

- **Script:** `scripts/06_create_pood_compatible_splits.py`. Historical first-pass pipeline; **only standard columns and length OOD remain applicable**. Its inline mmseqs2/Foldseek/metapredict sections are superseded (§4.3–4.4).
- **Input/output:** update `splits/kcat_{train,val,test}.csv` in place.
- Add historical standard columns (`unique_id`, `file_type`, `aa_seq`, `labels`), later migrated project-wide to current `struct_file`/`label` with `file_type` removed. Compute `Default` / `OOD_ExtremeLength` / `OOD_ExtremeLong` / `OOD_ExtremeShort`. Historical counts in `output/kcat_pood_summary.json` are for reference only.
- **Unique IDs:** initially `AF-<acc>-F1-model_v6`, matching protein structure filenames. Because kcat rows are enzyme–reactant pairs, append `-<md5(reactant_smiles)[:8]>` to repeated protein IDs, ensuring uniqueness within a task/split. `(uniprot, reactant_smiles)` pairs have been verified unique; unrepeated IDs remain unchanged. Before the 2026-09 fix, protein-only IDs duplicated 10,591 train /1,488 val /3,516 test rows. The source `add_standard_columns` now fixes this; backup: `output/backup_before_unique_id_fix/`.

**Final schema:** train/val have six columns: `unique_id, aa_seq, struct_file, ligand_smiles, ligand_ecfp4, label`. Test adds 31 OOD/annotation columns. Store `label` as a string, reactant `ligand_smiles` with molecules separated by `.`, and `ligand_ecfp4` as a string of active bit indices.

## 4.3 Sequence and structural homology OOD

- **Canonical recomputation scripts:**
  - `.skills/homology-ood-annotation/scripts/seq_homology_ood.py`: mmseqs2 easy-search, query=test, target=train+val. `seq_Redundancy_<T>` means maximum identity <T%, T=90…30. No significant hit or best e-value >1e-3 sets `OOD_Orphan=True`.
  - `.skills/homology-ood-annotation/scripts/struct_homology_ood.py`: maximum Foldseek alntmscore; `TM-score_<T>` means max TM <T, T=0.9…0.3. Supports `--m8-cache` replay.
- Deduplicate sequence searches by `aa_seq` and structure searches by UniProt, then map results back to all substrate rows for each protein.
- **Superseded:** inline mmseqs2/Foldseek in `06_create_pood_compatible_splits.py`; always use the skill scripts above for recomputation.
- **Unification audit, 2026-08-30:** rerun sequence easy-search (`-s 7.5`), cache `analysis/output/seq_ood_unify_assessment/kcat/mmseqs_easysearch.m8`, backup `output/backup_before_seq_easysearch_unify/`. Only 34/4,412 Default samples change at boundaries: 83 seq50/40/30 and Orphan cells, all True→False. InD remains 366; old annotations and easy-search definitions are substantially equivalent.

## 4.4 Intrinsic disorder OOD

- **Canonical script:** `.skills/ood-annotation-toolkit/scripts/idr_ood.py`, metapredict v3.0.2. `idr_ratio` is the fraction of residues with disorder score >0.5, a continuous annotation available for all test rows. `OOD_IDR` requires continuous disordered regions to occupy >0.3 of the sequence.
- **Superseded:** inline metapredict in `06_create_pood_compatible_splits.py` and archived `reformat_splits_to_pood_format.py`.

## 4.5 CATH+TED domains: fold/superfamily holdouts

- **Wrapper:** `scripts/07_ted_augment_cath_ood.py`.
- **Input:** `splits/kcat_{train,val,test}.csv`; shared project-level `data/ted/ted_api_cache.json`.
- **Output:** `output/kcat_ted_aug_acc_list.csv` (7,036 accession keys), `output/kcat_acc_cath_ted_labels.csv`, `output/kcat_cath_ted_holdout.csv`, and `output/backup_before_ted_aug/`.
- **Workflow:**
  0. Generate the accession list from splits; formerly a historical one-off output without a retained script.
  1. `.skills/cath-ted-domain-annotation/scripts/cath_ted_annotate.py --chains output/kcat_ted_aug_acc_list.csv --out output/kcat_acc_cath_ted_labels.csv --offline`.
  2. `.skills/cath-ted-domain-annotation/scripts/cath_ted_holdout.py --labels ... --split-col splits --out output/kcat_cath_ted_holdout.csv`.
  3. Map holdouts by accession to each test row; Default=False receives NaN. The regex `^AF-(.+?)-F1-model_v6(?:-[0-9a-f]{8})?$` supports reactant-hash suffixes (§4.2).
- **Definition, UniProt mode:** labels = experimental CATH v4.4 (union of SIFTS reverse-mapped PDB chains) ∪ putative full-length TED, without domain filtering. H supplies C.A.T.H, T only C.A.T. Reference = **training** CATH+TED union, with symmetric train/val/test augmentation; any semantics; unannotated →False.
- **Coverage:** both 2,700 /ted-only 4,103 /cath-only 42 /none 191. Topology coverage rises from CATH-only 39.0% to CATH+TED **97.3%**. All TED requests hit cache: offline, zero new queries.
- **Current results, Default=True n=4,412:** fold holdout **43**, superfamily holdout **160**, unannotated →False **109** rows.

## 4.6 EC function OOD and value-bin long tails

**NewEC: `OOD_NewEC_L4` / `OOD_NewEC_L3`**

- Script: `.skills/ec-function-ood-annotation/scripts/add_unified_newec.py`; summary `output/newec_unified_summary.json`, backup `output/backup_before_newec_unify/`.
- Reference: EC labels in train+val. **All-not-in:** annotated samples are True only if every EC at that level is absent from train+val. Unannotated →False; hierarchy parsing allows partial `-`.
- EC source: `data/kcat_with_reactant_ecfp4.csv`, original uniprot→ec annotations. Split semicolon-separated multiple ECs. Reuse this historical source; current UniProt annotations introduce source drift.
- Default n=4,412, per-row results: L4 **853**, L3 **3**.
- Replaces deleted legacy `OOD_NewFunction` (L4 any semantics, 1,224 rows); old `add_additional_ood_labels.py` is archived.

**LongTailEC: four columns**

- Script: `.skills/ec-function-ood-annotation/scripts/add_unified_longtail_ec.py`; summary `output/longtail_ec_unified_summary.json`; backups `output/backup_before_longtail_unify/`, `output/backup_before_longtail_allnotin/`.
- Reference: **training-only row frequencies**. **All-not-in:** all ECs at the selected level must have training frequency <5 /<10, **including zero**. Unannotated →False; Default=False →NaN.
- Default n=4,412, per-row results: `OOD_LongTail_EC_L3_le5` **28**, `L3_le10` **75**, `L4_le5` **1,524**, `L4_le10` **2,163**.
- Deleted legacy `OOD_LongTailFunction` (217 rows) and transitional `OOD_LongTail_EC_L3` (85 rows).

**Regression value-bin long tails: `OOD_LongTail_KcatBin_le50` / `le100`**

- Script: `scripts/08_add_value_bin_longtail_ood.py`. Canonical reusable implementation: `.skills/ood-annotation-toolkit/scripts/value_bin_longtail_ood.py`.
- Bin width 0.5, edges [-6.5, 7.0). **Separate clipped boundary bins:** label==-6.0 and ==6.0 have training counts 3 /92, avoiding censored accumulation hiding the right tail. Default=True test rows are True if their training-bin count ≤50 /100; missing/unbinnable →False; Default=False →NaN.
- Default n=4,412: le50 **7**, le100 **133**. Summary `output/kcat_longtail_kcatbin_summary.json`, backup `output/backup_before_kcatbin_clipfix/`.
- Training labels (n=15,479) are approximately normal: mean=0.88, std=1.67, skew=0.11. ≤50 bins cover (-6.5,-3.5]; ≤100 additionally covers (-4.0,-3.0] and (5.0,5.5], giving symmetric tails.

## 4.7 Remove very long chains and finalize NaN/InD

- **Length >2000, manual project-wide step:** remove **38** test rows /**10** UniProts. Test 5,040→**5,002**; extreme rows 628→590. Backup: `output/len2000_removal_backup/kcat_test.before_len2000.csv`. Current test length range: 16–1928.
- **NaN convention:** on Default=False rows, all `OOD_*`, `seq_Redundancy_*`, and `TM-score_*` columns except `OOD_ExtremeShort/Long` are empty/NaN, meaning not computed. Since 2026-08-28, generating scripts (`06` and skills) write NaN at the appropriate pipeline step. `finalize_ood_columns.py --mode nan` is a fallback check only (currently a no-op); historical backup `output/backup_before_nondefault_nan/`. Continuous `idr_ratio` is unaffected.
- **InD recomputation:** `.skills/ood-annotation-toolkit/scripts/finalize_ood_columns.py --mode ind`; True if all `OOD_*`, `seq_Redundancy_*`, `TM-score_*` flags are False, treating NaN as False. Backup `output/backup_before_ind/`. Current InD: **366** (7.3% of test).

## 4.8 Validate splits

Read-only `.skills/ood-annotation-toolkit/scripts/check_splits.py` checks unique IDs within each split; pairwise train/val/test `unique_id` and UniProt intersections are **zero**. Train/val have six core columns; test has 37. Non-length OOD columns are empty on Default=False rows.

```bash
python ../../.skills/ood-annotation-toolkit/scripts/check_splits.py \
    --train splits/kcat_train.csv --val splits/kcat_val.csv --test splits/kcat_test.csv
```

## 4.9 Current split OOD statistics

Test: **5,002** rows, **4,412** Default=True and **590** Default=False extreme-length rows. Counts are per row; all non-length OOD True values fall within Default=True.

| OOD column | True count | Definition |
|----------|--------:|------|
| `Default` | 4,412 | Normal-length test rows |
| `OOD_ExtremeShort` | 59 | Length <60 |
| `OOD_ExtremeLong` | 531 | Length >1000, after removal of >2000 |
| `seq_Redundancy_90` | 3,936 | Maximum identity to train+val <90% |
| `seq_Redundancy_80` | 3,415 | <80% |
| `seq_Redundancy_70` | 3,063 | <70% |
| `seq_Redundancy_60` | 2,645 | <60% |
| `seq_Redundancy_50` | 2,145 | <50% |
| `seq_Redundancy_40` | 1,218 | <40% |
| `seq_Redundancy_30` | 386 | <30% |
| `OOD_Orphan` | 173 | No significant homologous hit |
| `TM-score_0.9` | 956 | Maximum TM-score <0.9 |
| `TM-score_0.8` | 339 | <0.8 |
| `TM-score_0.7` | 106 | <0.7 |
| `TM-score_0.6` | 13 | <0.6 |
| `TM-score_0.5` | 6 | <0.5 |
| `TM-score_0.4` | 6 | <0.4 |
| `TM-score_0.3` | 1 | <0.3 |
| `idr_ratio` | — | Disordered-residue fraction, continuous for all test rows |
| `OOD_IDR` | 48 | Continuous disordered regions >30% |
| `OOD_FoldHoldout` | 43 | CATH+TED topology holdout |
| `OOD_SuperfamilyHoldout` | 160 | CATH+TED superfamily holdout |
| `OOD_NewEC_L4` | 853 | All EC L4 absent from train+val |
| `OOD_NewEC_L3` | 3 | All EC L3 absent from train+val |
| `OOD_LongTail_EC_L3_le5` | 28 | All EC L3 training frequencies <5 |
| `OOD_LongTail_EC_L3_le10` | 75 | <10 |
| `OOD_LongTail_EC_L4_le5` | 1,524 | All EC L4 training frequencies <5 |
| `OOD_LongTail_EC_L4_le10` | 2,163 | <10 |
| `OOD_LongTail_KcatBin_le50` | 7 | Training count in the 0.5-log bin ≤50 |
| `OOD_LongTail_KcatBin_le100` | 133 | ≤100 |
| `InD` | 366 | Pure InD: 7.3% of test |

**Current split distribution consistency:**

| Dimension | Train | Val | Test |
|------|-------|-----|------|
| Samples | 15,479 | 2,150 | 5,002 |
| UniProts | 4,888 | 662 | 1,486 |
| Mean length | 410.9 | 421.5 | 489.7 |
| Length range | 61–990 | 72–965 | 16–1928 |
| Length <60 | 0 | 0 | 59 |
| Length >1000 | 0 | 0 | 531 |

## 4.10 Tools and parameters

| Tool | Version | Purpose | Main parameters |
|------|------|------|----------|
| mmseqs2 | v18.8cc5c | Sequence homology | `easy-search`, `--format-output query,target,fident,evalue` |
| Foldseek | commit 718d421 | Structural similarity | `search -a`, `--format-output query,target,alntmscore` |
| metapredict | v3.0.2 | Intrinsic disorder | Default disorder prediction |
| RDKit | — | ECFP4 fingerprints | Morgan radius=2, 2048 bits |

---

# Reproduction sequence

```bash
cd datasets/enzyme_kinetics_prediction

# Stage 1: obtain and parse sources; download in 1.1 is manual.
python scripts/01_download_alphafold_missing.py        # 1.2 Download missing AlphaFold structures.
python scripts/02_process_kcat_with_pdbs.py            # 1.3 PDB filtering and sequence consistency.
python scripts/03_decompress_pdbs_and_truncate_q9y233.py  # 1.3 Decompression and Q9Y233 truncation.

# Stage 2: aggregate and compute labels.
python scripts/04_compute_reactant_ecfp4.py            # 2.1 Reactant filtering and ECFP4.

# Stage 4: splitting and OOD annotation.
python scripts/05_create_length_ood_split.py           # 4.1 Length-based OOD split.
python scripts/06_create_pood_compatible_splits.py     # 4.2 Standard/length columns; mmseqs/Foldseek/IDR superseded.
# 4.3 Sequence/structure: .skills/homology-ood-annotation/scripts/seq_homology_ood.py and struct_homology_ood.py.
# 4.4 IDR: .skills/ood-annotation-toolkit/scripts/idr_ood.py.
python scripts/07_ted_augment_cath_ood.py              # 4.5 CATH+TED: invokes two skill scripts.
# 4.6 EC: .skills/ec-function-ood-annotation/scripts/add_unified_newec.py and add_unified_longtail_ec.py.
python scripts/08_add_value_bin_longtail_ood.py        # 4.6 Value-bin long tails; or toolkit value_bin_longtail_ood.py.
# 4.7 Remove >2000 aa manually; finalize_ood_columns.py --mode nan / --mode ind.

# 4.8 Read-only validation.
python ../../.skills/ood-annotation-toolkit/scripts/check_splits.py \
    --train splits/kcat_train.csv --val splits/kcat_val.csv --test splits/kcat_test.csv
```

---

# Output files

| File/directory | Description |
|-----------|------|
| `splits/kcat_{train,val,test}.csv` | Final splits: 15,479 /2,150 /5,002 rows |
| `pdbs/` | 7,058 AlphaFold v6 .pdb files, including full-length Q9Y233 backup |
| `data/kcat_filtered_with_pdb.csv` / `data/kcat_with_reactant_ecfp4.csv` | Stage 1/2 outputs: 22,736 /22,690 rows |
| `output/kcat_split_summary.json` | §4.1 script-time summary, test=5,040 |
| `output/kcat_pood_summary.json` | Historical first-pass OOD counts, reference only |
| `output/kcat_ecfp4_summary.json` | §2.1 ECFP4 summary |
| `output/kcat_longtail_kcatbin_summary.json` | Value-bin long-tail summary |
| `output/longtail_ec_unified_summary.json` / `newec_unified_summary.json` | Unified EC summaries |
| `output/kcat_acc_cath_ted_labels.csv` / `kcat_ted_aug_acc_list.csv` / `kcat_cath_ted_holdout.csv` | CATH+TED annotations and holdouts |
| `output/kcat_test_ood_*.csv` | Historical OOD subsets with both old/new naming conventions; see notes |
| `output/backup_*/`, `output/len2000_removal_backup/` | Backups before split updates |
| `output/scripts_archive/` | Archived historical scripts |

---

# Archived scripts: `output/scripts_archive/` (unused)

| Script | Reason archived |
|------|----------|
| `add_additional_ood_labels.py` | Legacy `OOD_NewFunction` / `OOD_LongTailFunction` replaced by canonical all-not-in `OOD_NewEC_*` / `OOD_LongTail_EC_*`; old columns deleted |
| `reformat_splits_to_pood_format.py` | Core-column reorder, idr_ratio, and extreme-row NaN handled by project-wide schema migration/finalizer; historical `pdb_file`/`labels`/`file_type` schema obsolete |
| `analyze_catpred_data.py`, `analyze_distributions.py`, `compare_zenodo_github.py` | Exploratory statistics outside the main workflow |
| `split_kcat_random.py`, `create_ecfp4_splits.py`, `move_low_similarity_to_train.py`, `create_ood_test_labels.py`, `finalize_splits_and_replace_sequences.py` | Superseded random-split and writeback workflow |

# Notes

1. **Obsolete OOD columns:** `OOD_NewFunction`, `OOD_LongTailFunction`, and `OOD_LongTail_EC_L3` were replaced and removed. `output/kcat_pood_summary.json` and some `kcat_test_ood_*.csv` files retain old definitions; do not use them as current statistics.
2. **Two subset naming conventions:** `output/` contains newer `kcat_test_ood_ExtremeShort.csv` and older `kcat_test_ood_OOD_ExtremeShort.csv`, reflecting historical script conventions. These have not been cleaned up.
3. **Schema migration:** archived `reformat_splits_to_pood_format.py` wrote `pdb_file`/`labels`/`file_type`; a project-wide migration replaced them with `struct_file`/`label` and removed `file_type`. Current headers are authoritative.
4. **Manual steps:** the 38-row >2000-aa removal (§4.7) and schema migration were one-off applications of project-wide rules; only backups remain. Apply these rules when reproducing.
5. **Historical counting changes:** before 2026-08-27, OOD statistics (e.g. old seq_Redundancy_90=4,558) included extreme-length rows in numerators/denominators. Under the current NaN convention, non-length OOD counts refer only to 4,412 Default=True rows. §4.9 is authoritative.
6. **Rows versus accessions:** multiple substrates per UniProt make row-level EC OOD counts much larger than unique-accession counts. This document consistently uses rows, matching split writeback.
7. **Putative TED annotations:** prioritize experimental CATH and fill gaps using the union with TED. Rare comma-separated TED dual labels are parsed as unions; no comma-containing pseudo-labels remain.
8. **EC source drift:** retain original annotations from `data/kcat_with_reactant_ecfp4.csv`; do not substitute current UniProt annotations.
9. **Truncated Q9Y233:** `pdbs/AF-Q9Y233-F1-model_v6.pdb` is the 779-residue truncation; `AF-Q9Y233-F1-model_v6_full.pdb` is the full-length backup. Split `struct_file` references the former.
