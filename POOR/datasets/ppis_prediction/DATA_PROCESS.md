# PPIS Prediction Dataset Processing

> Complete workflow from raw DIPS-Plus v1.3.0 (Zenodo) to the final Protein–Protein Interaction Site Prediction (PPIS) dataset, reproducible in one pass from source data.
> Location: `datasets/ppis_prediction/`
> Source: https://zenodo.org/records/8140981
> Reference: Morehead et al., *DIPS-Plus: The enhanced database of interacting protein structures for interface prediction*, Scientific Data 2023.

**Task:** single-chain, residue-level binary classification. Given a chain's sequence and structure, predict whether each residue belongs to a protein–protein interface. The final **chain-centric** dataset identifies each sample by a unique chain sequence with one aggregated label set. The earlier pair-centric workflow is retained only as historical reference below.

**Current measured split sizes:**

| File | Rows | Description |
|------|------|------|
| `ppis_train.csv` | **6,043** | Training |
| `ppis_val.csv` | **672** | Validation |
| `ppis_test.csv` | **1,980** | 1,685 Default + 295 extreme-length samples, with OOD columns |

Core fields: `unique_id` (sequence MD5, or `seq_hash`), `aa_seq` (PDB-derived), `struct_file` (`{unique_id}.pdb` in `pdbs/`), `label` (residue interface labels: 1 interface, 0 non-interface, -1 missing; exactly aa_seq length), and `struct_label` (aligned to PDB residue order; §4.7).

---

## Stage 1: Obtain and parse source data

> Obtain DIPS-Plus v1.3.0, extract and index it, and parse metadata and pair-level interface labels.

### 1.1 Download: `scripts/01_download_dips_plus.py`

| Item | Details |
|------|------|
| Zenodo record | 8140981, DOI `10.5281/zenodo.8140981`, version 1.3.0 |
| File | `final_raw_dips.tar.gz`, 15.7 GB, MD5 `a4af4e14162aa88a59c05fef5d562088` verified |
| Downloader | `aria2c`, four connections, resumable; optional `--msa` archive was not downloaded (see notes) |

### 1.2 Extract and inspect: `scripts/02_extract_and_inspect.py`

Extract into `data/raw/` (~88 GB): official `pairs-postprocessed*.txt` split lists, `*_pdb_metadata.csv`, and pair files bucketed by the first two PDB-ID characters (`aa/`…`zz/`). Pair filenames are `{pdb_id}.pdb{model}_{pair_index}.{ext}`, using dill or HDF5.

### 1.3 Index files: `scripts/03_build_pair_index.py`

Scan all pairs: **42,112** unique complexes. Both dill and HDF5 exist; prefer HDF5, read with `hickle`.

Outputs: `output/dips_plus_pair_index.csv` (42,112 rows), `output/dips_plus_split_index.csv`, `output/extraction_summary.json`.

### 1.4 Parse metadata: `scripts/04_parse_pair_files.py`

Read each eight-element hickle structure: `complex_name`, `df0`, `df1`, `pos_idx`, `neg_idx`, `srcs`, `id`, `sequences`. Write 42,112 rows to `output/dips_plus_metadata.csv`, including `unique_id`, `pdb_id`, `file_path`, `file_format`, receptor/ligand sequences and lengths, chain IDs, atom/residue counts, `pos_pairs`/`neg_pairs`, and `has_idr_*`.

### 1.5 Pair-level interface labels: `scripts/05_generate_interface_labels.py`

Derive residue labels from presampled `pos_idx` positive interface atom pairs, based on a 6 Å atom-distance threshold. Label a residue 1 if any of its atoms appears in `pos_idx`.

Outputs: 42,112 `output/interface_labels/{unique_id}_labels.npz` files with `receptor_labels`/`ligand_labels`, `output/interface_label_distribution.csv`, and `output/interface_label_summary.json`. Overall pair-level interface residue fraction: 14.8%.

> This script was briefly archived with the pair-centric route and restored as an official step in 2026-08. Pair labels feed chain-centric aggregation and provide the reference implementation of pos_idx semantics.

---

## Stage 2: Aggregate and compute labels

> Union interface residues across all occurrences of each unique chain sequence.

### 2.1 Chain-centric aggregation: `scripts/06_build_chain_centric_dataset.py`

1. Identify each complete amino-acid sequence by MD5 (`seq_hash`).
2. Traverse all 42,112 pairs and extract receptor/ligand interface labels using §1.5 semantics.
3. **Union** interface positions across all occurrences of each sequence.

| Metric | Value |
|------|------|
| Original pairs | 42,112 |
| Unique chain sequences | **22,223** |
| Repeated unique chains | 12,422 (55.9%) |
| Mean interface fraction | 28.3%; union covers multiple binding regions, versus pair-level 14.8% |

Outputs: `output/chain_centric_metadata.csv` (22,223 rows), `output/chain_interface_labels/{seq_hash}_labels.npz`, `output/chain_centric_summary.json`.

---

## Stage 3: Redundancy reduction and refinement

> Deduplicate at 95% identity and extract the most complete correctly assigned chain structure.

### 3.1 Sequence deduplication: `scripts/07_dedup_chain_centric.py`

Run mmseqs2 `cluster` on 22,223 unique chains at 95% identity / 80% coverage / cov-mode 0. Keep the automatically selected representative per cluster.

| Metric | Value |
|------|------|
| Representatives | **8,850** |
| Reduction | 60.2% |
| Mean / maximum cluster size | 2.51 / 64 |

Outputs: `output/chain_centric_dedup_metadata.csv`, `output/chain_centric_dedup_summary.json`.

> Historical `dedup_chain_centric_and_link_structures.py` linked representatives to pair-centric structures in `output/chain_representative_structures/`. Obsolete code was removed; §3.2 fills structure-source fields directly.

### 3.2 Representative structures: `scripts/08_extract_representative_structures.py`

1. Scan all pair metadata; select the occurrence with **the most atoms** per unique sequence.
2. Load `df0`/`df1` and select **by sequence**: exact match, unique length match, then highest sequence similarity. These DataFrames can be swapped relative to receptor/ligand metadata; role alone is unsafe.
3. Write single-chain PDBs with chain ID `A`; fill `structure_source_unique_id` and `structure_source_role` (`{role}->{matched_df}`) in deduplicated metadata.

> Sequence matching incorporates the one-off `fix_representative_structures.py` patch; role-only selection originally extracted ~6,786 chains incorrectly. Its resname spacing bug, shifting chain/resseq columns right, is not retained. Current output uses standard columns, with coordinates starting in columns 31–38.

| Metric | Value |
|------|------|
| Representatives with structures | 8,850 / 8,850 (100%) |
| PDB/metadata sequences identical after correction | 8,741 |
| Equal length with residue substitutions, mostly modified residues | 63 |
| Remaining length mismatches | 46; aligned in §4.3 or cleaned in §4.4 |

Output: 8,850 `pdbs/{seq_hash}.pdb` files.

---

## Stage 4: Split and annotate OOD

> Construct 72:8:20 splits, annotate OOD, and apply canonical finalization.

### 4.1 Split: `scripts/09_split_chain_centric.py`

1. Reserve chains outside **[60, 1000]** as `test_extremelength`.
2. Randomly split remaining chains **72:8:20** (seed=42). Representatives have already been deduplicated at 95% identity to control sequence leakage.

At this stage, 8,850 samples comprise train 6,153 / val 683 / test 1,710 / extreme-length 304.

Outputs: `output/chain_centric_splits/chain_centric_{train,val,test,test_extremelength}.csv`, `output/chain_centric_split_summary.json`.

### 4.2 Initial sequence/structure/IDR OOD: `scripts/10_mark_ood_chain_centric.py`

Reference: Train + Val. Outputs: `output/chain_centric_splits/chain_centric_test_ood.csv`, `output/ood_chain_centric_summary.json`.

- **Sequence homology** (mmseqs2 search, `OOD_LowHomology_30..90`, `OOD_Orphan`) and **structural similarity** (foldseek, `OOD_NewFold_0.3..0.9`) are **superseded** by `.skills/homology-ood-annotation/seq_homology_ood.py` / `struct_homology_ood.py`, using easy-search. Both current annotations were recomputed: TM on 2026-08-28 using cropped structures (§4.4); sequence/Orphan on 2026-08-29. Legacy default mmseqs search sensitivity missed hits and overestimated OOD: Orphan 522→458, seq_Redundancy_90 1,570→1,457, InD 64→134. Backup: `output/backup_before_seq_easysearch_unify/`; audit: `analysis/output/ppis_seq_easysearch_audit/`; m8: `output/ood_mmseqs_test_vs_trainval.m8`.
- **OOD_IDR remains applicable:** use DIPS-Plus **flDPnn**, flagging when >30% of interface residues have disorder propensity >0.5. This differs from toolkit metapredict and must not be mixed with it.

### 4.3 Final files and label alignment: `scripts/11_build_final_ppis_splits.py`

1. Merge four intermediate splits; extract aa_seq from each `pdbs/{seq_hash}.pdb`.
2. Load `output/chain_interface_labels/{seq_hash}_labels.npz`. Use directly if PDB sequence and labels have equal length; otherwise map through **Biopython global alignment**: identity ≥0.85, or identity ≥0.75 and coverage ≥0.75. PDB insertions receive -1; skip failed alignments.
3. Combine normal and extreme-length test samples; generate Default/ExtremeLong/ExtremeShort flags. Compute regular OOD only for Default=True.
4. Initial CATH-name holdouts and SIFTS NewFunction are later replaced by skills (§4.5–4.6). Sparse-interface `OOD_LongTail` remains applicable: positive fraction <5% among valid labels.

Alignment: 8,849 retained samples, comprising 8,807 direct matches and 42 aligned mappings; one skipped (`8b00203a…`, PDB length 279 versus metadata 488, unalignable).

### 4.4 Crop and clean structures: `scripts/12_crop_ppis_structures.py`

Require a chain exactly matching aa_seq. Otherwise select the best chain, globally align, and **crop** insertions for one-to-one residue correspondence: structure length = aa_seq length = label length. Remove rows with identity <0.95 or missing structural residues as unrecoverable.

- Remove **154** samples: train −109 / val −11 / test −34; crop **61** structures. Final size: **8,695**.
- Quarantine original PDBs in `output/quarantine_structures/`; pre-removal splits in `output/quarantine_splits/`.
- Inline Biopython `protein_letters_3to1_extended` mapping replaces the former dependency on `analysis/scripts/verify_aa_seq_vs_structure.py`.

### 4.5 CATH+TED annotation and holdouts: task wrapper `scripts/13_ted_augment_cath_ood.py`

Canonical implementation: `.skills/cath-ted-domain-annotation/`, using `cath_ted_annotate.py` and `cath_ted_holdout.py`.

1. Map seq_hash → source pair/role → pdb/chain via `src_chain`, reusing archived `add_cath_holdout_ood.py` logic. Write 8,695 chains to `output/ppis_ted_aug_chain_list.csv`.
2. `cath_ted_annotate.py --offline`: union CATH v4.4 name matches with putative TED annotations, requiring overlap ≥50% of domain length. REST cache: `data/ted/ted_api_cache.json`. Write `output/ppis_chain_cath_ted_labels.csv` and `output/ppis_chain_domains.csv`.
3. `cath_ted_holdout.py`: training CATH+TED union reference; any semantics (any unseen topology/superfamily →True, unannotated →False). Update test `OOD_FoldHoldout` / `OOD_SuperfamilyHoldout`; Default=False →NaN (§4.7).

Current Default test n=1,685: fold holdout **65**, superfamily holdout **235**. Sources across all chains: both 7,125 / cath-only 985 / ted-only 323 / none 262. CATH-only backup: `output/backup_before_ted_aug/`.

### 4.6 Canonical EC OOD: `.skills/ec-function-ood-annotation/`

- **NewEC**, `add_unified_newec.py`: train+val EC reference, all-not-in semantics. L4 **144**, L3 **9**.
- **LongTailEC**, `add_unified_longtail_ec.py`: training-only frequencies, all ECs below 5 /10, including zero. L3_le5 **62**, L3_le10 **137**, L4_le5 **515**, L4_le10 **659**, among 1,685 Default test samples.
- Source: `data/sifts/sifts_chain_ec.tsv.gz`, mapped by §4.5 src_chain. Default test EC coverage: 873 samples (51.8%).
- EC long-tail flags coexist with sparse-interface `OOD_LongTail` (§4.3), which has different semantics.

### 4.7 Finalization and validation: skill chain

Run scripts under `.skills/` in order:

1. **homology-ood-annotation:** canonical sequence/Orphan/TM columns. TM recomputed 2026-08-28 using 1,980 test queries /6,715 train+val structures; cache `output/foldseek_easysearch_test_vs_trainval.m8`. Sequence/Orphan recomputed 2026-08-29 with easy-search -s 7.5, 1,980 test queries against deduplicated train+val; cache `output/ood_mmseqs_test_vs_trainval.m8`.
2. **ood-annotation-toolkit `idr_ood.py`: not applicable.** Retain DIPS-Plus flDPnn (§4.2).
3. **cath-ted-domain-annotation:** §4.5.
4. **ec-function-ood-annotation:** §4.6. Regression value-bin long-tail does not apply to residue classification.
5. **finalize**, `finalize_ood_columns.py`:
   - `--mode nan`: validate that 295 Default=False extreme-length rows have NaN in non-length OOD, seq, and TM columns. Since 2026-08-28, generating scripts 11/13 and skills write NaN directly; finalization is a no-op. That day, historical False seq/TM entries became NaN; backup `output/backup_nondefault_seqtm_nan/`.
   - `--mode ind`: all OOD/seq/TM flags False, treating NaN as False. Current **134/1,980 (6.8%)**.
6. **struct_label**, `struct-label-alignment/add_struct_label.py`: enumerate ATOM residues in file order by first (chain, resseq, icode), globally align labels to the PDB sequence, and assign -1 when unmappable. Added to all splits; backup `output/backup_before_struct_label/`.
   - **Resseq audit/fix, 2026-09-02**, `scripts/14_fix_struct_label_resseq.py`: audit 13 structure/aa_seq length mismatches with missing residues. Eleven mappings were already correct. Two training rows (`fc54381e1a5a…`, `16d04ef06e42…`) had C-terminal tags outside aa_seq coverage, incorrectly labeled 0/1; correct to -1 using resseq. Backup: `output/backup_20260902_structlabel_fix_ppis_train.csv`. After correction these two rows each contain one -1 residue.
   - Alignment inheritance can place gaps ambiguously when missing residues occur in repeats/low-complexity regions (SSP's 148-row misalignment event; `ss_prediction/DATA_PROCESS.md` §5.5 and skill notes). A comparable full LBS audit found no issue: 96,077 direct and 179 equal-length substitutions.
7. **check_splits.py:** all three pairwise unique_id intersections are **0**, label/aa_seq lengths agree everywhere, and no struct_file is missing from pdbs/.

Project-wide schema migration, without a standalone task script: `file_type→struct_file`, `labels→label`, `OOD_LowHomology_*→seq_Redundancy_*`, `OOD_NewFold_*→TM-score_*`, `OOD_CATH_FoldHoldout→OOD_FoldHoldout`, `OOD_NewFunction→OOD_NewEC_L4`, `ExtremeLong/Short→OOD_ExtremeLong/Short`.

### 4.8 Final measured split statistics

| Split | Samples | Fraction of normal 8,400 | Mean length | Mean interface fraction |
|------|--------|----------------------|---------|-------------|
| Train | 6,043 | 71.9% | 252.8 | 26.4% |
| Val | 672 | 8.0% | 262.3 | 26.3% |
| Test Default | 1,685 | 20.1% | 256.9 | 26.4% |
| Test extreme-length | 295 | — | 173.1 | 35.9%; 34 long chains 1,012–1,484 aa +261 short 7–59 aa |

True counts among 1,980 test rows; non-Default/Extreme/InD flags are evaluated within Default:

| OOD column | True | | OOD column | True |
|--------|-----:|---|--------|-----:|
| `Default` | 1,685 | | `TM-score_0.9` / `0.5` / `0.3` | 770 / 47 / 20 |
| `OOD_ExtremeLong` / `Short` | 34 / 261 | | `OOD_IDR` | 28 |
| `seq_Redundancy_90` / `70` / `50` / `30` | 1,457 / 1,298 / 1,066 / 602 | | `OOD_FoldHoldout` / `SuperfamilyHoldout` | 65 / 235 |
| `OOD_Orphan` | 458 | | `OOD_NewEC_L4` / `L3` | 144 / 9 |
| `OOD_LongTail` (sparse interface <5%) | 48 | | LongTailEC L3_le5/le10, L4_le5/le10 | 62 / 137 / 515 / 659 |
| `InD` | 134 | | | |

---

## Output files

```
datasets/ppis_prediction/
├── DATA_PROCESS.md                        # This document
├── splits/                                # Final three split files
│   ├── ppis_train.csv / ppis_val.csv / ppis_test.csv
├── pdbs/                                  # 8,833 representative PDBs; see note 7
├── data/raw/                              # DIPS-Plus pairs (~88 GB)
├── logs/                                  # Run logs
├── scripts/                               # Active scripts 01–13, described above
└── output/
    ├── dips_plus_pair_index.csv / dips_plus_metadata.csv
    ├── interface_labels/                  # 42,112 pair-label NPZs (§1.5)
    ├── chain_centric_metadata.csv / chain_interface_labels/  # 22,223 unique-chain NPZs
    ├── chain_centric_dedup_metadata.csv / chain_centric_splits/
    ├── ppis_chain_cath_ted_labels.csv / ppis_chain_domains.csv / ppis_ted_aug_chain_list.csv
    ├── newec_unified_summary.json / longtail_ec_unified_summary.json
    ├── backup_before_*/                   # Split backups before updates
    ├── quarantine_splits/ / quarantine_structures/  # §4.4
    └── scripts_archive/                   # Archived scripts below
```

## Quality checks

1. **Leakage:** measured pairwise unique_id intersections are zero; representatives were deduplicated at 95% identity.
2. **Consistency:** all 8,695 samples have label length equal to aa_seq and corresponding structural residues. Recorded -1 masked residues: 140, all in test. No missing struct_file.
3. **Distributions:** train/val/Default-test mean interface fractions are ~26.4%, with similar length distributions.

## Archived scripts: `output/scripts_archive/`

| Script | Reason |
|------|---------|
| `fix_representative_structures.py` | Sequence-based DataFrame selection merged into 08, excluding its resname-format bug |
| `18_fix_pdb_column_layout.py` | One-off PDB column/concatenated-chain repair applied to data; backup `output/backup_before_pdb_layout_fix/`; root cause fixed in 08 |
| `add_function_longtail_ood.py` | Initial EC/interface long-tail extension; EC replaced by §4.6, interface logic incorporated in 11 |
| `evaluate_cath_coverage.py` | Diagnostic CATH coverage output only |
| `18_add_longtail_ec_ood.py` | Legacy any-semantics single column replaced by four all-not-in columns |
| `add_cath_holdout_ood.py` | CATH-name-only holdouts replaced by §4.5 CATH+TED |
| `generate_interface_labels.py` | Restored as active 05 in 2026-08; identical archive copy deduplicated |
| `split_dataset.py`, `filter_by_pairwise_identity*.py`, `filter_by_sequence_identity.py`, `split_82_and_mark_lowhomology.py`, `mark_ood.py`, `mark_ood_82split.py`, `extract_chain_structures.py`, `fix_mismatched_structures.py` | Historical pair-centric route replaced by chain-centric workflow |

## Notes

1. **Storage:** core archive 15.7 GB; extracted data/raw ~88 GB; mmseqs2/foldseek intermediates require additional space.
2. **MSA archive not downloaded:** `interim_external_feats_dips_msas.tar.gz`, 12.0 GB, omitted for space.
3. **IDR:** use flDPnn propensity >0.5 in >30% of interface residues, not toolkit metapredict.
4. **Homology definitions:** both sequence and TM use canonical easy-search. TM recomputed 2026-08-28 on cropped structures; legacy search -a 1 overestimated OOD: TM0.9 1,600→770 and TM0.3 1,365→20. Backup `output/backup_before_foldseek_easysearch_unify/`; audit `analysis/output/foldseek_easysearch_audit/ppis/`. Sequence/Orphan recomputed 2026-08-29: seq90 1,570→1,457, seq30 686→602, Orphan 522→458, InD 64→134. Backup `output/backup_before_seq_easysearch_unify/`; audit `analysis/output/ppis_seq_easysearch_audit/`; m8 `output/ood_mmseqs_test_vs_trainval.m8`.
5. **PDB coordinate formatting:** existing pdbs/ already uses standard columns 31–38. The pre-cleanup step-08 template left three spaces after resseq, shifting replayed coordinates left and potentially losing minus signs. Four spaces now restore standard layout; 60 sampled chains reproduce existing PDBs byte for byte.
6. **Quarantine:** `output/quarantine_structures/` holds 62 pre-crop PDBs; `output/quarantine_splits/` is provenance only.
7. **PDB count:** 8,833 files versus 8,695 referenced samples; 138 retained structures belong to removed rows (§4.4), preserved under project constraints. Follow struct_file references.
8. **Test source:** DIPS-Plus has no standard test set; this project randomly splits representatives 72:8:20. For literature comparisons, official val may be used as test instead.
9. **Non-script step:** schema migration (§4.7) is project-wide, with no independent task script.
10. **History:** pair-centric→chain-centric; CATH-only→CATH+TED (experimental CATH prioritized, putative TED fills gaps); initial EC columns→canonical all-not-in. Backups in output/backup_before_*/.
11. **Citation:** cite Morehead et al., Scientific Data 2023 when using DIPS-Plus.

---

## Historical pair-centric route: archived

> Originally constructed at complex-pair level, retained for comparison. Scripts are in `output/scripts_archive/`; intermediates remain in output/, including interface_labels/, structures/, dips_plus_metadata_filtered_global_95id.csv, and ood_82split_summary.json.

1. **Labels:** generate_interface_labels.py, now restored as 05, derives pair receptor/ligand labels from pos_idx.
2. **Deduplication:** filter_by_pairwise_identity_global.py requires receptor–receptor and ligand–ligand identity ≥95%, coverage ≥80%; mmseqs2 self-search and greedy maximal independent set reduce 42,112→10,775 pairs.
3. **Splitting:** split_82_and_mark_lowhomology.py splits by PDB ID 8:2, then the training pool 9:1: train 7,821 / val 870 / test 2,084.
4. **OOD:** mark_ood_82split.py uses strict two-chain criteria for LowHomology/NewFold/Orphan/ExtremeLength/IDR.
5. **Structures:** extract_chain_structures.py writes 21,550 chain PDBs under output/structures/.

**Why chain-centric:** the final task predicts interfaces for one chain. Pair-level construction assigns different interfaces to identical sequences in different complexes, conflicting with single-chain model inputs. Unioning labels by unique sequence yields one chain and one label set per sample.
