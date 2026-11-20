# Ligand Binding Site Dataset Processing

> Complete raw-data reproduction workflow for a residue-level protein ligand binding site (LBS) dataset built from PDB structures.
> Location: `datasets/ligand_binding_site/`.
> Source: RCSB PDB (`https://www.rcsb.org/`); local project-root mirror in `data/pdb/`.

## Task definition and current scale

Binding sites are **protein residues within 4.5 Å of ligand heavy atoms**. Ligands include small-molecule cofactors, drugs, and metal-ion complexes, excluding water, nucleic acids, standard amino acids, common crystallization additives, and unresolved UNK/UNL (§Stage 2). Each chain–ligand pair has a binary array aligned to aa_seq: 1=binding, 0=nonbinding.

Current split sizes, measured from `splits/`:

| Metric | Train | Val | Test |
|------|------:|------:|------:|
| **Samples** | **71,093** | **7,938** | **17,225** |
| Unique protein chains | 65,462 | 6,718 | 14,328 |
| Unique PDBs | 60,136 | 5,174 | 11,379 |
| Unique ligand types | 31,940 | 2,826 | 5,125 |

Three mutually exclusive subsets cover the entire test set:

| Subset | Definition | Samples | Percentage |
|------|------|------:|-----:|
| `Default` | 60 ≤ length ≤ 1000 | **13,723** | 79.7% |
| `OOD_ExtremeShort` | length < 60 | **1,974** | 11.5% |
| `OOD_ExtremeLong` | length > 1000 | **1,528** | 8.9% |

Structures: 96,979 `{unique_id}.pdb` files in pdbs/, including 723 redundant files for discarded missing-SMILES samples (Notes).

---

## Stage 1: Obtain and parse sources (PDB →chain–ligand labels)

> Traverse local PDB/mmCIF files, identify protein–ligand complexes, compute residue labels, and aggregate into a chain–ligand dataset.

### Core scale

| Metric | Value |
|------|------|
| **Candidate PDB structures** | 255,227 |
| Successfully processed | 128,725 |
| Timed out | 1,099 |
| Parsing failures | 77 |
| Pairs before aggregation, including duplicates | 763,605 |
| **Pairs after unique_id deduplication** | **645,774** |
| Unique ligand types | 46,869 |

### Step 1: Batch-extract protein–ligand complexes

**Script:** scripts/01_extract_binding_sites.py; core parsing in the unnumbered imported library scripts/extract_binding_sites_core.py.

- Input: all `.pdb.gz` and `.cif.gz` files under project-root data/pdb/.
- Output: output/binding_sites_batch_{NNNNN}.csv.
- Parameters: heavy-atom distance 4.5 Å, using scipy.spatial.cKDTree; minimum ligand heavy atoms 5; minimum binding residues >3 (implementation accepts at least 3).
- Parse structures with BioPython PDBParser /MMCIFParser.

### Step 2: Aggregate batches

**Script:** scripts/02_aggregate_batches.py.

- Input: 221 output/binding_sites_batch_*.csv batches.
- Deduplicate unique_id, keeping the first: 763,605 →**645,774** pairs.
- Outputs: output/ligand_binding_sites_raw.csv and output/ligand_binding_sites_std.csv.

> **Archive note:** process-pool crashes interrupted the main run. One-off 03_diagnose_remaining_files.py classified files OK/Timeout/Error/Crash with 10-s limits; 05_reprocess_timeouts.py retried at 60 s, but the first 40 smallest files still timed out and produced nothing, so it was stopped. Both and the one-off 02_resume_extract_binding_sites.py are archived under output/scripts_archive/. Accept skipping 1,099 large timed-out files, approximately 0.4% of PDBs.

### Key distributions

Top ten ligands by chain–ligand pair count among 128,725 processed structures:

| Ligand | Pairs |
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

### Outliers and filtered samples

| Type | Count | Reason |
|------|------|------|
| Timeout | 1,099 | Parsing exceeded 10-s diagnostic /60-s retry limits |
| Parse error | 77 | Missing mmCIF columns (_atom_site.id, _atom_site.pdbx_PDB_ins_code) or invalid B-factors |

---

## Stage 2: Label aggregation and computation (semantics and ligand filters)

> Labels are computed during Stage-1 extraction in scripts/extract_binding_sites_core.py using the following definitions. No additional script runs here.

### Labels

- Binding residue: any protein heavy atom lies ≤4.5 Å from any ligand heavy atom.
- label: binary residue array aligned to aa_seq, serialized as bracketed space-separated integers, e.g. [0 0 1 0]. **Store as a string**, requiring len(label)==len(aa_seq).

### Ligand filters

| Item | Rule |
|--------|------|
| Water | Exclude HOH, DOD, WAT, etc. |
| Nucleic acids | Exclude A/T/G/C/U/DA/DT/... |
| Standard amino acids | Exclude all 20 standard amino acids |
| Modified amino acids | Map 34 types to standard residues, e.g. MSE→M, and treat as protein |
| Crystallization additives | Exclude 127 common additives; list in extract_binding_sites_core.py |
| Ligand size | At least five heavy atoms |
| Site size | Retain chain–ligand pairs with >3 binding residues |

### Validation and quality control

- All label arrays match aa_seq lengths.
- Record num_binding_residues and num_ligand_heavy_atoms per pair.
- Deduplicate unique_id={pdb_id}_{chain_id}_{ligand_name}_{ligand_chain}_{ligand_resnum} to avoid repeated pairs.

---

## Stage 3: Deduplication and refinement

> Deduplicate sequences within ligands, extract individual chains and validate sequence–structure consistency, annotate ligand SMILES/ECFP4, and discard unparseable samples.

### 3.1 Ligand-wise 95% sequence deduplication

Deduplicate ligand_binding_sites_std.csv to reduce homogeneous-sample overfitting and support fair evaluation. Group by **ligand name**; cluster protein chains at 95% identity and retain one representative per cluster. Similar sequences associated with different ligands are not mutually removed, preserving ligand-dependent binding diversity.

| Metric | Value |
|------|------|
| Pairs before deduplication | 645,774 |
| **Pairs after deduplication** | **109,954** |
| PDBs before | 128,725 |
| PDBs after | 88,208 |
| Unique ligand types | 46,869, unchanged |

**Step 1: iterative mmseqs2 clustering**, scripts/03_dedup_by_ligand_95.py.

- Input output/ligand_binding_sites_std.csv; output output/ligand_binding_sites_std_dedup95.csv.
- Add a 150-aa ligand-specific prefix per ligand_name so different ligands cannot cluster together.
- Iterate `mmseqs easy-cluster --min-seq-id 0.95 -c 0.95 --cov-mode 1 --cluster-mode 1` until no removal.
- Representative priority: most binding residues →best resolution →longest chain →most ligand heavy atoms.

**Step 2: exact cleanup**, scripts/04_exact_cleanup_dedup95.py.

- Ligands with ≤50 samples: construct a graph from exact edlib global alignments with identity≥0.95; retain one per connected component.
- Ligands with >50 samples: run separate mmseqs easy-cluster clustering.
- Invoke 03_dedup_by_ligand_95.py through subprocess for final convergence, ensuring no remaining ≥95%-similar pairs.

**Validation:** rerun easy-cluster with the same parameters separately for the five largest ligand groups. All sequences form singletons, max_cluster_size=1, redundant=0. Another deduplication iteration immediately converges, removing zero in round two.

Top ten ligands before/after deduplication:

| Ligand | Before | After |
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

### 3.2 Protein chain structure extraction

> Extraction consumes the split lists after Stage-4 temporal splitting and ligand filtering. Grouped here by refinement topic; in a linear run it follows Stage-4 Step 2.

**Script:** scripts/07_extract_protein_chains.py.

- Group by PDB ID and parse each source once. Use fast manual ATOM parsing for PDB gzip and parse the _atom_site loop for CIF gzip.
- Keep **protein residues only**, standard/common modified amino acids; skip nucleic acids, ligands, and additives. Retain modified MSE/SEP/TPO even as HETATM.
- Truncate output chain IDs to one character for PDB compatibility; renumber atoms per chain, avoiding coordinate corruption from serial numbers exceeding five columns.
- Keep only the first NMR model: parse_pdb_file stops at first ENDMDL; parse_cif_file retains pdbx_PDB_model_num=1. Strip quotes from mmCIF atom names such as "C2'". Root-cause fixes are integrated; 549 historical multimodel files were re-extracted by a one-off script (archive list).
- Accept uppercase/lowercase source filenames; 16 processes, approximately six minutes.
- Output pdbs/{unique_id}.pdb: all then-current 96,979 samples extracted successfully, zero missing sources.

### 3.3 Sequence–structure consistency QC

**Script:** scripts/13_verify_seq_structure_consistency.py, retained as final diagnostic/QC.

- Approximately 0.80% initially differed between aa_seq and extracted structure sequences due to modified residue names, repeated residue numbers, and non-protein atoms.
- Historical differences: output/seq_structure_mismatches.csv. Current split sequence–structure agreement: **100.00%**.

> **Archive note:** one-off 15_patch_inconsistent_from_structure.py recomputed aa_seq/label with core protein-identification rules; 16_reextract_inconsistent_biopython.py re-extracted and renumbered repeated residues. Agreement improved from 96,201/96,979 (99.20%) to 100%. Repairs are applied, root causes fixed in core/07, both scripts archived.

### 3.4 Ligand SMILES/ECFP4 and missing-sample removal

**Scripts:** scripts/08_add_ligand_smiles_ecfp4.py and scripts/09_drop_missing_smiles_and_rebuild.py.

- Derive ligand_name from unique_id; read SMILES/InChI from PDB CCD data/pdb/components.cif.gz, _pdbx_chem_comp_descriptor.
- RDKit computes canonical SMILES and Morgan ECFP4, radius=2,nBits=2048. For aromatic/charged SMILES failing kekulization, use sanitize=False and a fallback omitting kekulization.
- Store ECFP4 as sparse bit indices, e.g. [12,45,200].
- Script 09 drops empty ligand_smiles: 548 train, 20 test, zero val. Refresh NewFunction/LongTail incrementally by default; --recompute_full_ood invokes 11 for full recomputation.

| Metric | Train | Val | Test |
|------|------:|------:|------:|
| Before SMILES/ECFP4 | 71,641 | 7,938 | 17,400 |
| After missing-SMILES removal | **71,093** | **7,938** | **17,380** |

Validation: missing SMILES train/val/test=0/0/0; empty fingerprint lists=0; label/sequence length mismatches=0.

> Script 09's has_ec_annotation/OOD_NewFunction are historical pre-unification names, replaced by OOD_NewEC_L4, etc. (Stage 4). The script was adapted to label/struct_file/OOD_Extreme* in 2026-08.

---

## Stage 4: Splitting and OOD annotation

> Split train/val/test chronologically, filter chemical similarity, regroup extreme lengths, and add complete test OOD flags.
> Canonical definitions reside in .skills/; use their scripts for recomputation rather than alternative definitions.

### Step 1: Temporal train/val/test splitting

**Script:** scripts/05_create_time_based_splits.py.

- Input: output/ligand_binding_sites_std_dedup95.csv.
- Sort unique pdb_id values by deposition date and split 70%/10%/20%.
- All pairs from a PDB stay in one split, preventing leakage.
- Assign 37 PDBs with unparseable deposition dates, 49 pairs, to **test**.

| Split | PDBs | Chain–ligand pairs | Ligand types | Dates |
|------|--------|-------------------|-----------|---------|
| **Train** | 61,731 | 73,608 | 32,630 | 1974-06-06 ~ 2021-04-28 |
| **Val** | 8,818 | 12,487 | 6,480 | 2021-04-28 ~ 2022-12-22 |
| **Test** | 17,659 | 23,859 | 12,316 | 2022-12-23–2026-05-14, including 20 undated PDBs |

At initial splitting, all pairwise train/val/test pdb_id intersections are zero.

> Historical random 80/10/10 pdb_id splitting, seed42 (06_create_splits.py), is fully superseded by temporal splitting and archived.

### Step 2: Validation/test chemical similarity filtering

**Script:** scripts/06_filter_test_by_ligand_similarity.py.

- Parse CCD SMILES; generate RDKit Morgan ECFP4, radius=2,nBits=2048.
- For each val/test ligand, compute maximum Tanimoto to all training ligands. Keep samples with maximum **>0.6**; leave training unchanged.

| Split | Before | After | Missing SMILES | Fingerprint failures | Below threshold |
|------|--------|--------|-------------|----------|------------|
| Val  | 12,487 | **8,555** | 31 | 234 | 3,667 |
| Test | 23,859 | **14,816** | 35 | 475 | 8,533 |

### Step 3: Regroup extreme lengths and remove oversized chains

**Script:** scripts/10_rebuild_splits_with_extreme_subset.py.

- protein_length<60 →OOD_ExtremeShort; >1000 →OOD_ExtremeLong.
- Move all extreme samples from train/val/test into **test**; normal test rows receive Default=True.
- Recompute Default OOD columns by subprocess calling scripts/11_build_ood_splits.py.
- After regrouping and §3.4 SMILES removal: train71,093 /val7,938 /test17,380 =Default13,723 +short1,974 +long1,683.
- **>2000 removal, no standalone script:** globally exclude single chains >2000. All **155** affected samples are test ExtremeLong, lengths2,003–4,584. Test17,380→**17,225**; backup output/len2000_removal_backup/.

> Script 10's ExtremeShort/ExtremeLong/is_time_cutoff are legacy names, later unified as OOD_ExtremeShort/Long, etc. Core label/struct_file names and undefined ROOT path were fixed in 2026-08.

### Step 4: Sequence-homology OOD, Default only

**Canonical tool:** `.skills/homology-ood-annotation/scripts/seq_homology_ood.py`. Query=Default test, target=train+val; mmseqs2 **easy-search -s 7**, default significance e≤1e-3, output query,target,pident,evalue. Threshold max pident at90/80/…/30. Orphan means no significant hit, implemented as best e>1e-3. Deduplicate queries by chain_key={pdb}_{chain}.

- Inline mmseqs in scripts/11_build_ood_splits.py is **superseded** by the skill. Old search -e 100 -c 0.0 caused reversed Orphan/seq30 counts; archived one-off 24_fix_seq_homology_easysearch.py repaired them using easy-search -s 7.
- Cache: output/ood/mmseqs_max_identity_easysearch.csv, keyed by chain_key; verify that key space before replay, not unique_id.
- Default=False sequence/Orphan flags are NaN, not computed.
- Verified nesting: Orphan⊆seq30⊆…⊆seq90.

### Step 5: Structural-similarity OOD, Default only

**Script:** scripts/11_build_ood_splits.py embeds Foldseek calls; canonical `.skills/homology-ood-annotation/scripts/struct_homology_ood.py`.

- Cache uses **unified easy-search**, recomputed2026-08-28: `foldseek easy-search query target out tmp --format-output query,target,alntmscore`, default e≤10,sensitivity9.5. Query=14,328 unique test chain_keys, target=72,180 unique train+val chains, symlink directories avoid copies. Threshold max alntmscore at0.9/0.8/0.7/0.6/0.5/0.4/0.3.
- Old search -a -e inf --max-seqs 1000 +convertalis cache had 2,433/11,406 (21%) invalid scores>1, possibly wrong output columns; valid scores were systematically high, median inflation0.074 on937 chains. This severely undercounted OOD (TM0.5:4→241). Comparison scripts/reports in output/foldseek_easysearch_compare/ (run_compare.py,compare_report.txt,easysearch_max_tmscore.csv); old cache/test backup output/backup_before_foldseek_easysearch_unify/.
- Cache output/ood/foldseek_max_tmscore.csv:14,328 test chain_keys. Default=False flags receive NaN.

### Step 6: Intrinsic-disorder OOD, Default only

**Canonical tool:** `.skills/ood-annotation-toolkit/scripts/idr_ood.py`, metapredict v3 (measured3.0.2); local scripts/11_build_ood_splits.py.

- Clean sequences in **map mode**: B→N,Z→Q,J→L,U→C,O→K; other nonstandard symbols→A.
- idr_ratio is the **fraction of all sequence residues** with disorder score>0.5 (residue-ratio).
- **OOD_IDR =idr_ratio>0.3**.
- Equivalent command: idr_ood.py --mode residue-ratio --threshold 0.3 --clean map, matching ss/func/fold.
- Default cache output/ood/idr_predictions.csv, keyed by cleaned sequence, is not retained; only output/ood/idr_fastas/ remains. Default=False idr_ratio is empty and IDR flags NaN.

### Step 7: New EC function OOD, Default only

**Canonical scripts:** `.skills/ec-function-ood-annotation/scripts/`: ec_annotate.py +ec_ood.py, batch driver add_unified_newec.py.

- Source: SIFTS chain data/sifts/sifts_chain_ec.tsv.gz; first two unique_id components →PDB_CHAIN. Train+val has3,380 unique ECs.
- Reference=train+val ECs. **All-not-in:** annotated and all ECs at the selected level absent →True; unannotated →False. Partial `-` allowed:3.4.22.-→L3=3.4.22.
- OOD_NewEC_L4 recomputed in place by canonical driver; OOD_NewEC_L3 added during unification.
- EC-annotated Default samples: **6,307/13,723 (46.0%)**. Both columns NaN on Default=False.
- Legacy EC logic in22_add_newfunction_longtail_ood.py is replaced by the skill; script archived.

### Step 8: Long-tail EC OOD, Default only

**Canonical driver:** `.skills/ec-function-ood-annotation/scripts/add_unified_longtail_ec.py`. Original single-column26_add_longtail_ec_ood.py archived.

- Columns: OOD_LongTail_EC_L3_le5 /L3_le10 /L4_le5 /L4_le10.
- Annotated samples whose **all ECs** at that level have **training** count<5 /<10, **including zero**, →True (all-not-in); unannotated →False.
- Source/keys as Step7; all four columns NaN on Default=False.

### Step 9: Binding-site-ratio long tails, Default only

**Script:** scripts/11_build_ood_splits.py.

- Parse residue-aligned binary label; binding_site_ratio=binding residues/chain length.
- Ratio<0.01 sets OOD_LongTail=True. Default=False ratio empty, flag NaN.

### Step 10: CATH+TED domain holdouts, Default only

**Wrapper:** scripts/12_ted_augment_cath_ood.py calls `.skills/cath-ted-domain-annotation/scripts/cath_ted_annotate.py` and cath_ted_holdout.py for annotation/holdouts.

- Chain unit: first two unique_id components identify PDB/chain, e.g.10BL_A_ADP_A_503→10bl_A. Different ligands on one chain share annotations.
- Union experimental CATHv4.4 name-matched cath-domain-list.txt labels and putative TED REST labels cached in data/ted/ted_api_cache.json. Retain domains overlapping chain UniProt intervals by≥50% domain length; H provides C.A.T.H,T only C.A.T.
- Reference=**training** CATH+TED union, symmetric augmentation across splits. Any unseen topology/superfamily →True; unannotated →False; Default=False →NaN.

Scale and coverage:

- Chain list output/lbs_ted_aug_chain_list.csv:86,508 chains, historically generated once;25,336 unique UniProts. All TED cache hits, offline;3,902 chains lack SIFTS UniProt mappings.
- Sources: both41,953 /ted-only29,238 /cath-only4,674 /none10,643.
- Topology coverage: CATH-only53.9%→union**87.7%**.
- Outputs: output/lbs_chain_cath_ted_labels.csv (chain labels), output/lbs_chain_domains.csv (228,136 domain rows), output/lbs_chain_holdout_flags.csv (generated when rerunning wrapper12).

> Old CATH routes20_annotate_cath_domains.py /21_annotate_cath_domains_mmseqs2.py are superseded by CATH+TED and archived.

### Step 11: Structure-residue labels, struct_label

**Script:** `.skills/struct-label-alignment/add_struct_label.py`; dry-run first, then write.

- Last column struct_label exists in all splits. Label **every structure residue**, enumerated by first ATOM occurrence of (chain,resseq,icode) in file order, including residues without CA.
- Globally align PDB residue sequence to aa_seq; inherit label at aligned positions. Unmapped structure residues receive **-1**, ignored during training/evaluation.
- Validation: zero changes to other cells; all len(struct_label)==PDB residue counts; when count==len(aa_seq), struct_label exactly matches label. Rows/total residues with-1: train0/0,val0/0,test0/0.
- Backup: output/backup_before_struct_label/.

### Step 12: Validate NaN and recompute InD

**Tool:** `.skills/ood-annotation-toolkit/scripts/finalize_ood_columns.py --mode both`.

- **NaN:** on Default=False extreme rows, all non-length OOD dimensions (seq_Redundancy_*,TM-score_*,Orphan,IDR,NewEC_*,LongTail*,FoldHoldout,SuperfamilyHoldout) are empty/NaN, meaning not computed. idr_ratio/binding_site_ratio also exist only for Default. Since2026-08-28, generating scripts10/11/12 and skills write NaN directly; finalizer --mode nan is only a fallback no-op on current data.
- **InD:** all OOD_*,seq_Redundancy_*,TM-score_* False (NaN treated as False) →True, recomputed by finalizer --mode ind.
- Validator: `.skills/ood-annotation-toolkit/scripts/check_splits.py`.

### Measured current OOD statistics: percentages use Default n=13,723

| OOD column | Definition | True | Percentage |
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
| OOD_Orphan | No significant hit, e-value>1e-3 | 1,310 | 9.5% |
| OOD_IDR | Disordered fraction>30% | 140 | 1.0% |
| OOD_NewEC_L4 | All L4 ECs absent from train+val | 271 | 2.0% |
| OOD_NewEC_L3 | All L3 ECs absent from train+val | 19 | 0.1% |
| OOD_LongTail | Binding-site ratio<1% | 330 | 2.4% |
| OOD_FoldHoldout | Any topology absent from train | 15 | 0.1% |
| OOD_SuperfamilyHoldout | Any superfamily absent from train | 135 | 1.0% |
| OOD_LongTail_EC_L3_le5 | All L3 training counts<5 | 39 | 0.3% |
| OOD_LongTail_EC_L3_le10 | All L3 training counts<10 | 68 | 0.5% |
| OOD_LongTail_EC_L4_le5 | All L4 training counts<5 | 803 | 5.9% |
| OOD_LongTail_EC_L4_le10 | All L4 training counts<10 | 1,142 | 8.3% |

- InD, all OOD flags False: **5,209/17,225 (30.2% of test)**.
- Unannotated CATH+TED samples, both holdouts False:4,234 Default samples.
- All non-length OOD flags are computed only for Default;3,502 extreme rows receive NaN.

### Validation and quality control

| Check | Result |
|--------|------|
| Pairwise split unique_id intersections | 0 |
| Default pdb_id overlap with train/val | 0 |
| Extreme test rows sharing train pdb_id, Note5 | 516 rows /284 PDBs |
| Extreme test rows sharing val pdb_id | 284 rows /127 PDBs |
| label/aa_seq length mismatch | 0 |
| Default+ExtremeShort+ExtremeLong equals total test | 17,225 ==17,225 |
| Non-length OOD flags NaN on extreme rows | Yes,3,502 rows |
| Orphan⊆seq30⊆…⊆seq90 | Yes |
| TM0.3⊆…⊆TM0.9 on Default | Yes |
| Missing SMILES train/val/test | 0/0/0 |
| Split samples missing structures | 0 |
| aa_seq/structure sequence agreement | 100.00% |

---

## Output files

### Directory layout

```
datasets/ligand_binding_site/
├── DATA_PROCESS.md              # This document
├── data/                        # Task data, currently empty; shared sources/SIFTS in project data/
├── logs/                        # Run logs
├── pdbs/                        # 96,979 extracted chain PDBs
├── splits/                      # Final train/val/test CSVs only
│   ├── ligand_binding_site_train.csv
│   ├── ligand_binding_site_val.csv
│   └── ligand_binding_site_test.csv
├── scripts/                     # Numbered processing scripts
│   ├── 01_extract_binding_sites.py            # Stage1 batch extraction
│   ├── extract_binding_sites_core.py          # Imported parsing/binding-site library
│   ├── 02_aggregate_batches.py                # Stage1 aggregation
│   ├── 03_dedup_by_ligand_95.py               # Stage3 ligand-wise deduplication
│   ├── 04_exact_cleanup_dedup95.py            # Stage3 exact cleanup
│   ├── 05_create_time_based_splits.py         # Stage4 temporal split
│   ├── 06_filter_test_by_ligand_similarity.py # Stage4 chemical filtering
│   ├── 07_extract_protein_chains.py           # Stage3 fast structure extraction
│   ├── 08_add_ligand_smiles_ecfp4.py          # Stage3 ligand features
│   ├── 09_drop_missing_smiles_and_rebuild.py  # Stage3 missing removal/rebuild
│   ├── 10_rebuild_splits_with_extreme_subset.py # Stage4 extreme regrouping
│   ├── 11_build_ood_splits.py                 # Stage4 Foldseek/IDR/LongTail; mmseqs superseded
│   ├── 12_ted_augment_cath_ood.py             # Stage4 CATH+TED wrapper
│   └── 13_verify_seq_structure_consistency.py # Diagnostic/QC sequence–structure validation
├── output/                      # Intermediates, caches, backups
│   ├── ligand_binding_sites_raw.csv / _std.csv / _std_dedup95.csv
│   ├── binding_sites_batch_*.csv              # Stage1 batches
│   ├── ood/                                   # mmseqs/Foldseek/IDR caches
│   ├── intermediate_splits/                   # Historical intermediate splits
│   ├── lbs_chain_cath_ted_labels.csv / lbs_chain_domains.csv / lbs_chain_holdout_flags.csv
│   ├── backup_*/                              # Before split writes
│   ├── len2000_removal_backup/                # Before >2000 removal
│   └── scripts_archive/                       # One-off/superseded scripts
```

### Key outputs

**Final splits:** datasets/ligand_binding_site/splits/.

| File | Rows | Format | Description |
|------|------|------|------|
| ligand_binding_site_train.csv | 71,093 | CSV | Final train, seven columns, no OOD |
| ligand_binding_site_val.csv | 7,938 | CSV | Final val, seven columns, no OOD |
| ligand_binding_site_test.csv | 17,225 | CSV | Six core +29 OOD/annotation +struct_label =36 columns |

**Structures:** datasets/ligand_binding_site/pdbs/.

- 96,979 unique_id.pdb files, each containing ATOM records of its corresponding protein chain only.
- All96,256 current split samples have structures, with exact aa_seq/structure agreement.
- 723 files correspond to discarded missing-SMILES samples, absent from splits; retained for reference.

**Intermediates:**

- output/ligand_binding_sites_raw.csv /ligand_binding_sites_std.csv /ligand_binding_sites_std_dedup95.csv: main Stage1/3 data.
- output/intermediate_splits/: historical dedup95/ligsim60 intermediates, temporal mappings, etc.
- output/ood/: mmseqs_max_identity_easysearch.csv,foldseek_max_tmscore.csv,idr_predictions.csv caches, query keys=chain_key.

#### Train/val fields: seven columns

| Field | Type | Description |
|------|------|------|
| unique_id | string | {pdb_id}_{chain_id}_{ligand_name}_{ligand_chain}_{ligand_resnum} |
| aa_seq | string | Protein chain amino-acid sequence |
| struct_file | string | Filename only, {unique_id}.pdb in pdbs/ |
| ligand_smiles | string | Canonical ligand SMILES; missing samples removed |
| `ligand_ecfp4` | string | Sparse Morgan ECFP4 indices, radius=2, nBits=2048, e.g. `[12,45,200]` |
| label | string | **Required string** of residue-aligned binary binding labels, bracketed and single-space-separated, e.g.[0 0 1 0] |
| struct_label | string | Aligned to **PDB residue order**, Stage4 Step11; unmapped residues=-1 |

#### Additional test fields: current test CSV header

| Field | Type | Description |
|------|------|------|
| `Default` | bool | Normal-length test, 60 ≤ length ≤ 1000 |
| `OOD_ExtremeShort` / `OOD_ExtremeLong` | bool | Length <60 / >1000, mutually exclusive with Default |
| `seq_Redundancy_90` ~ `seq_Redundancy_30` | bool / NaN | Maximum pident below the corresponding threshold; Default only, extreme rows NaN |
| `TM-score_0.9` ~ `TM-score_0.3` | bool / NaN | Maximum TM-score below the corresponding threshold; Default only, extreme rows NaN |
| OOD_Orphan | bool/NaN | No significant mmseqs2 train+val hit |
| idr_ratio | float/empty | Disordered fraction of full sequence, Default only |
| `OOD_IDR` | bool / NaN | `idr_ratio` > 0.3 |
| `OOD_NewEC_L4` / `OOD_NewEC_L3` | bool / NaN | All ECs at the corresponding level absent from train+val, all-not-in |
| binding_site_ratio | float/empty | Binding residues/chain length, Default only |
| `OOD_LongTail` | bool / NaN | `binding_site_ratio` < 0.01 |
| `OOD_FoldHoldout` / `OOD_SuperfamilyHoldout` | bool / NaN | Any CATH+TED topology / superfamily absent from training |
| `OOD_LongTail_EC_L3_le5` / `_L3_le10` / `_L4_le5` / `_L4_le10` | bool / NaN | All EC training counts at the corresponding level <5 / <10, including zero |
| InD | bool | All OOD_*,seq_Redundancy_*,TM-score_* False, NaN treated as False; recomputed by finalizer --mode ind |

---

## Diagnostics and quality control

| Script | Purpose |
|------|------|
| scripts/13_verify_seq_structure_consistency.py | Final QC: compare split aa_seq against extracted pdbs/*.pdb sequences; current agreement100.00% |

---

## Archived scripts: output/scripts_archive/

| Script | Reason |
|------|----------|
| 02_resume_extract_binding_sites.py | One-off resume |
| 03_diagnose_remaining_files.py | Crash recovery, classify OK/Timeout/Error/Crash |
| 05_reprocess_timeouts.py | Retry timed-out files; first40 smallest still failed, stopped with zero output |
| 06_create_splits.py | Random80/10/10 replaced by temporal05 |
| 11_extract_protein_chains.py | Slow extraction replaced by fast07 |
| 13_move_intermediate_files.py | One-off intermediate migration |
| 15_patch_inconsistent_from_structure.py | Applied sequence–structure repairs; root fix in core/07 |
| 16_reextract_inconsistent_biopython.py | Same, Biopython re-extraction |
| 20_annotate_cath_domains.py /21_annotate_cath_domains_mmseqs2.py | Old CATH route replaced by CATH+TED12 and skill |
| 22_add_newfunction_longtail_ood.py | EC replaced by canonical skill; LongTail duplicates11 |
| 23_restore_idr_ood.py | One-off IDR disaster recovery |
| 24_fix_seq_homology_easysearch.py | One-off easy-search -s 7 unification; use homology skill for recomputation |
| 26_add_longtail_ec_ood.py | Original one-column version replaced by four all-not-in EC columns |
| 26_fix_multimodel_pdbs.py | Applied re-extraction of549 multimodel NMR files, root fix in07; internal import_module('11_extract_protein_chains_fast') uses the pre-rename module, now07_extract_protein_chains.py |
| fix_lbs_aaseq.py | Undocumented one-off applied aa_seq repair; stale pdb_file/labels columns |

---

## Notes

1. **PDB version:** local project-root data/pdb/ mirror; record download time for reproducibility.
2. **Sites:** residues within4.5 Å of ligand heavy atoms; exclude water,nucleic acids,standard amino acids,and127 additives (Stage2).
3. **IDR:** residue-ratio>0.3,clean=map. Script11 also computes longest contiguous disorder fraction but does not use it for flagging. Old descriptions of contiguous disorder>30% were inaccurate; use Stage4 Step6.
4. **TM unification2026-08-28:** all seven0.9–0.3 levels match the skill. Easy-search cache replaces search -a -e inf, whose21% invalid>1 scores severely undercounted flags. Details in Step5/output/foldseek_easysearch_compare/compare_report.txt.
5. **PDB overlap scope:** initial temporal splits have no pdb_id overlap. Moving extreme chains from train/val to test introduces516/284 extreme test rows sharing train/val PDBs (different chains or ligands). **Default remains PDB-disjoint from train/val**, preserving regular ID/OOD evaluation.
6. **NA chain IDs:** some chains literally use "NA". Read CSV with keep_default_na=False; writing scripts use pd.read_csv(dtype=str,keep_default_na=False), then to_csv(index=False,lineterminator="\n").
7. **ECFP4 storage:** sparse-list strings, e.g.[12,45,200]. Decode with ast.literal_eval or json.loads, then convert to NumPy.
8. **Search direction:** sequence/structure OOD uses **Default test versus train+val**. Validation participates in references, not OOD evaluation. Extreme samples are excluded from these calculations and receive NaN.
9. **Timeouts:** skip1,099 large files exceeding10/60 s, approximately0.4% of PDBs. Their potential samples are absent.
10. **Redundant pdbs/:**723 structures belong to removed missing-SMILES samples; filter by struct_file for exact split alignment.
11. **>2000 removal, no script:**155 test ExtremeLong samples,2,003–4,584 aa, removed project-wide; backup output/len2000_removal_backup/.
12. **Applied NMR repair2026-08-28:**549 PDBs concatenated approximately20–25 models without separators, duplicating residues and overflowing99,999 atom serials. Seven yielded phantom residues and -1 struct_labels. Root fixes in07 stop at first ENDMDL, retain CIF model1, and unquote atom names. Archived26_fix_multimodel_pdbs.py re-extracted affected files; originals in output/backup_before_multimodel_fix/. Recomputed struct_label has zero -1.
13. **Script adaptations2026-08:**09/10/11 migrated pdb_file/labels/file_type to struct_file/label; fixed undefined ROOT in10 and stale backup path ROOT/"ligand_binding_site" in09 to task output/. Script09 normalizes boolean strings on CSV readback. These historical one-off steps should be replayed using the Stage4 skill chain.
