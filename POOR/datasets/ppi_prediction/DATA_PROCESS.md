# PPI Prediction Dataset Processing

> Complete, reproducible workflow from raw PPI data to the final protein–protein interaction dataset.
> Location: `datasets/ppi_prediction/`
> Source: [PINDER](https://github.com/pinder-org/pinder), release `2024-02`
> Reference: Kovtun et al., *PINDER: The Protein INteraction Dataset and Evaluation Resource*, 2024.

## Task and current scale

**Binary protein–protein interaction prediction:** each sample is a protein A/B pair, labeled 1 for interaction or 0 for non-interaction.

Measured current splits:

| File | Rows | Composition: positive bio / negative XTAL / negative random |
|------|------|------|
| `ppi_train.csv` | **58,968** | 29,484 /10,301 /19,183 |
| `ppi_val.csv` | **8,394** | 4,197 /2,104 /2,093 |
| `ppi_test.csv` | **16,264** | 8,132 /3,947 /4,185; 37 columns including OOD |

Test: Default=True 14,939; Default=False 1,325, triggered by OOD_ExtremeLong when either chain exceeds 1,000 residues. Positive/negative balance is exactly 1:1.

> Run scripts from datasets/ppi_prediction/. Superseded scripts remain in output/scripts_archive/ (see archive list).

---

## Stage 1: Obtain and parse PINDER indices and sequences

> Obtain PPI annotations and all protein sequences without downloading the complete ~700-GB PINDER structure collection.

### 1.1 Sources

Location: `datasets/ppi_prediction/data/`.

| File | Purpose | Size | Description |
|------|------|------|------|
| `data/index.parquet` | Main entry point | ~119 MB | Dimer metadata and split information |
| `data/metadata.parquet` | Detailed annotations | ~100 MB | Methods, resolution, interface statistics, labels |

**Raw PINDER 2024-02:** 2,319,564 index records; ~62,706 unique PDB IDs; 47,098 unique UniProt IDs; 1,180,362 BIO positives; 873,090 XTAL negative candidates.

| Field | Type | Description |
|------|------|------|
| `id` | string | PINDER system ID, e.g. `7rzb__A1_A0A229LVN5--7rzb__A2_A0A229LVN5` |
| `pdb_id` | string | Four-character PDB code |
| `uniprot_R` / `uniprot_L` | string | Receptor / ligand UniProt IDs |
| `cluster_id` | string | Interface cluster ID for deduplication |
| `split` | string | Original train/val/test/invalid split; not retained here (Stage 4 temporal splitting) |
| `label` | string | BIO biological interaction; XTAL crystal contact candidate negative |
| `complex_type` | string | homomer / heteromer |
| `length 1` / `length 2` | int | Receptor / ligand lengths |

Manual download:

```bash
aria2c -x 8 -s 8 -o data/index.parquet     "https://storage.googleapis.com/pinder/2024-02/index.parquet"
aria2c -x 8 -s 8 -o data/metadata.parquet  "https://storage.googleapis.com/pinder/2024-02/metadata.parquet"
```

### 1.2 Fetch UniProt sequences

**Script:** `scripts/01_fetch_uniprot_sequences.py`.

- Optionally preload shared `<repo>/data/uniprot/uniprot_sprot.fasta`; skip automatically if absent (not currently supplied).
- Batch-download missing sequences through UniProt REST API `rest.uniprot.org/uniprotkb/stream`, supporting resumption.
- Output: `data/uniprot_sequences_all.fasta`.
- Initial API coverage: 45,368 /47,098 IDs (96.33%); recover the remaining 1,730 in the next step.

### 1.3 Recover missing sequences from monomer PDBs

**Script:** `scripts/05_fetch_missing_sequences_from_pdb.py`.

- Identify missing sequences from `output/ppi_pairs_corrected_seq95.csv`.
- Batch-download 3,566 corresponding PINDER holo monomers into `data/pdb_structures_missing/` with aria2c.
- Parse ATOM records with Bio.PDB and append recovered sequences to the FASTA.
- Missing means `p not in seqs or not seqs[p]`: empty records also count as missing.
- Cascade rerun 04_build_pairs_with_corrected_ranking.py and 09_rebuild_with_xtal_negatives.py.

Recovered 1,691 proteins. Only one final sample is filtered for missing sequence: **99.999% coverage**.

> Applied repair, 2026-08-28: training pair `4p5x__A1_O75746--4p5x__A2_O75746` had empty aa_seq1/aa_seq2 because the FASTA contained the only empty record, >O75746, from an empty API response. The old check tested key presence only. Fill both FASTAs with current 678-aa UniProt SLC25A12, update the two train cells, and check empty sequences in script 05. Backup: output/backup_before_o75746_fix/. All splits now have zero empty-chain samples.

---

## Stage 2: Pair labels, representatives, and EC cache

> Deduplicate unordered pairs, preserve label provenance, select representatives, and cache UniProt EC annotations. Execution order is 02→03 (Stage 3)→04; this section groups by topic.

### 2.1 Label sources

PINDER metadata.label supplies biological relevance:

| Label | Meaning | Use |
|------|------|------|
| BIO | PRODIGY-cryst predicts biological interaction | Positive, label=1 |
| XTAL | Predicted nonphysiological crystal contact | Hard negative, label=0 |

Construct negatives and balance 1:1 during Stage 4 Step 1.

### 2.2 Deduplicate pairs across all original splits

**Script:** `scripts/02_build_unique_protein_pairs.py`.

- Combine train/val/test/invalid.
- Deduplicate unordered {protein_A_id, protein_B_id} pairs.
- Prefer BIO if a pair has both BIO and XTAL annotations.
- Output: `output/ppi_unique_pairs.csv`, 74,779 rows.

### 2.3 Corrected representative selection

**Script:** `scripts/04_build_pairs_with_corrected_ranking.py`.

Priority order:

1. Monomer availability: high-quality apo >low-quality apo >predicted >holo only.
2. Better effective resolution: use holo complex resolution as a proxy because PINDER metadata does not directly supply apo resolution or predicted confidence (see notes).
3. Greater total length.
4. Latest date.

Outputs: `output/ppi_pairs_corrected_uniprot.csv` (74,779), `output/ppi_pairs_corrected_seq95.csv` (65,214; **main downstream input**).

### 2.4 Cache UniProt EC numbers

**Script:** `scripts/08_fetch_uniprot_ec.py`.

- Batch-fetch semicolon-separated EC numbers through UniProt REST.
- Input: `data/uniprot_sequences_all.fasta`; output: `output/uniprot_ec_mapping.csv`, 47,059 entries.
- Used with SIFTS in Stage 4 EC OOD.

---

## Stage 3: Sequence redundancy and structure refinement

### 3.1 Sequence-identity deduplication

**Script:** `scripts/03_dedup_by_sequence_identity.py`.

- Cluster UniProt sequences with mmseqs2 at 95% identity.
- Treat two pairs as redundant if corresponding chains fall in the same sequence clusters.
- Output: `output/ppi_unique_pairs_seq95_dedup.csv`, 63,189 rows.

### 3.2 Select and download monomers

**Script:** `scripts/06_download_structures_and_filter_test.py`.

- For both chains of each BIO/XTAL representative: **apo >predicted AlphaFold >holo**.
- Batch-download chosen PDBs using aria2c into `data/pdb_structures_all/`.
- Outputs:
  - `output/structure_selection.csv`, 65,214 representative-to-structure rows, including both_chains_apo_or_predicted.
  - `output/ppi_test_structure_filtered.csv`, diagnostic structure-filtered test data. **Final splits come from script 09**, which incorporates this filter.

Sources among 65,214 representatives:

| R source | L source | Count |
|--------|--------|------|
| apo | apo | 5,648 |
| apo | predicted | 3,864 |
| apo | holo | 405 |
| predicted | apo | 2,975 |
| predicted | predicted | 48,130 |
| predicted | holo | 642 |
| holo | apo | 216 |
| holo | predicted | 695 |
| holo | holo | 2,639 |

R: predicted 51,747 /apo 9,917 /holo 3,550. L: predicted 52,689 /apo 8,839 /holo 3,686. Both chains apo/predicted: 60,617 /65,214. All 50,306 selected structure files downloaded successfully; see output/structure_download_summary.json.

### 3.3 Stage structures into pdbs/

**Script:** `scripts/07_stage_pdbs.py`, added 2026-08 to formalize a previously unscripted step.

- Required files: union of struct_file1/struct_file2 in current splits (44,178; fields assigned by 09 using structure_selection.csv), plus 2,224 historical pre-XTAL references in output/splits_pre_xtal_backup/. Disable historical inclusion with --no-pre-xtal.
- Idempotent: existing destination→skip; existing source→copy; missing source→record.
- Audit, 2026-08: dry-run expected list (46,402) and pdbs/ agree in both directions; sampled copies match existing files byte for byte.
- Currently data/pdb_structures_all/ retains 10,690 files; others were moved out. Fresh reproduction requires rerunning 06 to download all 50,306. The 3,904 selected files never referenced by any split version are excluded from staging.

### 3.4 Structure integrity diagnostics

**Script:** `scripts/diag_validate_structures.py`, outside the main pipeline.

- Parse downloaded PDBs with Bio.PDB, checking zero-size files, ATOM residues, and parsability.
- Outputs: output/structure_validation_summary.json and output/structure_validation_failed.csv.
- All **50,306** downloaded files pass: empty=0, no_atoms=0, parse_error=0.

---

## Stage 4: Splitting and OOD annotation

### Step 1: Rebuild splits with XTAL negatives and temporal splitting

**Script:** `scripts/09_rebuild_with_xtal_negatives.py`.

Construct negatives independently **after splitting** within each train/val/test pool:

1. **XTAL hard negatives:** retain all PINDER XTAL pairs. They resemble interfaces but are nonbiological, from real complexes. Protein overlap with BIO is 65.2%, reducing the random-pair versus real-complex shortcut.
2. **Random completion:** when XTAL is insufficient, sample two proteins from the split's entire deduplicated BIO+XTAL pool. **Self-pairs are allowed**, simulating negative homodimers and removing the seq1==seq2→label1 shortcut. Exclude all known BIO/XTAL pairs in output/ppi_unique_pairs.csv.
3. **1:1 balance:** total XTAL+random negatives equal positives within each split.
4. **pair_type:** all splits annotate bio, xtal, or random.

| Rule | Description |
|------|------|
| Temporal split | Sort BIO by release_date; split at cumulative 70%/80% (70:10:20). Apply the same date boundaries to XTAL. |
| Global hard length filter | Remove any pair with a chain >2,000 before splitting. |
| Train/val length | Chains ≤1,000; move longer pairs into extreme-long test. |
| Test structure filter | Require apo or predicted structures for both chains; applied inline. |
| Balance | BIO : (XTAL+random) =1:1 per split. |

Inputs: output/ppi_pairs_corrected_seq95.csv, output/structure_selection.csv, output/ppi_unique_pairs.csv, data/uniprot_sequences_all.fasta. Random seeds: 42/142/242. Set test Default=False when either chain exceeds 1,000. Outputs: splits/ppi_{train,val,test}.csv, output/xtal_rebuild_summary.json.

> Historically, 21_resample_negative_pairs.py generated random-only negatives without self-pairs, creating homodimer and random/real-complex shortcuts and inflated AUROC. Iterations 21→23→09 introduced XTAL and self-pairs. Old scripts are archived; original splits remain in output/splits_pre_xtal_backup/.

### Step 2: Recompute local OOD logic

**Script:** `scripts/10_recompute_all_ood.py`, operating on final splits.

- ExtremeShort: either chain <60; 715 Default samples. Another 16 short+long mixed pairs are Default=False and this flag becomes NaN.
- ExtremeLong: either chain >1,000; 1,325 samples, all Default=False. Default=True: 14,939.
- Extreme lengths are flagged, not removed here; Stage 4 Step 1 already restricts train/val to ≤1,000.
- Default=False: preserve only ExtremeLong, the pre-reservation reason. All other OOD, sequence, TM, and ExtremeShort flags become NaN, meaning not computed. Backup: output/backup_before_nondefault_full_nan/.
- Inline mmseqs2/Foldseek/IDR/EC sections are **superseded** by Steps 3–6.

### Step 3: Sequence homology for Default samples

**Canonical tool:** `.skills/homology-ood-annotation/scripts/seq_homology_ood.py`.

- Reference: all unique train+val protein sequences; query: test proteins. mmseqs2 easy-search -s 7.5 with default significance e-value ≤1e-3.
- Cache: output/ood_mmseqs_test_vs_trainval_easysearch.m8, skill keys unique_id|aa_seq{1,2}; replay via --m8-cache.
- seq_Redundancy_XX: XX=90..30, step 10; either chain maximum identity <XX→True.
- Orphan: either chain lacks significant train+val hits→True.

> Legacy mmseqs search --max-seqs10 -e 100 -c 0.0 underestimated low-threshold OOD (seq30 515, Orphan 211). Recomputed with canonical easy-search 2026-08-27. Backup: output/backup_before_easysearch_unify/; old output/ood_mmseqs_test_vs_trainval.m8 is provenance only. Local script 10 section is superseded.

### Step 4: Structural similarity for Default samples

**Canonical tool:** `.skills/homology-ood-annotation/scripts/struct_homology_ood.py`.

- Query: test monomer structures in pdbs/; target: train+val monomers.
- Foldseek easy-search, default e-value ≤10, output query, target, alntmscore.
- Cache: output/ood_foldseek_test_vs_trainval_easysearch.m8, replay via --m8-cache.
- TM-score_X.X: thresholds 0.9..0.3, step 0.1; either chain maximum alntmscore below threshold→True.

> Legacy search -a -e inf retained weak hits, inflating maximum alntmscore and altering OOD counts: TM 0.9 4,827→8,172; TM 0.3 1,610→943 after unification on 2026-08-27. Backup: output/backup_before_foldseek_easysearch_unify/; old output/ood_foldseek_tmscore.csv is provenance only. Local script 10 section is superseded.

### Step 5: Intrinsic disorder for Default samples

**Tool:** `.skills/ood-annotation-toolkit/scripts/idr_ood.py`, metapredict v 3 with --mode region-ratio --clean map --threshold 0.3.

- Predict disorder for all unique test sequences; cache output/idr_predictions.csv (7,418 sequences).
- Map nonstandard residues B→N, Z→Q, J→L, U→C, O→K, other→A.
- A chain is IDR-containing when its longest contiguous disorder>0.5 region occupies >30% of full length (**region-ratio**).
- Either chain IDR-containing→OOD_IDR=True.
- Local script 10 IDR section is superseded.

### Step 6: EC OOD for Default samples

**Canonical drivers:** `.skills/ec-function-ood-annotation/scripts/add_unified_newec.py` and add_unified_longtail_ec.py.

Collect EC through two sources:

- UniProt REST cache output/uniprot_ec_mapping.csv (47,059 records), from script 08.
- SIFTS `<repo>/data/sifts/sifts_chain_ec.tsv.gz`, by PDB+chain.
- bio/xtal: parse unique_id into pdb/chain/uniprot, union cached UniProt and SIFTS ECs. Random negatives: UniProt cache only.
- L 3 uses first three components when ≥3 exist; L 4 requires four. Partial '-' is allowed: 3.4.22.-→L 3 3.4.22.

| Columns | Reference | Definition |
|------|--------|------|
| OOD_NewEC_L4/L 3 | Train+val EC union | Annotated and all ECs at the level absent→True; unannotated→False. |
| OOD_LongTail_EC_L3_le5/le10, L 4_le5/le10 | Training per-sample frequencies | Annotated and every EC frequency <5/<10, including zero→True; unannotated→False. |

> Legacy OOD_LongTail_EC_L4 (train+val frequency<10, any semantics, 655) was removed and replaced by four columns. NewEC_L4 changed any 245→all-not-in 203. Backups: output/backup_before_newec_unify/, backup_before_longtail_unify/, backup_before_longtail_allnotin/. Local script 10 EC section is superseded.

### Step 7: CATH+TED holdouts for Default samples

**Wrapper:** scripts/11_ted_augment_cath_ood.py: annotate→aggregate pairs→canonical cath_ted_holdout.py→write splits.

UniProt mode; **pair-level labels union both protein_A_id/protein_B_id**:

- Union CATH v 4.4 experimental labels reverse-mapped through SIFTS with putative TED full-length annotations, without overlap filtering. REST cache: `<repo>/data/ted/ted_api_cache.json`. H level provides C.A.T.H; T level C.A.T only.
- Reference: training CATH+TED union, with symmetric train/val/test augmentation: 1,262 topologies /4,711 superfamilies.
- Any unseen topology/superfamily in either partner→True; both unannotated→False; Default=False→NaN.

Workflow:

1. `cath_ted_annotate.py --chains output/ppi_ted_aug_acc_list.csv --out output/ppi_acc_cath_ted_labels.csv --offline`; all 39,455 accessions cached.
2. Union both partner labels into output/ppi_pair_cath_ted_labels.csv: unique_id, topos, sfs, split.
3. `cath_ted_holdout.py --labels <pair-label-table> --split-col split --out output/ppi_pair_holdout_flags.csv`.
4. Write splits/ppi_test.csv; Default=False→NaN.

Sources among 39,455 accessions: both 19,756 /ted-only 14,066 /cath-only 1,365 /none 4,268. Topology coverage rises from CATH-only 53.5% to **89.2%**.

Default test n=14,939: fold44, superfamily265, both unannotated 950 (False).

Outputs: accession labels output/ppi_acc_cath_ted_labels.csv; historical one-off key list output/ppi_ted_aug_acc_list.csv; pair labels/flags above; backup output/backup_before_ted_aug/.

### Step 8: NaN convention and InD

**Tool:** `.skills/ood-annotation-toolkit/scripts/finalize_ood_columns.py`.

- Default=False: every OOD field is NaN except ExtremeLong, distinguishing not-computed from False. This includes ExtremeShort, sequence, and TM. Backups: output/backup_before_nondefault_nan/ and backup_before_nondefault_full_nan/.
- InD: all OOD/sequence/TM flags False (NaN treated as False), recomputed with --mode ind after all annotations. **1,944/16,264 (12.0%)**; backup output/backup_before_ind/.
- Validate using check_splits.py in the same toolkit.

### Measured test OOD counts

16,264 samples: positive 8,132 /negative 8,132; Default 14,939 /nondefault1,325. Unless noted, flags below apply only to Default; nondefault values are NaN.

| Column | True | Notes |
|--------|--------:|------|
| seq_Redundancy_90 | 10,655 | Default only |
| seq_Redundancy_80 | 9,856 | Default only |
| seq_Redundancy_70 | 9,027 | Default only |
| seq_Redundancy_60 | 8,098 | Default only |
| seq_Redundancy_50 | 6,982 | Default only |
| seq_Redundancy_40 | 5,496 | Default only |
| seq_Redundancy_30 | 3,653 | Default only |
| OOD_Orphan | 3,017 | Default only |
| TM-score_0.9 | 8,172 | Default only |
| TM-score_0.8 | 5,345 | Default only |
| TM-score_0.7 | 3,771 | Default only |
| TM-score_0.6 | 2,843 | Default only |
| TM-score_0.5 | 2,186 | Default only |
| TM-score_0.4 | 1,510 | Default only |
| TM-score_0.3 | 943 | Default only |
| OOD_ExtremeShort | 715 | Default only; 16 mixed pairs set to NaN |
| OOD_ExtremeLong | 1,325 | Nondefault count; sole retained pre-reservation reason |
| OOD_IDR | 3,743 | Default only |
| OOD_NewEC_L4 | 203 | Default only |
| OOD_NewEC_L3 | 10 | Default only |
| OOD_LongTail_EC_L3_le5 | 17 | Default only |
| OOD_LongTail_EC_L3_le10 | 72 | Default only |
| OOD_LongTail_EC_L4_le5 | 444 | Default only |
| OOD_LongTail_EC_L4_le10 | 698 | Default only |
| OOD_FoldHoldout | 44 | Default only |
| OOD_SuperfamilyHoldout | 265 | Default only |
| InD | 1,944 | 12.0% |

---

## Output files

```
datasets/ppi_prediction/
├── DATA_PROCESS.md                      # This document
├── scripts/
│   ├── 01_fetch_uniprot_sequences.py    # Stage 1: UniProt API sequences
│   ├── 02_build_unique_protein_pairs.py # Stage 2: combine/deduplicate pairs across original splits
│   ├── 03_dedup_by_sequence_identity.py # Stage 3: mmseqs2 95% identity
│   ├── 04_build_pairs_with_corrected_ranking.py  # Stage 2: corrected ranking
│   ├── 05_fetch_missing_sequences_from_pdb.py    # Stage 1: PDB sequence recovery; reruns 04/09
│   ├── 06_download_structures_and_filter_test.py # Stage 3: select/download monomers, filter test
│   ├── 07_stage_pdbs.py                 # Stage 3: staging; formalized 2026-08
│   ├── 08_fetch_uniprot_ec.py           # Stage 2: EC cache
│   ├── 09_rebuild_with_xtal_negatives.py  # Stage 4: XTAL +random negatives, temporal splits
│   ├── 10_recompute_all_ood.py          # Stage 4: length/Default/NaN active; other stages superseded
│   ├── 11_ted_augment_cath_ood.py       # Stage 4: CATH+TED holdouts
│   └── diag_validate_structures.py      # Diagnostic only
├── data/                                # Raw data and task caches
│   ├── index.parquet / metadata.parquet # PINDER 2024-02
│   ├── uniprot_sequences_all.fasta      # All sequences
│   ├── pdb_structures_all/              # Downloaded monomers: historically 50,306, now10,690 (§3.3)
│   └── pdb_structures_missing/          # 3,566 holo monomers for sequence recovery
├── output/                              # Intermediates, caches, backups
│   ├── structure_selection.csv          # 65,214 representative-to-monomer mappings
│   ├── uniprot_ec_mapping.csv           # 47,059 EC entries
│   ├── idr_predictions.csv              # 7,418 sequence predictions
│   ├── ood_mmseqs_test_vs_trainval_easysearch.m8     # Replayable canonical mmseqs cache
│   ├── ood_foldseek_test_vs_trainval_easysearch.m8   # Replayable canonical Foldseek cache
│   ├── ood_mmseqs_test_vs_trainval.m8 / ood_foldseek_tmscore.csv  # Legacy provenance only
│   ├── ppi_acc_cath_ted_labels.csv      # Accession-level labels
│   ├── scripts_archive/                 # Historical scripts, original filenames retained
│   ├── backup_before_*/                 # Split backups before updates
│   └── splits_pre_xtal_backup/          # Pre-XTAL splits
├── pdbs/                                # 46,402 final structure files
└── splits/
    ├── ppi_train.csv                    # 58,968 rows, nine columns, no OOD
    ├── ppi_val.csv                      # 8,394 rows, nine columns, no OOD
    └── ppi_test.csv                     # 16,264 rows, 37 columns with all OOD
```

Shared columns:

| Column | Description |
|------|------|
| unique_id | PINDER representative_id for bio/xtal; neg_xxxx for random pairs |
| protein_A_id /protein_B_id | UniProt IDs |
| aa_seq1 /aa_seq2 | Amino-acid sequences |
| struct_file1 /struct_file2 | Monomer filenames only, under pdbs/ |
| label | Pair-level 1 interaction /0 non-interaction |
| pair_type | bio /xtal /random |

Full test header (37 columns; actual file authoritative):

```
unique_id, protein_A_id, protein_B_id, aa_seq1, aa_seq2,
struct_file1, struct_file2, label, Default, pair_type,
OOD_IDR, OOD_NewEC_L4, OOD_ExtremeShort, OOD_ExtremeLong,
seq_Redundancy_90, seq_Redundancy_80, seq_Redundancy_70, seq_Redundancy_60,
seq_Redundancy_50, seq_Redundancy_40, seq_Redundancy_30, OOD_Orphan,
TM-score_0.9, TM-score_0.8, TM-score_0.7, TM-score_0.6, TM-score_0.5,
TM-score_0.4, TM-score_0.3,
OOD_FoldHoldout, OOD_SuperfamilyHoldout,
OOD_LongTail_EC_L3_le5, OOD_LongTail_EC_L3_le10,
OOD_LongTail_EC_L4_le5, OOD_LongTail_EC_L4_le10,
OOD_NewEC_L3, InD
```

---

## Quality checks

1. Pairwise split unique_id intersections are zero; random negatives exclude all known BIO/XTAL pairs in output/ppi_unique_pairs.csv.
2. Each split is exactly1:1 positive/negative. Negative XTAL fractions: train 34.9% /val 50.1% /test 48.5%, limited by period-specific availability; random pairs fill deficits.
3. Remove chains >2,000 globally before splitting; train/val chains ≤1,000; longer test chains remain as ExtremeLong/Default=False.
4. All50,306 downloaded PDBs pass Bio.PDB parsing diagnostics.
5. Sequence coverage 99.999%; one sample removed.
6. OOD updates preserve row counts, column order, and unique_id relative to backups in output/backup_before_*/.
7. Staging audit, 2026-08: expected current 44,178 plus historical 2,224 =46,402 matches pdbs/ in both directions.

## Diagnostic tools

| Script | Purpose |
|------|------|
| scripts/diag_validate_structures.py | Check empty/no-ATOM/unparseable files in data/pdb_structures_all/; write summary JSON and failures. |

## Archived scripts: output/scripts_archive/

| Scripts | Reason |
|------|----------|
| 02_build_ppi_pairs.py /03_dedup_ppi_pairs.py | Initial pairing/deduplication replaced by 02/03 |
| 06_build_pairs_with_custom_ranking.py | Experimental ranking replaced by corrected04 |
| 07_build_negative_pairs.py /08_time_split_dataset.py | Replaced by 09; script 05 cascade updated to 04→09 |
| 15_mark_ood.py /17_build_ood_splits.py /18_compute_ood_idr.py /19_compute_ec_ood.py /22_recompute_homology_newfold_ood.py | Dispersed OOD implementations replaced by 10 plus canonical skills |
| 16_fix_negative_structure_labels.py /20_add_pdb_file_columns.py | One-off column repairs already applied |
| 21_resample_negative_pairs.py | Random-only negatives with homodimer shortcut replaced by XTAL+random 09 |

## Notes

1. **PPI versus PPIS:** PPIS predicts residue interfaces; PPI predicts interaction for a protein pair.
2. **Resolution proxy:** metadata supplies holo resolution only. Ranking uses this for apo/predicted quality; actual apo resolution and AlphaFold pLDDT remain future improvements.
3. **Default:** only ExtremeLong triggers Default=False (1,325). Extreme-length samples are flagged, not removed at this step. Preserve only ExtremeLong and set all other OOD/sequence/TM/ExtremeShort fields to NaN.
4. **Homology:** since 2026-08-27, mmseqs2 and Foldseek use canonical skill easy-search. Local10 sections are superseded; old caches are provenance only.
5. **EC:** all-not-in; NewEC uses train+val, LongTail uses train frequency including 0. Legacy single LongTail_L4 (655) removed; canonical skill drivers write current columns.
6. **TED dual labels:** comma-separated foldseek/foldclass labels are unioned; no comma-containing pseudo-labels remain. All TED cache entries hit; annotation reproduces offline.
7. **Current raw-structure cache, audited2026-08:** 10,690 PDBs remain: 6,786 split references plus3,904 selected files never referenced by any split version. Most historical files moved elsewhere. pdbs/ with 46,402 files is authoritative; rerun06 before fresh staging.
8. **Historical dependency:** the complete46,402 staging list includes2,224 pre-XTAL references. Without output/splits_pre_xtal_backup/, only 44,178 current-reference files reproduce.
