# func_prediction Dataset Processing (EC / GO-MF / GO-CC / GO-BP)

> Complete linear workflow for reproducing `datasets/func_prediction/` from raw data.
> Four subtasks: **EC** (`ec`) and three **GO** branches (`go_mf` / `go_cc` / `go_bp`).
> Final outputs: 12 split files in `splits/`, one train/val/test file per subtask.

## Task definition and current scale

Start from experimental PDB structures, map them to UniProt with SIFTS, and retain the best experimental chain per UniProt. Build multilabel datasets (EC L3 `x.x.x.-`; GO terms), split by deposition date, and add multidimensional OOD annotations.

Measured current `splits/` sizes:

| Subtask | train | val | test | Default test | Test columns |
|:---|---:|---:|---:|---:|---:|
| **EC** | 15,389 | 1,687 | 5,047 | 4,243 | 31 |
| **GO-MF** | 18,179 | 2,084 | 5,609 | 4,609 | 34 |
| **GO-CC** | 13,878 | 1,730 | 5,613 | 4,256 | 34 |
| **GO-BP** | 15,232 | 1,742 | 5,142 | 4,121 | 34 |

Four processing stages plus an optional extension:

- **Stage 1:** Obtain and parse sources (SIFTS → experimental PDB structures → UniProt mapping → structural metadata and functional annotations).
- **Stage 2:** Deduplicate and refine (chain extraction → UniProt representative selection → length filtering → sequence clustering).
- **Stage 3:** Aggregate labels and split data (subtask allocation → temporal split → test filtering → label cleaning).
- **Stage 4:** Annotate OOD (merge extreme-length chains into test → annotate multiple dimensions → finalize InD/NaN).
- **Stage 5, optional:** Collect AlphaFold DB predictions for training-set expansion only.

---

## Stage 1: Obtain and parse sources (SIFTS structure–sequence–function mapping)

> Establish precise **PDB-chain-level** UniProt mappings via SIFTS and obtain structural metadata and initial function annotations. This is the principal **experimental-structure entry point**, capturing chains absent from Swiss-Prot cross-references, including TrEMBL and newly deposited structures.

**Script:** `scripts/01_parse_sifts.py` (SIFTS parsing, batch metadata queries, and functional annotation export).

### 1.1 Download core SIFTS mappings (no script)

- **Source:** EBI SIFTS (https://www.ebi.ac.uk/pdbe/docs/sifts/).
- **Core TSV flat files**, stored in project-level `POOR/data/sifts/`:
  - `sifts_chain_uniprot.tsv.gz` / `pdb_chain_uniprot.tsv.gz`: **residue-level** PDB chain ID ↔ UniProt accession mappings.
  - `sifts_chain_go.tsv.gz`: PDB chain ID ↔ GO ID and evidence codes from UniProt-GOA.
  - `sifts_chain_ec.tsv.gz`: PDB chain ID ↔ EC numbers.
  - `pdb_chain_cath.tsv.gz`: structural classification, used for Stage 4.6 OOD_FoldHoldout.

### 1.2 Extract PDB chain metadata

- **Source:** RCSB PDB GraphQL API (`https://data.rcsb.org/graphql`).
- **Queries:** batches of 50 entries/request; resumable, with checkpoints every 500 entries.
- **GraphQL fields:** `exptl.method`, `refine.ls_d_res_high` (X-ray resolution), `refine.ls_R_factor_R_work` / `refine.ls_R_factor_R_free`, and `em_3d_reconstruction.resolution` (EM resolution).
- **Resolution:** prefer X-ray `refine.ls_d_res_high`, falling back to EM `em_3d_reconstruction.resolution`.
- **Metadata coverage**, 238,401 entries: resolution 94.9%, R-work 79.7%, R-free 78.6%.

### 1.3 Obtain UniProt reference sequences

- **Purpose:** full-length UniProt sequences mapped to PDB chains, for sequence–structure alignment validation and coverage.
- **Source:** UniProt REST API, FASTA batches of 100 accessions/request.
- **Coverage:** approximately 587,000 sequences, 100%; outputs in `output/uniprot/` and `output/uniprot_reference_sequences.json`.

### 1.4 Obtain functional annotations through SIFTS

- **GO:** use SIFTS GO IDs/evidence codes directly. SIFTS GO mappings cover **all UniProtKB**, Swiss-Prot and TrEMBL.
- **EC:** obtain via SIFTS EC mappings.
- **Outputs:** `data/sifts/sifts_go_annotations.json` (822M), `data/sifts/sifts_ec_annotations.json` (43M).

### 1.5 Initial filtering: length and annotation

Filter approximately 958,000 source chains. **Defer structural-quality restrictions** (resolution, method, R-factor), applying only:

| Criterion | Threshold | Purpose |
|:---|:---|:---|
| SIFTS mapped length | `60 <= mapped_length <= 1000` | Exclude short fragments and oversized proteins |
| Function annotation | EC **or** GO required | Ensure downstream supervision |

- **Current recomputation** from `data/sifts/annotated_chains.csv`: **769,892** chains, **56,059** UniProts. Historical run: 769,730 / 56,057; differences reflect regenerated SIFTS drift (see Notes).
- **Historical methods:** X-RAY DIFFRACTION 442,283; ELECTRON MICROSCOPY 319,654; SOLUTION NMR 6,670; SOLID-STATE NMR 650; ELECTRON CRYSTALLOGRAPHY 301; other 172.

> **Rationale:** early resolution/method filtering would remove many EM structures important for large-complex representations. Apply quality criteria during later refinement.

### 1.6 Download initially selected structures

- **Source:** RCSB PDB (`https://files.rcsb.org/download/`); prefer `{pdb_id}.pdb.gz`, falling back to `{pdb_id}.cif.gz` on 404.
- **Concurrency/recovery:** 10 threads, up to five retries/file, exponential backoff for HTTP 429, checkpoints every 100 files.
- **Storage:** project-shared `data/pdb/`, approximately 255K gzip files.

### 1.7 Outputs

| File | Description |
| :--- | :--- |
| `data/sifts/sifts_pdb_chains.csv` | All candidate chains and structural metadata (41M) |
| `data/sifts/sifts_uniprot_mapping.json` | Residue-level PDB chain → UniProt mappings (154M) |
| `data/sifts/sifts_go_annotations.json` | Chain GO annotations (822M) |
| `data/sifts/sifts_ec_annotations.json` | Chain EC annotations (43M) |
| `data/sifts/annotated_chains.csv` | Candidate metadata, 827,186 rows; apply §1.5 to obtain the candidate pool |
| `output/uniprot/` | UniProt reference sequences and Swiss-Prot metadata |

---

## Stage 2: Deduplication and refinement (PDB–UniProt alignment)

> Extract chains, deduplicate by UniProt, retain **one best experimental PDB chain per UniProt**, and cluster sequences. Stage 1 keeps all candidates, including homomer copies; Stage 2 refines them.

### 2.1 UniProt-level representative selection

**Script:** `scripts/02_select_representative_chains.py`. Reimplemented to fill a gap: the historical script and `representative_chains.csv` were not migrated from the old workspace. Uses the priorities below.

**Strategy:** group by UniProt ID, retaining **one best PDB chain** using successive priorities:

| Priority | Criterion | Ordering |
|:---|:---|:---|
| 1 | Extraction mode | Tier A > tier B |
| 2 | Source quality | Swiss-Prot > TrEMBL (constant within groups; not discriminative) |
| 3 | EC availability | Annotated > unannotated |
| 4 | Experimental GO richness | Descending |
| 5 | Total GO richness | Descending |
| 6 | Resolution | Ascending |
| 7 | R-free | Ascending |
| 8 | Mapped length | Descending |
| 9 | PDB ID + chain ID | Ascending lexicographic deterministic tiebreak |

- **Candidate sources:** `data/sifts/annotated_chains.csv` (§1.5), `sifts_uniprot_mapping.json` (valid pdb_beg/pdb_end segment implies A; otherwise B), `sifts_go_annotations.json` (ijson-streamed experimental GO counts), and `output/uniprot_metadata.json` (reviewed status).
- **Output:** `output/representative_chains.csv`; historical **55,124** representatives: A 31,770 / B 23,354; Swiss-Prot 23,374 / TrEMBL 31,750; EC-annotated 24,369; experimental GO 14,100.
- **Known deviation:** historical B-tier boundary inference filtered identity ≥0.90 before ranking: 439,427 attempted, 419,370 passed. `b_boundary_inference.jsonl` was not migrated. This script admits all B-tier chains unless `--b-keep-list` supplies passing IDs. Rerun on 2026-08-28 without that list: 769,892 candidates /56,059 UniProts; **56,057** representatives; **38,835 /47,852 (81.2%)** exactly match final nonredundant representatives. Different chains were selected for 9,017 UniProts, mainly from missing B-tier filtering and A/B definition drift (Notes).

**Optional GO/EC refinement, skipped:** originally intended to require EC or experimental GO evidence (EXP/IDA/IPI/IMP/IGI/IEP/TAS/IC/HTP/HDA) for TrEMBL chains. All §1.5 candidates instead entered ranking. Consequently, 42.2% (23,258) representatives lack experimentally validated function annotations; applying this filter would reduce representatives to approximately 39,000.

### 2.2 Chain extraction: three tiers

**Script:** `scripts/03_extract_representative_structures.py`
Extract chains from source `.pdb.gz` / `.cif.gz` in `data/pdb/` to `pdbs/`. Paths are relative to the working directory, which must contain both directories.

| Tier | Condition | Extraction | Priority |
|:---|:---|:---|:---|
| **A** | Valid SIFTS `pdb_beg` / `pdb_end` | Exact official residue bounds | Highest |
| **B** | Missing pdb_beg/end, successful inference | Infer PDB bounds by aligning sp_beg/sp_end subsequences | Medium |
| **C** | Failed B inference or identity <90% | Whole-chain fallback or discard | Lowest |

- **Inputs:** `output/representative_chains.csv` (§2.1), `data/sifts/sifts_uniprot_mapping.json`, and `output/b_boundary_inference.jsonl`. The last was **not migrated**; regenerate using boundary inference below before rerunning.
- **Outputs:** `pdbs/{pdb_id}_{chain_id}.pdb` for single-character chain IDs, `.cif` for multicharacter IDs exceeding PDB field widths.
- **Formats:** stream `.pdb.gz`; parse `.cif.gz` with BioPython MMCIFParser and retain mmCIF output. Batch by pdb_id so each CIF is parsed once. Stop NMR entries at the first ENDMDL, retaining model 1; the root-cause fix is built in.
- **Extraction:** **55,124 structures**, 100% representative coverage:
  - Exact full-boundary A: 20,828 (37.8%); partly inferred A: 10,942 (19.9%).
  - Subsequence-inferred B: 23,354 (42.4%); 12 unusual chains (insertion codes/decreasing ranges) repaired manually; C: 0.
  - PDB: 47,621; mmCIF: 7,503.

**B-tier boundary inference:** locally align `UniProt[sp_beg:sp_end]` against the whole PDB chain using edlib HW and map back to residue numbers. Identity ≥0.90 retains **419,370** B candidates; exclude 19,906 low-quality chains and 151 without inference.

> **Historical note:** the old run aligned full-length UniProt, systematically inflating identity/coverage (95.5% had identity ≥0.90). Subsequence alignment yields stricter distributions; see Notes.

### 2.3 Length filtering and sequence clustering

#### 2.3.1 Count residues and retain extreme chains separately

**Script:** `scripts/04_extract_seqs_and_count_residues.py`: extract sequences and count deduplicated residues in `pdbs/`, using the standard 20 amino acids plus MSE/SEC/PYL; write to `output/`.

- **Rule:** normal length **[60,1000]** residues. **Retain shorter/longer chains separately** for Stage 4.1 test integration. `scripts/08_recover_extreme_length.py` collects initial candidates into `output/extreme_length_chains.csv`: **3,834** chains, 3,146 short +688 long, range 3–1,997. **Globally discard >2000-aa chains** for memory/cost reasons.
- **Result:** **53,902** normal chains; 1,222/55,124 (2.2%) outside range. Short-chain rates: A 1.2%, B 3.6%.

#### 2.3.2 95% sequence-identity clustering (MMseqs2, no script)

- **Tool:** MMseqs2 `mmseqs cluster --min-seq-id 0.95 -c 0.8 --cov-mode 0`.
- **Input:** `output/filtered_sequences.fa`, 53,902 normal-length sequences.
- **Result:** **47,852** clusters; 44,052 singletons (92.1%); maximum size 28.

#### 2.3.3 Within-cluster representative selection (no standalone script)

Apply §2.1 priorities within each cluster, retaining **47,852 nonredundant representatives**:

- A 28,266 (59.1%) /B 19,586 (40.9%); Swiss-Prot 22,515 (47.1%) /TrEMBL 25,337 (52.9%).
- Methods: X-RAY 33,977 (71.0%), EM 11,868 (24.8%), NMR 1,978 (4.2%), other 29 (0.1%).
- Mean resolution 2.42 Å; mean mapped length 300.5 residues; clustering reduction **11.2%**.

#### 2.3.4 Outputs

| File | Description |
| :--- | :--- |
| `output/filtered_sequences.fa` | 53,902 length-filtered sequences, approximately 11 MB |
| `output/chain_residue_counts.json` | Per-chain residue counts from 04, regenerated on rerun |
| `output/final_nonredundant_representatives.csv` | **47,852 final representatives**, including extract_mode, cluster_rep, cluster_size |
| `output/extreme_length_chains.csv` | 3,834 separately retained extreme chains, with sequences, deposition dates, and task labels |
| `pdbs/` | 50,123 chain structures: 47,852 representatives +2,271 extreme chains integrated into test; 43,567 PDB +6,556 CIF |

---

## Stage 3: Aggregate labels and split data

> Organize 47,852 representatives into EC/GO-MF/CC/BP datasets, split by **PDB deposition date**, then clean labels. Intermediates: `output/datasets_split_raw/`.

### 3.1 Assign tasks

**Script:** `scripts/06_split_datasets.py`

| Dataset | Assignment | Chains |
|:---|:---|---:|
| **EC** | `ec_count > 0` | 22,095 |
| **GO-MF** | At least one Molecular Function GO annotation | 39,646 |
| **GO-CC** | At least one Cellular Component GO annotation | 27,972 |
| **GO-BP** | At least one Biological Process GO annotation | 34,607 |

> A chain may belong to multiple datasets; split each independently. Parse namespaces from `go-basic.obo`, downloaded from Gene Ontology (Notes).

### 3.2 Temporal splitting

**Scripts:** `scripts/05_parse_deposition_dates.py` parses HEADER / `_pdbx_database_status.recvd_initial_deposition_date` from `data/pdb/` into `output/deposition_dates.json`; `scripts/06_split_datasets.py` splits ascending dates at quantiles.

| Dataset | Train:val:test | Training cutoff | Validation interval | Test start |
|:---|:---|:---|:---|:---|
| **EC** | 72 : 8 : 20 | ≤ 2019-06-25 | 2019-06-25 ~ 2021-04-29 | ≥ 2021-04-29 |
| **GO-MF** | 63 : 7 : 30 | ≤ 2019-05-19 | 2019-05-20 ~ 2020-10-07 | ≥ 2020-10-07 |
| **GO-CC** | 63 : 7 : 30 | ≤ 2020-10-14 | 2020-10-14 ~ 2021-11-22 | ≥ 2021-11-22 |
| **GO-BP** | 63 : 7 : 30 | ≤ 2020-01-15 | 2020-01-15 ~ 2021-04-20 | ≥ 2021-04-21 |

> Independent splits ensure exact internal ratios, so cutoffs differ. Temporal gaps prevent test leakage into train/validation.

### 3.3 Test filtering

| Dataset | Rule | Before | After | Retention |
|:---|:---|---:|---:|---:|
| **EC** | **No filter**: EC reliability does not rely on GO-style evidence codes | 4,419 | **4,419** | 100% |
| **GO-MF** | Reviewed **OR** any experimental GO evidence | 11,894 | **6,236** | 52.4% |
| **GO-CC** | Same | 8,392 | **4,692** | 55.9% |
| **GO-BP** | Same | 10,383 | **5,653** | 54.4% |

Ten experimental GO evidence codes: EXP, IDA, IPI, IMP, IGI, IEP, TAS, IC, HTP, HDA. IEA electronic annotations (approximately 51.7%) are noisier and filtered for evaluation. EC lacks a separate evidence-code system and is relatively reliable even for TrEMBL sources.

### 3.4 Label cleaning

**Script:** `scripts/07_clean_datasets.py`. Both input/output are `output/datasets_split_raw/`; overwrite split files with cleaned results and save `{task}_labels.csv`, `*_meta.csv`, and `cleaning_stats.json`.

#### 3.4.1 EC cleaning: merge L4 into L3

1. **Levels:** retain L3 (`x.x.x.-`) and L4 (`x.x.x.y`) EC numbers only.
2. **Global merge:** map all L4 numbers to L3 prefixes (`1.1.1.1` → `1.1.1.-`), avoiding parent–child errors in multilabel evaluation.
3. **Training vocabulary:** remove ECs absent from training; discard chains left without labels.

| Dataset | Before | After | Training label count |
|:---|---:|---:|---:|
| EC train | 15,908 | **15,389** | **260**, all L3 |
| EC val | 1,768 | **1,687** | 175 |
| EC test | 4,419 | **4,243** | 211 |

> Single-label: 94.1%; multilabel: 5.9%.

#### 3.4.2 GO cleaning: frequency filtering

Process namespaces separately. Retain terms with training counts **[5,5000]**: the lower bound removes noisy rare terms, the upper bound overly broad terms. Discard chains with no terms remaining.

| Task | Clean train | Clean val | Clean test | Labels |
|:---|---:|---:|---:|---:|
| **GO-MF** | **18,179** | 2,084 | 4,609 | 447 |
| **GO-CC** | **13,878** | 1,730 | 4,256 | 264 |
| **GO-BP** | **15,232** | 1,742 | 4,121 | 884 |

> GO-BP removes the largest fraction because BP terms are sparse: 67.6% occur only 1–10 times in training.

### 3.5 Outputs

| File | Description |
|:---|:---|
| `output/datasets_split_raw/{task}_{train,val,test}.csv` | Clean labels: chain_key, labels, num_labels |
| `output/datasets_split_raw/{task}_{train,val,test}_meta.csv` | Complete matching sample metadata |
| `output/datasets_split_raw/split_summary.csv` | Split statistics including date cutoffs |
| `output/datasets_split_raw/cleaning_stats.json` | Cleaning report |

> Between Stages 3/4, **schema export has no standalone script**: write cleaned splits to `splits/` using the project-wide `unique_id, aa_seq, struct_file, label` schema. Historical `datasets_split/` means current `splits/`.

---

## Stage 4: OOD annotation

> Add multidimensional flags to **test only** after cleaning. Leave train/val unchanged. Each dimension is boolean; non-length flags on extreme-length rows are NaN.
>
> **Canonical definitions:** use these four skills for recomputation rather than defining alternatives:
> - `.skills/homology-ood-annotation`: seq_Redundancy_*, OOD_Orphan, TM-score_* (mmseqs2/Foldseek; --m8-cache replay).
> - `.skills/cath-ted-domain-annotation`: fold/superfamily holdouts, called by scripts/13_ted_augment_cath_ood.py.
> - `.skills/ec-function-ood-annotation`: OOD_NewEC_* / OOD_LongTail_EC_*, unified all-not-in semantics.
> - `.skills/ood-annotation-toolkit`: IDR, fallback Default=False→NaN validation (generating scripts write NaN directly), InD recomputation, and split checks.

### 4.1 Merge extreme-length chains into test

Integrate chains retained separately in Stage 2.3.1 into each task's test set:

1. **Collect candidates** with `scripts/08_recover_extreme_length.py`: select initial chains whose UniProts are absent from the 47,852 representatives; retain the best per UniProt using has_ec → resolution → mapped_length → lexicographic ID. Extract whole-chain sequences/deposition dates from `data/pdb/`. Length <60 →ExtremeShort; 1000<len≤2000 →ExtremeLong; >2000 discarded. Output `output/extreme_length_chains.csv`, 3,834 rows.
2. **Merge test** with `scripts/09_integrate_extreme_to_splits.py`: clean task labels (EC L4→L3 and training presence; GO frequency [5,5000] and training presence, plus reviewed/experimental evidence as §3.3). Append to `splits/{task}_test.csv`, setting Default=False and the appropriate length flag True; other flags await finalization. Automatically back up to `output/splits_pre_extreme_backup/` before writing.
3. **Extract structures** with `scripts/10_extract_extreme_pdbs.py` / `scripts/11_extract_extreme_cifs.py`: export whole chains to `pdbs/`, text-filtering PDB or filtering CIF by auth_asym_id. Stop multimodel entries at the first ENDMDL.

| Task | Previous test | Added short | Added long | Updated test |
|:---|---:|---:|---:|---:|
| **EC** | 4,243 | 478 | 326 | **5,047** |
| **GO-MF** | 4,609 | 753 | 247 | **5,609** |
| **GO-CC** | 4,256 | 1,090 | 267 | **5,613** |
| **GO-BP** | 4,121 | 757 | 264 | **5,142** |

> **Default/NaN convention:** Default=True denotes normal [60,1000]-residue test samples. Extreme rows have Default=False; non-length OOD_* columns are **not computed**, receiving **NaN** with `.skills/ood-annotation-toolkit/scripts/finalize_ood_columns.py --mode nan`. Empty NaN differs from computed non-OOD False. All four test CSVs follow this convention; backup `output/backup_before_nondefault_nan/`.

### 4.2 OOD-TimeCutoff

The dataset is already split by deposition date (§3.2). All test samples naturally form time-cutoff OOD, serving as the baseline without a separate column.

### 4.3 Sequence-homology OOD (seq_Redundancy_* + OOD_Orphan)

- **Method:** high-sensitivity mmseqs2 `easy-search -s 7`, **query=test, target=train+val**. Take max pident per test sample; seq_Redundancy_<T> means max pident<T, T=90…30. **No significant hit at e-value≤1e-3 sets OOD_Orphan=True.**
- **Canonical tool:** `.skills/homology-ood-annotation/scripts/seq_homology_ood.py`. **Current values recomputed 2026-08-30**; replay m8 cache `analysis/output/seq_ood_unify_assessment/func_{ec,go_bp,go_cc,go_mf}/mmseqs_easysearch.m8` via --m8-cache. Backup `output/backup_before_seq_easysearch_unify/`; report `analysis/output/seq_ood_unify_assessment/REPORT.md`. Historical `output/test_vs_trainval.tsv` (6/9, search+convertalis, fractional fident) is for provenance only: old Orphan matches it rowwise, but old sequence flags differ widely (1,018 EC rows), systematically overestimating OOD. Their source remains uncertain, possibly historical low-sensitivity/truncated search or UniRef50 overlays.
- **Superseded:** archived orphan_mark.py used UniRef50 references and fident≥0.30, incompatible with current definitions and not replayable with the current skill (Notes).

| Task | Default test | 90 | 80 | 70 | 60 | 50 | 40 | 30 | OOD_Orphan |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **EC** | 4,243 | 94.0% | 85.5% | 78.2% | 67.9% | 52.9% | 34.2% | 17.2% | 498 (11.7%) |
| **GO-MF** | 4,609 | 93.3% | 83.7% | 74.5% | 63.8% | 51.4% | 37.6% | 24.4% | 902 (19.6%) |
| **GO-CC** | 4,256 | 92.7% | 82.5% | 74.0% | 65.1% | 55.2% | 43.9% | 31.6% | 1,059 (24.9%) |
| **GO-BP** | 4,121 | 93.0% | 83.3% | 73.9% | 63.0% | 51.3% | 38.8% | 25.7% | 834 (20.2%) |

> Measurements after 2026-08-30 unification; old values in `output/backup_before_seq_easysearch_unify/`. Sequence flag counts decrease (previous OOD overestimation), Orphan increases (previous underestimation), and InD increases (§4.11). Three historical Orphan⊈seq30 violations are resolved.

### 4.4 Structural-homology OOD (TM-score_*)

- **Method:** Foldseek easy-search, test structures against train+val structures; take max alntmscore. TM-score_<T> means max score<T, T=0.9…0.3.
- **Tool:** `.skills/homology-ood-annotation/scripts/struct_homology_ood.py`, also supporting --m8-cache replay.

| Task | Default test | 0.9 | 0.8 | 0.7 | 0.6 | 0.5 | 0.4 | 0.3 |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| **EC** | 4,243 | 32.1% | 13.5% | 7.1% | 3.6% | 2.5% | 1.3% | 0.7% |
| **GO-MF** | 4,609 | 39.9% | 19.7% | 11.5% | 6.8% | 4.4% | 3.1% | 2.0% |
| **GO-CC** | 4,256 | 46.0% | 27.2% | 16.8% | 10.3% | 7.0% | 4.7% | 3.4% |
| **GO-BP** | 4,121 | 41.6% | 21.8% | 12.4% | 7.7% | 5.3% | 3.6% | 2.4% |

> Train+val hit coverage is 97.6%–99.5%. TM<0.5 samples represent 2–7%, forming strict new-fold subsets.

### 4.5 OOD-IDR: high intrinsic disorder

- **Tool:** metapredict v3.0.2 on all representative sequences; strip X/B/Z/J/U first from 41 sequences with nonstandard residues.
- **Rule:** **absolute threshold**: disordered residue fraction IDR_ratio>0.3 →True.
- **Tool:** `.skills/ood-annotation-toolkit/scripts/idr_ood.py --mode residue-ratio --threshold 0.3 --clean strip`, preserving the historical definition; cache output/idr_features.csv.

> **Why not relative thresholds?** 86.3% of PDB chains have IDR_ratio=0; test and train+val distributions nearly coincide, so P90/P95 thresholds produce no OOD samples.

| Task | Default test | OOD_IDR |
|:---|---:|---:|
| **EC** | 4,243 | 27 (0.6%) |
| **GO-MF** | 4,609 | 168 (3.6%) |
| **GO-CC** | 4,256 | 175 (4.1%) |
| **GO-BP** | 4,121 | 157 (3.8%) |

### 4.6 OOD_FoldHoldout / OOD_SuperfamilyHoldout: CATH+TED domains

- **Workflow:** scripts/13_ted_augment_cath_ood.py wraps all four subtasks:
  1. Generate output/func_{sub}_ted_aug_chain_list.csv from splits. This replaces a historical one-off output; regeneration on 2026-08-28 matched the existing file rowwise.
  2. Call `.skills/cath-ted-domain-annotation/scripts/cath_ted_annotate.py` for chain labels.
  3. Call cath_ted_holdout.py from the same skill and write flags to splits.
- **Definition:** chain mode, unique_id=PDB_chain. Union CATH v4.4 experimental labels (name matching and SIFTS domain-to-chain interval mapping) with putative TED domains overlapping the chain UniProt interval by ≥50% of domain length. H supplies C.A.T.H; T only C.A.T. Reference: each task's **training** CATH+TED union, symmetrically augmented across splits. Any unseen topology/superfamily →True; unannotated →False; Default=False →NaN (§4.1).
- **Chain coverage:** CATH-only 43.4%–59.8%; union **86.2%–92.3%**: EC 91.1%, GO-BP 90.9%, GO-CC 86.2%, GO-MF 92.3%.

| Task | Test rows | Default=True | OOD_FoldHoldout | OOD_SuperfamilyHoldout |
|:---|---:|---:|---:|---:|
| **EC** | 5,047 | 4,243 | 31 | 103 |
| **GO-BP** | 5,142 | 4,121 | 50 | 151 |
| **GO-CC** | 5,613 | 4,256 | 57 | 185 |
| **GO-MF** | 5,609 | 4,609 | 31 | 186 |

**Outputs**, one set per {sub} ∈{ec,go_bp,go_cc,go_mf}: `output/func_{sub}_chain_cath_ted_labels.csv` (chain labels), `output/func_{sub}_chain_domains.csv` (domain details), `output/func_{sub}_ted_aug_chain_list.csv` (chain list), `output/backup_before_ted_aug/{sub}_test.csv` (backup).

> TED labels are putative foldseek/foldclass predictions; prioritize experimental CATH and use TED for gaps. Fixed uppercase split PDB/lowercase label-key mismatches and parsing of comma-separated dual TED labels (skill documentation).

### 4.7 OOD_NewEC: new enzyme function

- **All-not-in:** reference ECs from train+val; True only when annotated and **all** ECs at that level are absent. Unannotated →False; Default=False →NaN. Source: **SIFTS chain annotations**, data/sifts/sifts_chain_ec.tsv.gz.
- **Tool:** `.skills/ec-function-ood-annotation/scripts/add_unified_newec.py`.
- **func_ec retains L4 only**, OOD_NewEC_L4. L3 x.x.x.- is the prediction target, so no OOD_NewEC_L3.
- **Result:** **381** True, 9.0% of 4,243 Default samples.
- GO tasks' OOD_NewEC_L3/L4 are listed in §4.9, generated with LongTailEC.

> **Superseded:** old any semantics, JSON source, requiring four non-`-` EC components yielded 441. Replaced by canonical definitions; newfunction_mark.py archived. Backup output/backup_before_newec_unify/.

### 4.8 OOD_LongTail: predicted-label long tails

- **Rule:** task **training** labels occurring **≤10** times form the tail set. Any test label in it →True. Target granularity: EC L3, GO terms.
- **Script:** scripts/12_longtail_mark.py now reads/writes splits/ and flags only Default=True. Historically used output/datasets_split_raw/{task}_test_ood.csv; see Notes. Read-only recomputation on 2026-08-28 exactly matched current splits.

| Task | Training labels | Tail labels (≤10) | OOD test samples | Percentage |
|:---|---:|---:|---:|---:|
| **EC** | 260 | 104 (40.0%) | 163 | 3.8% |
| **GO-MF** | 447 | 127 (28.4%) | 295 | 6.4% |
| **GO-CC** | 264 | 68 (25.8%) | 235 | 5.5% |
| **GO-BP** | 884 | 372 (42.1%) | 684 | 16.6% |

### 4.9 OOD_LongTail_EC: SIFTS L4/L3 function tails

Unlike §4.8 target-label tails, measure training frequencies of **SIFTS chain-level EC L4/L3 annotations**, **not** the label column:

- **All-not-in:** annotated and all ECs at the given level have **training** counts <5 /<10, **including zero**, →True. Unannotated →False; Default=False →NaN. Including zero nests NewEC⊆LongTailEC.
- **func_ec retains two L4 columns**, OOD_LongTail_EC_L4_le5/le10. L3 conflicts with target granularity and duplicates §4.8, so both L3 columns were removed. Tool: `.skills/ec-function-ood-annotation/scripts/add_unified_longtail_ec.py`.

| Column | True count (Default n=4,243) | Percentage |
|:---|---:|---:|
| `OOD_LongTail_EC_L4_le5` | 1,155 | 27.2% |
| `OOD_LongTail_EC_L4_le10` | 1,752 | 41.3% |

> Much higher fractions than target-label tails: SIFTS L4 is finer than the 260-label target space. Backups output/backup_before_longtail_unify/ and output/backup_before_longtail_allnotin/.

- **Six EC OOD columns for GO tasks**, NewEC L3/L4 and LongTail_EC_{L3,L4}_le{5,10}, are generated together by `.skills/ec-function-ood-annotation/scripts/add_go_ec_ood.py` using the same definitions:

| Task | Default test | No EC (→False) | NewEC_L3 | NewEC_L4 | LT_L3_le5 | LT_L3_le10 | LT_L4_le5 | LT_L4_le10 |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| **GO-MF** | 4,609 | 2,993 | 2 | 183 | 34 | 106 | 483 | 736 |
| **GO-BP** | 4,121 | 2,928 | 8 | 118 | 55 | 83 | 370 | 509 |
| **GO-CC** | 4,256 | 3,220 | 9 | 136 | 44 | 97 | 431 | 550 |

> **GO-based functional OOD is infeasible here:** the closed vocabulary makes all test labels a subset of training, so NewGO is zero. All-not-in on full SIFTS GO annotations also yields zero because every chain carries common training IEA terms. LongTailGO is unavailable too: target-level le5 is always zero after vocabulary filtering at ≥5. Thus expand EC OOD with six columns rather than define GO OOD.

### 4.10 OOD_Combinatorial: new label combinations (no retained script)

- **Rule:** enumerate unordered label pairs (A,B), A<B. True if **at least one test pair** is absent from train+val pairs.
- **Rationale:** unseen co-occurrences may fail even when individual labels appeared in training.
- **Reproduction:** historically computed once with no retained standalone script; apply this definition to label in splits/{task}_test.csv.

| Task | Default test | Multilabel samples | OOD_Combinatorial | Test percentage |
|:---|---:|---:|---:|---:|
| **EC** | 4,243 | 313 (7.4%) | 58 | 1.4% |
| **GO-MF** | 4,609 | 2,547 (55.3%) | 449 | 9.7% |
| **GO-CC** | 4,256 | 2,998 (70.4%) | 472 | 11.1% |
| **GO-BP** | 4,121 | 2,008 (48.7%) | 808 | 19.6% |

### 4.11 InD: pure in-distribution flag

- **Rule:** all OOD_*, seq_Redundancy_*, and TM-score_* flags False (NaN treated as False) →InD=True.
- **Tool:** `.skills/ood-annotation-toolkit/scripts/finalize_ood_columns.py --mode ind`, after all OOD columns exist; backup output/backup_before_ind/.

| Task | InD=True | Total test rows | Percentage |
|:---|---:|---:|---:|
| **EC** | 151 | 5,047 | 3.0% |
| **GO-MF** | 194 | 5,609 | 3.5% |
| **GO-CC** | 192 | 5,613 | 3.4% |
| **GO-BP** | 172 | 5,142 | 3.3% |

> Recomputed after easy-search unification on 2026-08-30. Previously 121 /101 /77 /82; backup output/backup_before_seq_easysearch_unify/.

### 4.12 Output format

Final datasets: `splits/`, 12 files, {task} ∈{ec,go_mf,go_cc,go_bp}:

| File | Columns |
|:---|:---|
| `{task}_train.csv` / `{task}_val.csv` | unique_id, aa_seq, struct_file, label: four core columns, no OOD |
| `ec_test.csv` (31 columns) | Four core columns + seq_Redundancy_90…30, TM-score_0.9…0.3, OOD_IDR, OOD_Orphan, OOD_NewEC_L4, OOD_LongTail, OOD_Combinatorial, Default, OOD_ExtremeShort, OOD_ExtremeLong, OOD_FoldHoldout, OOD_SuperfamilyHoldout, OOD_LongTail_EC_L4_le5, OOD_LongTail_EC_L4_le10, InD |
| `go_*_test.csv` (34 columns) | Four core columns + seq_Redundancy_*, TM-score_*, OOD_IDR, OOD_Orphan, OOD_LongTail, OOD_Combinatorial, Default, OOD_ExtremeShort/Long, OOD_FoldHoldout, OOD_SuperfamilyHoldout, OOD_NewEC_L3, OOD_NewEC_L4, four OOD_LongTail_EC_{L3,L4}_le{5,10} columns, InD |

- Store label as a string; separate multilabels with semicolons. EC labels are L3 x.x.x.-.
- struct_file is a filename only ({pdb_id}_{chain}.pdb or .cif), located in pdbs/.

### 4.13 OOD summary

| OOD type | Tool/skill | Metric | Threshold/rule | Difficulty |
|:---|:---|:---|:---|:---|
| TimeCutoff | Inherent in §3.2 | deposition date | Single baseline | All test is OOD |
| seq_Redundancy_* + Orphan | mmseqs2 (homology) | Sequence identity /significant hit | 90…30; no hit | Lower threshold is harder |
| TM-score_* | Foldseek (homology) | Structural TM-score | 0.9…0.3 | Lower threshold is harder |
| OOD_IDR | metapredict (toolkit) | Disordered residue fraction | Absolute >0.3 | High-disorder proteins |
| FoldHoldout /SuperfamilyHoldout | CATH+TED (cath-ted) | Topology /superfamily | Held out against train | Unseen folds/superfamilies |
| OOD_NewEC_L4 | SIFTS EC (ec) | EC L4 presence | All-not-in, train+val | New enzyme functions: EC 9.0% |
| OOD_LongTail | Label frequency (scripts/12_longtail_mark.py) | Training label count ≤10 | Any semantics | Rare predicted labels |
| OOD_LongTail_EC_* | SIFTS EC (ec) | Training EC count <5/10, including zero | All-not-in | All functions are long-tail |
| OOD_Combinatorial | Label co-occurrence, no retained script | Unseen train+val pairs | Any semantics | New label combinations |
| ExtremeShort /ExtremeLong | Length statistics | Residues <60 />1000 (≤2000) | Absolute; Default=False | Extreme length |

---

## Stage 5, optional: collect AlphaFold DB predictions

> **Train/validation expansion only**, excluded from strict test. Complements experimental Stages 1/2 with Swiss-Prot entries lacking PDB structures.
> stage2.py /stage2_uniprot.py /stage2_finish.py have obsolete paths and are decoupled from current splits; archived in output/scripts_archive/ for history. Outputs such as output/alphafold_* and output/representative_proteins.csv remain.

Workflow:

1. **Download:** full Swiss-Prot AlphaFold DB predictions and global pLDDT summaries, reviewed entries only.
2. **Quality:** residue pLDDT from B-factors →output/alphafold_plddt_profiles.json; global pLDDT<70 used for training only.
3. **Function mapping:** Swiss-Prot sequence/GO evidence/EC/metadata →output/alphafold_uniprot_annotations.json.
4. **Consistency:** globally align AlphaFold ATOM sequence with UniProt FASTA; require identity≥0.95 and coverage≥0.90.
5. **Cross-deduplication:** prefer experimental structures; exclude UniProts already represented in Stage 2. AlphaFold-only complement →output/source_overlap_map.json.

| Scenario | AlphaFold used? | Notes |
|:---|:---|:---|
| Strict benchmark | No | Experimental PDB test only |
| Training expansion /large-scale pretraining | Yes | Swiss-Prot entries without PDB |
| Predicted-structure evaluation | Report separately | Build a separate AlphaFold test set |

---

## Data quality checks

Read-only `.skills/ood-annotation-toolkit/scripts/check_splits.py` checks on all four subtasks (2026-08-27):

1. **Leakage:** all pairwise train/val/test unique_id intersections **empty** in all subtasks.
2. **Schema:** train/val have four core columns, no OOD; test includes complete OOD blocks (EC 31, GO 34). No deprecated hyphenated OOD-* or labels columns.
3. **NaN:** non-length OOD_* values on all Default=False extreme rows are **empty**, as required.
4. **Nesting fixed 2026-08-30:** one EC/GO-MF/GO-CC row each (7QS4_A,7WFF_c,8YWA_A) had Orphan=True but seq30=False. All three yield no train+val hits at e≤1e-3, validating Orphan; seq30=False came from historical relaxed-e-value hits (0.01–1, pident31–40). Full easy-search recomputation (§4.3) restores Orphan⊆seq30.
5. **Outliers:** globally exclude >2000-aa chains (§2.3.1/4.1) and low-identity B-tier inference <0.90 (19,906 chains, §2.2).

---

## Key files

| File/directory | Description | Stage |
|:---|:---|:---|
| data/sifts/sifts_*.json / *.tsv.gz | Project-shared chain mappings and GO/EC/CATH annotations | 1 |
| data/sifts/annotated_chains.csv | Initial candidate metadata, 827,186 rows | 1 |
| data/pdb/ | Project-shared source PDB/mmCIF, approximately 255K gzip files | 1 |
| pdbs/ | 50,123 chain structures: 43,567 PDB +6,556 CIF | 2 /4 |
| output/representative_chains.csv | UniProt representative list from 02; regenerated on rerun | 2 |
| output/filtered_sequences.fa | 53,902 length-filtered sequences | 2 |
| output/final_nonredundant_representatives.csv | 47,852 final representatives | 2 |
| output/extreme_length_chains.csv | 3,834 separately retained extreme chains | 2 /4 |
| output/deposition_dates.json | Deposition date per chain | 3 |
| output/datasets_split_raw/ | Split/cleaning intermediates: labels, metadata, statistics | 3 |
| output/test_vs_trainval.tsv, output/trainval_db* | mmseqs2 caches, replayable via --m8-cache | 4 |
| output/idr_features.csv | metapredict IDR cache | 4 |
| output/func_{sub}_chain_cath_ted_labels.csv and related files | CATH+TED outputs per subtask | 4 |
| output/backup_before_*/ | Backups before split writes | 4 |
| output/scripts_archive/ | Legacy scripts, see archive list | — |
| splits/ | **12 final split files** | 4 |
| output/alphafold_*, output/representative_proteins.csv, etc. | Optional AlphaFold outputs | 5 |

---

## Archived scripts: output/scripts_archive/, no longer used

| Script | Reason |
|------|----------|
| fix_multimodel_pdbs.py | One-off NMR repair applied across pdbs/: 1,940 truncated files, backup output/backup_before_multimodel_fix/. Root-cause first-ENDMDL stop incorporated into scripts/03_extract_representative_structures.py and scripts/10_extract_extreme_pdbs.py |
| stage2.py /stage2_uniprot.py /stage2_finish.py | Optional Stage-5 expansion; obsolete workspace paths, decoupled from current splits; retained in Stage-5 documentation |
| newfunction_mark.py | Legacy any-semantics JSON-based NewFunction requiring four non-`-` parts; replaced by canonical all-not-in EC skill, 441→381 |
| orphan_mark.py | UniRef50 reference and fident≥0.30, unlike current train+val reference |
| fetch_uniprot_metadata.py | Old temporal-split helper with invalid hard-coded paths |

## Notes

1. **func_ec: L4 function OOD only.** L3 x.x.x.- is the target, so omit OOD_NewEC_L3 /OOD_LongTail_EC_L3_*. Target-level OOD_LongTail (§4.8) already covers L3 rarity. GO targets do not conflict with EC levels; retain all six EC columns.
2. **Structure formats:** single-character chain IDs →PDB; multicharacter →CIF due to field widths. CIF covers **20.4%**, 784/3,834, of extreme chains. Dispatch parsers by struct_file extension.
3. **Representative reconstruction deviation (2026-08-28):** 02 is a gap-filling reimplementation. Including all B candidates selects 56,057 representatives; **38,835/47,852 (81.2%)** match final chains, 9,017 UniProts differ. Causes: (a) missing historical B identity≥0.90 inference file b_boundary_inference.jsonl, admitting 554,852 B chains versus 419,370 historically; (b) A/B classification drift, A=214,878 versus 330,303 historically, with undocumented partially inferred A rules; (c) regenerated SIFTS pool 769,892 versus 769,730. Historical-definition deviations did not modify any pdbs/ files. Exact reconstruction requires rerunning B inference and supplying --b-keep-list.
4. **Historically inflated B identity:** full-length UniProt alignment inflated identity/coverage. UniProt[sp_beg:sp_end] subsequence alignment is stricter; account for this when rerunning §2.2.
5. **Three nesting violations fixed 2026-08-30:** 7QS4_A (EC),7WFF_c (GO-MF),8YWA_A (GO-CC) had Orphan=True/seq30=False. No hits at e≤1e-3 with -s7/-s7.5 validate Orphan; relaxed e=0.01–1 yields weak pident31–40 hits explaining historical seq30. Orphan matches the 6/9 output/test_vs_trainval.tsv cache rowwise; sequence flags differ widely, with unconfirmed provenance. All four tasks were recomputed using skill easy-search -s7 (§4.3), restoring nesting.
6. **No standalone script:** MMseqs2 clustering (`mmseqs cluster --min-seq-id 0.95 -c 0.8 --cov-mode 0`), within-cluster ranking (§2.1), SIFTS flat-file/go-basic.obo downloads, Stage3→4 schema export, OOD_Combinatorial (§4.10), and whole-library multimodel PDB scanning (output/dup_multimodel_files.json).
7. **Download go-basic.obo separately:** scripts/06_split_datasets.py needs https://purl.obolibrary.org/obo/go/go-basic.obo under output/. Not distributed; output/go_id_to_namespace.json is the existing parsed cache used by scripts/07_clean_datasets.py.
8. **Long-tail I/O deviation:** historical scripts/12_longtail_mark.py ran before extreme integration on output/datasets_split_raw/{task}_train.csv and {task}_test_ood.csv, not retained. Current version uses splits/, Default=True only, with read-only recomputation matching counts 163/295/235/684.
9. **Old paths:** historical datasets_split/ means splits/. Stage-3 intermediates are output/datasets_split_raw/. Shared SIFTS/PDB/UniRef50/CATH/TED caches reside in POOR/data/. Old data/intermediate/ was not migrated; its outputs now use output/ (representative_chains.csv, b_boundary_inference.jsonl, filtered_sequences.fa, etc.).
10. **NMR repairs applied:** 1,940/43,567 PDB files concatenated models without MODEL separators. Truncation at the first nonconsecutive repeated residue fixed them, backup output/backup_before_multimodel_fix/; rescanning found zero residuals. Extraction now stops at the first ENDMDL; one-off fix_multimodel_pdbs.py is archived. All 6,556 CIFs were confirmed single-model.
11. **Default=False sequence/TM NaN unification, 2026-08-28:** 58,548 historical False cells in seq_Redundancy_*/TM-score_* on extreme rows changed to NaN (not computed); backup output/backup_nondefault_seqtm_nan/. InD unchanged, verified. Generating/skill scripts now write NaN directly; finalize_ood_columns.py --mode nan is a fallback no-op on current data.
