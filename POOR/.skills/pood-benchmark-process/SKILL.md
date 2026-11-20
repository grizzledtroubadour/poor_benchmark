# POOD-Benchmark Data Processing Standards

## Description

Standards for processing raw data in POOD-Benchmark (`raw_data/`):

- Internal task-directory organization.
- Dataset split formats and output paths.
- OOD test-set definitions and naming.
- DATA_PROCESS.md documentation standards.
- Quality control and reproducibility.

## When to use

- Process data in a POOD task directory: fold_classification, ppi_prediction, func_prediction, ss_prediction, ppis_prediction, ligand_binding_site, or antibody_antigen.
- Write, review, or revise DATA_PROCESS.md.
- Generate split files in `{task_dir}/splits/`.
- Design OOD test sets or adjust thresholds.
- Check leakage, distribution consistency, or outliers.

## 1. Task-directory organization

| Directory | Purpose | Contents |
|:---|:---|:---|
| `data/` | Task-specific sources | Raw/downloaded data and task caches; isolate tasks to avoid cross-task mixing |
| `output/` | Intermediate outputs | Processing results, summary tables/statistics, caches, and backups |
| `logs/` | Processing logs | Running/completed script logs for background monitoring, diagnosis, and reproduction |
| `pdbs/` | Final structures | PDB/mmCIF chain structures corresponding to split samples and used as model inputs |
| `splits/` | Final dataset splits | **Only** final `*_train.csv`, `*_val.csv`, `*_test.csv` |
| `scripts/` | Processing scripts | All task-processing scripts; filenames indicate stage and function |

> data/output/logs contain process artifacts; pdbs/splits contain final products ready for training/evaluation.

## 2. Split output paths

All final splits go in `{task_dir}/splits/`, e.g.:

- `fold_classification/splits/`
- `func_prediction/splits/`
- `ppi_prediction/splits/`

> func_prediction/splits/ contains both ec_*.csv and go_{bp,cc,mf}_*.csv.

### 2.1 Filenames

| Split | Filename |
|:---|:---|
| Train | `{task_name}_train.csv` |
| Validation | `{task_name}_val.csv` |
| Test | `{task_name}_test.csv` |

Only the three final split files belong in splits/, never intermediates.

## 3. Split format

### 3.1 Core columns

| Field | Type | Description |
|------|------|------|
| `unique_id` | string | Unique sample identifier. Single-chain tasks typically use `{pdb_id}_{chain}` or CATH domain ID; PPI uses a protein-pair ID |
| `aa_seq` | string | Protein/chain amino-acid sequence. Multichain proteins may separate chains with `\|`; remove separators for downstream tools. Pair tasks use aa_seq1/aa_seq2 |
| `struct_file` | string | **Filename only**, no directory, e.g. 7E3T_A.pdb in task pdbs/. Extension .pdb/.cif determines format. Pair tasks use struct_file1/struct_file2 |
| `label` | string | **Must be stored as a string**. Semicolon-separated multilabel targets; array strings for residue labels |

### 3.2 Task-specific extensions

#### Single-chain single-label tasks: e.g. fold_classification

```
unique_id, aa_seq, struct_file, label
```

#### Single-chain multilabel tasks: EC/GO

```
unique_id, aa_seq, struct_file, label
```

- Separate multiple labels with ASCII semicolon `;`.
- EC example: `2.7.6.-;3.1.7.-`.
- GO example: `GO:0006164;GO:0006189`.

#### Pair-level tasks: PPI

| Field | Type | Description |
|------|------|------|
| `unique_id` | string | Unique protein-pair ID |
| `protein_A_id` / `protein_B_id` | string | Protein A/B identifiers, e.g. UniProt |
| `aa_seq1` / `aa_seq2` | string | Protein A/B sequences |
| `struct_file1` / `struct_file2` | string | Structure filenames only, in pdbs/ |
| `label` | int/string | Interaction: 1=yes, 0=no |

#### Residue-level tasks: secondary structure prediction

```
unique_id, aa_seq, struct_file, label
```

- Label arrays must **align residue by residue** to aa_seq and be stored as strings.
- **Serialize bracketed, single-space-separated integers**, `[v1 v2 v3 ... vn]`, e.g. `[0 0 1 0 -1]`; **no commas**. Mask unresolved/missing residues with -1.
- Require `len(label) == len(aa_seq)`.

#### Ligand-related tasks: kcat/LBA/LBS

Add ligand_smiles and ligand_ecfp4 to core columns.

### 3.3 OOD flags: test only

Append boolean True/False OOD columns to test CSVs after core columns. Non-length OOD columns for extreme-length samples may be NaN (§3.4).

| Column | Description |
|------|------|
| `Default` | Regular test samples not preselected by extreme-length filtering |
| `seq_Redundancy_90`–`seq_Redundancy_30` | Sequence homology OOD, 10% steps |
| `TM-score_0.9`–`TM-score_0.3` | Structural similarity OOD, 0.1 steps |
| `OOD_Orphan` | No significant train+val homologous hit |
| `OOD_ExtremeShort` / `OOD_ExtremeLong` | Extreme chain length <60 / >1000 |
| `OOD_IDR` | High intrinsic-disorder fraction |
| `OOD_FoldHoldout` | Held-out fold/CATH topology |
| `OOD_SuperfamilyHoldout` | Held-out CATH superfamily |
| `OOD_NewEC_L4` / `OOD_NewEC_L3` / `OOD_NewEC_L1` | Unseen EC functions; suffix indicates full level 4, level 3, or top-level class |
| `OOD_LongTail` | Infrequent training **prediction-label** classes for the task itself: CATH architecture, GO/EC, secondary structure, interface/binding-site ratios, etc. |
| `OOD_LongTail_EC_L4` / `OOD_LongTail_EC_L3` | Infrequent training **EC function** classes, including as OOD dimensions for non-EC tasks; suffix indicates EC level |
| `OOD_Combinatorial` | Label combinations absent from train+val |
| `InD` | All OOD_*, seq_Redundancy_*, TM-score_* False, treating NaN as False. Formerly Pure_ID; recompute with toolkit finalize_ood_columns.py --mode ind after all OOD columns are ready |

> **OOD annotation skills:** canonical definitions and tools reside in homology-ood-annotation (seq/Orphan/TM), cath-ted-domain-annotation (fold/superfamily), ec-function-ood-annotation (NewEC/LongTailEC), and ood-annotation-toolkit (IDR, value-bin LongTail, NaN/InD finalization, split audit). Call the corresponding scripts when adding/recomputing columns rather than inventing separate definitions.

> **Threshold variants:** use suffixes for multiple thresholds, e.g. optimal_ph OOD_LongTail_EC_L3_le5/le10 (EC L3 frequencies ≤5/≤10) and OOD_LongTail_pHbin_le50/le100 (pH-bin counts ≤50/≤100, prediction-label tails).

### 3.4 Default semantics

Default identifies **regular test samples**, not preselected for extreme length, and serves as the baseline evaluation subset.

| Task type | Default=True | Default=False | Notes |
|:---|:---|:---|:---|
| **LBA** | All test samples | None | Train/val retain extreme lengths; test does not exclude them |
| **PPI** | Both chains ≤1000 | Either chain >1000 | Preselect pairs containing >1000-aa chains into test as OOD_ExtremeLong; train/val chains ≤1000 |
| **All other tasks** | Length in [60,1000] | Length <60 or >1000 | Extreme samples moved from train/val to test before/after splitting, or reserved as separate OOD subsets |

#### Extreme-length requirements: tasks with Default=False rows; not LBA

1. **Other OOD flags become NaN.** These samples do not participate in sequence/structure/IDR or other OOD computations. Preserve only the triggering OOD_ExtremeShort/OOD_ExtremeLong reason flag; use NaN rather than False elsewhere.
2. **Exclude from ID-Test statistics.** Default=False samples are not counted as in-distribution test data.
3. **Evaluate separately** as OOD_ExtremeLength/Long/Short subsets.

> PPI preselects pairs with a >1000-aa chain into test: Default=False, ExtremeLong=True, 1,325 rows. Other flags, including seq/TM/ExtremeShort, are NaN since 2026-08-27; historically seq/TM were computed for all rows. LBA uses Default=True throughout.

#### Hard exclusion: single chains >2000 aa

In addition to flagging/moving extreme samples, **remove every sample containing a chain longer than 2000 aa from all splits**, rather than retaining it as an OOD subset. Very long sequences/structures impose memory/processing costs and lie in the extreme distribution tail, unsuitable for routine evaluation.

- Unit: individual chain length; for paired tasks **either chain** >2000 removes the pair.
- Applies to all tasks.

### 3.5 Constraints

1. **Store label as string**, avoiding floating-point interpretation of 1.10 or 3.2.1.-.
2. **Complete test annotations:** every OOD flag has an explicit value, not an unspecified missing value.
3. **Identifiable ID-Test:** test samples with no OOD flags form the in-distribution subset.
4. **No unique_id overlap** between train/val/test.

## 4. Standard OOD definitions

| Type | Definition | Tool | Threshold/rule | Scientific purpose |
|:---|:---|:---|:---|:---|
| **LowHomology** | Maximum test-to-train+val sequence identity below threshold | mmseqs2 search | **<90**, <80, <70, <60, <50, <40, <30 | Generalization with decreasing sequence similarity |
| **Orphan** | No significant train+val homolog | mmseqs2 search | Max identity=0 or e-value >1e-3 | Generalization without homologous training examples |
| **FamilyHoldout** | Held-out protein families | CATH/SCOPe/Pfam/ECOD | Group holdout | Unseen-family generalization |
| **SuperfamilyHoldout** | Held-out superfamilies | Same as above | Group holdout | Remote homology recognition |
| **NewFold** | Maximum structural similarity below threshold | Foldseek | **<0.9**, <0.8, <0.7, <0.6, <0.5, <0.4, <0.3 | Fold-space generalization |
| **FoldHoldout** | Held-out fold types | CATH/SCOPe/ECOD | Group holdout | Unseen structural folds |
| **IDR** | Substantial intrinsic disorder | metapredict v3/flDPnn/IUPred | Regions with disorder >0.5 cover >30% of sequence | Disorder-rich proteins |
| **ExtremeLength** | Protein lengths in distribution tails | Length statistics | <60 or >1000 | Extreme-length robustness |
| **NewFunction** | Held-out function classes | EC/GO/InterPro | Group holdout | Unseen functions |
| **LongTailFunction** | Infrequent training **functional** labels | EC/GO/InterPro frequencies | EC count **<10**, adjustable by task | Long-tail functions |
| **LongTail** | Infrequent training **prediction labels** | Label frequencies | Class count **<10** | Prediction-label tail generalization |
| **CombinatorialFunction** | Unseen test label combinations | Multilabel combination statistics | Pair (A,B) never co-occurs in train+val | Compositional generalization |

## 5. DATA_PROCESS.md standards

### 5.1 Top-level structure

```markdown
# Task name + Dataset Processing
> Purpose and processing location

---

## Stage N: Name and objective
> Objective statement

### N.M Substep
#### N.M.K Finer-grained operation, if needed

---

### X. Output files
### X+1. Dimension statistics
### X+2. Notes
```

### 5.2 Required stage content

1. **Overview:** core task and scientific rationale.
2. **Scale:** required at first mention; bold key numbers.
3. **Sources:** input-file table with paths, purposes, and notes.
4. **Steps:** purpose, method, script path, and results per step.
5. **Validation:** methods, metric distributions, and anomaly explanations.

### 5.3 Required tables

Use tables for dataset scale, file inventories, filtering conditions, comparisons, and OOD split results.

### 5.4 Script/tool records

Record each script's path (scripts/script_name.py), one-sentence purpose, key parameters (concurrency/threshold/mode), success/failure counts, and runtime.

For external tools record name/version, core parameters, search direction, and decision metric.

### 5.5 Quality control

Include a separate **Data quality checks** section:

1. **Leakage:** whether pairwise train/val/test intersections are empty.
2. **Distribution consistency:** compare length, resolution, class distributions, and other key features.
3. **Outliers:** list removed samples and reasons.

### 5.6 Output records

At each stage's end list output location, filename, row count/size, format, and description.

## 6. Text and number style

| Item | Requirement |
|:---|:---|
| Language | Match the project's main language |
| Numbers | Thousands separators; one decimal for percentages |
| Bold | Key conclusions, central numbers, recommended thresholds, and important warnings |
| Code/paths | Backticks around paths, scripts, columns, and field values |
| Scientific names | Preserve mmseqs2, Foldseek, PDB, UniProt, etc. |

## 7. Checklist

Before finalizing DATA_PROCESS.md:

- [ ] Header states processing location.
- [ ] Every stage opens with an objective.
- [ ] Tables record all scale changes, input→output.
- [ ] Script paths and parameters are documented.
- [ ] External tool versions and key parameters are recorded.
- [ ] OOD splits state definitions, tools, thresholds, scale, and scientific purpose.
- [ ] Leakage checks and results are recorded.
- [ ] Outliers/filtered samples are listed.
- [ ] Complete output inventory includes paths, formats, and sizes.
- [ ] Key fields/columns have types and meanings.
- [ ] Split files conform to `{task_dir}/splits/`.
- [ ] splits/ contains only final train/val/test files.
- [ ] data/output/logs/pdbs/splits/scripts follow directory standards.
- [ ] Train/val have no OOD columns; test has all flags.
- [ ] Default matches task type: LBA all True; PPI preselected >1000-aa chains False; others by extreme length.
- [ ] Other OOD flags on Default=False extreme rows are NaN or False.
- [ ] label is stored as string.
- [ ] Notes cover known issues and special handling.
- [ ] Non-obvious strategies have design explanations.
