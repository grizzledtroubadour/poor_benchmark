# Protein Design Archive (PDA): Additional Test-Set Construction (v2)

> Extract designed protein chains from PDA to build **additional test sets** for POOR.
> Current tasks: **fold prediction** (`splits/cath_design.csv`) and **EC in func_prediction** (`splits/ec_design.csv`).
>
> Source: `POOR/data/protein-design-archive/`, a PDA website repository snapshot (Nature Biotechnology 2025). Core file: `backend/scripts/data.json`, with 1,879 entries /2,389 chains.
>
> **v2 design, rebuilt 2026-08, replacing v1:**
> 1. Pool: all **chain_type D+U** chains (2,207), no longer restricted to synthetic-construct sources. **N (natural-source) chains are excluded entirely.**
> 2. **Exact sequence deduplication only**, without 95% identity clustering: small sequence changes may encode different functions.
> 3. Exclude train/val overlap by **PDB ID match OR maximum mmseqs pident ≥95**.
> 4. EC labels come **only from the chain's own annotations** (PDB depositor ∪ SIFTS), without within-cluster transfer.
>
> v1 (synthetic-only +95% deduplication) outputs are archived in `output/backup_v1/`.

**Current measured splits:**

| File | Rows | Description |
|------|------|------|
| `splits/cath_design.csv` | **122** | Additional fold test set, 29 topology labels; schema matches `fold_classification/splits/cath_test.csv` |
| `splits/ec_design.csv` | **31** | Additional EC test set, 15 L3 labels, one multilabel sample; schema matches `func_prediction/splits/ec_test.csv` exactly |

---

## Stage 1: Obtain and parse sources

### 1.1 Extract designed chains — `scripts/01_extract_chains.py`

- Input: `POOR/data/protein-design-archive/backend/scripts/data.json`.
- Extract all `chain_type ∈ {D, U}` chains: **2,207** (D 1,223 /U 984) across 1,840 PDB entries. Skip 182 N chains.
- For comma-separated `chain_id` copies of the same entity (1,017 records), select the first chain.
- `unique_id = {pdb}_{chain}`; use `chain_seq_nat` as sequence.
- Retain entry-level PDA CATH fields `cath_codes` / `cath_cat_set` / `label_candidate`. Only entries with one distinct C.A.T are unambiguous: 250 records.

**Output:** `output/design_chains_all.csv`.

### 1.2 Extract single-chain structures — `scripts/02_extract_structures.py`

- Prefer `POOR/data/pdb/{PDB}.pdb.gz`, otherwise `{PDB}.cif.gz`, with uppercase-filename fallback.
- PDB: filter ATOM/HETATM by chain column and write `.pdb`. CIF: filter text rows by `_atom_site.auth_asym_id` and write `.cif`, avoiding BioPython stalls on large CIFs.
- For NMR/multimodel entries, stop at the first `ENDMDL`, retaining model 1 only.
- Success: **2,207 /2,207**.

**Output:** temporary pre-deduplication structures in `output/chains_all/`.

---

## Stage 2: Aggregate and compute labels

### 2.1 Fold (CATH) labels

Use the **entry-level** PDA CATH annotations retained in §1.1.

> **Important granularity limit:** PDA `cath_full` contains entry-level code sets from Foldseek→CATH searches of designed structures (upstream `merizo_to_json.py`), **without chain/domain mappings**. The fold test set therefore keeps only chains whose entries have **one distinct C.A.T** (`label_candidate` nonempty). Multilabel C.A.T entries are not disambiguated.

### 2.2 EC labels: strict chain-level annotations without transfer

Implemented in the first half of `scripts/06_build_ec_testset.py`. Each chain uses annotations from **its own entry**:

- PDB: map `MOLECULE:/CHAIN:/EC:` subfields in COMPND by chain.
- CIF: map entity-level `_entity.pdbx_ec` to chains through `_atom_site.(auth_asym_id, label_entity_id)`.
- Union with SIFTS `data/sifts/sifts_ec_annotations.json`, using uppercase `{PDB}_{chain}` keys.

When **identical sequences** recur across entries, annotated members must have identical L3 label sets. Otherwise remove the whole group. Current conflicts: zero; write conflicts to `output/ec_dropped_conflict.csv`.

EC annotations match **69 D+U chains** (`output/design_ec_annotations.csv`), leaving **43** after exact deduplication, with zero conflicts.

---

## Stage 3: Redundancy reduction and refinement

### 3.1 Exact deduplication — `scripts/03_dedup_exact.py`

- Merge **only identical `aa_seq` values**: 2,207→**1,726**.
- Prefer a representative with successful structure extraction; otherwise use the first member.
- **Copy**, rather than symlink, representative structures into `pdbs/`: 1,726 files for reuse by downstream tasks.

**Output:** `output/design_chains_dedup_exact.csv` and `pdbs/`.

---

## Stage 4: Select test sets and annotate OOD

> These are **additional test sets**, with no train/val split. Select labels, remove overlap against each task's existing train/val, annotate OOD, and finalize using canonical skills.

### 4.1 Select fold test set — `scripts/04_build_fold_testset.py`

| Filter | Rule | Count |
|------|------|------|
| Unambiguous label | One distinct entry-level C.A.T (`label_candidate`) | 213 |
| Supported label | C.A.T belongs to fold training's 1,163 topology labels | 211 |
| PDB overlap exclusion | pdb_id in fold train/val | −32 |
| Sequence overlap exclusion | mmseqs2 `easy-search -s 7.5` against train+val; max pident ≥95 | −88; union −89 |
| **Final** | | **122**, 29 topology labels |

**Output:** `output/cath_design_base.csv`, `output/removed_overlap.csv`, and `output/design_mmseqs_vs_trainval.m8`, reused in §4.2.

### 4.2 Annotate fold OOD — `scripts/05_annotate_ood.py`

Schema matches `fold_classification/splits/cath_test.csv`. Definitions follow canonical skills:

- `seq_Redundancy_*` / `OOD_Orphan`: mmseqs2 easy-search versus fold train+val; `TM-score_*`: Foldseek easy-search. Aligned with `.skills/homology-ood-annotation/` (`seq_homology_ood.py` / `struct_homology_ood.py`).
- `OOD_IDR`: metapredict v3 residue-ratio at 0.3, aligned with the fold definition in `.skills/ood-annotation-toolkit/idr_ood.py`.
- `OOD_SuperfamilyHoldout`: all candidate C.A.T.H codes from PDA `cath_codes` whose C.A.T matches the label are absent from train/val.
- `OOD_NewEC`: SIFTS EC7 prior, the fold convention; largely inapplicable because designed proteins usually lack EC7 annotations.
- `OOD_ExtremeShort`: <60 aa; `OOD_LongTail`: training label frequency ≤10.

> **Default is always True:** no samples were reserved before splitting during construction, matching the fold convention that extreme-short (ES) rows are also Default=True. Changed from all False on 2026-08-27; backup `output/backup_before_default_true/`. Other dependent OOD columns remain uncomputed/empty for 45 ES rows. **Known lag:** the main fold pipeline now annotates dependent OOD on ES rows, but this additional set has not been updated.

| OOD axis | True count | Notes |
|--------|--------:|------|
| OOD_ExtremeShort | 45 | 37%; designs are enriched in miniproteins/short peptides |
| OOD_Orphan | 37 | 48% of 77 dependency-eligible rows |
| seq_Redundancy_30 | 39 | 51% |
| TM-score_0.9 /0.5 | 44 /1 | Few novel folds |
| OOD_SuperfamilyHoldout /OOD_NewEC /OOD_LongTail | 0 | — |
| OOD_IDR | 2 | — |
| **InD** | **8** | 6.6% |

**Checks:** zero PDB overlap, maximum pident 94.5 <95, no identical sequences, labels ⊆ training labels, all structure files present, bidirectional InD consistency, all Default=True, unique IDs.

**Output:** `splits/cath_design.csv`, `output/design_foldseek_vs_trainval.m8`.

### 4.3 Select EC test set — second half of `scripts/06_build_ec_testset.py`

- Normalize L4→L3 and join multiple labels with `;`. Require every L3 label to belong to the 260 EC training labels: 43→**43**.
- Exclude overlap with EC train/val by PDB ID OR mmseqs2 `easy-search -s 7` max pident ≥95: PDB exclusions 2, identity exclusions 12, union 12. **Final: 31**, 15 L3 labels.

**Output:** `output/ec_design_base.csv`, `output/ec_removed_overlap.csv`, `output/design_ec_mmseqs_vs_trainval.m8`.

### 4.4 Annotate EC OOD — `scripts/07_annotate_ec_ood.py`

Schema matches `func_prediction/splits/ec_test.csv` **exactly**, without idr_ratio:

- `seq_Redundancy_*` / `OOD_Orphan`: mmseqs2 `-s 7` against EC train+val, reusing §4.3's .m8. `TM-score_*`: Foldseek against EC train+val structures. Aligned with `.skills/homology-ood-annotation/`.
- `OOD_IDR`: metapredict v3 residue-ratio at 0.3, aligned with the function-task definition in `.skills/ood-annotation-toolkit/idr_ood.py`.
- `OOD_NewEC_L4`: any level-4 EC absent from train+val L4; False if no L4 annotation. `OOD_LongTail`: training frequency ≤10. `OOD_Combinatorial`: label pair never co-occurs in train+val.
- All Default=True (§4.2). No extreme-length samples in this batch; all OOD columns computed for every row.

| OOD axis | True count | Fraction of 31 |
|--------|--------:|------:|
| seq_Redundancy_90 /30 | 26 /2 | 84% /6% |
| OOD_Orphan | 2 | 6% |
| TM-score_0.9 /0.5 | 7 /4 | 23% /13% |
| OOD_IDR | 1 | 3% |
| NewEC_L4 /LongTail /Combinatorial /ES /EL | 0 | — |
| **InD** | **4** | 12.9% |

**All checks pass:** zero PDB overlap or identical sequences with EC train/val, max pident 94.6 <95, labels ⊆ training labels (15 L3), all structure files present, exact schema match to ec_test.csv, no missing Boolean values, unique IDs, all sequence-redundancy flags True for orphan rows, D+U-only input pool.

**Output:** `splits/ec_design.csv`, `output/design_ec_foldseek_vs_trainval.m8`.

### 4.5 Canonical skill finalization: scripts under `.skills/`

1. **Four LongTailEC columns:** batch driver `ec-function-ood-annotation/add_unified_longtail_ec.py`, applied to `cath_design.csv`. Annotated design chains are True if **all** L3/L4 ECs have frequency <5 /<10 in fold `cath_train.csv`, including zero (all-not-in). ES rows receive NaN; non-ES rows, including NewEC-prior rows, are evaluated, as for OOD_LongTail. SIFTS EC coverage is sparse: only three of 77 non-ES rows are annotated. **Only `OOD_LongTail_EC_L4_le10` has one True**; the other three columns have zero. Summary: `output/longtail_ec_unified_summary.json`. Training frequencies match fold exact keys: 3,053 L4 classes.
2. **InD recomputation:** `ood-annotation-toolkit/finalize_ood_columns.py --mode ind` updates both split files. Pure InD requires all `OOD_*` / `seq_Redundancy_*` / `TM-score_*` flags False. `--mode nan` has no effect because it targets Default=False and all these rows are True.
3. Including `OOD_LongTail_EC_*` reduces cath_design InD from legacy Pure_ID=9 to **8**.

---

## Directory layout

```
protein-design-archive/
├── DATA_PROCESS.md                      # This document, v2
├── scripts/
│   ├── 01_extract_chains.py             # Complete D+U chain list
│   ├── 02_extract_structures.py         # Single-chain structure extraction
│   ├── 03_dedup_exact.py                # Exact sequence deduplication only
│   ├── 04_build_fold_testset.py         # Fold filtering
│   ├── 05_annotate_ood.py               # Fold OOD annotation
│   ├── 06_build_ec_testset.py           # EC annotation extraction and filtering
│   └── 07_annotate_ec_ood.py            # EC OOD annotation
├── output/
│   ├── design_chains_all.csv            # 2,207 D+U chains
│   ├── chains_all/                      # 2,207 pre-deduplication single-chain structures
│   ├── design_chains_dedup_exact.csv    # 1,726 exact-deduplicated chains reused downstream
│   ├── cath_design_base.csv / removed_overlap.csv
│   ├── design_mmseqs_vs_trainval.m8 / design_foldseek_vs_trainval.m8
│   ├── design_ec_annotations.csv        # 69 chain-level EC matches
│   ├── ec_design_base.csv / ec_removed_overlap.csv
│   ├── design_ec_mmseqs_vs_trainval.m8 / design_ec_foldseek_vs_trainval.m8
│   ├── cath_design_longtail_ec_summary.json   # Initial single-column LongTailEC summary
│   ├── longtail_ec_unified_summary.json       # Unified four-column LongTailEC summary
│   ├── dup_multimodel_files.json        # Files affected by multimodel repair
│   ├── s04_tmp/ ... s07_tmp/            # Step 4–7 caches, including Foldseek symlinks
│   ├── scripts_archive/                 # Archived scripts listed below
│   ├── backup_before_longtail_aug/      # Split backups before LongTailEC iterations
│   ├── backup_before_longtail_unify/    # and InD recomputation
│   ├── backup_before_longtail_allnotin/
│   ├── backup_before_ind/
│   ├── backup_before_multimodel_fix/    # Original PDBs before multimodel repair
│   └── backup_v1/                       # Archived synthetic-only +95% deduplication v1
├── pdbs/                                # 1,726 exact-deduplicated design-chain structures
└── splits/
    ├── cath_design.csv                  # Additional fold test set: 122
    └── ec_design.csv                    # Additional EC test set: 31
```

## Archived scripts: `output/scripts_archive/`

| Script | Reason archived |
|------|---------|
| `fix_multimodel_pdbs.py` | One-off repair of 203 NMR PDBs that concatenated models without MODEL separators: truncate at the first non-contiguous repeated residue. Backups: `output/backup_before_multimodel_fix/`. Data are repaired; `02_extract_structures.py` now stops at the first ENDMDL. Validation found zero residual duplicates or malformed lines; all 203 files parse in Bio.PDB as single models, as do all 52 CIFs. |
| `08_add_longtail_ec_ood.py` | Initial any-semantics single-column `OOD_LongTail_EC_L4` replaced by the EC skill's canonical four-column all-not-in implementation (§4.5). |

## Notes

1. **Entry-level CATH:** PDA lacks chain mappings; entries with multiple C.A.T labels (the majority) are excluded from the fold test set.
2. **Computed/depositor labels:** fold CATH labels derive from PDA Foldseek→CATH searches; EC labels from PDB depositors ∪ SIFTS. They need not have independent experimental validation.
3. **U chains:** unknown-source chains are treated as designs because PDA includes them; a small number of ambiguously annotated natural chains may remain.
4. **Overlap removal:** exclude PDB ID overlap OR maximum train/val pident ≥95%. Do not deduplicate against fold/EC **test** sets, since both are evaluation sets.
5. **Exact-only internal deduplication:** ≥95%-homologous design variants, such as successive designs in a series, may coexist. These highly correlated samples should not be treated as independent during evaluation.
6. **OOD_NewEC is largely inapplicable** to fold design proteins, which often lack UniProt/EC7 annotations.
7. **Known cath_design ES annotation lag:** dependent OOD columns on 45 ES rows remain empty, unlike current fold annotations. `finalize_ood_columns.py --mode nan` does not change these Default=True rows.
8. **Default always True:** neither additional set reserves samples before splitting. Changed from all False on 2026-08-27; backup `output/backup_before_default_true/`.
9. **History:** the initial single-column any-semantics LongTailEC script is archived. Since 2026-08-27, use project-wide four-column all-not-in (§4.5) and recompute InD through `finalize_ood_columns.py --mode ind`.
