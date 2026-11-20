---
name: ec-function-ood-annotation
description: Generate standardized EC annotations for arbitrary tasks and compute POOD-Benchmark OOD_NewEC_L3/L4 (all-not-in, train+val reference) and OOD_LongTail_EC_L3/L4_le5/le10 (all-not-in, training frequencies including zero). Use to add or recompute functional novelty and long-tail EC OOD columns.
---

# EC annotation and NewEC / LongTailEC OOD

## Description

This skill implements unified POOD-Benchmark EC OOD annotation: first generate sample EC sets using SIFTS chain annotations or UniProt reverse mapping plus API completion, then compute six OOD columns. It complements structural classification OOD in `cath-ted-domain-annotation`.

## When to use

- Add/recompute `OOD_NewEC_L3` / `OOD_NewEC_L4`.
- Add/recompute `OOD_LongTail_EC_{L3,L4}_le{5,10}`.
- Generate sample EC sets from standardized sources.

## Canonical definitions: consistent across tasks

| Decision | Convention | Rationale |
|:---|:---|:---|
| Hierarchy `levels()` | ≥3 EC components → L3=first three; exactly four → L4=full number; **partial "-" allowed**, e.g. `3.4.22.-` → L3=`3.4.22` | Matches unified drivers; discarding partial numbers misses annotations |
| **NewEC** reference | ECs in **train+val** | Finalized 2026-08-27 |
| **NewEC** criterion | **all-not-in:** annotated and every EC at the given level absent from reference → True; unannotated → False | Multi-EC samples must have entirely unseen functions; any semantics are too permissive (ppis previously 142→162) |
| **LongTailEC** reference | **Training-only**, per-sample counts: training samples carrying each class | Matches func/optimal_ph/kcat LongTail definitions |
| **LongTailEC** criterion | **all-not-in:** annotated and every EC at the given level has training frequency <5 / <10, **including zero** → True; unannotated → False | Including zero makes NewEC ⊆ LongTailEC: nested scenarios |
| Default/extreme rows | **Default=False preselected extreme-length rows always receive NaN**: project-wide since 2026-08-27, meaning not computed rather than False. Since 2026-08-28, skill/task scripts write NaN at generation; `.skills/ood-annotation-toolkit/scripts/finalize_ood_columns.py --mode nan` is fallback validation only | Both PDA design sets have Default=True throughout (additional test sets without preselection), so this rule does not affect them |

> Both NewEC and LongTailEC use all-not-in (finalized 2026-08-27): NewEC measures entirely unfamiliar functions; LongTailEC requires every function to be long-tail. Early LongTailEC used any semantics, later unified to all-not-in. func_ec `OOD_LongTail` concerns prediction labels at L3, distinct from LongTailEC at L4.

## Sources

| File/API | Purpose |
|:---|:---|
| `data/sifts/sifts_chain_ec.tsv.gz` | PDB+chain → EC; sole chain-mode source |
| `data/sifts/pdb_chain_uniprot.tsv.gz` | Accession → PDB chain reverse mapping; UniProt source ① |
| `data/uniprot_ec_cache.json` | **Project-wide shared UniProt EC cache**; UniProt source ②, reused across tasks |
| UniProt REST `rest.uniprot.org/uniprotkb/search?fields=accession,ec` | Query cache misses in batches of 100, with 0.5s intervals |

> **EC source drift:** current UniProt annotations may differ from original dataset annotations. In kcat, UniProt removed `1.11.1.15` for Q71SQ4/Q9BP39 and changed P91883 from `1.11.1.15` to `1.11.1.24`. Recompute existing tasks using their historical sources below; new tasks may use default sources.

### Historical task sources: preserve on recomputation

| Task | EC source |
|:---|:---|
| kcat | `datasets/enzyme_kinetics_prediction/data/kcat_with_reactant_ecfp4.csv`: original uniprot→ec annotations; split multi-EC strings by `;` |
| optimal_ph | `datasets/enzyme_optimal_ph/data/metadata.csv`, ec_id; **782 UniProts have multi-EC strings separated by `; ` and must be split**, fixed 2026-08-27 |
| ppi | UniProt API cache (`output/uniprot_ec_mapping.csv`) ∪ SIFTS chain ECs; union both for bio/xtal rows, UniProt only for random negatives |
| ligand_affinity | SIFTS chain ECs, unioned across all PDB chains |
| ss / lbs | SIFTS chain ECs |
| ppis | SIFTS chain ECs via seq_hash→pair_meta→(pdb,chain) mapping |
| fold / PDA | `data/sifts/sifts_ec_annotations.json`; NewEC is a presplit EC7 L1 prior, **outside this skill's definition** |
| func_ec | SIFTS chain ECs; **L4 only** for NewEC/LongTailEC. L3 is the prediction level and `OOD_LongTail` already covers its tail. Do not compute LongTailEC from label, which contains L3 values `x.x.x.-` |
| func go_mf/go_bp/go_cc | SIFTS chain ECs: six columns, NewEC L3+L4 and four LongTailEC columns, added 2026-08-27 using `scripts/add_go_ec_ood.py`. GO OOD is not viable: closed vocabulary makes NewGO always zero, and all-not-in on full GO annotations is also always zero |

## Usage

### 1. EC annotation: ec_annotate.py

```bash
# Chain mode: input contains pdb and chain columns.
python .skills/ec-function-ood-annotation/scripts/ec_annotate.py \
    --input my_chain_list.csv --out my_ec_labels.csv

# UniProt mode: input contains uniprot; --offline uses cache only.
python .skills/ec-function-ood-annotation/scripts/ec_annotate.py \
    --input my_acc_list.csv --out my_ec_labels.csv --offline
```

Append `ec_numbers` (semicolon-separated full ECs) and `has_ec` after input columns.

### 2. OOD computation: ec_ood.py

```bash
# Single file: train/val/test identified by split, with ';' for multiple memberships.
python .skills/ec-function-ood-annotation/scripts/ec_ood.py \
    --labels all_ec_labels.csv --split-col split --out test_flagged.csv

# Separate test and reference files; ref also requires split to distinguish train/val.
python .skills/ec-function-ood-annotation/scripts/ec_ood.py \
    --labels test_ec.csv --ref-labels trainval_ec.csv --out test_flagged.csv
```

Output test rows only, appending six columns: `OOD_NewEC_L3/L4`, `OOD_LongTail_EC_L3_le5/le10`, `OOD_LongTail_EC_L4_le5/le10`.

### 3. Write splits: caller responsibility

1. **Back up** test CSV in the task's `output/backup_*` directory.
2. Mask Default=False extreme-length rows to NaN by project convention.
3. Recompute L4 in place; append L3/LongTail columns, matching the 2026-08-27 unification.
4. Check unchanged row count, column order, unique_id, and non-target column values.

### 4. Batch reruns of existing tasks

Three drivers from the 2026-08-27 unification were moved from `analysis/scripts/` into this skill's `scripts/`. Repository root resolves through `Path(__file__).resolve().parents[3]`.

| Script | Action | Tasks |
|:---|:---|:---|
| `add_unified_newec.py` | Write NewEC L3/L4 using all-not-in and train+val | kcat, optimal_ph, ligand_affinity, ppi, ss, lbs, ppis, func_ec (L4 only) |
| `add_unified_longtail_ec.py` | Write four LongTailEC columns using all-not-in and training frequencies including zero | Above plus fold and PDA cath_design non-ES rows; fold/PDA use `sifts_ec_annotations.json` |
| `add_go_ec_ood.py` | Extend six NewEC/LongTailEC columns to GO tasks | func go_mf/go_bp/go_cc |

Drivers encode task-source adapters listed above, write `splits/` directly, and generate `output/*_summary.json`. **Warning:** they do not create backups or perform NaN finalization. Back up before rerunning; afterward run toolkit `finalize_ood_columns.py --mode both` and `check_splits.py`. For detailed definitions or a new single task, prefer `ec_annotate.py` + `ec_ood.py`.

## Validation

- Chain-mode regression baseline: ss Default subset, measured 2026-08-27. NewEC_L3=1, NewEC_L4=246, LongTail_L3=37/77, LongTail_L4=739/1113 (all-not-in). Recomputed flags agree row by row.
- Inspect ≥5 NewEC=True samples and confirm all their ECs at the selected level are absent from train+val.
- Report **row-level** counts explicitly. kcat has multiple substrate rows per accession; row and unique-accession counts differ greatly, e.g. L4<10: 2,355 versus 788.

## Cautions

1. **Duplicate unique_id:** do not write via `flags.loc[uids]`, which expands duplicate labels. Use a dict mapping, e.g. `groupby(level=0).first().to_dict()`.
2. **Pandas float rewriting:** CSV roundtrips with float columns such as `idr_ratio` can reduce 17-digit precision to 16. Diff against backups and restore original text.
3. **Fold:** `OOD_NewEC` is a presplit EC7 L1 prior, outside this skill. func_ec keeps L4 only because L3 is its prediction target.
4. **optimal_ph:** training has only ~8.6k samples with UniProt-level frequency counts; L4 LongTail covers 64–78% and has limited discrimination. L3 is more informative.
5. Failed UniProt API batches cache empty sets; manually remove those keys before retrying.
6. Without a reference `split` column, all rows count as training, inflating LongTail frequencies. Always include split in two-file mode.
