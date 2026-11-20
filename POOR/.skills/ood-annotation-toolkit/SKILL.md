---
name: ood-annotation-toolkit
description: Generic POOD-Benchmark computational OOD and finalization tools: metapredict IDR OOD (idr_ood.py, two historical modes), regression value-bin LongTail OOD (value_bin_longtail_ood.py), Default=False NaN masking and InD recomputation (finalize_ood_columns.py), and read-only split consistency/leakage checks (check_splits.py). Use for IDR/value-bin annotation, NaN/InD finalization, or split audits.
---

# OOD annotation toolkit: computational OOD and finalization

## Description

Four cross-task tools implement IDR OOD, value-bin LongTail OOD, NaN/InD finalization, and split consistency checks. Complements `homology-ood-annotation` (search homology), `cath-ted-domain-annotation` (structural classification), and `ec-function-ood-annotation` (EC functions).

## When to use

- Compute `idr_ratio` and `OOD_IDR` with metapredict v3.
- Add regression `OOD_LongTail_*Bin_le*` columns for kcat/affinity/optimal_ph-like tasks.
- Apply Default=False→NaN conventions or recompute test `InD`.
- Audit train/val/test leakage, schemas, and Default/InD statistics without writing.

## Definitions

### idr_ood.py: task-specific historical definitions require explicit mode selection

| Decision | Convention |
|:---|:---|
| Predictor | metapredict v3 (tested 3.0.2); residue disorder >0.5 means disordered |
| `--mode residue-ratio` | Fraction of disordered residues: ss/func/fold/kcat/affinity/optimal_ph. Affinity's sum of region lengths divided by sequence length is mathematically equivalent |
| `--mode region-ratio` | **Longest continuous disordered region** divided by sequence length: ppi, from archived `18_compute_ood_idr.py` |
| `--threshold` | Default 0.3 (ss/func/fold/ppi); **affinity/optimal_ph historically use 0.1** |
| `--clean strip` (default) | Remove nonstandard amino-acid characters: fold/ss |
| `--clean map` | B→N, Z→Q, J→L, U→C, O→K, others→A: ppi |
| Paired chains | `--seq-col aa_seq1,aa_seq2`; either chain exceeding threshold sets OOD_IDR (any semantics, ppi) |
| Cache | Incremental by cleaned sequence; CSV stores sequence, residue_ratio, max_region_ratio. Save each batch; changing mode/threshold requires no predictions |

Historical task parameters to preserve on recomputation:

| Task | Mode | Threshold | Cleaning | Source |
|:---|:---|:---|:---|:---|
| ss_prediction | residue-ratio | 0.3 | strip (historically removes XBZJU only, approximately strip) | `scripts/09_run_idr_prediction.py` |
| func_prediction | residue-ratio | 0.3 | strip | DATA_PROCESS.md IDR section |
| fold_classification | residue-ratio | 0.3 | strip | `scripts/05_annotate_test_ood.py`; cache `output/idr_ratio_full.csv` |
| ppi_prediction | **region-ratio** | 0.3 | **map** | `output/scripts_archive/18_compute_ood_idr.py` / `scripts/10_recompute_all_ood.py` |
| ligand_binding_affinity | residue-ratio | **0.1** | strip, also removing `\|` chain separators | `scripts/10_compute_structure_idr_ood.py` |
| enzyme_optimal_ph | residue-ratio | **0.1** | strip | `scripts/04_compute_ood_markers.py`; former compute_ood_idr.py archived and incorporated |
| enzyme_kinetics | residue-ratio, continuous `idr_ratio` only | — | strip | DATA_PROCESS.md |

### value_bin_longtail_ood.py: regression value-bin long tails

| Decision | Convention |
|:---|:---|
| Criterion | Default=True test row's bin has **training** frequency <=threshold → True; missing/out-of-range labels → False |
| Separate clipped bins | `--clip lo,hi`: values equal to lo/hi form separate bins, preventing censored kcat boundary mass from hiding tails |
| Default=False | NaN by default (`--nondefault nan`, generation-time project convention); false retained only for compatibility |
| Binning | `pd.cut` uses left-open/right-closed `(a,b]`; edges = arange(start, end, width) |

Preserve historical task parameters:

| Task | Label | Bins | Clip | Thresholds | Prefix |
|:---|:---|:---|:---|:---|:---|
| kcat | log10(kcat) | [-6.5, 7.0), width 0.5 | -6,6 | 50,100 | `OOD_LongTail_KcatBin` |
| ligand_binding_affinity | -log10(K) | [-0.5, 16.5), width 0.5 | None | 50,100 | `OOD_LongTail_AffBin` |
| enzyme_optimal_ph | pH | [1.5, 13.5), width 0.5 | None | 50,100 | `OOD_LongTail_pHbin` |

### finalize_ood_columns.py: NaN and InD finalization

| Decision | Convention |
|:---|:---|
| `--mode nan` | Set **OOD_-prefixed** columns except `OOD_ExtremeShort`/`OOD_ExtremeLong` to NaN on Default=False rows, meaning not computed. **Only OOD_ columns** are included; seq_Redundancy_*/TM-score_*/idr_ratio are outside scope, matching original `nondefault_ood_nan.py` |
| Exemption | `--exempt-nondefault` skips nan for intentionally all-Default=False special files. No current task needs it: both PDA design sets have Default=True throughout |
| `--mode ind` | InD means all `OOD_*`, `seq_Redundancy_*`, `TM-score_*` flags False, treating NaN as False. Overwrite existing InD; rename/recompute legacy `Pure_ID` in place (fold/PDA); otherwise append |
| Byte protection | Read all columns with `dtype=str, keep_default_na=False`; assert non-target cells unchanged before writing |

### check_splits.py: read-only audit

Check pairwise train/val/test unique_id intersections; required columns; deprecated columns; missing labels; absence of OOD columns on train/val; column-set agreement (train==val, test adds only OOD and valid annotations such as idr_ratio); Default/InD counts; Default=False NaN violations; and `OOD_Orphan ⊆ seq_Redundancy_30`. Any issue exits with status 1.

## Usage

```bash
S=.skills/ood-annotation-toolkit/scripts

# 1. IDR using ss/func/fold definitions.
python $S/idr_ood.py --test ssp_test.csv --out ssp_idr.csv \
    --mode residue-ratio --threshold 0.3 --cache output/idr_cache.csv
# Paired-chain ppi IDR.
python $S/idr_ood.py --test ppi_test.csv --seq-col aa_seq1,aa_seq2 \
    --mode region-ratio --clean map --out ppi_idr.csv --cache output/idr_cache.csv

# 2. kcat value-bin LongTail.
python $S/value_bin_longtail_ood.py --train kcat_train.csv --test kcat_test.csv \
    --out kcat_test.csv --bin-start -6.5 --bin-end 7.0 --bin-width 0.5 \
    --clip -6,6 --thresholds 50,100 --col-prefix OOD_LongTail_KcatBin

# 3. Finalize NaN then InD; --dry-run reports without writing.
python $S/finalize_ood_columns.py --test-csv cath_test.csv --mode both --dry-run
python $S/finalize_ood_columns.py --test-csv cath_test.csv --mode both

# 4. Read-only audit; specify paired-chain core columns for ppi.
python $S/check_splits.py --train ppi_train.csv --val ppi_val.csv --test ppi_test.csv \
    --seq-cols aa_seq1,aa_seq2 --struct-cols struct_file1,struct_file2
```

## Validation

- IDR baseline: fold residue-ratio @0.3, strip, must reproduce `splits/cath_test.csv` idr_ratio/OOD_IDR. Cached reruns finish in seconds.
- On current fold/ppi, finalizer --dry-run should report unchanged InD and no NaN clearing; both compliant as measured 2026-08-27.
- check_splits.py passes fold and ppi, measured 2026-08-27.
- Writing scripts (value_bin/finalize): back up, update only target columns, and diff after writing.

## Cautions

1. **Do not mix IDR definitions:** ppi region-ratio and other tasks' residue-ratio are not interchangeable. One cache stores both and can be shared, but flags must use the task's historical mode.
2. **metapredict environment:** not available in all Python environments; `/root/miniconda3/bin/python` has 3.0.2 on the original machine. Missing libraries produce an error suggesting cached replay.
3. **Duplicate unique_id**, e.g. multi-substrate kcat: this toolkit writes by row without ID remapping. See the EC skill for groupby cautions.
4. **Excluded analysis scripts:** migrate_column_schema.py performed nine completed one-off renames and had no reuse value. survey_columns.py's cross-task column-frequency survey was not incorporated (check_splits.py audits three files for one task only). Neither belongs to this toolkit; originals were deleted.
5. Use `--nondefault false` for historical optimal_ph LongTail-bin compatibility; new tasks use default nan.
