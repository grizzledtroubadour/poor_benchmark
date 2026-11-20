---
name: homology-ood-annotation
description: Compute standardized POOD-Benchmark search-based homology OOD columns for arbitrary tasks: mmseqs2 sequence homology (seq_Redundancy_90..30 + OOD_Orphan) and foldseek structural homology (TM-score_0.9..0.3). Search test→train+val, supporting m8 reuse and paired chains (aa_seq1/aa_seq2). Use to add or recompute sequence/structure homology OOD.
---

# Search-based homology OOD: mmseqs2 sequences and foldseek structures

## Description

Unified search-based homology workflows: seq_homology_ood.py uses mmseqs2 easy-search; struct_homology_ood.py uses foldseek easy-search. Both search test against train+val, take maximum similarity, and apply threshold ladders. Complements classification holdouts in cath-ted-domain-annotation.

## When to use

- Add/recompute `seq_Redundancy_90`–`seq_Redundancy_30` and `OOD_Orphan`.
- Add/recompute `TM-score_0.9`–`TM-score_0.3`.
- Replay existing mmseqs2/foldseek m8 results for reruns or threshold changes.

## Canonical definitions

| Decision | Convention | Rationale |
|:---|:---|:---|
| Direction | **query=test, target=train+val** | Matches fold annotate_test_ood.py; reversing loses the per-test multiple-homolog interpretation |
| Sequence search | mmseqs2 `easy-search -s 7.5 --format-output query,target,pident`, default e-value ≤1e-3 | Matches fold; hits are significant. **-s 7 (lbs) versus 7.5 has negligible observed effect:** among 8,864 lbs chains, seq_Redundancy ≥50 unchanged, seq40 Δ=3, seq30 Δ=9, Orphan Δ=18/932. Differences concern weak pident 22.7–42.0 hits near detection limits; see lbs `output/verify/sens_s{7,75}.m8` |
| Sequence flag | Maximum pident per test sample, no hit=0.0; `seq_Redundancy_<T>` means max pident <T, T=90,80,...,30 | 10% threshold ladder |
| Orphan | **No significant hit →OOD_Orphan=True** | Default e-value ≤1e-3 makes no hit equivalent to repaired lbs best e-value >1e-3. Ensures **OOD_Orphan ⊆ seq_Redundancy_30 ⊆ ... ⊆ seq_Redundancy_90**, asserted internally |
| Structure search | foldseek `easy-search --format-output query,target,alntmscore`; normalize query IDs by removing paths and `.pdb` | Matches fold normalization |
| Structure flag | Maximum alntmscore per test structure, no hit=0.0; `TM-score_<T>` means max <T, T=0.9,...,0.3 | 0.1 threshold ladder |
| Paired-chain merge | Search chains separately; **any semantics**: (A<T) or (B<T); Orphan means A or B has no hit; analogous for TM | Matches ppi 10_recompute_all_ood.py `(A_max_id<thr)\|(B_max_id<thr)` and `(A==0)\|(B==0)`. **Do not take the maximum across chains first**: either chain without homology makes the pair OOD. Optional `--merge all` |
| Self-alignment | Skip q==t while parsing m8 | Matches repaired lbs; avoids fixed 100 pident for overlapping test/reference IDs |
| Cache | Existing `--m8-cache` →replay without external tools; absent →search and save | Search is expensive; classification replay is cheap |

## Historical implementations and differences

| Task | Implementation | Difference |
|:---|:---|:---|
| fold_classification | `scripts/05_annotate_test_ood.py`, **baseline** | None; replay `output/ood_mmseqs_test_vs_trainval.m8` / `ood_foldseek_test_vs_trainval.m8` directly |
| ligand_binding_site | Archived `output/scripts_archive/24_fix_seq_homology_easysearch.py` | Uses -s 7, adds evalue output, Orphan best e-value >1e-3 (equivalent to no hit); deduplicates queries by chain_key `{pdb}_{chain}`, not unique_id |
| ppi_prediction | `scripts/10_recompute_all_ood.py` | Both tools unified to **easy-search** on 2026-08-27 using skill scripts. Replay caches `output/ood_mmseqs_test_vs_trainval_easysearch.m8`, `ood_foldseek_test_vs_trainval_easysearch.m8`. Earlier mmseqs used `search --max-seqs 10 -e 100 -c 0.0` (seq30 515→3,653; Orphan 211→3,017). Earlier foldseek used `search -a -e inf`; a gap with no max_tmscore in [0.3,0.4) made TM0.3 and TM0.4 identical. After recomputation: TM0.9 4,827→8,172; TM0.3 1,610→943. Old caches are provenance only. Paired-chain any merge is the default here |
| func_prediction | Archived `scripts/orphan_mark.py`; current easy-search -s 7 in task DATA_PROCESS.md §4.3 | Historical reference was UniRef50, not train+val; Orphan also required fident ≥0.30. Different semantics: do not replay old OOD_Orphan with this tool |
| Legacy lbs | `scripts/11_build_ood_splits.py`, now easy-search -s 7 | Earlier `mmseqs search -e 100 -c 0.0` weak alignments inflated maximum pident, yielding Orphan 1,285 ≫ SeqID<30% 134. Fixed by archived 24_fix |

> Prefer task m8-cache replay when recomputing existing tasks. New searches use this skill's defaults; document parameter differences from historical implementations in DATA_PROCESS.md.

## Usage

### 1. Sequence homology

```bash
# Single-chain task: fold definition.
python .skills/homology-ood-annotation/scripts/seq_homology_ood.py \
    --test cath_test.csv --ref cath_train.csv cath_val.csv \
    --out test_flagged.csv --m8-cache output/ood_mmseqs_test_vs_trainval.m8 \
    --threads 16

# Paired-chain task: ppi any merge.
python .skills/homology-ood-annotation/scripts/seq_homology_ood.py \
    --test ppi_test.csv --ref ppi_train.csv ppi_val.csv \
    --seq-col aa_seq1,aa_seq2 --out ppi_flagged.csv --threads 16

# Replay without mmseqs2; m8 query keys must match unique_id.
python .skills/homology-ood-annotation/scripts/seq_homology_ood.py \
    --test cath_test.csv --ref cath_train.csv \
    --out replay.csv --m8-cache existing.m8
```

Output retains all test rows and adds max_pident (paired: one `max_pident_<seq-col>` per chain), seven True/False-string seq_Redundancy flags, and OOD_Orphan.

Options: `--seq-col` (aa_seq), `--thresholds` (90,80,70,60,50,40,30), `--merge {any,all}`, `--sensitivity` (7.5), `--id-col`, `--tmp-dir`.

### 2. Structural homology

```bash
# Directory mode: all *.pdb/*.cif.
python .skills/homology-ood-annotation/scripts/struct_homology_ood.py \
    --test-structs test_pdbs/ --ref-structs trainval_pdbs/ \
    --out tm_flags.csv --m8-cache output/ood_foldseek_test_vs_trainval.m8

# Manifest mode: one structure path per line.
python .skills/homology-ood-annotation/scripts/struct_homology_ood.py \
    --test-structs test_pdbs.txt --ref-structs trainval_pdbs.txt --out tm_flags.csv
```

One output row per test structure: `query_id, max_alntmscore, TM-score_0.9 ... TM-score_0.3`.

### 3. Write splits: caller responsibility

As in the EC skill: back up, read all columns with dtype=str/keep_default_na=False, overwrite target columns only, and mask Default=False by task convention (toolkit finalizer --mode nan). Diff all cells after writing.

## Validation

- Internal monotonicity assertions: Orphan ⊆ seq30 ⊆ ... ⊆ seq90, and analogous TM columns; violations fail immediately.
- Fold regression baseline: 12,911 test rows, m8 replay measured 2026-08-27, agrees row by row with all seven seq flags +Orphan and seven TM columns in splits/cath_test.csv.
- For new tasks inspect ≥5 Orphan=True samples, confirming no significant sequence/structure hits against train+val.
- Number of queries with hits should be comparable in order of magnitude to m8 row count. If anomalously low, check key matching: paired keys use `unique_id|<seq-col>`; historical protein_id/chain_key caches are not interchangeable.

## Cautions

1. **Cache keys:** generated query keys are unique_id, or `unique_id|<seq-col>` for paired chains. Fold caches match directly; historical ppi protein_id and lbs chain_key differ. Verify key spaces before replay or all queries may become Orphan.
2. **Threshold names:** remove trailing zeros in columns: 90→seq_Redundancy_90; 0.9→TM-score_0.9. Align custom thresholds with existing names.
3. foldseek requires directories; symlink structures as in fold to avoid ~10 GB copying. `--tmp-dir` retains intermediate FASTA/m8 files.
4. Missing mmseqs2/foldseek on PATH yields explicit errors suggesting m8 replay.
5. This skill computes all rows rather than masking Default=False extreme-length rows. Toolkit finalize_ood_columns.py applies the NaN convention.
