---
name: struct-label-alignment
description: Generate struct_label columns for residue-level task split CSVs (ss_prediction, ligand_binding_site, ppis_prediction, etc.). Assign a residue label to each PDB residue in ATOM file order, allowing structure models to train/evaluate from pdb plus struct_label. Use to map sequence-position labels to structure residues.
---

# struct-label-alignment

Map CSV residue labels (`label`, length equal to `aa_seq`) from sequence positions to **PDB structure residue order** in `struct_label`.

## Semantics

- Enumerate PDB residues in ATOM file order by first occurrence of `(chain, resseq, icode)`, **including residues without CA atoms**, consistent with parsers such as Biopython.
- Globally align the PDB residue sequence (three-letter to one-letter conversion; modified residues use script mapping `RES3`, unknown residues become `X`) to `aa_seq` using Biopython PairwiseAligner: match=2, mismatch=-1, open=-2, extend=-0.5. Each structure residue inherits its aligned sequence position's label.
- Equal-length sequences use a positional fast path; substituted residues still count as aligned.
- **-1 = ignore index:** residues without a sequence mapping, such as expression tags present only in the structure, receive -1 and must be skipped during training/evaluation.
- Output uses the `label` format: `[ 8 8 3 ... ]` (brackets and space-separated values).
- Assert `len(struct_label) == number of PDB ATOM residues`.

## Usage

```bash
python .skills/struct-label-alignment/add_struct_label.py \
    --csv datasets/<task>/splits/<split>.csv \
    --pdb-dir datasets/<task>/pdbs \
    [--seq-col aa_seq] [--label-col label] [--struct-col struct_file] \
    [--out-col struct_label] [--workers 32] [--dry-run]
```

- Inspect `--dry-run` statistics (rows_with_-1 / total_-1 / rows_with_nonX_mismatch) before writing.
- CSV writes modify the file **in place**, appending the column. Back up first to `<task>/output/backup_before_struct_label/`.
- CSV I/O: `pd.read_csv(dtype=str, keep_default_na=False)` and `to_csv(index=False, lineterminator="\n")`, already implemented.

## Applied datasets

- Train/val/test for ss_prediction, ligand_binding_site, and ppis_prediction: nine CSVs, 2026-08.
- On 2026-08-28, ss_prediction `struct_label` was replaced with direct DSSP annotations transferred from the struct_native view. The script is archived at `ss_prediction/output/scripts_archive/patch_struct_label_from_native.py`; the native view at `ss_prediction/output/splits_struct_native_deprecated/`. This skill remains applicable to alignment-based generation for lbs/ppis and subsequent tasks.

## Caution

- **Do not recompute rows where `aa_seq` differs from the structure sequence using `--overwrite`.** The inherited mapping (structure sequence → `aa_seq` → main `label`) can shift alignment frames in low-complexity/repeated regions and transfer labels from incorrect positions. In ss_prediction on 2026-09-01, this changed 12,199 residue labels across 148 rows, with only 71% agreement against direct DSSP labels; see recovery in `ss_prediction/DATA_PROCESS.md` §5.5. The positional fast path is safe only when `aa_seq` is the structure-extracted sequence, as for short-chain SSP rows restored in 2026-09. Freeze columns with more authoritative sources, such as transferred direct DSSP annotations, rather than overwriting them.
