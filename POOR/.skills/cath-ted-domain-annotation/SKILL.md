---
name: cath-ted-domain-annotation
description: Generate unified CATH+TED topology/superfamily domain annotations for arbitrary PDB-chain tasks and compute POOD-Benchmark OOD_FoldHoldout / OOD_SuperfamilyHoldout. Use to add or recompute structural classification and CATH holdout columns for ppis, ligand_binding_affinity, func_prediction, ss_prediction, or other datasets.
---

# CATH+TED domain annotation and holdout OOD

## Description

Unified domain annotation: **prioritize CATH v4.4 experimental labels and fill gaps with putative TED labels (union semantics)**. Produce chain-level topology (C.A.T) and superfamily (C.A.T.H) labels, then optionally compute `OOD_FoldHoldout` / `OOD_SuperfamilyHoldout`.

## When to use

- Add topology/superfamily annotations to PDB-chain samples.
- Compute/recompute CATH holdout columns.
- Improve CATH-only coverage for chains absent from CATH, such as newer PDB structures.

## Canonical definitions: consistent across tasks

| Decision | Convention | Rationale |
|:---|:---|:---|
| Label priority | CATH v4.4 experimental labels; TED fills gaps by union | TED uses computational foldseek/foldclass predictions and misses ~6–7% of existing CATH labels; it cannot replace CATH |
| TED domain filtering | Retain domains whose chopping overlaps the chain's UniProt interval by **≥50% of domain length** | TED domains use full-length UniProt; exclude domains outside the chain. Sensitivity at 25%/75% was robust |
| TED levels | `H` → C.A.T + C.A.T.H; `T` → C.A.T only | T lacks superfamily information and cannot support SuperfamilyHoldout |
| Holdout reference | **Training only**, with training labels also augmented by CATH+TED | Augmenting only test makes TED-only labels falsely unseen and inflates holdouts |
| Criterion | Any unseen label → True; **unannotated → False** | Matches existing task scripts |
| Sample unit | Input chain list determines granularity: one row per chain, or union all chains per PDB for complex tasks such as ligand | Matches archived `ligand_binding_affinity/output/scripts_archive/16_cath_ec_ood.py` and add_cath_holdout_ood.py |

## Sources

| File | Purpose |
|:---|:---|
| `data/CATHv44/cath-classification-data/cath-domain-list.txt` | Official CATH v4.4 domain classification; experimental baseline |
| `data/CATHv44/cath-classification-data/cath-domain-boundaries-seqreschopping.txt` | SEQRES domain boundaries, converted to UniProt coordinates |
| `data/sifts/pdb_chain_cath.tsv.gz` | PDB+chain → UniProt SP_PRIMARY; preferred source |
| `data/sifts/sifts_chain_uniprot.tsv.gz` | PDB+chain → UniProt and SP_BEG/SP_END intervals for TED overlap filtering |
| `data/ted/ted_api_cache.json` | **Project-wide shared TED API cache**; already contains 9,889 accessions |
| `https://ted.cathdb.info/api/v1/uniprot/summary/<acc>` | Query cache misses; ~10 req/s is safe with concurrency |

> TED bulk files on Zenodo record 13908086 are frequently rate-limited (403) on the original network. API plus cache is the verified primary route.

## Usage

### 1. Annotate: cath_ted_annotate.py

Input mode is detected by columns:

- **Chain:** `pdb`, `chain`; CATH name matching plus 50% TED domain-overlap filtering.
- **UniProt:** `uniprot`, for full-length samples such as AFDB-keyed datasets. Reverse-map SP_PRIMARY from the union of `pdb_chain_cath` and `sifts_chain_uniprot` to chains, match cath-domain-list names, and union all PDB-chain labels for that accession. TED uses all annotated accession domains **without overlap filtering**, because the sample is full-length.

```bash
python .skills/cath-ted-domain-annotation/scripts/cath_ted_annotate.py \
    --chains my_chain_list.csv --out output/my_labels.csv          # Chain mode
python .skills/cath-ted-domain-annotation/scripts/cath_ted_annotate.py \
    --chains my_acc_list.csv --out output/my_labels.csv --offline  # UniProt mode
```

Chain-mode output appends:

| Column | Description |
|:---|:---|
| `cath_topos` / `cath_sfs` | CATH experimental labels, semicolon-separated sets |
| `ted_topos` / `ted_sfs` | TED additions after overlap filtering |
| `topos` / `sfs` | Union; **use these for downstream holdout computation** |
| `uniprot` | Chain's UniProt accession |
| `annotation_source` | `cath` / `ted` / `both` / `none` |
| `cath_domains` | CATH domains with UniProt intervals, e.g. `1qqfA00:5-120,200-260`; empty after domain ID if unmapped |
| `ted_domains` | Retained TED domains and chopping in UniProt coordinates |
| `cath_chain_cov` / `aug_chain_cov` | Fraction of chain UniProt interval covered by CATH / CATH+TED domains |

**CATH boundaries** from `data/CATHv44/cath-classification-data/cath-domain-boundaries-seqreschopping.txt` use SEQRES numbering and are converted through SIFTS segment mappings `RES_BEG/RES_END ↔ SP_BEG/SP_END`. CATH labels still use name matching regardless of boundaries; boundaries support interval/coverage analysis. `--domains-out domains.csv` adds domain details (`mapped_frac`, `in_chain_frac`, `kept`). CATH domains with `in_chain_frac <0.5` generally lack a SIFTS chain-to-UniProt mapping.

Key options: `--overlap` (0.5), `--threads` (8), `--offline` (cache only), `--cache` (project shared cache).

### 2. Holdouts: cath_ted_holdout.py

```bash
# Single file: train/val/test identified by split; ';' separates multiple memberships.
python .skills/cath-ted-domain-annotation/scripts/cath_ted_holdout.py \
    --labels output/my_labels.csv --split-col split --out output/test_flagged.csv

# Separate test and reference files:
python .skills/cath-ted-domain-annotation/scripts/cath_ted_holdout.py \
    --labels test_labels.csv --ref-labels trainval_labels.csv --out test_flagged.csv
```

Output test rows only, appending True/False `OOD_FoldHoldout` and `OOD_SuperfamilyHoldout`.

> **Required:** task wrappers **must call `cath_ted_holdout.py`** when recomputing holdouts for split files. Do not independently reimplement reference scope, any semantics, or unannotated→False logic; this script is authoritative.

### 3. Complex-level tasks: PDB union

First union chain `topos`/`sfs` by sample, then run holdouts. See groupby in `datasets/ligand_binding_affinity/output/scripts_archive/16_cath_ec_ood.py`.

## Validation

- Inspect ≥20 chains: retained TED labels must be a subset of raw accession API labels after overlap filtering.
- Coverage/holdout counts should be comparable to measurements on 2026-08-26:

| Task | CATH-only topology coverage | CATH+TED | FoldHoldout | SuperfamilyHoldout |
|:---|:---|:---|:---|:---|
| ppis test (Default) | 92.8% | 97.4% | 64 →65 | 231 →235 |
| ligand test | 70.8% | 90.9% | 24 →38 | 45 →73 |

> Ligand values are final after the 2026-08-26 split update, including chain-seeding fixes (+1/+1) and corrected comma-separated dual TED labels across 82 domains, unioning both prediction methods.

- CATH/TED agreement for doubly annotated samples: exact agreement 88–91%; nonempty intersection >99%. Investigate mapping errors if substantially different.

## Cautions

1. **Symmetric augmentation:** annotate train/val/test together so the reference also receives TED additions.
2. TED uses AFDB **v4**; map local v6 structures by accession without version handling.
3. TED labels are putative. State in DATA_PROCESS.md that CATH experimental labels take priority and TED fills gaps.
4. About 8% of accessions are absent from TED (`not_found`, mostly new UniProt entries); these remain unannotated.
5. TED and CATH domain boundaries differ; extra TED labels on CATH-annotated chains (333 ppis chains) partly reflect boundary differences. Interpret complementarity accordingly.
6. Occasionally TED `cath_label` contains comma-separated dual labels from foldseek/foldclass: 82 domains. Union them. **Direct `split(".")` creates false labels such as `3.30.160,3`**, previously causing false holdout flips; protection is implemented.
7. **UniProt-mode caveat:** PDB constructs for one accession may cover different regions or mutants; reverse-mapping unions can slightly overcover. Unioning both partners in pair tasks such as ppi amplifies sensitivity of any semantics.
