# POOD-Benchmark — Protein Out-of-Distribution Benchmark

A benchmark suite for training and evaluating deep-learning models on protein tasks
**under distribution shift**. It bundles 9 tasks (12 sub-datasets) spanning
protein-level regression/classification, protein–ligand tasks, protein–protein
interaction, and residue-level prediction. Every test set is annotated with
multiple **out-of-distribution (OOD)** axes so models can be evaluated both
in-distribution (ID) and across controlled distribution shifts. An additional
evaluation-only set of **designed proteins** (`protein-design-archive`) extends
the fold-classification and EC tasks.

Each task ships two things:

- `splits/` — `train` / `val` / `test` CSV files (the training/evaluation data).
- `pdbs/` — the protein structure files (`.pdb` / `.cif`) referenced by each row.

---

## 1. Repository layout

This repository contains the **data-processing documentation and code** — not the
data files themselves (`splits/`, `pdbs/`, and intermediate `output/` are not
tracked in git; see §2 for how to obtain them).

```
├── benchmark_analysis.md            # dataset sizes, OOD coverage matrix, data-quality reports
├── datasets/
│   ├── <task>/
│   │   ├── DATA_PROCESS.md          # full processing pipeline: sources, filtering, splits, OOD thresholds
│   │   └── scripts/                 # numbered pipeline scripts (01_..., 02_..., ...)
│   └── ...                          # one directory per task (see §3)
└── .skills/                         # shared pipeline modules reused across tasks
    ├── cath-ted-domain-annotation/  # CATH+TED domain annotation → topology/superfamily holdout
    ├── ec-function-ood-annotation/  # EC annotation → NewEC / LongTailEC OOD axes
    ├── homology-ood-annotation/     # sequence (MMseqs2) & structure (Foldseek) homology OOD
    ├── ood-annotation-toolkit/      # IDR, value-bin long-tail, column finalization, QA checks
    ├── struct-label-alignment/      # per-residue struct_label alignment
    └── pood-benchmark-process/      # overall pipeline conventions
```

**To reproduce a dataset from scratch**, follow the task's `DATA_PROCESS.md` and
run its numbered scripts in order; shared steps call the modules under `.skills/`.
Raw data sources must be downloaded separately (links and versions are given in
each `DATA_PROCESS.md`):

- RCSB PDB structures and SIFTS residue-level mappings
- CATH v4.4 domain classification and TED domain annotations
- UniProtKB metadata / sequences, UniRef50
- AlphaFold predicted structures (for entries without experimental structures)
- PDBbind 2020 (LBA), DIPS-Plus (PPIS), Protein Design Archive (designed proteins)

---

## 2. Download

Pre-built `splits/` + `pdbs/` for all tasks: **<https://huggingface.co/datasets/mqwang1221/poor_benchmark>**

---

## 3. Tasks and statistics

Sizes are `train / val / test`.

### 3.1 Protein-level regression

| Task | Dir | Input | `label` | Sizes | Suggested metric |
|------|-----|-------|---------|-------|------------------|
| Enzyme turnover (kcat) | `enzyme_kinetics_prediction` | protein `aa_seq` + **ligand** (reactant) | `log10(kcat)` (float) | 15,479 / 2,150 / 5,002 | RMSE, Pearson/Spearman r |
| Enzyme optimal pH | `enzyme_optimal_ph` | protein `aa_seq` | optimal pH, float (≈1.5–12.5) | 7,794 / 866 / 2,955 | RMSE/MAE, Pearson r |
| Ligand binding affinity (LBA) | `ligand_binding_affinity` | protein `aa_seq` + **ligand** | binding affinity −log₁₀(Kd/Ki), i.e. pKd/pKi (float) | 8,366 / 930 / 2,367 | RMSE, Pearson r |

**How to train:** encode the protein (sequence and/or structure); for kcat and LBA
also encode the ligand (§5) and concatenate; regress the scalar `label` with an
MSE/Huber loss.

### 3.2 Protein-level classification

| Task | Dir | `label` | Classes | Sizes | Suggested metric |
|------|-----|---------|---------|-------|------------------|
| Fold classification | `fold_classification` | CATH `Class.Architecture.Topology` string (e.g. `3.40.50`) — **load as `str`** | 1,163 | 43,473 / 6,132 / 12,911 | accuracy / macro-F1 |
| Enzyme Commission (EC) | `func_prediction` (`ec_*`) | `;`-separated EC numbers | multi-label | 15,389 / 1,687 / 5,047 | Fmax, AUPR |
| GO Biological Process | `func_prediction` (`go_bp_*`) | `;`-separated GO terms | multi-label | 15,232 / 1,742 / 5,142 | Fmax, AUPR |
| GO Cellular Component | `func_prediction` (`go_cc_*`) | `;`-separated GO terms | multi-label | 13,878 / 1,730 / 5,613 | Fmax, AUPR |
| GO Molecular Function | `func_prediction` (`go_mf_*`) | `;`-separated GO terms | multi-label | 18,179 / 2,084 / 5,609 | Fmax, AUPR |

**How to train:** encode the protein; for fold use a softmax over 1,163 classes
(cross-entropy); for EC/GO build the label vocabulary from the **training** labels
and train a multi-label classifier (binary cross-entropy per label).

### 3.3 Protein–protein interaction (pair classification)

| Task | Dir | Input | `label` | Sizes | Suggested metric |
|------|-----|-------|---------|-------|------------------|
| PPI | `ppi_prediction` | two chains `aa_seq1` + `aa_seq2` (+ `struct_file1/2`) | binary interaction (0/1) | 58,968 / 8,394 / 16,264 | AUROC, AUPR, accuracy |

**How to train:** encode each chain (shared encoder), combine the two embeddings
(concat / difference / attention), predict interaction with binary cross-entropy.

### 3.4 Residue-level (token) prediction

`label` is a per-residue array `[v v v ...]` with `len(label) == len(aa_seq)`;
`-1` marks masked residues to ignore in loss and metrics.

| Task | Dir | Per-residue `label` | Sizes | Suggested metric |
|------|-----|---------------------|-------|------------------|
| Secondary structure (8-state) | `ss_prediction` | DSSP state `0–8` (`H0 G1 I2 E3 B4 T5 S6 P7 C8`), `-1` unresolved | 29,622 / 3,317 / 8,372 | Q8 (per-residue accuracy) |
| PPI interface site | `ppis_prediction` | interface residue 1 / non-interface 0 (`-1` unlabeled) | 6,043 / 672 / 1,980 | AUPR, AUROC, F1 |
| Ligand binding site (LBS) | `ligand_binding_site` | binding residue 1 / 0; **+ ligand** (`ligand_smiles` / `ligand_ecfp4`) | 71,093 / 7,938 / 17,225 | AUPR, AUROC, F1 |

**How to train:** produce a per-residue representation (sequence and/or structure
encoder) and a per-residue head. For LBS, condition the residue head on the ligand
embedding (§5). Interface/binding-site labels are highly imbalanced — prefer
AUPR/F1 and mask `-1` positions.

### 3.5 Designed-protein test sets (evaluation only)

`protein-design-archive/splits/` provides extra test sets built from designed
proteins — no train/val, use the corresponding task's training data:

| File | Task | Rows | Schema |
|------|------|-----:|--------|
| `cath_design.csv` | fold classification | 122 | same columns as `cath_test.csv` |
| `ec_design.csv` | EC | 31 | same columns as `ec_test.csv` (without the CATH-holdout / LongTail-EC columns) |

---

## 4. Common CSV schema

| Column | Present in | Meaning |
|--------|-----------|---------|
| `unique_id` | all | sample id (unique within a task; no overlap across train/val/test) |
| `aa_seq` | single-chain tasks | amino-acid sequence (model input) |
| `aa_seq1` / `aa_seq2` | PPI | the two chains of a protein pair |
| `struct_file` | single-chain tasks | structure **file name** (no directory) located under the task's `pdbs/`; extension = file type |
| `struct_file1` / `struct_file2` | PPI | structure file names of the two chains |
| `label` | all | task target (see per-task sections; parsing rules below) |
| `ligand_smiles`, `ligand_ecfp4` | kcat, LBA, LBS | ligand representations (see §5) |
| `pair_type` | PPI (train/val/test) | pair provenance: `bio` (positive, label 1), `xtal` (crystal-packing negative), `random` (random-pairing negative); negatives are balanced 1:1 with positives per split |
| `struct_label` | SSP, PPIS, LBS (train/val/test) | per-residue label aligned to the `.pdb` ATOM records (see §6) |
| `Default` + `InD` + `OOD_*`, `seq_Redundancy_*`, `TM-score_*` | **test only** | OOD annotations (see §7) |
| `idr_ratio`, `binding_site_ratio` | some test sets | numeric auxiliaries behind `OOD_IDR` (fraction of disordered residues) and LBS `OOD_LongTail` (fraction of binding residues) |

**Loading `label` correctly** — the format depends on the task:

```python
import pandas as pd, numpy as np

# (a) Regression labels (kcat, LBA, pH): plain float
#     -> read normally, label is a float.

# (b) Single-label classification (fold/CATH): a CATH "Class.Architecture.Topology"
#     code e.g. "3.40.50". MUST be read as string, otherwise pandas may mangle
#     numeric-looking codes.
df = pd.read_csv("fold_classification/splits/cath_test.csv", dtype={"label": str})

# (c) Multi-label classification (EC, GO): ";"-separated tokens
labels = df["label"].str.split(";")          # e.g. "3.2.1.-;2.7.7.-" or "GO:0006412"

# (d) Binary (PPI): int 0/1.

# (e) Per-residue arrays (SSP, PPIS, LBS; also struct_label): "[v v v ...]" —
#     bracketed, space-separated ints, len(label) == len(aa_seq). -1 marks masked
#     (unresolved / unlabeled) residues to be ignored in the loss/metrics.
def parse_residue_labels(s):
    return np.array(s.strip("[]").split(), dtype=int)
```

---

## 5. Ligand featurization (kcat, LBA, LBS)

The three ligand-aware tasks provide **both** a SMILES string and a precomputed
ECFP4 fingerprint, so you can pick either input modality:

- **`ligand_smiles`** — canonical SMILES. Feed directly to a molecular GNN /
  SMILES transformer, or featurize with your own scheme (RDKit descriptors,
  learned encoders, etc.). Multi-molecule reactants (kcat) are joined with `.`.
- **`ligand_ecfp4`** — an **ECFP4 / Morgan fingerprint (radius = 2, 2048 bits)**
  stored as the **list of active (ON) bit indices**, e.g. `[31, 52, 55, ...]`.
  Reconstruct a fixed-length binary vector before feeding an MLP:

```python
def ecfp4_to_vector(s, n_bits=2048):
    v = np.zeros(n_bits, dtype=np.float32)
    idx = [int(x) for x in s.strip("[]").replace(",", " ").split()]  # handles both delimiters
    v[idx] = 1.0
    return v
```

A typical ligand-aware model concatenates a **protein embedding** (from `aa_seq`
and/or the structure) with a **ligand embedding** (from `ligand_smiles` via a GNN,
or from `ecfp4_to_vector(ligand_ecfp4)` via an MLP) and predicts the target.

---

## 6. Structure input

`struct_file` (or `struct_file1/2`) names a file inside the task's `pdbs/`
directory; the extension (`.pdb` / `.cif`) is the format. Sequence-only models
can ignore these; structure-based models (GNN / GVP / structure encoders) load
the coordinates from `pdbs/<struct_file>`.

**Residue-level tasks, structure input:** the residue-level tasks
(`ss_prediction`, `ligand_binding_site`, `ppis_prediction`) carry a
`struct_label` column in **train/val/test alike**: one label per residue of the
`.pdb` file, in ATOM-record order, so it aligns 1:1 with a structure encoder's
per-residue outputs (`-1` = ignore index). Structure models only need
`pdbs/<struct_file>` + `struct_label`.

---

## 7. Training vs. OOD evaluation

- **Train on** `*_train.csv`, tune on `*_val.csv`. These contain only the core
  columns — **no OOD annotations** — and are your in-distribution training data.
- **Evaluate on** `*_test.csv`. Besides the label, the test set carries boolean
  OOD-annotation columns; report ID performance and per-axis OOD performance:

| Column(s) | OOD axis |
|-----------|----------|
| `Default` | `True` = sample produced by the normal split process; `False` only for samples deliberately set aside **before** splitting (e.g. a pre-split extreme-length pool) |
| `InD` | `True` = no OOD flag set on any axis (the purest ID subset) |
| `seq_Redundancy_90 … seq_Redundancy_30` | max sequence identity to train+val below 90%…30% (sequence-similarity shift; nested: `<30` ⊆ … ⊆ `<90`) |
| `TM-score_0.9 … TM-score_0.3` | max structural TM-score to train+val below 0.9…0.3 (fold-similarity shift; nested) |
| `OOD_Orphan` | no significant homolog in train+val (a subset of `seq_Redundancy_30`) |
| `OOD_ExtremeShort` / `OOD_ExtremeLong` | chain < 60 aa / > 1000 aa |
| `OOD_IDR` | high intrinsic-disorder content |
| `OOD_FoldHoldout` / `OOD_SuperfamilyHoldout` | held-out CATH topology / superfamily (based on CATH+TED domain annotation; the fold task itself uses only `OOD_SuperfamilyHoldout`, since topology is its label) |
| `OOD_NewEC_L4` / `OOD_NewEC_L3` | **all** EC classes (level 4 / 3) of the protein are absent from the reference set (all-not-in). The fold task uses `OOD_NewEC` instead — an EC7 (translocase) L1 holdout |
| `OOD_LongTail_EC_L3_le5/le10`, `OOD_LongTail_EC_L4_le5/le10` | all EC classes (level 3 / 4) are rare in the reference set (train frequency ≤ 5 / ≤ 10, including 0) |
| `OOD_LongTail` | rare **predicted-label** category (task-specific: rare CATH topology, GO/EC label, SS state, sparse binding/interface, …) |
| `OOD_LongTail_KcatBin_le50/le100` (kcat), `OOD_LongTail_pHbin_le50/le100` (pH), `OOD_LongTail_AffBin_le50/le100` (LBA) | regression value falls in a bin (histogram over train values) with ≤ 50 / ≤ 100 training samples |
| `OOD_Combinatorial` | label combination unseen in train+val (multi-label tasks: EC, GO) |

**Axis availability by dataset** (✓ = column present in the test set):

| OOD axis | kcat | pH | CATH | EC | GO-BP | GO-CC | GO-MF | LBA | LBS | PPI | PPIS | SSP |
|---|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| `OOD_ExtremeShort` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| `OOD_ExtremeLong` | ✓ | ✓ | – | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| `seq_Redundancy_90…30` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| `TM-score_0.9…0.3` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| `OOD_Orphan` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| `OOD_IDR` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| `OOD_FoldHoldout` | ✓ | ✓ | – | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| `OOD_SuperfamilyHoldout` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| `OOD_NewEC` (EC7 holdout) | – | – | ✓ | – | – | – | – | – | – | – | – | – |
| `OOD_NewEC_L4` | ✓ | ✓ | – | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| `OOD_NewEC_L3` | ✓ | ✓ | – | – | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| `OOD_LongTail` (label) | – | – | ✓ | ✓ | ✓ | ✓ | ✓ | – | ✓ | – | ✓ | ✓ |
| `OOD_LongTail_*Bin_le50/le100` | ✓ | ✓ | – | – | – | – | – | ✓ | – | – | – | – |
| `OOD_LongTail_EC_L3_le5/le10` | ✓ | ✓ | ✓ | – | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| `OOD_LongTail_EC_L4_le5/le10` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| `OOD_Combinatorial` | – | – | – | ✓ | ✓ | ✓ | ✓ | – | – | – | – | – |

`Default` and `InD` are present in every test set. Samples set aside before
splitting (`Default=False`, e.g. extreme-length pools) do not participate in the
other OOD axes — those columns are `NaN` for them (load with
`pd.read_csv(..., dtype=...)` treating empty as NaN, or map empty strings to
`False` for counting). Exact thresholds and per-task definitions are in each
task's `DATA_PROCESS.md`; per-axis sample counts are aggregated in
`benchmark_analysis.md`.

---

## 8. Data notes

- **Sequences** use the 20 standard amino acids plus IUPAC ambiguity codes
  `X B Z J U O` (rare); make sure your tokenizer has an unknown/`X` fallback.
  In LBA, `|` separates the chains of a multi-chain complex — split on it if
  your model consumes single chains.
- **Reproducibility:** every task directory has a `DATA_PROCESS.md` documenting
  sources, filtering, splits, and exact OOD thresholds. `benchmark_analysis.md`
  aggregates dataset sizes, the OOD coverage matrix, and data-quality reports.

## 9. Citation

If you use this benchmark, please cite the accompanying paper. *(citation to be added)*
