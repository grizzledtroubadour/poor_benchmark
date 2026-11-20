# Fold Classification: Dataset Processing

> Complete linear workflow from the CATH S95 list and local PDB structures to final OOD splits, organized for one-pass reproduction from raw data.
> Location: `datasets/fold_classification/`.
>
> **Task:** classify nonredundant CATH v4.4 S95 domains at C.A.T (topology) level. String labels such as `3.40.50` cover **1,163** classes; unique_id is the CATH domain ID, e.g. `2da1A01`.

**Current split sizes, measured from files:**

| File | Rows | Description |
|------|-----:|------|
| `splits/cath_train.csv` | **43,473** | unique_id,aa_seq,struct_file,label; all 1,163 topologies |
| `splits/cath_val.csv` | **6,132** | Same schema |
| `splits/cath_test.csv` | **12,911** | Four core plus 27 annotation columns: 31 total, including OOD and InD |

---

## Stage 1: Obtain sources and extract CATH S95 domains

> Extract ATOM coordinates for 62,915 CATH v4.4 S95 domains from experimental PDB structures, writing individual domain PDBs and metadata.

### 1.1 Sources

| Data | Location | Purpose |
|------|------|------|
| `cath-classification-data/cath-domain-list-S95.txt` | `data/CATHv44/` | **Entry point:** 62,915 nonredundant S95 domains |
| `cath-classification-data/cath-domain-boundaries.txt` | `data/CATHv44/` | PDB-numbered domain boundaries, Domall CDF 2.0 |
| `cath-classification-data/cath-domain-boundaries-seqreschopping.txt` | `data/CATHv44/` | Alternative contiguous SEQRES boundaries |
| `sequence-data/cath-domain-seqs-S95.fa` | `data/CATHv44/` | Authoritative domain sequences for aa_seq |
| Experimental PDBs, mainly `.pdb.gz` with `.cif.gz` fallback | `data/pdb/` | RCSB structures from https://files.rcsb.org/download/ |

Each S95 row has 12 fields: seven-character domain_id (PDBID+Chain+DomainNum), four C.A.T.H levels, S35/S60/S95/S100 cluster numbers, count, length, and resolution (999=NMR, 1000=obsolete). **S95 numbers are not globally unique**; they are unique only within S35/S60 context. Full cluster key: `(C,A,T,H,S35,S60,S95)`.

### 1.2 Download missing PDBs — scripts/01_download_missing_for_s95.py

- Input: S95 list and local data/pdb/ index. Automatically derive output/missing_for_s95.txt; no manual step.
- Output: missing PDB/CIF gzip files in data/pdb/; output/download_missing_s95.log; output/download_failed_s95.txt if needed.
- Ten download threads; prefer PDB, fall back to CIF after 404. Of 41,403 unique S95 PDB IDs, 246 were missing; **246/246 downloads succeeded**, reaching 100% coverage.

### 1.3 Extract all domains and rebuild summary — scripts/02_extract_cath_s95_domains.py

- Input: S95 list, domain boundaries, data/pdb/.
- Output: output/cath_s95_domains/{domain_id}.pdb, later moved to pdbs/; **62,901-row** output/cath_s95_summary.csv with domain_id,pdb_id,chain,dom_num,type,C,A,T,H,S35,S60,S95,S100,length_expected,n_residues,n_atoms,resolution; output/extract_checkpoint.pkl.
- ProcessPoolExecutor with cpu_count()-1 workers. Whole-chain domains (dom_num=00) extract the entire chain; chopped domains filter ATOM records by segments in cath-domain-boundaries.txt.
- Convert extracted CIF domains to PDB and renumber atoms. Stop NMR entries at the first ENDMDL, retaining model 1.
- Checkpoint every 500 domains for resumability.
- **Rebuild the full summary after extraction** by scanning output n_residues/n_atoms and joining S95 metadata. Earlier checkpoint-loop CSV appending omitted batches; archived fix_summary.py repaired it once. The reconstruction logic is now integrated.
- Results: **62,901** successes, 14 empty domains without coordinates; approximately 15–20 minutes.

### 1.4 Final structure library

Move all 62,901 domain PDBs into pdbs/. struct_file={domain_id}.pdb resolves here; Foldseek and downstream models use this library directly.

### 1.5 CATH hierarchy: 62,901 extracted domains

| Level | Classes | Description |
|------|------|------|
| Class (C) | 5 | Broad secondary-structure composition |
| Architecture (C.A) | 43 | Large-scale arrangement |
| **Topology (C.A.T)** | **1,472** | Connectivity; **task classification level** |
| Superfamily (C.A.T.H) | 6,631 | Evolutionary homologous superfamilies; too fine, with 3,207 singletons |

Classes: 1 Mainly Alpha 13,417 (21.3%); 2 Mainly Beta 16,227 (25.8%); 3 Alpha Beta 31,613 (**50.3%**); 4 Few Secondary Structures 675 (1.1%); 6 Special 969 (1.5%). Largest superfamily: 2.60.40.10, 4,542 domains. Distribution is strongly long-tailed.

---

## Stage 2: Aggregate labels and build the sample pool

> **Script:** scripts/03_build_domain_pool.py.
> **Label:** C.A.T string, e.g. 3.40.50. Superfamily is too fine (6,631 classes; 81% too small), architecture too coarse (43); neither is used.

- Input: output/cath_s95_summary.csv and data/CATHv44/sequence-data/cath-domain-seqs-S95.fa.
- Output: **62,896-row** output/domain_pool.csv, columns unique_id,aa_seq,struct_file.
- Use authoritative CATH domain reference aa_seq, not ATOM-derived sequences that omit unresolved loops. Assert each {domain_id}.pdb exists.
- Exclude five domains with length_expected >1000: 1ej6B00,1u6gC00,3a6pA00,3egwA02,3gjxD00, based on CATH-listed lengths.
- **Additional removal, 2026-09-01:** 3w3uA01 has length_expected=999 but actual reference length 1,047. Remove it from train (43,474→43,473) and pdbs/. Its topology 1.25.10 retains 91 training representatives.
- The pool is **row-equivalent**, verified by diff, to the union of v1 C.A-level backups in output/splits_ca_backup/. Keep those for historical comparison/rollback.

---

## Stage 3: Redundancy reduction and refinement

- **Sequence redundancy:** already reduced by CATH's 95%-identity S95 representative selection; no extra clustering.
- **Pre-split filtering**, implemented in 04_build_topology_split.py:

| Operation | Count | Description |
|------|------|------|
| Remove singleton topologies | **302** classes /domains | Cannot populate both train and test |
| Reserve EC7 prior | **413** domains | All EC7-translocase carriers removed from train/val pool and assigned to test; any EC7 on multi-EC domains triggers reservation |
| Remove all-EC7 topologies | Six topologies, **77** domains | No training-pool members remain |
| Retained EC7 prior | **336** domains | 413−77, all test, **2.6%** of test |

List 379 excluded domains in output/dropped_domains.csv, with topology, superfamily, and is_newec.

**Why EC7?** Compare effects of EC-class holdouts on test composition:

| Holdout | Domains | Topologies | Test fraction |
|------------|-----:|--------------:|:----------:|
| **EC7, adopted** | **413** | **40** | **3.2%** |
| EC6 | 1,584 | 93 | 11.5% |
| EC5 | 1,818 | 125 | 13.0% |
| EC1/2/3 | 5,364–9,595 | 215–389 | 31.9–47.5% |

EC1/2/3 contain most enzymes and would move almost half the dataset to test. EC7 has the smallest impact; introduced in 2018, its mechanisms differ from classical metabolic enzymes, motivating new-function evaluation.

> **Definition:** fold OOD_NewEC is a **pre-split EC7 (L1) prior**, not the EC skill's train+val all-not-in L3/L4 definition. **Do not recompute with that skill**; fold is explicitly excluded there.

**Archived NMR repair:** a 2026-08-28 scan found 4,076/62,901 PDBs concatenating models without MODEL separators. Fix extraction in 02 to stop at the first ENDMDL. Existing files were truncated to model 1, originals backed up at output/backup_before_multimodel_fix/. The one-off repair script is archived; rescanning found zero residual duplicates.

---

## Stage 4: Split and annotate OOD

### 4.1 Topology-stratified 7:1:2 split — scripts/04_build_topology_split.py, seed=42

- Input: output/domain_pool.csv from 03, falling back to row-equivalent v1 output/splits_ca_backup/ if absent; output/cath_s95_summary.csv for C.A.T.H; data/sifts/sifts_ec_annotations.json for chain EC7 priors.
- Output: splits/cath_{train,val,test}.csv; test initially contains OOD_NewEC, with other flags added in §4.2; output/dropped_domains.csv.
- Apply Stage 3 filtering. Randomly split non-EC7 samples **7:1:2 by topology**, seed=42. If only one pool member remains, force it into train; reserved priors supply test coverage. Add all EC7 priors to test.

| Split | Domains | Notes |
|------|----------:|------|
| **Train** | **43,473** | All **1,163** topologies; class frequency min/median/max=1/6/5,584 |
| **Val** | **6,132** | — |
| **Test** | **12,911** | 336 EC7 priors, OOD_NewEC=True |
| Total | **62,516** | 62,896 pool −379 filtered −one additional 3w3uA01 |

Built-in assertions: pairwise unique_id disjointness; all 1,163 topologies present in train and test; C.A.T string labels; zero EC7 in train/val.

> Store labels as strings. Load with `dtype={'label': str}` to avoid parsing 1.10 as float and collapsing classes.

### 4.2 Test OOD annotation — scripts/05_annotate_test_ood.py

> **All 12,911 rows participate in every dependent OOD computation**, including 336 NewEC priors and 767 ExtremeShort rows. Default=True throughout: short chains are not separately reserved before splitting, but follow the normal stratified split. The Default=False NaN convention therefore has no applicable rows.

All statistics below use current splits/cath_test.csv.

#### Sequence homology: seq_Redundancy_90…30 and OOD_Orphan

- mmseqs2 `easy-search -s 7.5 --format-output query,target,pident`; default e-value significance threshold ≤1e-3.
- Query: full test; target: train+val, 49,606 domains.
- Maximum pident per query, zero for no hit; seq_Redundancy_T means max pident <T, T=90,80,…,30. Orphan means no significant hit.
- Cache: output/ood_mmseqs_test_vs_trainval.m8, reused automatically.
- Aligned with .skills/homology-ood-annotation/seq_homology_ood.py; this script serves as that skill's baseline. Replay with --m8-cache.

| Threshold | Domains | Test fraction |
|------|----------:|--------:|
| <90% | 10,257 | 79.4% |
| <80% | 9,118 | 70.6% |
| <70% | 8,114 | 62.8% |
| <60% | 7,055 | 54.6% |
| <50% | 5,586 | 43.3% |
| <40% | 3,876 | 30.0% |
| **<30%** | **2,353** | **18.2%** |

**Orphan: 1,978** domains, 15.3%. check_splits verifies Orphan ⊆seq30 ⊆…⊆seq90.

#### Structural homology: TM-score_0.9…0.3

- Foldseek `easy-search --format-output query,target,alntmscore`. Symlink query/target directories to avoid copying ~10 GB. Strip paths and .pdb suffixes from query IDs.
- Maximum alntmscore per test structure, zero for no hit. Flag max <T, T=0.9…0.3.
- Cache: output/ood_foldseek_test_vs_trainval.m8.
- Aligned with .skills/homology-ood-annotation/struct_homology_ood.py.

| Threshold | Domains | Test fraction |
|------|----------:|--------:|
| <0.9 | 4,022 | 31.2% |
| <0.8 | 1,755 | 13.6% |
| <0.7 | 839 | 6.5% |
| <0.6 | 495 | 3.8% |
| <0.5 | 348 | 2.7% |
| <0.4 | 273 | 2.1% |
| **<0.3** | **237** | **1.8%** |

12,681/12,911 test structures have Foldseek hits.

#### OOD_SuperfamilyHoldout: post-split unseen superfamilies

After the normal 7:1:2 split, flag test C.A.T.H superfamilies from cath_s95_summary.csv absent from train/val. No prior reservation or distribution modification; these ordinary split samples participate in all dependent OOD axes.

**761** domains, 5.9%, across **697** superfamilies and **223** topologies.

> Naturally held-out superfamilies are almost all tiny: the probability that all members fall into 20% test decays exponentially with size. This axis measures rare/orphan superfamilies never seen in training.
>
> Fold uses native domain CATH labels, **not the CATH+TED skill**. No OOD_FoldHoldout column: topology is the class label and every topology must appear in train, making that holdout incompatible with the design.

#### OOD_ExtremeShort and Default

- Flag len(aa_seq)<60 on **all rows**, including 24 short NewEC priors; length and Default semantics are independent.
- **767** domains, 5.9%.
- Default=True everywhere: short chains follow the standard 7:1:2 split. Default=False is reserved for separately held-out pre-split samples, absent here, matching PPI-style conventions.
- ES rows have all dependent OOD annotations and idr_ratio. No test sample exceeds 1,000 aa after pool filtering; no OOD_ExtremeLong column.

#### OOD_LongTail: predicted-label rarity

Flag test topologies with **training frequency ≤10**, across all test rows. **753** rare topologies, **1,081** test domains, 8.4%.

#### OOD_IDR: intrinsic disorder

- metapredict v3; residue propensity >0.5 denotes disorder. Strip nonstandard characters before prediction.
- idr_ratio: disordered-residue fraction; >**0.3** sets OOD_IDR=True.
- Incremental sequence cache output/idr_ratio_full.csv; also inherit historical output/splits_ca_backup/cath_test.csv results to avoid recomputation.
- Canonical equivalent: toolkit idr_ood.py `--mode residue-ratio --threshold 0.3 --clean strip`.

**401** domains, 3.1%. All 12,911 idr_ratio values present, range 0.0–1.0.

### 4.3 Skill finalization: four EC tails and InD

Execute .skills/ scripts in order:

1. **Four EC tails:** .skills/ec-function-ood-annotation/scripts/add_unified_longtail_ec.py, project-wide batch driver with fold EC-source adaptation, writing splits directly. All-not-in: annotated samples True only when **all** ECs at the selected level have **training** frequency <5 /<10, **including zero**; unannotated False. Use full ECs from data/sifts/sifts_ec_annotations.json, the same source as the NewEC prior.

   | Column | True | Test fraction |
   |------|--------:|--------:|
   | OOD_LongTail_EC_L3_le5 | 368 | 2.9% |
   | OOD_LongTail_EC_L3_le10 | 415 | 3.2% |
   | OOD_LongTail_EC_L4_le5 | 1,459 | 11.3% |
   | OOD_LongTail_EC_L4_le10 | 2,199 | 17.0% |

   > All 336 NewEC priors necessarily meet LongTailEC: EC7 has zero training frequency, included by definition. This nesting is expected.

2. **InD:** toolkit finalize_ood_columns.py --mode ind. True if all OOD_*, seq_Redundancy_*, TM-score_* are False, treating NaN as False. Include four LongTailEC columns, recomputing after all annotations. Rename old Pure_ID in place. --mode nan has no applicable rows because all Default=True. **1,805** domains, **14.0%**.
3. **Validate:** toolkit check_splits.py; all checks passed, 2026-08-27.

### 4.4 Final files

Location: datasets/fold_classification/splits/.

| File | Rows | Columns |
|------|------|------|
| cath_train.csv | 43,473 | unique_id,aa_seq,struct_file,label |
| cath_val.csv | 6,132 | Same |
| cath_test.csv | 12,911 | Four core plus 27 annotations, 31 total |

Test annotation schema:

| Column | Type | Description |
|------|------|------|
| Default | bool | All True; no pre-split length holdout |
| InD | bool | 1,805 pure InD, 14.0% |
| seq_Redundancy_90…30 | bool | Seven sequence thresholds |
| TM-score_0.9…0.3 | bool | Seven structural thresholds |
| OOD_Orphan | bool | No significant train+val homolog |
| OOD_SuperfamilyHoldout | bool | Superfamily absent from train/val, post-split |
| OOD_NewEC | bool | Pre-split EC7 translocase prior |
| OOD_ExtremeShort | bool | <60 residues, including 24 NewEC priors |
| OOD_LongTail | bool | Topology training frequency ≤10 |
| idr_ratio | float | metapredict disordered-residue fraction |
| OOD_IDR | bool | idr_ratio >0.3 |
| OOD_LongTail_EC_L3_le5/le10, OOD_LongTail_EC_L4_le5/le10 | bool | EC L3/L4 training frequency <5 /<10, including zero |

### 4.5 OOD summary: all 12,911 test rows

| Dimension | Definition | Domains | Test fraction |
|----------|----------|----------:|--------:|
| **NewEC prior** | Pre-split EC7 holdout | **336** | **2.6%** |
| seq_Redundancy_90…30 | Max pident <90%…30% | 10,257 /9,118 /8,114 /7,055 /5,586 /3,876 /2,353 | 79.4%…18.2% |
| TM-score_0.9…0.3 | Max alntmscore <0.9…0.3 | 4,022 /1,755 /839 /495 /348 /273 /237 | 31.2%…1.8% |
| OOD_Orphan | No significant mmseqs2 hit | 1,978 | 15.3% |
| OOD_SuperfamilyHoldout | Post-split unseen superfamily | 761 | 5.9% |
| OOD_ExtremeShort | <60 aa, including 24 priors | 767 | 5.9% |
| OOD_LongTail | Topology training frequency ≤10 | 1,081 | 8.4% |
| OOD_LongTail_EC_L3_le5 /le10 | Training frequency <5 /<10, including zero | 368 /415 | 2.9% /3.2% |
| OOD_LongTail_EC_L4_le5 /le10 | Same, L4 | 1,459 /2,199 | 11.3% /17.0% |
| OOD_IDR | idr_ratio >0.3 | 401 | 3.1% |
| **InD** | No OOD flags | **1,805** | **14.0%** |

**Subset counts for reference:**

**ExtremeShort, 767 rows:** seq90…30=534/491/440/397/323/270/258; TM0.9…0.3=410/284/224/189/166/151/140; Orphan 258; superfamily holdout 151; LongTail 110; IDR 122; LongTailEC L3le5/L3le10/L4le5/L4le10=24/25/54/72; NewEC 24 short EC7 priors.

**NewEC prior, 336 rows:** seq90…30=334/330/322/313/292/249/163; TM0.9…0.3=220/115/56/33/24/21/21; Orphan 132; superfamily holdout 36; LongTail 43; IDR 18.

### 4.6 Quality checks: all pass

| Check | Result |
|--------|------|
| Pairwise train/val/test unique_id disjointness | ✅ |
| All 1,163 topologies in both train and test | ✅ |
| C.A.T string labels | ✅ |
| Zero train/val EC7, checked against SIFTS | ✅ |
| Every superfamily-holdout label absent from train/val | ✅ |
| Orphan ⊆seq30 and monotone sequence ladder | ✅ |
| ES Default=True; real dependent OOD values and idr_ratio retained | ✅ |
| No missing OOD values | ✅ |
| InD equivalent to no OOD flag in both directions | ✅ |
| Every struct_file exists in pdbs/ | ✅ |
| Toolkit check_splits.py | ✅ All pass, 2026-08-27 |

---

## Diagnostic tools

- Toolkit check_splits.py checks complete columns, monotone nesting, InD consistency, and other §4.6 rules.
- 04_build_topology_split.py /05_annotate_test_ood.py include assertions and summary prints (§4.1/4.2).

## Archived scripts: output/scripts_archive/

| Script | Reason archived |
|------|----------|
| fix_summary.py | One-off repair of early 02 CSV-appending defect; full post-extraction reconstruction now integrated |
| fix_multimodel_pdbs.py | Truncate 4,076 NMR files to model 1; applied to existing data and root cause fixed in 02 at first ENDMDL |
| extract_cath_domains_test.py | Development-only 50-domain extraction validation; mean identity 0.969 confirms boundary parsing |
| fill_es_ood.py, fill_prior_ood.py | One-off dependent-OOD completion during 2026-08-27 patching; integrated into current 05 |
| add_longtail_ec_ood.py | Initial single-column OOD_LongTail_EC_L4 replaced by canonical four-column EC skill |

## Notes

1. **Missing residues:** unresolved loops often yield n_residues <length_expected; this is normal.
2. **Expression tags:** some reference sequences include His-tags/vector residues absent from structures.
3. **Chain 0:** CATH chain='0' means no PDB chain field; SIFTS fallback tries ['0','A',' '].
4. **CIF conversion:** domains extracted from .cif.gz are written as atom-renumbered PDBs.
5. **S95 not globally unique:** use full (C,A,T,H,S35,S60,S95) cluster key.
6. **OOD overlap:** domains may satisfy multiple dimensions, standard in OOD evaluation.
7. **EC coverage:** only some chains have SIFTS ECs; unannotated LongTailEC is False, not NaN.
8. **Fold NewEC exception:** pre-split, single-column EC7 prior differs from the EC skill's all-not-in L3/L4. Do not recompute it with that skill; four LongTailEC columns do use the skill (§4.3).
9. **No CATH+TED skill:** native C.A.T.H domain labels require no TED supplementation; fold holdout is incompatible with topology-stratified training coverage (§4.2).
10. **Former manual steps:** 03 now removes five >1000-expected-length domains. It also reconstructs the v1 C.A→v2 C.A.T pool with verified row-equivalence; no residual manual step.
11. **Historical deviation audit:** new 03 pool and v1 backup are exactly row-equivalent; no deviations found.
12. **History:** v1, 2026-06, used 43 C.A classes with insufficient discrimination; backups in output/splits_ca_backup/. v2, 2026-08, uses C.A.T 7:1:2 stratification, narrows EC5+EC6 to EC7, annotates superfamily holdout post-split, and removes fold holdout. Since 2026-08-27 all ES/prior rows have dependent OOD annotations; Pure_ID renamed InD and recomputed over all flags. Backups: output/backup_*/.
