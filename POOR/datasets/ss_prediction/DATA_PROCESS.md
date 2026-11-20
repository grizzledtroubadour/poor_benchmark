# Secondary Structure Prediction (SSP) Dataset Processing

> Complete, reproducible workflow from experimental PDB structures to SSP data.
> Location: `datasets/ss_prediction/`.
> **Task:** secondary-structure prediction for high-quality PDB chains. Residue labels are DSSP states mapped as H 0 G1 I 2 E 3 B 4 T5 S6 P 7 C 8; -1 is the ignore index. unique_id uses PDB_chain, e.g. 10GW_A.

**Measured current splits:**

| File | Count | Fraction | Deposition range | Length range | Description |
|--------|------|------|-------------|---------|------|
| splits/ssp_train.csv | **29,622** | 71.7% | 1994-03–2019-06 | 31–998 | In-distribution training, including 834 short chains of 31–59 aa |
| splits/ssp_val.csv | **3,317** | 8.0% | 2019-06–2021-02 | 31–931 | Validation, including 118 short chains |
| splits/ssp_test.csv | **8,372** | 20.3% | 1999-07–2026-05 | 31–1,836 | 7,997 temporal samples +293 short +82 long; 34 columns including OOD, InD, struct_label |

> Current short-chain definition:31–59 aa (§1.3, §4.1). Total 41,311 =39,984 main-pool chains ≥60 +1,245 short +82 extreme-long.

## Design principles

1. **Direct PDB sources:** obtain RCSB metadata and chain information directly, without requiring SIFTS→UniProt mappings.
2. **Quality first:** prioritize X-ray crystallography with resolution≤2.5 Å and R-free≤0.25 for reliable DSSP labels.
3. **Traceability:** retain raw data, intermediates, statistics, checkpoints, and reproducible outputs.

---

## Stage 1: Obtain and parse source data

> Obtain entry metadata for the local PDB mirror, filter quality, and extract chain sequences.

### 1.1 Entry metadata: scripts/01_fetch_pdb_metadata_api.py

Source: RCSB GraphQL API, https://data.rcsb.org/graphql.

| GraphQL field | Property | Description |
|-------------|---------|------|
| exptl.method | Experimental method | X-RAY DIFFRACTION /SOLUTION NMR /ELECTRON MICROSCOPY, etc. |
| refine.ls_d_res_high | X-ray resolution | Å; core coordinate-quality metric |
| refine.ls_R_factor_R_work | R-work | Working-set R-factor |
| refine.ls_R_factor_R_free | R-free | Independent free-set validation |
| em_3d_reconstruction.resolution | EM resolution | Cryo-EM fallback when X-ray resolution is absent |
| rcsb_accession_info.deposit_date | Deposition date | Temporal splitting |

Query50 entries/request; checkpoint every 500 entries in output/api_checkpoint.json. Retry each failed batch up to 5 times with exponential backoff. Delay0.3s between requests; back off on 429.

Input: data/pdb/ listing,255,227 entries with .pdb.gz/.cif.gz.

| Output | Description | Fields |
|------|------|------|
| output/pdb_metadata_api.csv | 254,975 metadata rows | pdb_id, method, resolution, r_work, r_free, deposition_date |
| output/api_checkpoint.json | Resume checkpoint | processed |

Actual run:5,105 batches,~2.5 hours;254,975/255,227 successful (**99.9%**).

### 1.2 Entry quality: scripts/02_filter_entries.py

| Filter | Threshold | Rationale |
|------|------|------|
| Method | X-RAY DIFFRACTION | High coordinate precision and DSSP reliability |
| Resolution | ≤**2.5 Å** | Empirical hydrogen-bond assignment boundary |
| R-free | ≤**0.25** | High-quality independent validation |

Input: output/pdb_metadata_api.csv. Outputs: output/stage1_filtered_entries.csv (121,160) and output/stage1_filtering_stats.json.

| Step | Remaining | Removed | Description |
|------|---------|--------|------|
| All successful API entries | **254,975** | — | API results |
| X-ray | **204,781** |19.7%|80.3% of entries |
| Resolution≤2.5 Å | **160,940** |21.4%|Mean 1.80 Å |
| R-free≤0.25 | **121,160** |**24.7%**|Mean R-free**0.211** |

> Final mean resolution1.80 Å, mean R-free0.211. Early entries without R-free, e.g.1UBQ, are excluded.

### 1.3 Chain sequences: scripts/03_extract_chain_sequences.py

Extract chain sequences for deduplication from quality-filtered PDBs; do not write coordinates here.

- Read .pdb.gz/.cif.gz for 121,160 entries.
- Extract protein chains from ATOM/HETATM, deduplicating CA atoms.
- Map three-letter to one-letter residues, including MSE→MET and SEC→CYS.
- For NMR, use model 1 only.
- **MIN_RESIDUES=31:** include31–59-aa short chains, exclude≤30 at source.

Inputs: output/stage1_filtered_entries.csv and data/pdb/.

| Output | Rows | Description |
|------|------|------|
| output/stage2_chain_sequences.fa |255,471|All chain FASTA sequences |
| output/stage2_chain_metadata.csv |255,471|pdb_id, chain_id, num_residues, has_mse, file_type, seq |

> Existing disk files are v 1 MIN_RESIDUES=60 outputs with 249,018 rows. Current 03 produces255,471=249,018+6,453 short chains. Historical short extraction is output/short_chains.fa and output/short_chains_meta.csv:22,231 chains 1–59 aa, including 6,453 of 31–59 aa (see history).

Measured main-pool statistics, supplemented by short-chain counts (31–59 aa, median~32):

| Metric | Value |
|------|------|
| Entries processed |121,160|
| Entries with extracted chains |118,182 (97.5%)|
| Total chains |**255,471** =249,018≥60 +6,453 short|
| Single-chain entries |56,045|
| Multichain entries |~62,137|
| Maximum chains/entry |88, macromolecular complex|
| Mean main-pool length |271.3 residues|
| Length range |31–2,082|

---

## Stage 2: Global sequence deduplication

### 2.1 MMseqs2: scripts/04_run_mmseqs2_clustering.py

Cluster 255,471 chains to 41,337 representatives:40,067 main +1,270 short.

Use easy-linclust with `--min-seq-id 0.95 -c 0.8 --cov-mode 0 --threads 32`, **one stage over all chains≥31**. Supply FASTA directly: mmseqs18.8cc5c accepts it; prebuilt database inputs cause createdb errors. select_representatives chooses longest chain, then lexicographic chain_key.

> easy-linclust is input-order sensitive, permitting small boundary differences. Published data used a historical two-stage route: separately cluster≥60 and<60 pools, cross-deduplicate, then remove≤30 representatives. Under longest-first selection, short chains absorbed into≥60 clusters cannot become representatives under≥95% identity/0.8 bidirectional coverage, making this equivalent to cross-deduplication apart from small tie-break differences. Historical1–59-aa clustering:22,231 inputs→9,229 clusters→150 cross-deduplicated→9,079 survivors, including 1,270 of 31–59 aa.

Inputs: stage 2 FASTA and metadata.

| Output | Rows | Description |
|------|------|------|
| output/stage4_nonredundant_sequences.fa |41,337|Representative FASTA |
| output/stage4_cluster_representatives.csv |41,337|chain_key, file_type, seq;40,067 main+1,270 short |
| output/stage4_mmseqs_clusters.tsv |—|Clusters |
| output/stage4_clustering_stats.json |—|Statistics |

> Existing stage 4_* files are v 1 main-pool 40,067 outputs. Historical output/short_chain_representatives.csv has 9,079 rows including≤30. Preserve output/short_mmseqs_clusters.tsv and output/short_cross_dedup.m8 for provenance. Current 04 emits combined 41,337 rows.

Historical≥60 main-pool statistics:

| Metric | Value |
|------|------|
| Input chains |249,018|
| Clusters |**40,067**|
| Singletons |11,729|
| Multimember clusters |28,338|
| Reduction |**83.91%**|

---

## Stage 3: Coordinates, DSSP labels, and aggregation

> Extract single-chain PDBs and DSSP labels for 41,337 representatives (40,067 main+1,270 short), then aggregate.

### 3.1 Extract coordinates: scripts/05_extract_representative_chains.py

Extract representative chains from compressed PDB/mmCIF in parallel. Save only model 1 for NMR (`io.set_structure(structure[0])`).

DSSP v 4 requires correctly placed HEADER/EXPDTA/CRYST1 records; malformed ANISOU can fail parsing. BioPython PDBParser/MMCIFParser plus PDBIO filters problematic records; append minimal HEADER, EXPDTA, CRYST1, END records.

About 367 chains (0.9%) have multicharacter IDs such as AAA,1D, B 2. Truncate IDs and modify BioPython chain objects. Fifty-five from large viral entries (4Y4O,4YBB, etc.) suffered residue mismatches from truncated-ID collisions; independent Structure objects per chain fixed this.

Inputs: output/stage4_cluster_representatives.csv and data/pdb/. Recorded output: **41,327** PDBs=40,067 main+1,260 short. Of 1,270 short representatives,23 fail because multicharacter mmCIF IDs exceed PDB limits, analogous to historical main-pool failures; one additional short chain fails DSSP (§3.2).

### 3.2 DSSP: scripts/06_run_dssp_parallel.py

mkdssp v 4.6.1;24 workers,~1,669 files/worker.

| Original DSSP | Three-state | Structure |
|--------------|-------------|---------|
| H α-helix |H|Helix|
| G 3₁₀-helix |H|Helix|
| I π-helix |H|Helix|
| P Polyproline II |H|Helix|
| E β-sheet |E|Sheet|
| B β-bridge |E|Sheet|
| T turn |C|Coil|
| S bend |C|Coil|
| space/coil |C|Coil|

Input: pdbs/.

| Output | Description |
|------|------|
| output/dssp/ |40,067 main-pool raw DSSP files|
| output/dssp_short/ |1,260 short-chain DSSP files|
| output/dssp_labels.csv |40,067 rows: chain_key, seq, ss 8, ss 3|

Main pool: **40,067/40,067** success. Recorded short success: **1,246/1,260**; one parse_fail removed, with remaining1,246 including subsequently excluded7UYL_K (§4.1). Main residues:**10,407,197**, lengths 55–2,082; short chains 31–59.

Three-state distribution: H 3,890,548 (**37.4%**), E 2,510,675 (**24.1%**), C 4,005,974 (**38.5%**).

| State | Structure | Residues | Fraction |
|------|---------|--------|------|
|H|α-helix|3,199,406|30.7%|
|G|3₁₀-helix|419,045|4.0%|
|I|π-helix|56,498|0.5%|
|P|Polyproline II, DSSPv4|215,599|2.1%|
|E|β-sheet|2,386,784|22.9%|
|B|β-bridge|123,891|1.2%|
|T|Turn|1,181,359|11.4%|
|S|Bend|863,036|8.3%|
|C|Coil|1,961,579|18.8%|

### 3.3 Aggregate: scripts/07_aggregate_ssp_full_data.py

Inputs: stage 4 representatives for row membership/aa_seq/file_type; output/dssp_labels.csv ss 8 for labels. Output: output/ssp_full_data.csv, columns unique_id, aa_seq, file_type, labels. Labels are space-separated strings such as [8 8 3 ...], exactly aa_seq length.

Map H 0 G1 I 2 E 3 B 4 T5 S6 P 7 C 8. Full aa_seq can exceed DSSP because DSSP drops incomplete-backbone residues. Map exact sequences directly; otherwise (5,737 main-pool rows), use difflib.SequenceMatcher with default autojunk=True: transfer equal blocks, map replace blocks positionally, set uncovered residues to -1. Row-level diff verifies agreement with existing data.

> Existing ssp_full_data.csv is v 1,40,067 main-pool rows including later-removed4O9X_A. Short labels use the same logic, with 235 difflib alignments. Current 07 emits **41,313**=40,067+1,246 short, including later-removed7UYL_K.

---

## Stage 4: Splitting and OOD annotation

### 4.1 Temporal splits: scripts/08_create_ssp_datasets.py

Inputs: ssp_full_data.csv, stage 4 representatives for lengths, pdb_metadata_api.csv for dates. Outputs: splits/ssp_{train,val,test}.csv. Train/val: unique_id, aa_seq, struct_file, label; test adds Default, OOD_ExtremeLong.

- Remove chains>2,000: only 4O9X_A (2,082); assign to no split.
- QC_EXCLUDE={7UYL_K}:71/120 structural residues are UNK (59%). AA_MAP skips them, yielding aa_seq/label49 while structure/struct_label120, breaking alignment (QC below).
- **82 extreme-long chains**,1,000<length≤2,000: append to test by unique_id, without temporal splitting; Default=False, ExtremeLong=True.
- Remaining **41,229**,31≤length≤1,000: stable deposition_date sort; quantile split **72:8:20**, n_train=int(n*.72), n_val=round(n*.08), remainder test. Cutoffs are dynamic. Published cutoffs 2019-06-06/2021-02-09 originate from historical≥60 quantiles followed by short-chain insertion. Replay cutoffs may differ with pool composition.

Final train**29,622**/val**3,317**/test**8,372**=7,997 main temporal+293 short+82 long; total**41,311**.

Recorded properties: date sorting places an entry's chains in one split; unique_id intersections are all zero. Reproduction may differ by linclust boundary ties, dynamic cutoffs recalculated on the full pool, and row order (historically short chains appended; replay inserts chronologically). These are documented permitted differences, not quality changes.

### 4.2 OOD overview

Reference: Train+Val 32,939, including 952 short chains. Annotate8,372 test rows in ssp_test.csv.

> **NaN convention:**82 Default=False long chains do not enter homology or other OOD calculations. Preserve ExtremeLong=True, set other OOD fields NaN and InD=False. The293 short chains participate normally (Default=True, ExtremeShort=True). Canonical finalization: .skills/ood-annotation-toolkit/scripts/finalize_ood_columns.py --mode nan. Counts below use **8,290 Default samples**.

| Column | Definition |
|------|------|
|Default|8,290 temporal test samples True, including 293 short;82 long False|
|ExtremeShort|293 chains 31–59 aa, participating in OOD|
|ExtremeLong|82 chains>1,000|
|InD|All OOD/seq/TM flags False, NaN treated False;**1,123/8,290 (13.5%)**, after 2026-09-02 short-chain recomputation|

#### Sequence OOD: seq_Redundancy_90..30 and Orphan

Tool:.skills/homology-ood-annotation/scripts/seq_homology_ood.py. MMseqs2 easy-search, test queries versus train+val; maximum pident<T flags True, no hit counts 0. Orphan means no hits, ensuring Orphan⊆seq30.

**Current recomputation:2026-09-02**, after short-chain inclusion; easy-search -s 7.5, e≤1e-3, no coverage filtering. Replayable cache: output/ood_search/mmseqs_easysearch_test_vs_trainval_v3.m8. Superseded caches:2026-08-30 analysis/output/seq_ood_unify_assessment/ssp/mmseqs_easysearch.m8 and output/ood_search/mmseqs_easysearch_test_vs_trainval_v2.m8. Historical ood_search intermediates miss hits and overestimate OOD; provenance only.

| Threshold | Count | Fraction |
|------|---------|---------|
|<90%|6,894|83.2%|
|<80%|6,139|74.1%|
|<70%|5,548|66.9%|
|<60%|5,035|60.7%|
|<50%|4,327|52.2%|
|<40%|3,283|39.6%|
|<30%|2,016|24.3%|
|Orphan|**1,566**|**18.9%**|

Measured2026-09-02; denominator8,290.

#### Structural OOD: TM-score_0.9..0.3

Tool:.skills/homology-ood-annotation/scripts/struct_homology_ood.py. Foldseek easy-search, test structures versus train+val, including 31–59-aa targets; maximum alntmscore<T, no hit0. Current 2026-09-02 cache: output/ood_search/foldseek_easysearch_test_vs_trainval_v3.m8, replayable via --m8-cache. Short targets can raise some maximum scores above 1 (v 2:305 queries, max 1.072), documented as expected short-target normalization behavior.

| Threshold | Count | Fraction |
|------|---------|---------|
|<0.9|2,863|34.5%|
|<0.8|1,356|16.4%|
|<0.7|696|8.4%|
|<0.6|347|4.2%|
|<0.5|188|2.3%|
|<0.4|127|1.5%|
|<0.3|102|1.2%|

> Correction2026-08-28: old 0.9–0.6 columns came from foldseek search -e inf --max-seqs1000 with the entire pdbs/ as target, including test self-hits. The m8 had778 invalid scores>1, max 1.5, and underestimated counts (TM 0.6:19→252; TM 0.9:2,194→2,713). Existing scripts could not reproduce their provenance. Recomputed with canonical easy-search and added0.5/0.4/0.3; the old claim of≤5 samples below 0.5 is invalid. Comparison: output/foldseek_easysearch_compare/; backup: output/backup_before_foldseek_easysearch_unify/. Old output/ood_search/foldseek_result*.m8 is provenance only.

#### IDR: scripts/09_run_idr_prediction.py

metapredict v 3, aligned with toolkit idr_ood.py --mode residue-ratio --threshold 0.3 --clean strip.

1. Strip nonstandard X/B/Z/J/U before prediction (none present here).
2. Disorder residue: score>0.5.
3. IDR_ratio=number of such residues/full length.
4. Ratio>0.3→True.

Output: output/idr_predictions.csv (40,067 main rows) and test flags. Short-chain canonical cache: output/idr_cache_short.csv; flags: output/short_idr_flags.csv.

| Scope | Count | Fraction |
|------|-------------|------|
|Main 40,067|474|**1.18%**|
|Default test 8,290|159|**1.9%**|

Mean main IDR ratio**2.05%**, well below 0.3. Highly disordered proteins are scarce in experimental structures because crystallography struggles with disorder. Short chains have higher ratios, but only 40/293 test short chains are positive; most short chains entered train/val by date.

#### LongTail: scripts/10_run_longtail_prediction.py

Rare DSSP states: I π-helix 0.54%, B β-bridge 1.19%, P polyproline-II 2.07%. LT_ratio=(I+B+P residues)/valid residues; flag **top 5%**. Reference: all 41,311 including short chains. Current 95th-percentile threshold**0.082914**, output/longtail_threshold_v3.json. History:0.081181 at 40,067→0.093333 at 49,007→0.082914.

Output: historical output/longtail_predictions.csv (40,067 main) and test flags.

| Scope | Count | Fraction | Threshold |
|------|------------------|------|--------------|
|All41,311|2,066|**5.0%**|≥0.082914|
|Default test 8,290|436|**5.3%**|—|

Minimum flagged LT_ratio8.29%,~2.3× global mean 3.65%. Almost all top proteins are P-rich; I/B rarely dominate one protein.

#### New EC: OOD_NewEC_L3/L 4

Tool:.skills/ec-function-ood-annotation/scripts/add_unified_newec.py. Source: SIFTS data/sifts/sifts_chain_ec.tsv.gz. Reference: train+val ECs at the level. All-not-in: annotated and every EC absent→True; unannotated→False. L 3: first3 components when≥3; L 4: all 4, allowing '-'. Default=False→NaN; backup: output/backup_before_newec_unify/.

| Column | True among 8,290 Default,2,722 annotated |
|------|------:|
|NewEC_L4|**245**|
|NewEC_L3|**1**|

EC coverage**32.8%**. L 3 too small; prioritizeL4. Legacy run_ec_newcategory.py used any semantics/date-derived membership; archived. Current SSP-only recomputation2026-09-02 follows expanded short-chain references.

#### Long-tail EC

Tool:.skills/ec-function-ood-annotation/scripts/add_unified_longtail_ec.py. Training-only per-sample EC frequencies; all ECs below 5/10, including 0→True; unannotated→False. Four columns OOD_LongTail_EC_{L 3, L 4}_le{5,10}; Default=False→NaN. Backup: output/backup_before_longtail_unify/.

| Column | True among 8,290 Default,2,722 annotated |
|------|------:|
|L 3_le5|37|
|L 3_le10|77|
|L 4_le5|738|
|L 4_le10|1,118|

L 4 fraction 13.5% reflects temporal splitting: ECs in new entries are naturally rare in pre2019-06 training. Current SSP-only recomputation2026-09-02 includes short chains.

#### CATH+TED holdouts: scripts/11_ted_augment_cath_ood.py

Canonical annotation/holdout scripts in .skills/cath-ted-domain-annotation:

0. Build lowercase pdb/chain plus splits from unique_id→output/ss_ted_aug_chain_list.csv,**41,311** chains excluding4O9X_A/7UYL_K.
1. cath_ted_annotate.py --offline; most TED accessions hit data/ted/ted_api_cache.json;309 uncached receive none. Outputs: output/ss_chain_cath_ted_labels.csv, output/ss_chain_domains.csv (94,998 domain rows).
2. cath_ted_holdout.py --split-col splits; new train 29,622 gives1,189 topologies/4,214 superfamilies. Update test holdouts; Default=False→NaN.

Chain mode unique_id=PDB_chain. Labels union experimental CATHv4.4 name matches from cath-domain-list.txt, with intervals mapped through SIFTS, and putative cached TED domains overlapping chain UniProt ranges by≥50% of domain length. H providesC.A.T.H; T providesC.A.T only. Training CATH+TED union reference, symmetric across splits; any unseen topology/superfamily→True, unannotated→False.

Across41,311 filtered/deduplicated PDB chains: CATH-only coverage 58.5%; TED raises it to**82.8%**. Sources: both 19,994 /ted-only 10,366 /cath-only 3,852 /none 7,099.

| Default test 8,290 | True |
|------|--------:|
|FoldHoldout|20|
|SuperfamilyHoldout|96|
|Unannotated→False|3,511|

#### Finalization

1. Default=False long rows: all OOD/sequence flags NaN except ExtremeLong; InD=False. Since 2026-08-28 scripts 09/10/11 and skills write NaN directly; finalize_ood_columns.py --mode nan only validates (no-op currently). That day574 historical True sequence cells across82 rows became NaN; backup: output/backup_nondefault_seqtm_nan/.
2. --mode ind: recompute using all completed OOD columns;**1,123/8,290,13.5%**, after 2026-09-02 short-chain update.
3. Validate with .skills/ood-annotation-toolkit/scripts/check_splits.py.

Current test header,34 columns:

```
unique_id, aa_seq, struct_file, label, Default,
seq_Redundancy_90, seq_Redundancy_80, ..., seq_Redundancy_30,
TM-score_0.9, TM-score_0.8, TM-score_0.7, TM-score_0.6, TM-score_0.5, TM-score_0.4, TM-score_0.3,
OOD_Orphan, OOD_ExtremeShort, OOD_ExtremeLong, OOD_IDR, OOD_LongTail,
OOD_NewEC_L4, OOD_NewEC_L3,
OOD_FoldHoldout, OOD_SuperfamilyHoldout,
OOD_LongTail_EC_L3_le5, OOD_LongTail_EC_L3_le10,
OOD_LongTail_EC_L4_le5, OOD_LongTail_EC_L4_le10,
InD, struct_label
```

### 4.3 Structure-residue labels

All splits end with struct_label: one DSSP label per structural residue, enumerating first ATOM occurrences by(chain, resseq, icode) in file order, including residues without CA. Structure models parse PDB and read these labels directly, no alignment needed.

Canonical writer:.skills/struct-label-alignment/scripts/add_struct_label.py. Authoritative labels come directly from DSSP in the archived struct_native view, output/splits_struct_native_deprecated/. Equal-length matching enumerations (~98.7%) copy values directly. For rows lacking CA (train 384 /val 12 /test 44), align full-ATOM sequence to native CA sequence; missing-CA residues→**-1**, ignored in training/evaluation. Backup: output/backup_before_struct_label_native/.

Validation: all other cells unchanged; equal-length rows match native labels100%; sampled structural-label lengths match residue counts. -1 counts: train 871 /val 62 /test 176. Use main splits+struct_label; struct_native view is redundant and archived (note11).

> Historical patch_struct_label_from_native.py applied the migration and is archived. Current recomputation goes through canonical add_struct_label.py.

### 4.4 Quality checks

| Check | Result |
|--------|------|
|Pairwise split unique_id intersections|All0, verified|
|label length=aa_seq length|Generated per chain from ssp_full_data.csv; asserted rowwise in07|
|Test OOD completeness|8,290 Default including 293 short have True/False;82 long have NaN except ExtremeLong|
|InD recomputation|**1,123**,2026-09-02, matching files; history1,135 (2026-08-30 sequence unification)→1,179 (short chains by cutoff)→1,123 (≤30 removal)|
|Nested flags|Orphan⊆seq30⊆…⊆seq90; monotonic counts verified|

Applied QC events:

- **7UYL_K removal,2026-09-02:**71/120 UNK residues led to 49-residue aa_seq/label but120-residue structure/struct_label and 71 incorrect positions versus direct DSSP. Remove test row and pdbs/7UYL_K.pdb; backup: output/backup_20260902/ssp_test.before_7uylk_removal.csv. Current 08 pins QC_EXCLUDE. Seven other short rows differ by 1–3 missing-CA residues; their struct_label matches DSSP and is unaffected.
- **struct_label restoration,2026-09-01/02:**short-chain recovery with add_struct_label.py --overwrite changed 148 rows with aa_seq≠structure sequence. Chained alignment via aa_seq inherited labels incorrectly in low-complexity regions. Direct DSSP checks on 125 rows found old 99.996% agreement (2 differences/46,206 residues), new only 71%. **Restore 148 rows from output/backup_20260901/ssp_test.csv to native-derived labels**; pre-restoration backup: output/backup_20260901/ssp_test.before_struct_label_restore.csv. New short rows use CA-derived aa_seq and direct mapping;500 sampled agree 99.97%. **Do not recompute struct_label with --overwrite when aa_seq differs from the structural sequence.**

---

## Diagnostics and quality control

- `.skills/ood-annotation-toolkit/scripts/check_splits.py`: complete split validation.
- `output/verify/`: trial outputs from `08_create_ssp_datasets.py --out-dir output/verify`, historically verifying v 1 72:8:20 quantile cutoffs and rowwise agreement with splits. Current 08 restores dynamic 72:8:20 quantiles; this directory records the historical state only.
- Archived `ss_struct_alignment_probe.py` diagnoses alignment between structure residues and non-`-1` label positions; see the archive list.

## Archived scripts (`output/scripts_archive/`)

| Script | Reason |
|------|----------|
| `extract_subset.py` | Redundant subset extraction during resume, replaced by `05_extract_representative_chains.py`. Contains an **unfixed multimodel bug**: `io.set_structure(structure)` writes all models. Do not reuse; use 05 for correct handling |
| `run_ec_newcategory.py` | Initial NewEC with any semantics and date-derived splits; superseded by `.skills/ec-function-ood-annotation/scripts/add_unified_newec.py` |
| `build_ss_struct_native.py` | Builds deprecated struct_native views, archived in `output/splits_struct_native_deprecated/`; replaced by main-split `struct_label` |
| `patch_struct_label_from_native.py` | One-off migration patch applied to splits; recompute with `.skills/struct-label-alignment/scripts/add_struct_label.py` |
| `ss_struct_alignment_probe.py` | Diagnostic-only structure–label alignment probe |
| `fix_multimodel_pdbs.py` | Applied one-off repair retaining model 1 in 26 NMR PDBs; root cause fixed in `05_extract_representative_chains.py` using `structure[0]` |
| `add_longtail_ec_ood.py` | Initial single `OOD_LongTail_EC_L4` patch, superseded by four canonical EC columns |
| `extract_pdb_metadata.py`, `extract_representative_structures.py`, `run_dssp_batch.py` | Early pipeline versions superseded by 01/05/06 |
| `12_extract_short_chains.py`, `13_cluster_short_chains.py`, `14_extract_short_structures.py`, `15_build_short_test_rows.py` | Historical short-chain recovery (2026-09-01): <60 extraction, clustering/cross-deduplication, structures, and DSSP/test append. Incorporated into 03 (MIN_RESIDUES=31), 04 (single-stage all-pool clustering), and generic representative-list 05/06; see History |
| `16_redistribute_short_chains.py`, `17_apply_ood_updates.py` | Normal splitting of short chains at original cutoffs and OOD recomputation (2026-09-02); incorporated into 08 (72:8:20 quantiles +QC_EXCLUDE) |
| `18_drop_short_chains_le30.py`, `19_apply_ood_updates_v3.py` | Remove ≤30-aa chains and recompute OOD (2026-09-02); fixed in 03 as MIN_RESIDUES=31 |

## Notes

1. **API rate limits:** RCSB GraphQL limits requests; scripts include backoff. Increase `REQUEST_DELAY` if HTTP429 occurs frequently.
2. **Resume:** `01_fetch_pdb_metadata_api.py` automatically resumes from checkpoints after termination.
3. **Missing R-free:** some older entries, such as 1UBQ, lack R-free and are excluded during filtering.
4. **Formats:** prefer `.pdb.gz`, falling back to `.cif.gz`. Truncate multicharacter mmCIF chain IDs to one character during extraction.
5. **DSSP v 4:** PDBs require HEADER, EXPDTA, CRYST1, and END records, with correctly positioned ANISOU records. BioPython PDBIO extraction avoids these issues.
6. **Long chains:83→82.** Stages1/2 contain 83 representatives >1000 residues, but `4O9X_A` (2,082 residues) is removed during splitting by the global >2000 rule. Final long test count is **82**, total test **8,372**. Old 83/8,080 counts preceded this removal; historical short-chain paths produced test totals 17,020 /10,456 /10,455 (History).
7. **Removed `is_time_cutoff`:** an early helper column absent from current headers; equivalent information is in Default.
8. **NewEC history:** initial `run_ec_newcategory.py`, using any semantics/date-derived splits, yielded L4=287/L3=5. Canonical all-not-in with actual split membership yields **L4=246/L3=1**; backup `output/backup_before_newec_unify/`.
9. **LongTailEC history:** initial single `OOD_LongTail_EC_L4` column replaced by four `OOD_LongTail_EC_{L3,L4}_le{5,10}` columns; original script archived, backup `output/backup_before_longtail_unify/`.
10. **TM direction fix:** early threshold direction was wrong, corrected with backup `output/ssp_test.before_tm_inversion.csv`. On2026-08-28, four levels were found to use legacy search -e inf with invalid values>1. Recomputed all with canonical easy-search and added0.5/0.4/0.3 levels; see Stage 4 structural OOD and `output/backup_before_foldseek_easysearch_unify/`.
11. **Deprecated struct_native views (2026-08-28):** retain `output/splits_struct_native_deprecated/` because analysis scripts still reference it. Its flags stop at OOD_NewEC_L4/L3 and InD, lacking later fold/superfamily/LongTailEC columns. Its InD=871 was not updated, unlike main splits799. Historical only; **do not use for OOD evaluation**. Replaced by main `struct_label` (§4.3).
12. **NMR repair:** scanning on 2026-08-28 found26 multimodel PDBs:23 with MODEL/ENDMDL,3 concatenated without markers. Root fix in05 retains model 1; one-off applied repairs backed up to `output/backup_before_multimodel_fix/`, repair script archived. Recompute struct_label with --overwrite afterward; residue-key deduplication leaves annotation values unaffected.
13. **Historical file_type deviation found during consolidation:**787 rows in `output/ssp_full_data.csv` incorrectly say pdb. Sources exist only as `.cif.gz` in data/pdb/, correctly recorded as cif in stage 2_chain_metadata.csv. Downstream08 ignores this column, so splits are unaffected. Rebuilding with 07_aggregate_ssp_full_data.py corrects it; these787 cells are the sole differences, with unique_id/aa_seq/labels matching all 40,067 rows.
14. **Previously unscripted stages:** no additional manual step precedes MMseqs2 stage 04. Former unscripted ssp_full_data.csv and ss_ted_aug_chain_list.csv now have07 and integration into11, respectively; no remaining pipeline gaps.

## Key file index

| Path | Stage | Description |
|---------|------|------|
| `scripts/01_fetch_pdb_metadata_api.py` | 1 | GraphQL metadata retrieval |
| `scripts/02_filter_entries.py` | 1 | Quality filtering |
| `scripts/03_extract_chain_sequences.py` | 1 | Chain sequence extraction |
| `scripts/04_run_mmseqs2_clustering.py` | 2 | MMseqs2 deduplication |
| `scripts/05_extract_representative_chains.py` | 3 | Single-chain PDB extraction |
| `scripts/06_run_dssp_parallel.py` | 3 | Batch DSSP annotation |
| `scripts/07_aggregate_ssp_full_data.py` | 3 | Aggregate ssp_full_data.csv |
| `scripts/08_create_ssp_datasets.py` | 4 | Dataset splitting |
| `scripts/09_run_idr_prediction.py` | 4 | metapredict IDR and OOD flags |
| `scripts/10_run_longtail_prediction.py` | 4 | Eight-state long-tail computation and flags |
| `scripts/11_ted_augment_cath_ood.py` | 4 | Chain-list generation, CATH+TED annotation, fold/superfamily holdouts |
| `.skills/homology-ood-annotation/scripts/seq_homology_ood.py` | 4 | Canonical seq_Redundancy_*/Orphan |
| `.skills/homology-ood-annotation/scripts/struct_homology_ood.py` | 4 | Canonical TM-score_* |
| `.skills/ec-function-ood-annotation/scripts/add_unified_newec.py` | 4 | Canonical NewEC L3/L4 |
| `.skills/ec-function-ood-annotation/scripts/add_unified_longtail_ec.py` | 4 | Four canonical LongTailEC columns |
| `.skills/cath-ted-domain-annotation/scripts/cath_ted_annotate.py` | 4 | Chain CATH+TED annotation |
| `.skills/cath-ted-domain-annotation/scripts/cath_ted_holdout.py` | 4 | Recompute fold/superfamily holdouts |
| `.skills/ood-annotation-toolkit/scripts/idr_ood.py` | 4 | Canonical residue-ratio IDR, matching 09 |
| `.skills/ood-annotation-toolkit/scripts/finalize_ood_columns.py` | 4 | Default=False NaN finalization and InD |
| `.skills/struct-label-alignment/scripts/add_struct_label.py` | 4 | Recompute/write struct_label |
| `output/pdb_metadata_api.csv` | 1 | Full metadata,254,975 rows |
| `output/stage1_filtered_entries.csv` | 1 | Filtered entries,121,160 |
| `output/stage2_chain_sequences.fa` | 1 | Existing v 1:249,018 chains≥60; rerun current03:255,471, including6,453 short 31–59-aa chains |
| `output/stage4_cluster_representatives.csv` | 2 | Existing v 1 main pool 40,067; current04 rerun41,337, including1,270 short chains |
| `pdbs/` | 3 | 41,327 single-chain PDBs:40,067 main +1,260 short 31–59 |
| `output/dssp_labels.csv` | 3 | Main-pool 40,067 labels; short DSSP in output/dssp_short/,1,260 files |
| `output/ssp_full_data.csv` | 3 | Existing v 1:40,067, including removed4O9X_A; rerun41,313 |
| `output/idr_predictions.csv` | 4 | Main-pool 40,067; short cache output/idr_cache_short.csv |
| `output/longtail_predictions.csv` | 4 | Historical main-pool 40,067; current thresholds output/longtail_threshold_v3.json |
| `output/ss_chain_cath_ted_labels.csv` | 4 | 41,311 chain annotations |
| `output/ss_chain_domains.csv` | 4 | 94,998 domain rows |
| `output/ss_ted_aug_chain_list.csv` | 4 | 41,311 chain rows, generated by 11 |
| `splits/ssp_train.csv` | 4 | 29,622=28,788 main +834 short 31–59-aa |
| `splits/ssp_val.csv` | 4 | 3,317=3,199 main +118 short |
| `splits/ssp_test.csv` | 4 | 8,372=7,997 main +293 short +82 long; complete OOD/InD/struct_label |
| `output/short_chains.fa` / `output/short_chains_meta.csv` | History | Recovery sequences/metadata:22,231 chains 1–59 aa;6,453 at 31–59 entered the pool |
| `output/short_chain_representatives.csv` | History | 9,079 survivors, including≤30;1,270 at 31–59 enter representatives |
| `output/short_clustering_stats.json` | History | 1–59 pool:9,229 clusters→150 removed→9,079 |
| `output/short_dssp_failures.csv` | History | 138 short-chain structure/DSSP failures |
| `output/backup_20260901/ssp_test.csv` | History | Before ExtremeShort append,8,079 rows |
| `output/backup_20260902/` | History | Three split backups before short redistribution,test 17,020 |
| `output/backup_20260902_min30/` | History | Before≤30 removal,test 10,455 |
| `output/ood_search/mmseqs_easysearch_test_vs_trainval_v2.m8` | History | Sequence search after short splitting; current cache v 3 |
| `output/ood_search/foldseek_easysearch_test_vs_trainval_v2.m8` | History | Structure search after short splitting; current cache v 3 |
| `output/ood_search/mmseqs_easysearch_test_vs_trainval_v3.m8` | 4 | Current sequence cache, short reference included, after min 30 |
| `output/ood_search/foldseek_easysearch_test_vs_trainval_v3.m8` | 4 | Current structure cache, short reference included, after min 30 |
| `output/idr_cache_short.csv` / `output/short_idr_flags.csv` | History | 8,941 short IDR predictions/flags, including subsequently removed≤30 |
| `output/longtail_threshold_v2.json` / `output/longtail_threshold_v3.json` | History/4 | Thresholds 0.081181→0.093333→0.082914, current v 3 |
| `output/splits_struct_native_deprecated/` | History | Archived native-structure splits, replaced by main struct_label |

## History

Actual evolution of short-chain handling is now linearized into Stages1/2/4. All scripts below are archived in `output/scripts_archive/`.

1. **2026-09-01 recovery (12–15):**03 originally used MIN_RESIDUES=60. Re-extract22,231 chains 1–59 aa from source PDBs →easy-linclust9,229 clusters/representatives →cross-deduplicate against existing splits, removing150 →9,079 survivors →123 structure failures +15 DSSP failures →8,941 appended entirely to test with Default=False/ExtremeShort=True, giving intermediate test 17,020.
2. **2026-09-02 normal splitting at original cutoffs (16/17):** reconstruct08 cutoffs: train≤2019-06-06,val≤2021-02-09. Allocate short chains train 5,672/val892/test2,377 and recompute all OOD against expanded references. Remove7UYL_K by QC that day (§4.4).
3. **2026-09-02 ≤30 removal (18/19):** restrict short chains to31–59; remove 7,695≤30-aa chains (train 4,838/val774/test2,083) and orphan structures. Recompute OOD, yielding current41,311=29,622/3,317/8,372.
4. **2026-09-02 documentation linearization:** integrate short handling into03 (MIN_RESIDUES=31),04 (single-stage full-pool clustering),08 (72:8:20 quantiles+QC_EXCLUDE); archive12–19.

**Equivalence:** from-scratch reproduction differs from current data in three ways: (a) easy-linclust boundary ties, since input order can slightly change cluster membership; (b) **cutoffs move with pool composition**: requantiling the pool including short chains yields dates different from released2019-06-06/2021-02-09, changing a few near-boundary assignments (§4.1); (c) row order: historical short chains were appended to train/val, while reruns insert them by date. These are accepted deviations and do not affect data quality. Some disk intermediates (stage 2/stage4/ssp_full_data) remain v 1; section-specific counts apply and reruns update them.
