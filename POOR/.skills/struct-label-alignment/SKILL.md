---
name: struct-label-alignment
description: 为残基级任务（ss_prediction、ligand_binding_site、ppis_prediction 等）的 split csv 生成 struct_label 列——对 pdb 结构文件中的每个残基（ATOM 记录文件顺序）给出对应的残基级标注，使结构模型只需 pdb + struct_label 即可训练/评估。当需要将序列上的逐残基标注映射到结构残基时使用。
---

# struct-label-alignment

将 csv 中按序列位置组织的残基级标注（`label`，长度 == `aa_seq`）映射为按 **pdb 结构残基顺序** 组织的 `struct_label`。

## 语义

- pdb 残基枚举：ATOM 记录中按文件顺序、以 `(chain, resseq, icode)` 首次出现为准（**含无 CA 原子的残基**，与 Biopython 等解析器的读取顺序一致）。
- 映射方式：pdb 残基序列（三字母→一字母，修饰残基映射见脚本 `RES3`，未知为 `X`）与 `aa_seq` 做全局比对（Biopython PairwiseAligner，match=2/mismatch=-1/open=-2/extend=-0.5）；每个结构残基继承比对到的序列位置的标注。
- 等长时走位置直映射快速路径（替换残基仍算对齐）。
- **-1 = ignore index**：结构残基无法映射回序列（如仅存在于结构中的表达标签）时记 -1，训练/评估时应跳过。
- 输出格式与 `label` 一致：`[ 8 8 3 ... ]`（方括号+空格分隔）。
- 断言：`len(struct_label) == pdb ATOM 残基数`。

## 用法

```bash
python .skills/struct-label-alignment/add_struct_label.py \
    --csv datasets/<task>/splits/<split>.csv \
    --pdb-dir datasets/<task>/pdbs \
    [--seq-col aa_seq] [--label-col label] [--struct-col struct_file] \
    [--out-col struct_label] [--workers 32] [--dry-run]
```

- 先 `--dry-run` 看统计（rows_with_-1 / total_-1 / rows_with_nonX_mismatch），确认无异常再正式写入。
- 写入是**就地**修改 csv（列追加到末列）；写入前务必备份到 `<task>/output/backup_before_struct_label/`。
- csv 读写约定：`pd.read_csv(dtype=str, keep_default_na=False)` + `to_csv(index=False, lineterminator="\n")`（脚本已内置）。

## 已应用

- ss_prediction / ligand_binding_site / ppis_prediction 的 train/val/test（9 个 csv，2026-08）。
- ss_prediction 的 `struct_label` 后于 2026-08-28 改为直接移植 struct_native 视图的 DSSP 直达标注（脚本已归档至 `ss_prediction/output/scripts_archive/patch_struct_label_from_native.py`，native 视图已归档至 `ss_prediction/output/splits_struct_native_deprecated/`）；本 skill 脚本仍用于 lbs/ppis 及后续任务的比对式生成。

## 注意事项

- **`aa_seq ≠ 结构序列` 的行不要用本脚本 `--overwrite` 重算**：本脚本的继承式比对（结构序列→`aa_seq`→主 `label`）在低复杂度/重复区会发生错帧，导致残基继承错位位置的标注（ss_prediction 2026-09-01 实证：148 行 12,199 个残基标签翻转，对 DSSP 直接标注一致率仅 71%，恢复后见 `ss_prediction/DATA_PROCESS.md` §5.5）。仅当 `aa_seq` 即结构提取序列（如 SSP 2026-09 恢复的短链行）时直通路径才是安全的；已有更权威标注来源（如 DSSP 直达标注移植）的列应冻结，不要覆盖。
