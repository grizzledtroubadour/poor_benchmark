#!/usr/bin/env python3
"""统一 LongTailEC OOD 口径并落盘到 10 个任务的 test csv。

统一定义（2026-08-27 确定）：
- 参考集：train only（逐样本计数：携带某 EC 类别的 train 样本数）
- 层级：EC L3 与 L4
- 阈值：train 频次 < 5 / < 10（**含 0**，即 train 中未出现的类别也算）
- 语义：样本所有 EC（对应层级）都满足频次条件 → True（all-not-in）；无 EC 注释 → False
- 新列（4 列，追加到末尾）：
  OOD_LongTail_EC_L3_le5 / OOD_LongTail_EC_L3_le10 /
  OOD_LongTail_EC_L4_le5 / OOD_LongTail_EC_L4_le10
- Default=False（极端长度预挑出）行一律置 NaN（"未计算"语义，全项目统一
  惯例，流水线位置直接写入；finalize --mode nan 仅作兜底校验），仅 Default/非ES 子集参与计算
- 旧 LongTailEC 列被新列替换（kcat 的 OOD_LongTail_EC_L3、optimal_ph 的
  OOD_LongTail_EC_L3_le5/le10、ppi/ligand_affinity 及本次新增的 OOD_LongTail_EC_L4）

用法：python3 .skills/ec-function-ood-annotation/scripts/add_unified_longtail_ec.py
"""
import csv
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

# 仓库根：.skills/ec-function-ood-annotation/scripts/ 上溯 3 级
POOR = Path(__file__).resolve().parents[3]
SIFTS_EC_TSV = POOR / "data" / "sifts" / "sifts_chain_ec.tsv.gz"
SIFTS_EC_JSON = POOR / "data" / "sifts" / "sifts_ec_annotations.json"

THRS = (5, 10)
NEW_COLS = ["OOD_LongTail_EC_L3_le5", "OOD_LongTail_EC_L3_le10",
            "OOD_LongTail_EC_L4_le5", "OOD_LongTail_EC_L4_le10"]


def levels(ec):
    parts = ec.split(".")
    return (".".join(parts[:3]) if len(parts) >= 3 else None,
            ec if len(parts) == 4 else None)


# ---------------------------------------------------------------- EC 源
def load_sifts_chain_ec():
    m = defaultdict(set)
    with gzip.open(SIFTS_EC_TSV, "rt") as f:
        next(f)
        for line in f:
            p = line.strip().split("\t")
            if len(p) >= 4 and p[3] and p[3] != "-":
                m[(p[0].lower(), p[1])].add(p[3])
    return m


def load_ec_json():
    return json.load(open(SIFTS_EC_JSON))


def ec_json_exact(j, pdb, ch):
    ann = j.get(f"{pdb.upper()}_{ch}")
    if not ann:
        return set()
    return {a["ec_number"] for a in ann.get("annotations", [])
            if a.get("ec_number") and a["ec_number"] != "-" and a["ec_number"].count(".") == 3}


def ec_json_fallback(j, pdb, ch):
    for k in (f"{pdb.upper()}_{ch}", f"{pdb.upper()}_A", f"{pdb.upper()}_0"):
        if k in j:
            return {a["ec_number"] for a in j[k].get("annotations", [])
                    if a.get("ec_number") and a["ec_number"] != "-" and a["ec_number"].count(".") == 3}
    return set()


# ---------------------------------------------------------------- 核心计算
def compute(train_ec_sets, test_ids, test_ec_of):
    """返回 DataFrame（index=test_ids，列=NEW_COLS，bool）。"""
    f3, f4 = Counter(), Counter()
    for s in train_ec_sets:
        for x in {levels(e)[0] for e in s} - {None}:
            f3[x] += 1
        for x in {levels(e)[1] for e in s} - {None}:
            f4[x] += 1
    out = {}
    for uid in test_ids:
        s = test_ec_of(uid)
        l3s = {levels(e)[0] for e in s} - {None}
        l4s = {levels(e)[1] for e in s} - {None}
        out[uid] = {
            "OOD_LongTail_EC_L3_le5": bool(l3s) and all(f3.get(x, 0) < 5 for x in l3s),
            "OOD_LongTail_EC_L3_le10": bool(l3s) and all(f3.get(x, 0) < 10 for x in l3s),
            "OOD_LongTail_EC_L4_le5": bool(l4s) and all(f4.get(x, 0) < 5 for x in l4s),
            "OOD_LongTail_EC_L4_le10": bool(l4s) and all(f4.get(x, 0) < 10 for x in l4s),
        }
    return pd.DataFrame.from_dict(out, orient="index").reindex(test_ids)


def apply(test, flags, mask, nan_for_excluded):
    """把 flags 写入 test：mask 行用计算值，其余行 False 或 NaN。返回 (test, counts)。"""
    for col in NEW_COLS:
        if nan_for_excluded:
            test[col] = pd.Series(pd.NA, index=test.index, dtype=object)
        else:
            test[col] = False
        # unique_id 可能重复（如 kcat 同酶多底物），同一 uid 的 flags 相同，用 dict 映射
        lut = flags[col].groupby(level=0).first().to_dict()
        test.loc[mask, col] = test.loc[mask, "unique_id"].map(lut).to_numpy()
    counts = {col: int((test[col] == True).sum()) for col in NEW_COLS}  # noqa: E712
    return test, counts


def run_task(name, task_dir, test_name, train_name, ec_of, train_ids_of,
             mask_of, nan_for_excluded, drop_cols, test_ec_items=None):
    test_path = task_dir / "splits" / test_name
    train_path = task_dir / "splits" / train_name
    test = pd.read_csv(test_path)
    train_ids = train_ids_of(pd.read_csv(train_path))
    test = test.drop(columns=[c for c in drop_cols if c in test.columns and c not in NEW_COLS])
    mask = mask_of(test)
    if test_ec_items is None:
        flags = compute([ec_of(u) for u in train_ids],
                        test.loc[mask, "unique_id"].tolist(), ec_of)
    else:
        items = test_ec_items(test.loc[mask])
        flags = compute([ec_of(u) for u in train_ids],
                        test.loc[mask, "unique_id"].tolist(),
                        dict(zip(test.loc[mask, "unique_id"], items)).get)
    test, counts = apply(test, flags, mask, nan_for_excluded)
    test.to_csv(test_path, index=False)
    summary = {"test_rows": int(len(test)), "computed_rows": int(mask.sum()),
               "dropped_old_cols": drop_cols, **counts}
    out_dir = task_dir / "output"
    out_dir.mkdir(exist_ok=True)
    (out_dir / "longtail_ec_unified_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"{name}: computed={int(mask.sum())}, " +
          " ".join(f"{c.split('EC_')[1]}={n}" for c, n in counts.items()))
    return counts


def main():
    sifts = load_sifts_chain_ec()
    ecj = load_ec_json()
    pdb_ec = defaultdict(set)
    for (p, c), es in sifts.items():
        pdb_ec[p] |= es

    # ---- kcat（UniProt accession → 原始数据 EC） ----
    km = pd.read_csv(POOR / "datasets/enzyme_kinetics_prediction/data/kcat_with_reactant_ecfp4.csv",
                     usecols=["uniprot", "ec"], dtype=str)
    kcat_ec = km.dropna(subset=["ec"]).groupby("uniprot")["ec"].apply(lambda x: {e.strip() for v in x for e in str(v).split(";") if e.strip() and e.strip() != "-"}).to_dict()
    acc_of = lambda uid: uid[3:].rsplit("-F1-model", 1)[0]
    run_task("kcat", POOR / "datasets/enzyme_kinetics_prediction",
             "kcat_test.csv", "kcat_train.csv",
             lambda u: kcat_ec.get(acc_of(u), set()),
             lambda df: df["unique_id"],
             lambda t: t["Default"] == True, nan_for_excluded=True,  # noqa: E712
             drop_cols=["OOD_LongTail_EC_L3"])

    # ---- optimal_ph（metadata ec_id） ----
    meta = pd.read_csv(POOR / "datasets/enzyme_optimal_ph/data/metadata.csv",
                       usecols=["uniprot_id", "ec_id"], dtype=str)
    ph_ec = meta.dropna(subset=["ec_id"]).groupby("uniprot_id")["ec_id"].apply(
        lambda x: {e.strip() for v in x for e in str(v).split(";") if e.strip() and e.strip() != "-"}).to_dict()
    run_task("optimal_ph", POOR / "datasets/enzyme_optimal_ph",
             "optimal_ph_prediction_test.csv", "optimal_ph_prediction_train.csv",
             lambda u: ph_ec.get(u, set()),
             lambda df: df["unique_id"],
             lambda t: t["Default"] == True, nan_for_excluded=True,  # noqa: E712
             drop_cols=["OOD_LongTail_EC_L3_le5", "OOD_LongTail_EC_L3_le10"])

    # ---- ligand_affinity（pdb → 所有链 EC） ----
    run_task("ligand_affinity", POOR / "datasets/ligand_binding_affinity",
             "ligand_binding_affinity_test.csv", "ligand_binding_affinity_train.csv",
             lambda u: pdb_ec.get(u.lower(), set()),
             lambda df: df["unique_id"],
             lambda t: t["Default"] == True, nan_for_excluded=True,  # noqa: E712
             drop_cols=["OOD_LongTail_EC_L4"])

    # ---- ppi（UniProt 缓存 + SIFTS 双路，bio/xtal 解析 unique_id，random 走 protein_*_id） ----
    ec_csv = pd.read_csv(POOR / "datasets/ppi_prediction/output/uniprot_ec_mapping.csv")
    ec_csv["ec_numbers"] = ec_csv["ec_numbers"].fillna("")
    uniprot_ec = {r["uniprot_id"]: {x.strip() for x in str(r["ec_numbers"]).split(";") if x.strip()}
                  for _, r in ec_csv.iterrows()}

    def ppi_row_ec(row):
        uid = str(row["unique_id"])
        ecs = set()
        if row.get("pair_type", "bio") in ("bio", "xtal"):
            for part in uid.split("--"):
                sub = part.split("_")
                if len(sub) < 4:
                    continue
                ecs |= uniprot_ec.get(sub[-1], set())
                ecs |= sifts.get((sub[0].lower(), sub[-2][0]), set())
        else:
            for side in ("A", "B"):
                pid = row.get(f"protein_{side}_id")
                if isinstance(pid, str) and pid:
                    ecs |= uniprot_ec.get(pid, set())
        return ecs

    # ppi 按行（而非 unique_id）取 EC，单独处理
    ppi_dir = POOR / "datasets/ppi_prediction"
    ppi_test = pd.read_csv(ppi_dir / "splits/ppi_test.csv")
    ppi_train = pd.read_csv(ppi_dir / "splits/ppi_train.csv")
    ppi_test = ppi_test.drop(columns=["OOD_LongTail_EC_L4"], errors="ignore")
    f3, f4 = Counter(), Counter()
    for _, row in ppi_train.iterrows():
        s = ppi_row_ec(row)
        for x in {levels(e)[0] for e in s} - {None}:
            f3[x] += 1
        for x in {levels(e)[1] for e in s} - {None}:
            f4[x] += 1
    mask = ppi_test["Default"] == True  # noqa: E712
    cols_data = {c: [] for c in NEW_COLS}
    for _, row in ppi_test.iterrows():
        s = ppi_row_ec(row)
        l3s = {levels(e)[0] for e in s} - {None}
        l4s = {levels(e)[1] for e in s} - {None}
        vals = [bool(l3s) and all(f3.get(x, 0) < 5 for x in l3s),
                bool(l3s) and all(f3.get(x, 0) < 10 for x in l3s),
                bool(l4s) and all(f4.get(x, 0) < 5 for x in l4s),
                bool(l4s) and all(f4.get(x, 0) < 10 for x in l4s)]
        for c, v in zip(NEW_COLS, vals):
            cols_data[c].append(v)
    for c in NEW_COLS:
        # Default=False（极端长度）行置 NaN（全项目统一惯例，流水线位置直接写入）
        ppi_test[c] = [v if m else pd.NA for v, m in zip(cols_data[c], mask)]
    counts = {c: int((ppi_test[c] == True).sum()) for c in NEW_COLS}  # noqa: E712
    ppi_test.to_csv(ppi_dir / "splits/ppi_test.csv", index=False)
    (ppi_dir / "output/longtail_ec_unified_summary.json").write_text(json.dumps(
        {"test_rows": len(ppi_test), "computed_rows": int(mask.sum()),
         "dropped_old_cols": ["OOD_LongTail_EC_L4"], **counts}, indent=2))
    print(f"ppi: computed={int(mask.sum())}, " +
          " ".join(f"{c.split('EC_')[1]}={n}" for c, n in counts.items()))

    # ---- ss ----
    key_ss = lambda u: (u[:4].lower(), u[5:])
    run_task("ss", POOR / "datasets/ss_prediction",
             "ssp_test.csv", "ssp_train.csv",
             lambda u: sifts.get(key_ss(u), set()),
             lambda df: df["unique_id"],
             lambda t: t["Default"] == True, nan_for_excluded=True,  # noqa: E712
             drop_cols=["OOD_LongTail_EC_L4"])

    # ---- lbs ----
    key_lbs = lambda u: (u.split("_")[0].lower(), u.split("_")[1])
    run_task("lbs", POOR / "datasets/ligand_binding_site",
             "ligand_binding_site_test.csv", "ligand_binding_site_train.csv",
             lambda u: sifts.get(key_lbs(u), set()),
             lambda df: df["unique_id"],
             lambda t: t["Default"] == True, nan_for_excluded=True,  # noqa: E712
             drop_cols=["OOD_LongTail_EC_L4"])

    # ---- ppis（seq_hash 映射链） ----
    with open(POOR / "datasets/ppis_prediction/output/dips_plus_metadata.csv") as f:
        pair_meta = {r["unique_id"]: r for r in csv.DictReader(f)}
    with open(POOR / "datasets/ppis_prediction/output/chain_centric_dedup_metadata.csv") as f:
        meta_by_hash = {r["seq_hash"]: r for r in csv.DictReader(f)}

    def ppis_ec(h):
        m = meta_by_hash.get(h, {})
        uid, role = m.get("structure_source_unique_id", ""), m.get("structure_source_role", "")
        if not uid or not role:
            return set()
        role = role.split("->")[0].strip()
        pm = pair_meta.get(uid)
        pdb = uid.split(".")[0].lower()
        ch = pm.get("receptor_chain", "") if role == "receptor" else (
            pm.get("ligand_chain", "") if role == "ligand" else "")
        return sifts.get((pdb, ch), set()) if ch else set()

    run_task("ppis", POOR / "datasets/ppis_prediction",
             "ppis_test.csv", "ppis_train.csv",
             ppis_ec,
             lambda df: df["unique_id"],
             lambda t: t["Default"] == True, nan_for_excluded=True,  # noqa: E712
             drop_cols=["OOD_LongTail_EC_L4"])

    # ---- fold（非 ExtremeShort 行计算，ES 行 NaN） ----
    run_task("fold", POOR / "datasets/fold_classification",
             "cath_test.csv", "cath_train.csv",
             lambda u: ec_json_exact(ecj, u[:4], u[4]),
             lambda df: df["unique_id"],
             lambda t: t["OOD_ExtremeShort"] != True, nan_for_excluded=True,  # noqa: E712
             drop_cols=["OOD_LongTail_EC_L4"])

    # ---- PDA cath_design（非 ES 行计算，ES 行 NaN；train 参考 = fold train） ----
    pda_dir = POOR / "datasets/protein-design-archive"
    fold_train_ids = pd.read_csv(
        POOR / "datasets/fold_classification/splits/cath_train.csv",
        usecols=["unique_id"])["unique_id"].tolist()
    pda_test = pd.read_csv(pda_dir / "splits/cath_design.csv")
    pda_test = pda_test.drop(columns=["OOD_LongTail_EC_L4"], errors="ignore")
    mask = pda_test["OOD_ExtremeShort"] != True  # noqa: E712
    flags = compute([ec_json_exact(ecj, u[:4], u[4]) for u in fold_train_ids],
                    pda_test.loc[mask, "unique_id"].tolist(),
                    lambda u: ec_json_fallback(ecj, u.split("_")[0], u.split("_", 1)[1]))
    pda_test, counts = apply(pda_test, flags, mask, nan_for_excluded=True)
    pda_test.to_csv(pda_dir / "splits/cath_design.csv", index=False)
    (pda_dir / "output/longtail_ec_unified_summary.json").write_text(json.dumps(
        {"test_rows": len(pda_test), "computed_rows": int(mask.sum()),
         "dropped_old_cols": ["OOD_LongTail_EC_L4"], **counts}, indent=2))
    print(f"PDA: computed={int(mask.sum())}, " +
          " ".join(f"{c.split('EC_')[1]}={n}" for c, n in counts.items()))

    # ---- func_ec（仅 L4 两列；L3 级 LongTailEC 与任务预测层级冲突，不加） ----
    # func ec 的 label 为 L3 级（x.x.x.- 形式），其 OOD_LongTail 已是 L3 级 EC 长尾；
    # LongTailEC 使用 SIFTS 链级真实 L4 注释（与 NewEC_L4 同源），all-not-in。
    func_dir = POOR / "datasets/func_prediction"
    key_f = lambda u: (u[:4].lower(), u[5:])
    func_test = pd.read_csv(func_dir / "splits/ec_test.csv")
    func_train_ids = pd.read_csv(func_dir / "splits/ec_train.csv",
                                 usecols=["unique_id"])["unique_id"].tolist()
    func_test = func_test.drop(
        columns=["OOD_LongTail_EC_L3_le5", "OOD_LongTail_EC_L3_le10"],
        errors="ignore")
    mask = func_test["Default"] == True  # noqa: E712
    flags = compute([sifts.get(key_f(u), set()) for u in func_train_ids],
                    func_test.loc[mask, "unique_id"].tolist(),
                    lambda u: sifts.get(key_f(u), set()))
    for col in ["OOD_LongTail_EC_L4_le5", "OOD_LongTail_EC_L4_le10"]:
        lut = flags[col].groupby(level=0).first().to_dict()
        # Default=False（极端长度）行置 NaN（全项目统一惯例，流水线位置直接写入）
        func_test[col] = pd.Series(pd.NA, index=func_test.index, dtype=object)
        func_test.loc[mask, col] = func_test.loc[mask, "unique_id"].map(lut).to_numpy()
    counts = {c: int((func_test[c] == True).sum())  # noqa: E712
              for c in ["OOD_LongTail_EC_L4_le5", "OOD_LongTail_EC_L4_le10"]}
    func_test.to_csv(func_dir / "splits/ec_test.csv", index=False)
    (func_dir / "output/longtail_ec_unified_summary.json").write_text(json.dumps(
        {"test_rows": len(func_test), "computed_rows": int(mask.sum()),
         "note": "L4 only, SIFTS chain source, all-not-in; L3 cols removed",
         "dropped_old_cols": ["OOD_LongTail_EC_L3_le5", "OOD_LongTail_EC_L3_le10"],
         **counts}, indent=2))
    print(f"func_ec: computed={int(mask.sum())}, " +
          " ".join(f"{c.split('EC_')[1]}={n}" for c, n in counts.items()))


if __name__ == "__main__":
    main()
