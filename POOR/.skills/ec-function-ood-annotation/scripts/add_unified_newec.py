#!/usr/bin/env python3
"""统一 NewEC OOD 口径并落盘到 7 个任务的 test csv（2026-08-27）。

统一定义：
- 参考集：train+val 中出现的 EC 集合（对应层级）
- 语义：**all-not-in** —— 样本有 EC 注释、且其所有 EC（对应层级）都不在
  train+val 集合中 → True；无 EC 注释 → False
- 层级解析与 LongTailEC 统一脚本一致（levels()：段数 ≥3 → L3=前三段；
  段数 ==4 → L4=完整编号；允许部分 "-"）
- L4 列（OOD_NewEC_L4）就地重算，L3 列（OOD_NewEC_L3）追加到末尾
- 各任务 Default=False（极端长度预挑出）行一律置 NaN（"未计算"语义，
  全项目统一惯例，流水线位置直接写入；finalize --mode nan 仅作兜底校验）；
  affinity 全部 Default=True，不触发该惯例
- 排除：fold（EC7 prior）、PDA、func go_*；func_ec 仅 L4（L3 是预测目标层级）

用法：python3 .skills/ec-function-ood-annotation/scripts/add_unified_newec.py
"""
import csv
import gzip
import json
from collections import defaultdict
from pathlib import Path

import pandas as pd

# 仓库根：.skills/ec-function-ood-annotation/scripts/ 上溯 3 级
POOR = Path(__file__).resolve().parents[3]
SIFTS_EC_TSV = POOR / "data" / "sifts" / "sifts_chain_ec.tsv.gz"

L3_COL, L4_COL = "OOD_NewEC_L3", "OOD_NewEC_L4"


def levels(ec):
    parts = ec.split(".")
    return (".".join(parts[:3]) if len(parts) >= 3 else None,
            ec if len(parts) == 4 else None)


def load_sifts_chain_ec():
    m = defaultdict(set)
    with gzip.open(SIFTS_EC_TSV, "rt") as f:
        next(f)
        for line in f:
            p = line.strip().split("\t")
            if len(p) >= 4 and p[3] and p[3] != "-":
                m[(p[0].lower(), p[1])].add(p[3])
    return m


def trainval_sets(train_ec_sets, val_ec_sets):
    tv3, tv4 = set(), set()
    for s in list(train_ec_sets) + list(val_ec_sets):
        tv3 |= {levels(e)[0] for e in s} - {None}
        tv4 |= {levels(e)[1] for e in s} - {None}
    return tv3, tv4


def newec_flags(ec_set, tv3, tv4):
    l3s = {levels(e)[0] for e in ec_set} - {None}
    l4s = {levels(e)[1] for e in ec_set} - {None}
    return (bool(l3s) and all(x not in tv3 for x in l3s),
            bool(l4s) and all(x not in tv4 for x in l4s))


def apply(test, flags_by_uid, mask, nan_for_excluded, add_l3=True):
    """flags_by_uid: {unique_id: (l3_flag, l4_flag)}。L4 就地重算，L3 追加（add_l3=False 时不加）。"""
    lut3 = {u: f[0] for u, f in flags_by_uid.items()}
    lut4 = {u: f[1] for u, f in flags_by_uid.items()}
    if nan_for_excluded:
        new_l3 = pd.Series(pd.NA, index=test.index, dtype=object)
        l4_vals = pd.Series(pd.NA, index=test.index, dtype=object)
    else:
        new_l3 = pd.Series(False, index=test.index)
        l4_vals = pd.Series(False, index=test.index)
    l4_vals[mask] = test.loc[mask, "unique_id"].map(lut4).to_numpy()
    test[L4_COL] = l4_vals          # 就地重算（保持列位置）
    counts = {L4_COL: int((test[L4_COL] == True).sum())}  # noqa: E712
    if add_l3:
        new_l3[mask] = test.loc[mask, "unique_id"].map(lut3).to_numpy()
        test[L3_COL] = new_l3       # 追加到末尾
        counts[L3_COL] = int((test[L3_COL] == True).sum())  # noqa: E712
    return test, counts


def run_task(name, task_dir, test_name, train_name, val_name, ec_of,
             mask_of, nan_for_excluded, add_l3=True):
    test_path = task_dir / "splits" / test_name
    test = pd.read_csv(test_path)
    train_ids = pd.read_csv(task_dir / "splits" / train_name, usecols=["unique_id"])["unique_id"]
    val_ids = pd.read_csv(task_dir / "splits" / val_name, usecols=["unique_id"])["unique_id"]
    tv3, tv4 = trainval_sets([ec_of(u) for u in train_ids],
                             [ec_of(u) for u in val_ids])
    mask = mask_of(test)
    flags = {u: newec_flags(ec_of(u), tv3, tv4) for u in test.loc[mask, "unique_id"]}
    test, counts = apply(test, flags, mask, nan_for_excluded, add_l3=add_l3)
    test.to_csv(test_path, index=False)
    summary = {"test_rows": int(len(test)), "computed_rows": int(mask.sum()), **counts}
    (task_dir / "output" / "newec_unified_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"{name}: computed={int(mask.sum())}, " +
          " ".join(f"{c.replace('OOD_NewEC_', '')}={n}" for c, n in counts.items()))
    return counts


def main():
    sifts = load_sifts_chain_ec()
    pdb_ec = defaultdict(set)
    for (p, c), es in sifts.items():
        pdb_ec[p] |= es

    # ---- kcat ----
    km = pd.read_csv(POOR / "datasets/enzyme_kinetics_prediction/data/kcat_with_reactant_ecfp4.csv",
                     usecols=["uniprot", "ec"], dtype=str)
    kcat_ec = km.dropna(subset=["ec"]).groupby("uniprot")["ec"].apply(lambda x: {e.strip() for v in x for e in str(v).split(";") if e.strip() and e.strip() != "-"}).to_dict()
    acc_of = lambda uid: uid[3:].rsplit("-F1-model", 1)[0]
    run_task("kcat", POOR / "datasets/enzyme_kinetics_prediction",
             "kcat_test.csv", "kcat_train.csv", "kcat_val.csv",
             lambda u: kcat_ec.get(acc_of(u), set()),
             lambda t: t["Default"] == True, nan_for_excluded=True)  # noqa: E712

    # ---- optimal_ph ----
    meta = pd.read_csv(POOR / "datasets/enzyme_optimal_ph/data/metadata.csv",
                       usecols=["uniprot_id", "ec_id"], dtype=str)
    ph_ec = meta.dropna(subset=["ec_id"]).groupby("uniprot_id")["ec_id"].apply(
        lambda x: {e.strip() for v in x for e in str(v).split(";") if e.strip() and e.strip() != "-"}).to_dict()
    run_task("optimal_ph", POOR / "datasets/enzyme_optimal_ph",
             "optimal_ph_prediction_test.csv", "optimal_ph_prediction_train.csv",
             "optimal_ph_prediction_val.csv",
             lambda u: ph_ec.get(u, set()),
             lambda t: t["Default"] == True, nan_for_excluded=True)  # noqa: E712

    # ---- ligand_affinity ----
    run_task("ligand_affinity", POOR / "datasets/ligand_binding_affinity",
             "ligand_binding_affinity_test.csv", "ligand_binding_affinity_train.csv",
             "ligand_binding_affinity_val.csv",
             lambda u: pdb_ec.get(u.lower(), set()),
             lambda t: t["Default"] == True, nan_for_excluded=True)  # noqa: E712

    # ---- ppi（行级 EC，双路并集） ----
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

    ppi_dir = POOR / "datasets/ppi_prediction"
    ppi_test = pd.read_csv(ppi_dir / "splits/ppi_test.csv")
    ppi_train = pd.read_csv(ppi_dir / "splits/ppi_train.csv")
    ppi_val = pd.read_csv(ppi_dir / "splits/ppi_val.csv")
    tv3, tv4 = set(), set()
    for d in (ppi_train, ppi_val):
        for _, row in d.iterrows():
            s = ppi_row_ec(row)
            tv3 |= {levels(e)[0] for e in s} - {None}
            tv4 |= {levels(e)[1] for e in s} - {None}
    mask = ppi_test["Default"] == True  # noqa: E712
    # Default=False（极端长度）行置 NaN（全项目统一惯例，流水线位置直接写入）
    l3_vals = pd.Series(pd.NA, index=ppi_test.index, dtype=object)
    l4_vals = pd.Series(pd.NA, index=ppi_test.index, dtype=object)
    l3_list, l4_list = [], []
    for _, row in ppi_test.iterrows():
        f3, f4 = newec_flags(ppi_row_ec(row), tv3, tv4)
        l3_list.append(f3)
        l4_list.append(f4)
    l3_vals[mask] = [v for v, m in zip(l3_list, mask) if m]
    l4_vals[mask] = [v for v, m in zip(l4_list, mask) if m]
    ppi_test[L4_COL] = l4_vals
    ppi_test[L3_COL] = l3_vals
    counts = {L3_COL: int((ppi_test[L3_COL] == True).sum()),  # noqa: E712
              L4_COL: int((ppi_test[L4_COL] == True).sum())}  # noqa: E712
    ppi_test.to_csv(ppi_dir / "splits/ppi_test.csv", index=False)
    (ppi_dir / "output/newec_unified_summary.json").write_text(json.dumps(
        {"test_rows": len(ppi_test), "computed_rows": int(mask.sum()), **counts}, indent=2))
    print(f"ppi: computed={int(mask.sum())}, L3={counts[L3_COL]}, L4={counts[L4_COL]}")

    # ---- ss ----
    run_task("ss", POOR / "datasets/ss_prediction",
             "ssp_test.csv", "ssp_train.csv", "ssp_val.csv",
             lambda u: sifts.get((u[:4].lower(), u[5:]), set()),
             lambda t: t["Default"] == True, nan_for_excluded=True)  # noqa: E712

    # ---- lbs ----
    run_task("lbs", POOR / "datasets/ligand_binding_site",
             "ligand_binding_site_test.csv", "ligand_binding_site_train.csv",
             "ligand_binding_site_val.csv",
             lambda u: sifts.get((u.split("_")[0].lower(), u.split("_")[1]), set()),
             lambda t: t["Default"] == True, nan_for_excluded=True)  # noqa: E712

    # ---- ppis ----
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
             "ppis_test.csv", "ppis_train.csv", "ppis_val.csv",
             ppis_ec,
             lambda t: t["Default"] == True, nan_for_excluded=True)  # noqa: E712

    # ---- func_ec（仅 L4；L3 是预测目标层级，不加） ----
    # EC 源用 SIFTS 链级 tsv（与 LongTailEC_L4 同源），all-not-in，train+val 参考
    run_task("func_ec", POOR / "datasets/func_prediction",
             "ec_test.csv", "ec_train.csv", "ec_val.csv",
             lambda u: sifts.get((u[:4].lower(), u[5:]), set()),
             lambda t: t["Default"] == True, nan_for_excluded=True,  # noqa: E712
             add_l3=False)


if __name__ == "__main__":
    main()
