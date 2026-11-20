#!/usr/bin/env python3
"""Add six EC function OOD columns to func_prediction go_mf / go_bp / go_cc test sets.

Background: GO OOD (NewGO/LongTailGO) is not viable: the closed label vocabulary makes NewGO always zero;
all-not-in over full GO annotations is also always zero. The selected alternative extends unified EC OOD
(NewEC / LongTailEC) to the three go_* test sets.

Definition (identical to add_unified_newec.py / add_unified_longtail_ec.py):
- EC source: SIFTS chain-level EC（data/sifts/sifts_chain_ec.tsv.gz），
  uid mapping (u[:4].lower(), u[5:])，as in func_ec
- levels()：at least three components -> L3=first three; four components -> L4 (partial "-" allowed)
- OOD_NewEC_L3/L4：reference train+val: annotated and all ECs at the given level absent from reference -> True
  (all-not-in); no EC -> False
- OOD_LongTail_EC_{L3,L4}_le{5,10}：reference training-only per-sample frequencies, including zero;
  annotated and all EC training frequencies at the given level < threshold -> True (all-not-in); no EC -> False
- Default=False rows (preselected extreme-length samples): all six columns are NaN (not computed, distinct from False)
- Append six columns to the CSV.

Usage: python3 .skills/ec-function-ood-annotation/scripts/add_go_ec_ood.py
"""
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

# Repository root: three levels above .skills/ec-function-ood-annotation/scripts/.
POOR = Path(__file__).resolve().parents[3]
SIFTS_EC_TSV = POOR / "data" / "sifts" / "sifts_chain_ec.tsv.gz"
FUNC = POOR / "datasets" / "func_prediction"

L3_COL, L4_COL = "OOD_NewEC_L3", "OOD_NewEC_L4"
LT_COLS = ["OOD_LongTail_EC_L3_le5", "OOD_LongTail_EC_L3_le10",
           "OOD_LongTail_EC_L4_le5", "OOD_LongTail_EC_L4_le10"]


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


def main():
    sifts = load_sifts_chain_ec()
    ec_of = lambda u: sifts.get((u[:4].lower(), u[5:]), set())  # noqa: E731

    summary = {}
    for task in ["go_mf", "go_bp", "go_cc"]:
        test_path = FUNC / "splits" / f"{task}_test.csv"
        test = pd.read_csv(test_path)
        train_ids = pd.read_csv(FUNC / "splits" / f"{task}_train.csv",
                                usecols=["unique_id"])["unique_id"]
        val_ids = pd.read_csv(FUNC / "splits" / f"{task}_val.csv",
                              usecols=["unique_id"])["unique_id"]

        # References: train+val for NewEC; training-only per-sample frequencies for LongTailEC.
        tv3, tv4 = set(), set()
        f3, f4 = Counter(), Counter()
        for u in train_ids:
            s = ec_of(u)
            l3s = {levels(e)[0] for e in s} - {None}
            l4s = {levels(e)[1] for e in s} - {None}
            tv3 |= l3s
            tv4 |= l4s
            for x in l3s:
                f3[x] += 1
            for x in l4s:
                f4[x] += 1
        for u in val_ids:
            s = ec_of(u)
            tv3 |= {levels(e)[0] for e in s} - {None}
            tv4 |= {levels(e)[1] for e in s} - {None}

        def flags(u):
            s = ec_of(u)
            l3s = {levels(e)[0] for e in s} - {None}
            l4s = {levels(e)[1] for e in s} - {None}
            return {
                L3_COL: bool(l3s) and all(x not in tv3 for x in l3s),
                L4_COL: bool(l4s) and all(x not in tv4 for x in l4s),
                LT_COLS[0]: bool(l3s) and all(f3.get(x, 0) < 5 for x in l3s),
                LT_COLS[1]: bool(l3s) and all(f3.get(x, 0) < 10 for x in l3s),
                LT_COLS[2]: bool(l4s) and all(f4.get(x, 0) < 5 for x in l4s),
                LT_COLS[3]: bool(l4s) and all(f4.get(x, 0) < 10 for x in l4s),
            }

        mask = test["Default"] == True  # noqa: E712
        cols = [L3_COL, L4_COL] + LT_COLS
        lut = {u: flags(u) for u in test.loc[mask, "unique_id"]}
        for col in cols:
            ser = pd.Series(pd.NA, index=test.index, dtype=object)
            ser[mask] = test.loc[mask, "unique_id"].map(
                {u: f[col] for u, f in lut.items()}).to_numpy()
            test[col] = ser
        test.to_csv(test_path, index=False)

        counts = {col: int((test[col] == True).sum()) for col in cols}  # noqa: E712
        summary[task] = {"test_rows": int(len(test)),
                         "computed_rows": int(mask.sum()),
                         "no_ec_rows": int(sum(1 for u in test.loc[mask, "unique_id"]
                                               if not ec_of(u))),
                         **counts}
        print(f"{task}: computed={int(mask.sum())}, " +
              " ".join(f"{c.replace('OOD_', '')}={n}" for c, n in counts.items()))

    (FUNC / "output" / "go_ec_ood_summary.json").write_text(
        json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
