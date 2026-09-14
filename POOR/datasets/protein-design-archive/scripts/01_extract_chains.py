#!/usr/bin/env python3
"""Step 1 (v2): Extract ALL non-natural (chain_type D+U) chain records from PDA.

v2 changes vs v1:
  - pool = all chains with chain_type in {D, U} (not just synthetic-construct
    source); chain_type N (natural) chains are excluded entirely.
- Comma-merged chain_ids (e.g. "A,B,C") are identical copies of one entity;
  take the FIRST chain id as the representative.
- unique_id = {pdb}_{chain}
- label_candidate: entry-level distinct C.A.T set from cath_full; only
  unambiguous when exactly ONE distinct C.A.T exists.

Output: output/design_chains_all.csv
"""

import json
from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parent.parent
POOR = BASE.parent.parent
DATA_JSON = POOR / "data" / "protein-design-archive" / "backend" / "scripts" / "data.json"
OUT_DIR = BASE / "output"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def main():
    with open(DATA_JSON) as f:
        data = json.load(f)

    rows = []
    n_comma = 0
    for e in data:
        pdb = e["pdb"].lower()
        cath_codes = [c["code"] for c in (e.get("cath_full") or []) if c.get("code")]
        cat_set = sorted({".".join(code.split(".")[:3]) for code in cath_codes})
        label_cand = cat_set[0] if len(cat_set) == 1 else ""
        for c in e.get("chains", []):
            ctype = c.get("chain_type") or "U"
            if ctype == "N":  # natural-source chains are not used at all
                continue
            chain_id = c["chain_id"]
            if "," in chain_id:
                n_comma += 1
                chain_id = chain_id.split(",")[0]
            seq = (c.get("chain_seq_nat") or "").strip()
            if not seq:
                continue
            rows.append({
                "unique_id": f"{pdb}_{chain_id}",
                "pdb_id": pdb,
                "chain_id": chain_id,
                "chain_type": ctype,
                "chain_source": c.get("chain_source") or "",
                "aa_seq": seq,
                "seq_len": len(seq),
                "cath_codes": ";".join(cath_codes),
                "cath_cat_set": ";".join(cat_set),
                "label_candidate": label_cand,
                "release_date": e.get("release_date", ""),
                "exptl_method": ";".join(e.get("exptl_method") or []),
                "tags": ";".join(e.get("tags") or []),
            })

    df = pd.DataFrame(rows)
    assert df["unique_id"].is_unique, "unique_id collision!"
    df.to_csv(OUT_DIR / "design_chains_all.csv", index=False)

    print(f"chain records (D+U): {len(df)}")
    print(f"  chain_type: {df['chain_type'].value_counts().to_dict()}")
    print(f"  comma-merged chain_ids (took first): {n_comma}")
    print(f"  unique PDB entries: {df['pdb_id'].nunique()}")
    print(f"  unique sequences: {df['aa_seq'].nunique()}")
    print(f"  unambiguous label_candidate: {(df['label_candidate'] != '').sum()}")
    print(f"  seq_len < 60: {(df['seq_len'] < 60).sum()}, "
          f"> 1000: {(df['seq_len'] > 1000).sum()}")


if __name__ == "__main__":
    main()
