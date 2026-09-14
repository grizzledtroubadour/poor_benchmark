#!/usr/bin/env python3
"""按 resseq 连续段锚定修复"结构残基数 != aa_seq 长度"行的 struct_label。

这些行的结构缺失若干残基，add_struct_label.py 的序列比对在低复杂度区
缺口放置任意，导致界面标签错帧（详见 DATA_PROCESS.md 5.7 注记）。
本脚本利用 pdb 保留的 resseq 编号做确定性映射：
  - 结构残基按 (chain,resseq,icode) 首次出现枚举，resseq 连续段独立处理；
  - 每段在 aa_seq 中找**唯一精确匹配**锚定（序列一致率本就 100%）；
    唯一则段内按 resseq 顺序同步映射；不唯一/无匹配则该段残基置 -1；
  - 跨段校验：映射的 aa_seq 位置必须单调不重叠，否则整行报错不改。
消除了比对缺口放置的任意性。

Usage: python3 14_fix_struct_label_resseq.py [--apply]
"""
import os, sys
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
TD = os.path.join(ROOT, "datasets", "ppis_prediction")
SPLITS = ["ppis_train.csv", "ppis_val.csv", "ppis_test.csv"]

RES3 = {
    'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I',
    'LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S',
    'THR':'T','VAL':'V','TRP':'W','TYR':'Y','MSE':'M','SEC':'C','PYL':'K',
}

def parse_pdb(path):
    res=[]; seen=set()
    for l in open(path, errors="ignore"):
        if not l.startswith(("ATOM","HETATM")): continue
        key=(l[21], l[22:26], l[26])
        if key in seen: continue
        seen.add(key)
        res.append((int(l[22:26]), RES3.get(l[17:20].strip().upper(), 'X')))
    return res

def fix_row(residues, aa_seq, labels):
    """返回 (new_struct_label, stats) 或 (None, reason)。"""
    n = len(residues)
    # 切连续 resseq 段
    runs=[]; cur=[0]
    for i in range(1, n):
        if residues[i][0] == residues[i-1][0] + 1:
            cur.append(i)
        else:
            runs.append(cur); cur=[i]
    runs.append(cur)
    pos_map = {}
    unanchored = 0
    for run in runs:
        sub = "".join(residues[i][1] for i in run)
        hits=[]; start=0
        while True:
            j = aa_seq.find(sub, start)
            if j==-1: break
            hits.append(j); start=j+1
        if len(hits)==1:
            j0=hits[0]
            for k,i in enumerate(run):
                pos_map[i]=j0+k
        else:
            unanchored += len(run)
    # 单调性校验
    mapped = sorted(pos_map.values())
    if len(mapped)!=len(set(mapped)):
        return None, "overlap"
    new=[labels[pos_map[i]] if i in pos_map else -1 for i in range(n)]
    return new, {"unanchored": unanchored, "runs": len(runs)}

def main():
    apply = "--apply" in sys.argv
    for f in SPLITS:
        df = pd.read_csv(os.path.join(TD,"splits",f), dtype=str, keep_default_na=False)
        fixed=0
        for idx, r in df.iterrows():
            residues = parse_pdb(os.path.join(TD,"pdbs",r["struct_file"]))
            n=len(residues)
            cur=[int(x) for x in r["struct_label"].strip("[]").split()]
            if len(cur)!=n or n==len(r["aa_seq"]):
                continue  # 直通/等长路径无风险
            labels=[int(x) for x in r["label"].strip("[]").split()]
            new, st = fix_row(residues, r["aa_seq"], labels)
            if new is None:
                print(f"  [FAIL] {r['unique_id'][:12]} ({f}): {st}"); continue
            if new!=cur:
                d = np.array(new)!=np.array(cur)
                d1 = int((d & ((np.array(new)==1)|(np.array(cur)==1))).sum())
                print(f"  [FIX] {r['unique_id'][:12]} ({f}): {int(d.sum())} pos changed, {d1} at interface, unanchored={st['unanchored']}")
                if apply:
                    df.at[idx,"struct_label"]="[ "+" ".join(map(str,new))+" ]"
                fixed+=1
        if apply and fixed:
            df.to_csv(os.path.join(TD,"splits",f), index=False, lineterminator="\n")
        print(f"{f}: {fixed} rows fixed")
if __name__=="__main__":
    main()
