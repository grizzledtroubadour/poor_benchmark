#!/usr/bin/env python3
"""
生成 PPIS 最终数据划分文件。

输出三个文件到 datasets/ppis_prediction/splits/：
- ppis_train.csv
- ppis_val.csv
- ppis_test.csv（含 OOD 标记）

每行包含 4 个核心字段：
- unique_id: 序列 MD5 哈希，对应 PDB 文件 {unique_id}.pdb
- aa_seq: 从 PDB 文件中提取的氨基酸序列
- file_type: 文件类型，固定为 "pdb"
- labels: 逐残基界面标签列表（JSON 字符串），与 aa_seq 严格等长

测试集额外附加 OOD 列，布局如下：
- Default: 非 ExtremeLength 样本（长度在 [60, 1000]）
- OOD_LowHomology_90/80/.../30
- OOD_Orphan
- OOD_NewFold_0.9/0.8/.../0.3
- OOD_IDR
- OOD_CATH_FoldHoldout: CATH Topology (C.A.T) 未在 train 出现
- OOD_CATH_SuperfamilyHoldout: CATH Homologous Superfamily (C.A.T.H) 未在 train 出现
- OOD_NewFunction: 所有 EC 编号均未在 train/val 出现
- OOD_LongTail: 界面正样本比例低于阈值（默认 5%）
- ExtremeLong: 长度 > 1000
- ExtremeShort: 长度 < 60

除 ExtremeLong/ExtremeShort 外的所有 OOD 列均只对 Default 样本计算；
极端长度（Default=False）样本的上述 OOD 列直接写 NaN（空字符串，
"未计算"语义，全项目统一惯例），finalize_ood_columns.py --mode nan 仅作兜底校验。

当 PDB 序列长度与原始 labels 长度不一致时，使用全局序列对齐将 labels 映射到 PDB 序列。
"""

import argparse
import csv
import gzip
import json
import logging
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from Bio.Align import PairwiseAligner

Path("logs").mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("logs/build_final_ppis_splits.log"),
        logging.StreamHandler(sys.stdout),
    ],
)
logger = logging.getLogger(__name__)

RESIDUE_3TO1 = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
}

# 测试集 OOD 列顺序：Default + 常规 OOD + CATH OOD + Function OOD + ExtremeLong/ExtremeShort
OOD_COLUMNS = [
    "Default",
    "OOD_LowHomology_90", "OOD_LowHomology_80", "OOD_LowHomology_70",
    "OOD_LowHomology_60", "OOD_LowHomology_50", "OOD_LowHomology_40",
    "OOD_LowHomology_30", "OOD_Orphan",
    "OOD_NewFold_0.9", "OOD_NewFold_0.8", "OOD_NewFold_0.7",
    "OOD_NewFold_0.6", "OOD_NewFold_0.5", "OOD_NewFold_0.4",
    "OOD_NewFold_0.3", "OOD_IDR",
    "OOD_CATH_FoldHoldout", "OOD_CATH_SuperfamilyHoldout",
    "OOD_NewFunction", "OOD_LongTail",
    "ExtremeLong", "ExtremeShort",
]

# CATH / EC 文件路径（项目根 = POOR/）
CATH_DOMAIN_LIST = (
    Path(__file__).resolve().parents[3]
    / "data/CATHv44/cath-classification-data/cath-domain-list.txt"
)
EC_TSV_PATH = (
    Path(__file__).resolve().parents[3]
    / "data/sifts/sifts_chain_ec.tsv.gz"
)

# 全局对齐器：metadata 序列作为 target，pdb 序列作为 query
ALIGNER = PairwiseAligner()
ALIGNER.mode = "global"
ALIGNER.match_score = 2
ALIGNER.mismatch_score = -1
ALIGNER.open_gap_score = -0.5
ALIGNER.extend_gap_score = -0.1


def load_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def load_cath_domain_list(path):
    """Load CATH domain list: (pdb_id, chain_id) -> list of domain level dicts."""
    chain_domains = defaultdict(list)
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 9:
                continue
            domain_name = parts[0]
            pdb_id = domain_name[:4].lower()
            chain_id = domain_name[4]
            C, A, T, H = parts[1], parts[2], parts[3], parts[4]
            chain_domains[(pdb_id, chain_id)].append({"C": C, "A": A, "T": T, "H": H})
    return chain_domains


def get_source_chain_from_metadata(seq_hash, pair_meta, metadata_by_hash):
    """根据 chain metadata 还原来源 PDB ID 与 chain ID。"""
    m = metadata_by_hash.get(seq_hash, {})
    uid = m.get("structure_source_unique_id", "")
    role_info = m.get("structure_source_role", "")
    if not uid or not role_info:
        return None, None
    role = role_info.split("->")[0].strip()
    pm = pair_meta.get(uid)
    pdb_id = uid.split(".")[0].lower()
    if role == "receptor":
        chain = pm.get("receptor_chain", "") if pm else ""
    elif role == "ligand":
        chain = pm.get("ligand_chain", "") if pm else ""
    else:
        chain = ""
    return pdb_id, chain


def build_train_cath_sets(train_records, chain_domains, pair_meta, metadata_by_hash):
    """Build train sets of CATH Fold (C.A.T) and Superfamily (C.A.T.H) IDs."""
    fold_ids = set()
    superfamily_ids = set()
    for r in train_records:
        pdb_id, chain = get_source_chain_from_metadata(
            r["seq_hash"], pair_meta, metadata_by_hash
        )
        if not pdb_id or not chain:
            continue
        for d in chain_domains.get((pdb_id, chain), []):
            fold_ids.add(f"{d['C']}.{d['A']}.{d['T']}")
            superfamily_ids.add(f"{d['C']}.{d['A']}.{d['T']}.{d['H']}")
    return fold_ids, superfamily_ids


def compute_cath_ood(seq_hash, chain_domains, pair_meta, metadata_by_hash,
                     train_folds, train_superfamilies):
    """Return (fold_ood, superfamily_ood) for a chain. Missing CATH -> (False, False)."""
    pdb_id, chain = get_source_chain_from_metadata(seq_hash, pair_meta, metadata_by_hash)
    domains = chain_domains.get((pdb_id, chain), []) if pdb_id and chain else []
    if not domains:
        return False, False
    fold_ood = any(f"{d['C']}.{d['A']}.{d['T']}" not in train_folds for d in domains)
    superfamily_ood = any(
        f"{d['C']}.{d['A']}.{d['T']}.{d['H']}" not in train_superfamilies for d in domains
    )
    return fold_ood, superfamily_ood


def load_ec_annotations(path):
    """Load SIFTS PDB-chain to EC annotation mapping."""
    ec_map = defaultdict(set)
    with gzip.open(path, "rt") as f:
        next(f)  # header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            pdb_id = parts[0].lower()
            chain_id = parts[1]
            ec = parts[3]
            ec_map[(pdb_id, chain_id)].add(ec)
    return ec_map


def build_trainval_ec_set(train_records, val_records, ec_map, pair_meta, metadata_by_hash):
    """Aggregate EC numbers observed in train and val sets."""
    trainval_ec = set()
    for r in train_records + val_records:
        pdb_id, chain = get_source_chain_from_metadata(r["seq_hash"], pair_meta, metadata_by_hash)
        if pdb_id and chain:
            trainval_ec.update(ec_map.get((pdb_id, chain), set()))
    return trainval_ec


def compute_new_function(seq_hash, ec_map, pair_meta, metadata_by_hash, trainval_ec):
    """True iff chain has EC annotations and all are unseen in train+val."""
    pdb_id, chain = get_source_chain_from_metadata(seq_hash, pair_meta, metadata_by_hash)
    ecs = ec_map.get((pdb_id, chain), set()) if pdb_id and chain else set()
    if not ecs:
        return False
    return all(ec not in trainval_ec for ec in ecs)


def compute_longtail(labels, threshold=0.05):
    """True iff positive interface ratio among valid labels is below threshold."""
    valid = [x for x in labels if x >= 0]
    if not valid:
        return False
    pos_ratio = sum(1 for x in valid if x == 1) / len(valid)
    return pos_ratio < threshold


def extract_sequence_from_pdb(pdb_path):
    """从 PDB 文件中提取氨基酸序列。"""
    seen = set()
    seq = []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            if len(line) < 28:
                continue
            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue
            res_name = line[17:20].strip()
            chain_id = line[22].strip()
            res_seq = line[23:27].strip()
            ins_code = line[27].strip() if len(line) > 27 else ""
            res_id = (chain_id, res_seq, ins_code)
            if res_id in seen:
                continue
            seen.add(res_id)
            aa = RESIDUE_3TO1.get(res_name, "X")
            seq.append(aa)
    return "".join(seq)


def load_labels(label_path):
    data = np.load(label_path)
    return data["labels"].astype(int).tolist()


def align_labels_to_pdb(metadata_seq, labels, pdb_seq, min_identity=0.85):
    """
    使用全局序列对齐将基于 metadata 序列的 labels 映射到 PDB 序列。

    返回：
        aligned_labels: 长度等于 len(pdb_seq) 的列表
        identity: 对齐 identity
    """
    if not metadata_seq or not pdb_seq:
        return None, 0.0

    alignments = ALIGNER.align(metadata_seq, pdb_seq)
    alignment = next(iter(alignments), None)
    if alignment is None:
        return None, 0.0
    aligned_meta, aligned_pdb = alignment[0], alignment[1]

    matches = sum(1 for a, b in zip(aligned_meta, aligned_pdb) if a == b and a != "-")
    aligned_len = sum(1 for a, b in zip(aligned_meta, aligned_pdb) if a != "-" or b != "-")
    identity = matches / aligned_len if aligned_len > 0 else 0.0

    coverage = min(len(metadata_seq), len(pdb_seq)) / max(len(metadata_seq), len(pdb_seq))
    if identity < min_identity and not (identity >= 0.75 and coverage >= 0.75):
        return None, identity

    # 将 labels 映射到 pdb 序列位置
    pdb_labels = []
    label_idx = 0
    for a, b in zip(aligned_meta, aligned_pdb):
        if b == "-":
            # PDB 中缺失，跳过 metadata 的 label
            if a != "-":
                label_idx += 1
            continue
        if a == "-":
            # PDB 中插入，设为 -1 表示标注缺失
            pdb_labels.append(-1)
        else:
            # 对齐位置，传递 label
            if label_idx < len(labels):
                pdb_labels.append(int(labels[label_idx]))
            else:
                pdb_labels.append(0)
            label_idx += 1

    if len(pdb_labels) != len(pdb_seq):
        return None, identity

    return pdb_labels, identity


def build_rows(
    records,
    metadata_by_hash,
    struct_dir,
    label_dir,
    include_ood=False,
    pair_meta=None,
    chain_domains=None,
    train_folds=None,
    train_superfamilies=None,
    ec_map=None,
    trainval_ec=None,
    longtail_threshold=0.05,
):
    rows = []
    skipped = 0
    aligned = 0
    direct = 0
    missing_struct = 0
    missing_label = 0

    for r in records:
        seq_hash = r["seq_hash"]
        pdb_path = struct_dir / f"{seq_hash}.pdb"
        label_path = label_dir / f"{seq_hash}_labels.npz"

        if not pdb_path.exists():
            missing_struct += 1
            logger.warning("Missing structure: %s", pdb_path)
            continue
        if not label_path.exists():
            missing_label += 1
            logger.warning("Missing labels: %s", label_path)
            continue

        aa_seq = extract_sequence_from_pdb(pdb_path)
        labels = load_labels(label_path)
        meta = metadata_by_hash.get(seq_hash, {})
        meta_seq = meta.get("sequence", "")

        if len(aa_seq) == len(labels) and len(aa_seq) == len(meta_seq):
            # 完全一致
            final_labels = labels
            direct += 1
        else:
            # 尝试对齐
            aligned_labels, identity = align_labels_to_pdb(meta_seq, labels, aa_seq)
            if aligned_labels is None:
                skipped += 1
                logger.warning(
                    "Skipping %s: alignment identity %.3f too low or failed (meta=%d pdb=%d)",
                    seq_hash, identity, len(meta_seq), len(aa_seq),
                )
                continue
            final_labels = aligned_labels
            aligned += 1

        if len(aa_seq) != len(final_labels):
            logger.warning(
                "Length mismatch after alignment for %s: seq=%d labels=%d",
                seq_hash, len(aa_seq), len(final_labels),
            )
            skipped += 1
            continue

        row = {
            "unique_id": seq_hash,
            "aa_seq": aa_seq,
            "file_type": "pdb",
            "labels": json.dumps(final_labels, separators=(",", ":")),
        }

        if include_ood:
            length = int(meta.get("length", len(aa_seq)))
            is_default = 60 <= length <= 1000
            is_long = length > 1000
            is_short = length < 60

            row["Default"] = str(is_default)
            row["ExtremeLong"] = str(is_long)
            row["ExtremeShort"] = str(is_short)

            if is_default:
                # 常规 OOD 列直接复用原计算结果
                for col in [
                    "OOD_LowHomology_90", "OOD_LowHomology_80", "OOD_LowHomology_70",
                    "OOD_LowHomology_60", "OOD_LowHomology_50", "OOD_LowHomology_40",
                    "OOD_LowHomology_30", "OOD_Orphan",
                    "OOD_NewFold_0.9", "OOD_NewFold_0.8", "OOD_NewFold_0.7",
                    "OOD_NewFold_0.6", "OOD_NewFold_0.5", "OOD_NewFold_0.4",
                    "OOD_NewFold_0.3", "OOD_IDR",
                ]:
                    row[col] = r.get(col, "False")

                # CATH Fold / Superfamily holdout
                fold_ood, sf_ood = compute_cath_ood(
                    seq_hash, chain_domains, pair_meta, metadata_by_hash,
                    train_folds, train_superfamilies
                )
                row["OOD_CATH_FoldHoldout"] = str(fold_ood)
                row["OOD_CATH_SuperfamilyHoldout"] = str(sf_ood)

                # NewFunction：所有 EC 类别均未见于 train+val
                nf = compute_new_function(
                    seq_hash, ec_map, pair_meta, metadata_by_hash, trainval_ec
                )
                row["OOD_NewFunction"] = str(nf)

                # LongTail：界面正样本比例低于阈值
                lt = compute_longtail(final_labels, threshold=longtail_threshold)
                row["OOD_LongTail"] = str(lt)
            else:
                # 极端长度样本的常规 OOD 列直接写 NaN（空字符串，"未计算"语义）
                for col in [
                    "OOD_LowHomology_90", "OOD_LowHomology_80", "OOD_LowHomology_70",
                    "OOD_LowHomology_60", "OOD_LowHomology_50", "OOD_LowHomology_40",
                    "OOD_LowHomology_30", "OOD_Orphan",
                    "OOD_NewFold_0.9", "OOD_NewFold_0.8", "OOD_NewFold_0.7",
                    "OOD_NewFold_0.6", "OOD_NewFold_0.5", "OOD_NewFold_0.4",
                    "OOD_NewFold_0.3", "OOD_IDR",
                    "OOD_CATH_FoldHoldout", "OOD_CATH_SuperfamilyHoldout",
                    "OOD_NewFunction", "OOD_LongTail",
                ]:
                    row[col] = ""

        rows.append(row)

    logger.info(
        "Processed %d rows (direct=%d, aligned=%d, skipped=%d), missing_struct=%d, missing_label=%d",
        len(rows), direct, aligned, skipped, missing_struct, missing_label,
    )
    return rows


def write_csv(rows, path, include_ood=False):
    base_fields = ["unique_id", "aa_seq", "file_type", "labels"]
    fields = base_fields + (OOD_COLUMNS if include_ood else [])
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    logger.info("Wrote %d rows to %s", len(rows), path)


def main():
    parser = argparse.ArgumentParser(description="Build final PPIS split files")
    parser.add_argument("--train", default="output/chain_centric_splits/chain_centric_train.csv")
    parser.add_argument("--val", default="output/chain_centric_splits/chain_centric_val.csv")
    parser.add_argument("--test", default="output/chain_centric_splits/chain_centric_test_ood.csv")
    parser.add_argument("--test_extreme", default="output/chain_centric_splits/chain_centric_test_extremelength.csv")
    parser.add_argument("--metadata", default="output/chain_centric_dedup_metadata.csv")
    parser.add_argument("--pair_metadata", default="output/dips_plus_metadata.csv")
    parser.add_argument("--cath_list", default=str(CATH_DOMAIN_LIST))
    parser.add_argument("--ec_tsv", default=str(EC_TSV_PATH))
    parser.add_argument("--longtail_threshold", type=float, default=0.05,
                        help="Interface positive ratio threshold for OOD_LongTail (default: 0.05)")
    parser.add_argument("--out_dir", default="splits")
    parser.add_argument("--struct_dir", default="pdbs")
    parser.add_argument("--label_dir", default="output/chain_interface_labels")
    args = parser.parse_args()

    root = Path(__file__).resolve().parent.parent
    struct_dir = (root / args.struct_dir).resolve()
    label_dir = (root / args.label_dir).resolve()
    out_dir = (root / args.out_dir).resolve()

    train_records = load_csv((root / args.train).resolve())
    val_records = load_csv((root / args.val).resolve())
    test_records = load_csv((root / args.test).resolve())
    test_extreme_records = load_csv((root / args.test_extreme).resolve())

    metadata = load_csv((root / args.metadata).resolve())
    metadata_by_hash = {r["seq_hash"]: r for r in metadata}

    pair_metadata = load_csv((root / args.pair_metadata).resolve())
    pair_meta = {r["unique_id"]: r for r in pair_metadata}

    chain_domains = load_cath_domain_list(Path(args.cath_list))
    train_folds, train_superfamilies = build_train_cath_sets(
        train_records, chain_domains, pair_meta, metadata_by_hash
    )
    logger.info(
        "Loaded CATH domains: %d chains; train Fold IDs=%d Superfamily IDs=%d",
        len(chain_domains), len(train_folds), len(train_superfamilies),
    )

    ec_map = load_ec_annotations(Path(args.ec_tsv))
    trainval_ec = build_trainval_ec_set(train_records, val_records, ec_map, pair_meta, metadata_by_hash)
    logger.info(
        "Loaded EC annotations: %d chains; train+val unique EC numbers=%d",
        len(ec_map), len(trainval_ec),
    )

    logger.info(
        "Loaded train=%d val=%d test=%d test_extreme=%d",
        len(train_records), len(val_records), len(test_records), len(test_extreme_records),
    )

    train_rows = build_rows(train_records, metadata_by_hash, struct_dir, label_dir, include_ood=False)
    val_rows = build_rows(val_records, metadata_by_hash, struct_dir, label_dir, include_ood=False)
    test_rows = build_rows(
        test_records, metadata_by_hash, struct_dir, label_dir, include_ood=True,
        pair_meta=pair_meta, chain_domains=chain_domains,
        train_folds=train_folds, train_superfamilies=train_superfamilies,
        ec_map=ec_map, trainval_ec=trainval_ec,
        longtail_threshold=args.longtail_threshold,
    )
    test_extreme_rows = build_rows(
        test_extreme_records, metadata_by_hash, struct_dir, label_dir, include_ood=True,
        pair_meta=pair_meta, chain_domains=chain_domains,
        train_folds=train_folds, train_superfamilies=train_superfamilies,
        ec_map=ec_map, trainval_ec=trainval_ec,
        longtail_threshold=args.longtail_threshold,
    )

    write_csv(train_rows, out_dir / "ppis_train.csv", include_ood=False)
    write_csv(val_rows, out_dir / "ppis_val.csv", include_ood=False)
    write_csv(test_rows + test_extreme_rows, out_dir / "ppis_test.csv", include_ood=True)

    logger.info(
        "Final splits: train=%d val=%d test=%d (normal=%d + extreme=%d)",
        len(train_rows), len(val_rows), len(test_rows) + len(test_extreme_rows),
        len(test_rows), len(test_extreme_rows),
    )


if __name__ == "__main__":
    main()
