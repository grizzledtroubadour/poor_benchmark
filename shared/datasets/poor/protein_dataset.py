"""
Poor benchmark dataset reader for KPLM training pipeline.

适配 poor_benchmark 数据格式：
  - CSV 列名 struct_file（仅文件名）而非 pdb_path（完整路径），需通过 pdb_dir 拼接
  - 结构文件支持 .pdb 和 .cif 两种格式
  - EC 多标签分类的 label 为分号分隔的 EC 编号字符串（如 "3.2.1.-;3.1.7.-"），
    需通过 ec_vocab.txt 映射为 multi-hot 向量
"""

import os
import ast
import numpy as np
import torch
import torch.nn.functional as F
import pandas as pd
from torch.utils.data import Dataset
from sklearn.preprocessing import MultiLabelBinarizer

from shared.utils.common import pmap_multi
from shared.esm.sdk.api import ESMProtein
from shared.esm.utils.structure.protein_chain import ProteinChain
import biotite.structure as bs


# ---------------------------------------------------------------------------
# 结构文件加载
# ---------------------------------------------------------------------------

def _load_structure(pdb_dir, struct_file):
    """从 pdb_dir 目录加载 .pdb 或 .cif 结构文件，返回 ESMProtein 对象。

    默认读取全部多肽链：序列按链以 '|'（CHAIN_BREAK_STR）拼接，atom37
    坐标在链边界插入 NaN 行，与序列逐位对齐（官方多聚体格式）。单链文件
    行为与此前一致（无 '|'）。全链解析失败时回退到首链读取。
    """
    full_path = os.path.join(pdb_dir, struct_file)
    ext = os.path.splitext(struct_file)[1].lower()

    if ext == ".cif":
        # 优先使用移植自官方 esm 仓库的原生 mmCIF 解析（MmcifWrapper），
        # 不经过 PDB 格式转换，天然支持多字符 chain_id 和超宽 b_factor 列。
        from shared.esm.utils.structure.protein_chain import ProteinChain
        from shared.esm.utils.structure.protein_complex import ProteinComplex
        try:
            chains = list(ProteinChain.chain_iterable_from_mmcif(full_path))
            if len(chains) > 1:
                return ESMProtein.from_protein_complex(
                    ProteinComplex.from_chains(chains))
            return ESMProtein.from_protein_chain(chains[0])
        except Exception:
            pass  # 无 SEQRES 实体记录的预测结构 cif 等，走 biotite 回退路径
        # biotite 回退路径：从 ATOM 记录构建序列。chain_id 超 1 字符时
        # ProteinChain 内部转 PDB 会抛 BadStructureError，故逐链重映射为
        # 单字符 chain_id；b_factor 超 PDB 3 位整数列上限时截断。
        from biotite.structure.io.pdbx import CIFFile, get_structure
        cif_file = CIFFile.read(full_path)
        atom_array = get_structure(cif_file, model=1, extra_fields=["b_factor"])
        atom_array = atom_array[bs.filter_amino_acids(atom_array)]
        import string
        unique_cids = list(dict.fromkeys(atom_array.chain_id))
        remap = {cid: string.ascii_uppercase[i % 26]
                 for i, cid in enumerate(unique_cids)}
        chains = []
        for cid in unique_cids:
            sub = atom_array[atom_array.chain_id == cid].copy()
            sub.chain_id[:] = remap[cid]
            if hasattr(sub, "b_factor"):
                sub.b_factor = np.clip(sub.b_factor, -99.99, 999.99)
            chains.append(ProteinChain.from_atomarray(sub))
        if len(chains) > 1:
            return ESMProtein.from_protein_complex(
                ProteinComplex.from_chains(chains))
        return ESMProtein.from_protein_chain(chains[0])
    else:
        # .pdb 及其他格式：优先全链读取，失败回退首链
        try:
            return ESMProtein.from_pdb(full_path, chain_id="all")
        except Exception:
            return ESMProtein.from_pdb(full_path)


def _load_plddt(pdb_dir, struct_file):
    """读取结构文件全部链 CA 原子的 b_factor（AF2/ESMFold 结构中为 pLDDT）。

    供 SaProt 的 pLDDT 掩码使用（官方建议 AF2 结构中 pLDDT<70 的区域
    用 '#' 掩码）。多链按链拼接，链边界插入 100.0 占位（与 _load_structure
    的全链序列逐位等长；100.0 使 pLDDT<70 掩码不会误伤边界位）。
    读取失败或不含 b_factor 时返回 None。
    """
    full_path = os.path.join(pdb_dir, struct_file)
    ext = os.path.splitext(struct_file)[1].lower()
    try:
        if ext == ".cif":
            from biotite.structure.io.pdbx import CIFFile, get_structure
            arr = get_structure(CIFFile.read(full_path), model=1,
                                extra_fields=["b_factor"])
        else:
            from biotite.structure.io.pdb import PDBFile
            arr = PDBFile.read(full_path).get_structure(
                model=1, extra_fields=["b_factor"])
        arr = arr[arr.atom_name == "CA"]
        if arr.array_length() == 0:
            return None
        # 按链出现顺序拼接，链边界插入 100.0（与全链序列的 '|' 位对齐）
        chain_ids = list(dict.fromkeys(arr.chain_id))
        parts = []
        for cid in chain_ids:
            sub = arr[arr.chain_id == cid]
            parts.append(np.asarray(sub.b_factor, dtype=np.float32))
        sep = np.array([100.0], dtype=np.float32)
        joined = parts[0]
        for part in parts[1:]:
            joined = np.concatenate([joined, sep, part])
        return torch.tensor(joined, dtype=torch.float32)
    except Exception:
        return None


# ---------------------------------------------------------------------------
# ECFP4 指纹解析
# ---------------------------------------------------------------------------

def ecfp4_to_vector(ecfp4_str, n_bits=2048):
    """把 ECFP4 稀疏 bit 索引列表字符串（如 "[80, 294, 389]"）转成 0/1 向量。"""
    vec = torch.zeros(n_bits)
    try:
        indices = [int(i) for i in ast.literal_eval(str(ecfp4_str))]
        indices = [i for i in indices if 0 <= i < n_bits]
        if indices:
            vec[torch.tensor(indices, dtype=torch.long)] = 1.0
    except Exception:
        pass
    return vec


# ---------------------------------------------------------------------------
# 标签词表加载
# ---------------------------------------------------------------------------

def load_label_vocab(vocab_path):
    """加载 EC 词表文件，返回 {EC编号字符串: 整数索引} 映射。

    词表文件格式：每行 ``index\\tEC_number``，例如 ``0\\t1.1.1.-``
    """
    vocab = {}
    with open(vocab_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) == 2:
                idx, ec = int(parts[0]), parts[1]
                vocab[ec] = idx
    return vocab


# ---------------------------------------------------------------------------
# 标签解析（顶层函数，供 read_data_poor 与 PPI embedding bank 复用）
# ---------------------------------------------------------------------------

def parse_label_poor(label, task_type, num_classes, label_vocab):
    """解析单条样本的标签，无法解析时返回 None。"""
    if task_type == "multi_labels_classification":
        if label_vocab is not None:
            # poor_benchmark: 分号分隔的 EC 编号字符串 → multi-hot
            ec_labels = str(label).split(";")
            indices = [label_vocab[ec] for ec in ec_labels if ec in label_vocab]
            label = torch.zeros(int(num_classes))
            if indices:
                label[indices] = 1.0
        else:
            # 回退到原始逻辑：逗号分隔的整数索引
            mlb = MultiLabelBinarizer(classes=range(int(num_classes)))
            label = str(label)
            label = torch.tensor(
                mlb.fit_transform(
                    [[int(ele) for ele in label.split(",")]]
                ).flatten().tolist()
            )
    elif task_type == "contact":
        label = torch.load(label, weights_only=True)
    elif task_type in ("residual_classification", "residual_binary_classification"):
        label = torch.tensor(
            list(map(int, label.strip("[]").replace("\n", " ").split()))
        )
    elif task_type == "classification" and label_vocab is not None:
        # 字符串标签（如 CATH fold "1.10"）通过 label_vocab 映射为整数索引
        label_str = str(label)
        if label_str not in label_vocab:
            return None
        label = torch.tensor(label_vocab[label_str])
    elif task_type in ("regression", "pair_regression"):
        # 回归 label 必须走 float：pandas 会把连续 label 列读成 float64，
        # 先 int() 会把 4.31 静默截断成 4（poor_lba 96% 的样本会被改值）
        try:
            label = torch.tensor(float(label), dtype=torch.float32)
        except (TypeError, ValueError):
            return None
    else:
        try:
            label = int(label)
            label = torch.tensor(label)
        except Exception:
            try:
                label = float(label)
                label = torch.tensor(label)
            except Exception:
                return None
    return label


# ---------------------------------------------------------------------------
# 并行数据读取 worker（必须是顶层函数以便 pickle）
# ---------------------------------------------------------------------------

def read_data_poor(aa_seq, struct_file, label, unique_id,
                   task_type, num_classes, smiles, csv_name,
                   pdb_dir, label_vocab, sequence_only=False, ecfp4=None,
                   parse_structure=True):
    """读取单条数据，返回与 KPLM ProteinDataset 相同格式的 dict 或 None。"""
    try:
        if unique_id is None:
            unique_id = str(hash(aa_seq))

        # ---- 标签解析 ----
        label = parse_label_poor(label, task_type, num_classes, label_vocab)
        if label is None:
            return None

        # ---- 残基级任务：label 长度必须与序列一致（不依赖结构文件）----
        # sequence_only 或无结构文件的分支不会走到下方结构解析处的长度过滤，
        # 错位样本会在 __getitem__ 里被补 0（0 是合法类、参与 loss）或被
        # 静默截断，这里直接丢弃
        if task_type in ("contact", "residual_classification",
                         "residual_binary_classification"):
            if label.shape[0] != len(aa_seq):
                return None

        # ---- 结构文件解析 ----
        name = str(hash(struct_file)) if struct_file is not None else str(hash(aa_seq))

        # sequence_only 模式：跳过结构文件加载，直接使用 CSV 中的 aa_seq
        # 适用于 纯序列蛋白质语言模型，避免加载大量 PDB 文件
        if sequence_only or struct_file is None:
            return {
                "name": name,
                "seq": aa_seq,
                "X": None,
                "label": label,
                "unique_id": unique_id,
                "pdb_path": None,
                "smiles": smiles,
                "ecfp4": ecfp4,
                "struct_file": struct_file,
            }

        # parse_structure=False：模型只消费 pdb_path/struct_file（3Di/SST
        # 预计算缓存或自行重解析，如 saprot/prost_t5/protrek_struct/prosst/
        # mulan），跳过结构解析与 plddt 读取，数据集构建从几十分钟降到秒级。
        # 缓存未命中时模型侧回退 fallback 3Di（与结构缺失行为一致）。
        # 注意：残基级任务由调用方保证 parse_structure=True（加载期需要
        # 结构序列长度过滤错位样本）。
        if not parse_structure:
            if "|" in struct_file:
                full_path = "|".join(
                    os.path.join(pdb_dir, sf) for sf in struct_file.split("|"))
            else:
                full_path = os.path.join(pdb_dir, struct_file)
            return {
                "name": name,
                "seq": aa_seq,
                "X": None,
                "plddt": None,
                "label": label,
                "unique_id": unique_id,
                "pdb_path": full_path,
                "smiles": smiles,
                "ecfp4": ecfp4,
                "struct_file": struct_file,
            }

        if struct_file is not None:
            if "|" not in struct_file:
                structure = _load_structure(pdb_dir, struct_file)
                # 残基级任务：结构可能缺失残基，若结构序列或标签长度与 aa_seq
                # 不一致，embedding 与逐残基 label 会错位，直接丢弃该样本
                if task_type in ("contact", "residual_classification",
                                 "residual_binary_classification"):
                    if len(structure.sequence) != len(aa_seq) or label.shape[0] != len(aa_seq):
                        return None
                full_path = os.path.join(pdb_dir, struct_file)
                return {
                    "name": name,
                    "seq": aa_seq if "flip" in csv_name.lower() else structure.sequence,
                    "X": structure.coordinates,
                    "plddt": _load_plddt(pdb_dir, struct_file),
                    "label": label,
                    "unique_id": unique_id,
                    "pdb_path": full_path,
                    "smiles": smiles,
                    "ecfp4": ecfp4,
                    "struct_file": struct_file,
                }
            else:
                # PPI：多个结构文件用 | 分隔
                structures, sequences = [], []
                for _struct_file in struct_file.split("|"):
                    structure = _load_structure(pdb_dir, _struct_file)
                    structures.append(structure.coordinates)
                    sequences.append(structure.sequence)
                full_path = "|".join(
                    os.path.join(pdb_dir, sf) for sf in struct_file.split("|")
                )
                return {
                    "name": name,
                    "seq": "|".join(sequences),
                    "X": structures,
                    "label": label,
                    "unique_id": unique_id,
                    "pdb_path": full_path,
                    "smiles": smiles,
                    "ecfp4": ecfp4,
                    "struct_file": struct_file,
                }
        else:
            return {
                "name": name,
                "seq": aa_seq,
                "X": None,
                "label": label,
                "unique_id": unique_id,
                "pdb_path": None,
                "smiles": smiles,
                "ecfp4": ecfp4,
                "struct_file": struct_file,
            }
    except Exception:
        return None


def read_bank_sample(key, aa_seq, struct_file, pdb_dir, sequence_only=False,
                     parse_structure=True):
    """加载去重后的单条蛋白质（PPI embedding bank 用），返回 dict 或 None。"""
    try:
        if sequence_only or struct_file is None:
            return {
                "key": key,
                "name": str(hash(aa_seq)),
                "seq": aa_seq,
                "X": None,
                "label": torch.tensor(0),
                "unique_id": key,
                "pdb_path": None,
                "smiles": None,
                "ecfp4": None,
            }
        # 见 read_data_poor：只消费 pdb_path 的模型跳过结构解析
        if not parse_structure:
            return {
                "key": key,
                "name": str(hash(struct_file)),
                "seq": aa_seq,
                "X": None,
                "plddt": None,
                "label": torch.tensor(0),
                "unique_id": key,
                "pdb_path": os.path.join(pdb_dir, struct_file),
                "smiles": None,
                "ecfp4": None,
                "struct_file": struct_file,
            }
        structure = _load_structure(pdb_dir, struct_file)
        return {
            "key": key,
            "name": str(hash(struct_file)),
            "seq": structure.sequence,
            "X": structure.coordinates,
            "plddt": _load_plddt(pdb_dir, struct_file),
            "label": torch.tensor(0),
            "unique_id": key,
            "pdb_path": os.path.join(pdb_dir, struct_file),
            "smiles": None,
            "ecfp4": None,
        }
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Dataset 类
# ---------------------------------------------------------------------------

class PoorProteinDataset(Dataset):
    """Poor benchmark 数据集，与 KPLM 的 ProteinDataset 接口兼容。

    相比 ProteinDataset 的主要差异：
      - CSV 使用 ``struct_file`` 列（文件名），通过 ``pdb_dir`` 拼接完整路径
      - 支持 .pdb 和 .cif 结构文件
      - ``multi_labels_classification`` 的 label 为分号分隔的 EC 编号字符串，
        通过 ``label_vocab_path`` 指定的词表文件映射为 multi-hot 向量
      - ``parse_structure=False`` 时跳过结构文件解析（只保留 pdb_path），
        供只消费 3Di/SST 缓存的模型（saprot/prost_t5/protrek_struct/
        prosst/mulan）使用；消费坐标 X 的模型（esm3/esm_if1）与残基级
        任务必须保持 True
    """

    def __init__(self, csv_file, pretrain_model_name="esm2_650m",
                 max_length=1024, pretrain_model_interface=None,
                 task_name="pretrain", task_type="classification",
                 num_classes=None, pdb_dir=None, label_vocab_path=None,
                 sequence_only=False, parse_structure=True):
        self.max_length = max_length
        self.pretrain_model_name = pretrain_model_name
        self.task_name = task_name
        self.task_type = task_type
        self.num_classes = num_classes
        self.embedding_bank = None

        # 加载标签词表（EC 编号 → 整数索引）
        self.label_vocab = None
        if label_vocab_path is not None:
            self.label_vocab = load_label_vocab(label_vocab_path)
            print(f"PoorProteinDataset: loaded {len(self.label_vocab)} labels from {label_vocab_path}")

        # 读取 CSV
        csv_data = pd.read_csv(csv_file)

        # 检测 PPI 双序列格式（aa_seq1/aa_seq2 + struct_file1/struct_file2）
        is_ppi = "aa_seq1" in csv_data.columns and "aa_seq2" in csv_data.columns

        # PPI + 预提取 embedding 模式：对唯一蛋白质去重提取 embedding 存入 bank，
        # 每条 pair 样本只保存两个引用，避免按 pair 存储拼接 embedding 导致内存爆炸
        if is_ppi and pretrain_model_interface is not None:
            self.pretrain_model_interface = pretrain_model_interface
            self._build_ppi_embedding_bank(csv_data, pdb_dir, sequence_only,
                                           parse_structure)
            print(f"PoorProteinDataset(PPI bank): {len(self.data)} pairs, "
                  f"{len(self.embedding_bank)} unique proteins from {csv_file}")
            return

        # SMILES 列名适配：kcat/lba/lbs 用 ligand_smiles，其他用 smiles
        smiles_col = None
        for col in ("smiles", "ligand_smiles"):
            if col in csv_data.columns:
                smiles_col = col
                break

        # ECFP4 预计算指纹列优先（稀疏 bit 索引），省去 RDKit 现算并规避解析失败
        ecfp4_col = "ligand_ecfp4" if "ligand_ecfp4" in csv_data.columns else None

        path_list = []
        for i in range(len(csv_data)):
            row = csv_data.iloc[i]
            if is_ppi:
                # PPI: 两个蛋白质序列/结构用 | 连接，与 read_data_poor 的 PPI 分支匹配
                aa_seq = str(row["aa_seq1"]) + "|" + str(row["aa_seq2"])
                sf1 = row.get("struct_file1")
                sf2 = row.get("struct_file2")
                if pd.notna(sf1) and pd.notna(sf2):
                    struct_file = str(sf1) + "|" + str(sf2)
                else:
                    struct_file = None
            else:
                aa_seq = row.get("aa_seq")
                struct_file = row.get("struct_file")
                if pd.isna(struct_file):
                    struct_file = None

            smiles_val = row.get(smiles_col) if smiles_col else None
            if smiles_val is None or pd.isna(smiles_val):
                smiles_val = None
            if ecfp4_col is not None and pd.notna(row.get(ecfp4_col)):
                ecfp4_val = ecfp4_to_vector(row[ecfp4_col])
            else:
                ecfp4_val = None

            path_list.append((
                aa_seq,
                struct_file,
                row["label"],
                None if pd.isna(row.get("unique_id")) else row.get("unique_id"),
                task_type,
                num_classes,
                smiles_val,
                csv_file,
                pdb_dir,
                self.label_vocab,
                sequence_only,
                ecfp4_val,
                parse_structure,
            ))

        # 默认上限 16 线程：8 个训练进程并发时 -1（全部物理核）会造成
        # 数千线程的上下文切换开销（PPI 卡死事故的根因），可用
        # POOR_PMAP_NJOBS 环境变量覆盖
        self.data = pmap_multi(read_data_poor, path_list,
                               n_jobs=int(os.environ.get(
                                   "POOR_PMAP_NJOBS", min(16, os.cpu_count() or 1))))
        self.data = [d for d in self.data if d is not None]
        n_dropped = len(path_list) - len(self.data)
        if n_dropped > 0:
            print(f"PoorProteinDataset: dropped {n_dropped}/{len(path_list)} samples "
                  f"(parse failure or residual-task length mismatch)")
        # Preserve identifiers before inference_datasets replaces the dicts.
        self._unique_ids = [d.get("unique_id") for d in self.data]
        self.max_length = min(self.max_length, max([len(d["seq"]) for d in self.data]) + 2)
        self.pretrain_model_interface = pretrain_model_interface

        if pretrain_model_interface is not None:
            self.data = pretrain_model_interface.inference_datasets(
                self.data, task_name=self.task_name
            )
            # Restore identifiers so downstream evaluation can align with CSV rows.
            for d, uid in zip(self.data, self._unique_ids):
                d["unique_id"] = uid

        print(f"PoorProteinDataset: {len(self.data)} samples loaded from {csv_file}")

    def __len__(self):
        return len(self.data)

    def _build_ppi_embedding_bank(self, csv_data, pdb_dir, sequence_only,
                                  parse_structure=True):
        """PPI 去重 embedding bank：每条唯一蛋白质只提取一次 embedding，
        每个 pair 样本只保存两个引用（key_A/key_B），在 __getitem__ 时即时拼接。"""
        bank_args = {}
        pair_rows = []
        for i in range(len(csv_data)):
            row = csv_data.iloc[i]
            seq1, seq2 = str(row["aa_seq1"]), str(row["aa_seq2"])
            sf1, sf2 = row.get("struct_file1"), row.get("struct_file2")
            sf1 = str(sf1) if pd.notna(sf1) else None
            sf2 = str(sf2) if pd.notna(sf2) else None

            try:
                label = parse_label_poor(row["label"], self.task_type,
                                         self.num_classes, self.label_vocab)
            except Exception:
                continue
            if label is None:
                continue

            # 结构模型按 struct_file 去重（结构决定坐标），纯序列模型按序列去重
            if sequence_only:
                key1, key2 = seq1, seq2
            else:
                key1 = sf1 if sf1 is not None else seq1
                key2 = sf2 if sf2 is not None else seq2
            if key1 not in bank_args:
                bank_args[key1] = (key1, seq1, sf1, pdb_dir, sequence_only,
                                   parse_structure)
            if key2 not in bank_args:
                bank_args[key2] = (key2, seq2, sf2, pdb_dir, sequence_only,
                                   parse_structure)

            pair_rows.append({
                "name": f"{key1}|{key2}",
                "label": label,
                "key_A": key1,
                "key_B": key2,
                "unique_id": None if pd.isna(row.get("unique_id")) else row.get("unique_id"),
            })

        bank_samples = pmap_multi(read_bank_sample, list(bank_args.values()),
                                  n_jobs=int(os.environ.get(
                                      "POOR_PMAP_NJOBS", min(16, os.cpu_count() or 1))))
        bank_samples = [s for s in bank_samples if s is not None]
        valid_keys = {s["key"] for s in bank_samples}
        self.data = [d for d in pair_rows
                     if d["key_A"] in valid_keys and d["key_B"] in valid_keys]

        results = self.pretrain_model_interface.inference_datasets(
            bank_samples, task_name=self.task_name
        )
        self.embedding_bank = {
            s["key"]: {"embedding": r["embedding"],
                       "attention_mask": r["attention_mask"]}
            for s, r in zip(bank_samples, results)
        }

        # 每条蛋白的 embedding 分别裁剪到 max_length（config seq_len），再拼接；
        # 不做拼接后的整体裁剪，避免 B 蛋白尾部被丢
        per_protein_max = self.max_length
        for v in self.embedding_bank.values():
            v["embedding"] = v["embedding"][:per_protein_max]
            v["attention_mask"] = v["attention_mask"][:per_protein_max]

        # max_length 取裁剪后的真实最大拼接长度（无分隔位，≤ 2*per_protein_max），
        # __getitem__ 的 pad_data 只做 padding，不再截断
        if self.data:
            emb_lens = {k: v["embedding"].shape[0]
                        for k, v in self.embedding_bank.items()}
            self.max_length = max(
                emb_lens[d["key_A"]] + emb_lens[d["key_B"]] for d in self.data
            )

    def pad_data(self, data, dim=0, pad_value=0, max_length=1024):
        if data.shape[dim] < max_length:
            data = dynamic_pad(data, [0, max_length - data.shape[dim]],
                                dim=dim, pad_value=pad_value)
        else:
            start = 0
            data = data[start:start + max_length]
        return data

    def __getitem__(self, idx):
        if self.embedding_bank is not None:
            # PPI bank 模式：A/B embedding 直接拼接（无分隔位），再 pad 到 max_length
            max_length_batch = self.max_length
            d = self.data[idx]
            entry_A = self.embedding_bank[d["key_A"]]
            entry_B = self.embedding_bank[d["key_B"]]
            emb_A, emb_B = entry_A["embedding"], entry_B["embedding"]
            mask_A, mask_B = entry_A["attention_mask"], entry_B["attention_mask"]

            embedding = torch.cat([emb_A, emb_B], dim=0)
            attention_mask = torch.cat([mask_A, mask_B], dim=0) == 1

            embedding = self.pad_data(embedding, dim=0, pad_value=0,
                                      max_length=max_length_batch)
            attention_mask = self.pad_data(attention_mask, dim=0, pad_value=0,
                                           max_length=max_length_batch)
            label = d["label"]
            if self.task_type == "binary_classification":
                label = label[None].float()

            return {
                "name": d["name"],
                "embedding": embedding,
                "attention_mask": attention_mask,
                "label": label,
            }

        if self.pretrain_model_interface is not None:
            max_length_batch = self.max_length
            name = self.data[idx]["name"]
            embedding = self.pad_data(self.data[idx]["embedding"], dim=0,
                                     pad_value=0, max_length=max_length_batch)
            attention_mask = self.pad_data(self.data[idx]["attention_mask"], dim=0,
                                           pad_value=0, max_length=max_length_batch)
            label = self.data[idx]["label"]

            if self.task_type == "binary_classification":
                label = label[None].float()
            if self.task_type == "contact":
                label = (label == 0).int()
                label = F.pad(label, [0, max_length_batch - label.shape[0],
                                      0, max_length_batch - label.shape[0]])
            if self.task_type in ("residual_classification", "residual_binary_classification"):
                label = F.pad(label, [0, max_length_batch - label.shape[0]])

            result = {
                "name": name,
                "embedding": embedding,
                "attention_mask": attention_mask,
                "label": label,
            }

            if self.data[idx].get("smiles") is not None:
                result["smiles"] = self.data[idx]["smiles"]
            if self.data[idx].get("ecfp4") is not None:
                result["ecfp4"] = self.data[idx]["ecfp4"]

            return result
        else:
            max_length_batch = self.max_length
            label = self.data[idx]["label"]
            if self.task_type == "binary_classification":
                label = label[None].float()
            if self.task_type == "contact":
                label = (label == 0).int()
                label = F.pad(label, [0, max_length_batch - label.shape[0],
                                      0, max_length_batch - label.shape[0]])
            if self.task_type in ("residual_classification", "residual_binary_classification"):
                label = F.pad(label, [0, max_length_batch - label.shape[0]])

            data = {
                "name": self.data[idx]["name"],
                "seq": self.data[idx]["seq"],
                "X": self.data[idx]["X"],
                "label": label,
                "unique_id": self.data[idx]["unique_id"],
                "pdb_path": self.data[idx]["pdb_path"],
                "smiles": self.data[idx]["smiles"],
                "ecfp4": self.data[idx].get("ecfp4"),
            }
            return data


def dynamic_pad(tensor, pad_size, dim=0, pad_value=0):
    shape = list(tensor.shape)
    num_dims = len(shape)

    pad = [0] * (2 * num_dims)
    prev_pad_size, post_pad_size = pad_size
    pad_index = 2 * (num_dims - dim - 1)
    pad[pad_index] = prev_pad_size
    pad[pad_index + 1] = post_pad_size

    padded_tensor = F.pad(tensor, pad, mode="constant", value=pad_value)
    return padded_tensor
