"""
Poor benchmark dataset reader for KPLM training pipeline.

Support for the poor_benchmark data format:
  - CSVs use struct_file (filename), not pdb_path (full path); resolve it relative to pdb_dir.
  - Both .pdb and .cif structures are supported.
  - EC multilabel targets are semicolon-separated EC numbers (e.g. "3.2.1.-;3.1.7.-").
    Map them to multi-hot vectors using ec_vocab.txt.
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
# Structure loading
# ---------------------------------------------------------------------------

def _load_structure(pdb_dir, struct_file):
    """Load a .pdb or .cif structure from pdb_dir and return an ESMProtein.

    Read all polypeptide chains by default. Join sequences with '|' (CHAIN_BREAK_STR),
    inserting NaN atom37 coordinate rows at chain boundaries to align with the sequence
    (the official multimer format). Single-chain behavior is unchanged (no '|'). Fall back to the first chain if all-chain parsing fails.
    """
    full_path = os.path.join(pdb_dir, struct_file)
    ext = os.path.splitext(struct_file)[1].lower()

    if ext == ".cif":
        # Prefer native mmCIF parsing (MmcifWrapper) ported from the official esm repository.
        # Avoid PDB conversion to support multicharacter chain IDs and wide b_factor fields.
        from shared.esm.utils.structure.protein_chain import ProteinChain
        from shared.esm.utils.structure.protein_complex import ProteinComplex
        try:
            chains = list(ProteinChain.chain_iterable_from_mmcif(full_path))
            if len(chains) > 1:
                return ESMProtein.from_protein_complex(
                    ProteinComplex.from_chains(chains))
            return ESMProtein.from_protein_chain(chains[0])
        except Exception:
            pass  # Predicted CIFs without SEQRES entities use the biotite fallback.
        # biotite fallback: build the sequence from ATOM records. Multicharacter chain IDs
        # trigger BadStructureError during internal PDB conversion in ProteinChain, so remap
        # each chain to a single character and clip b_factor values to the PDB field limit.
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
        # .pdb and other formats: try all chains first, then fall back to the first chain.
        try:
            return ESMProtein.from_pdb(full_path, chain_id="all")
        except Exception:
            return ESMProtein.from_pdb(full_path)


def _load_plddt(pdb_dir, struct_file):
    """Read CA b_factor values across all chains (pLDDT for AF2/ESMFold structures).

    Used for SaProt pLDDT masking (the official recommendation masks AF2 regions
    with pLDDT <70 using '#'). Concatenate chains with a 100.0 boundary placeholder,
    matching the all-chain sequence in _load_structure without masking chain boundaries.
    Return None if reading fails or b_factor is unavailable.
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
        # Concatenate in chain order, inserting 100.0 at each '|' sequence boundary.
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
# ECFP4 fingerprint parsing
# ---------------------------------------------------------------------------

def ecfp4_to_vector(ecfp4_str, n_bits=2048):
    """Convert a sparse ECFP4 bit-index string (e.g. "[80, 294, 389]") to a binary vector."""
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
# Label vocabulary loading
# ---------------------------------------------------------------------------

def load_label_vocab(vocab_path):
    """Load the EC vocabulary as a mapping from EC number strings to integer indices.

    File format: one ``index\\tEC_number`` per line, e.g. ``0\\t1.1.1.-``.
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
# Label parsing (top-level function shared by read_data_poor and the PPI embedding bank)
# ---------------------------------------------------------------------------

def parse_label_poor(label, task_type, num_classes, label_vocab):
    """Parse one sample label; return None if parsing fails."""
    if task_type == "multi_labels_classification":
        if label_vocab is not None:
            # poor_benchmark: semicolon-separated EC numbers -> multi-hot vector.
            ec_labels = str(label).split(";")
            indices = [label_vocab[ec] for ec in ec_labels if ec in label_vocab]
            label = torch.zeros(int(num_classes))
            if indices:
                label[indices] = 1.0
        else:
            # Fall back to the original comma-separated integer-index format.
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
        # Map string labels (e.g. CATH fold "1.10") to integer indices with label_vocab.
        label_str = str(label)
        if label_str not in label_vocab:
            return None
        label = torch.tensor(label_vocab[label_str])
    elif task_type in ("regression", "pair_regression"):
        # Regression labels must use float: pandas reads continuous labels as float64.
        # Converting to int first silently truncates 4.31 to 4 (affecting 96% of poor_lba).
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
# Parallel data-loading worker (must be top-level to support pickling)
# ---------------------------------------------------------------------------

def read_data_poor(aa_seq, struct_file, label, unique_id,
                   task_type, num_classes, smiles, csv_name,
                   pdb_dir, label_vocab, sequence_only=False, ecfp4=None,
                   parse_structure=True):
    """Read one sample and return a KPLM ProteinDataset-compatible dict or None."""
    try:
        if unique_id is None:
            unique_id = str(hash(aa_seq))

        # ---- Label parsing ----
        label = parse_label_poor(label, task_type, num_classes, label_vocab)
        if label is None:
            return None

        # ---- Residue tasks: require label/sequence length agreement without structure parsing. ----
        # sequence_only and missing-structure branches bypass the structural length filter below.
        # Otherwise __getitem__ may pad misaligned labels with valid class 0 (included in loss)
        # or silently truncate them. Reject these samples here.
        if task_type in ("contact", "residual_classification",
                         "residual_binary_classification"):
            if label.shape[0] != len(aa_seq):
                return None

        # ---- Structure parsing ----
        name = str(hash(struct_file)) if struct_file is not None else str(hash(aa_seq))

        # sequence_only: skip structure loading and use aa_seq from the CSV.
        # Intended for sequence-only protein language models; avoids loading many PDB files.
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

        # parse_structure=False: models consume only pdb_path/struct_file via cached 3Di/SST
        # or their own parsers (saprot/prost_t5/protrek_struct/prosst/mulan).
        # Skipping structure and pLDDT parsing reduces dataset construction from minutes to seconds.
        # On cache misses the model uses fallback 3Di, as for missing structures.
        # Callers must set parse_structure=True for residue-level tasks to filter samples
        # whose structure sequence lengths do not align with their labels.
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
                # Structures may omit residues. If structure sequence or label lengths differ
                # from aa_seq, embeddings and residue labels are misaligned; reject the sample.
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
                # PPI: multiple structure filenames are separated by |.
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
    """Load a deduplicated protein for the PPI embedding bank; return a dict or None."""
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
        # As in read_data_poor: models consuming only pdb_path skip structure parsing.
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
# Dataset class
# ---------------------------------------------------------------------------

class PoorProteinDataset(Dataset):
    """Poor benchmark dataset, compatible with the KPLM ProteinDataset interface.

    Main differences from ProteinDataset:
      - CSVs use ``struct_file`` filenames, resolved relative to ``pdb_dir``.
      - Both .pdb and .cif structures are supported.
      - ``multi_labels_classification`` labels contain semicolon-separated EC numbers,
        mapped to multi-hot vectors using the vocabulary at ``label_vocab_path``.
      - ``parse_structure=False`` skips parsing and retains only pdb_path for models
        consuming cached 3Di/SST (saprot/prost_t5/protrek_struct/prosst/mulan).
        Models consuming coordinates X (esm3/esm_if1) and residue-level tasks
        must keep this option True.
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

        # Load the label vocabulary (EC numbers -> integer indices).
        self.label_vocab = None
        if label_vocab_path is not None:
            self.label_vocab = load_label_vocab(label_vocab_path)
            print(f"PoorProteinDataset: loaded {len(self.label_vocab)} labels from {label_vocab_path}")

        # Read the CSV.
        csv_data = pd.read_csv(csv_file)

        # Detect PPI paired-sequence columns (aa_seq1/aa_seq2 + struct_file1/struct_file2).
        is_ppi = "aa_seq1" in csv_data.columns and "aa_seq2" in csv_data.columns

        # PPI with precomputed embeddings: deduplicate proteins and store embeddings in a bank.
        # Each pair stores two references, avoiding memory blowup from duplicated pair embeddings.
        if is_ppi and pretrain_model_interface is not None:
            self.pretrain_model_interface = pretrain_model_interface
            self._build_ppi_embedding_bank(csv_data, pdb_dir, sequence_only,
                                           parse_structure)
            print(f"PoorProteinDataset(PPI bank): {len(self.data)} pairs, "
                  f"{len(self.embedding_bank)} unique proteins from {csv_file}")
            return

        # SMILES column: kcat/lba/lbs use ligand_smiles; other tasks use smiles.
        smiles_col = None
        for col in ("smiles", "ligand_smiles"):
            if col in csv_data.columns:
                smiles_col = col
                break

        # Prefer precomputed ECFP4 sparse bit indices to avoid RDKit work and parsing failures.
        ecfp4_col = "ligand_ecfp4" if "ligand_ecfp4" in csv_data.columns else None

        path_list = []
        for i in range(len(csv_data)):
            row = csv_data.iloc[i]
            if is_ppi:
                # PPI: join the two sequences/structures with |, matching read_data_poor.
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

        # Default to at most 16 workers: -1 uses all physical cores and eight concurrent
        # training processes cause excessive thread switching (the cause of PPI hangs).
        # Override with the POOR_PMAP_NJOBS environment variable.
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
        """PPI embedding bank: embed each unique protein once; each pair stores two references
        (key_A/key_B), concatenated on demand in __getitem__."""
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

            # Structure models deduplicate by struct_file; sequence-only models by sequence.
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

        # Clip each protein embedding to max_length (config seq_len), then concatenate.
        # Do not clip the combined pair, which would discard the tail of protein B.
        per_protein_max = self.max_length
        for v in self.embedding_bank.values():
            v["embedding"] = v["embedding"][:per_protein_max]
            v["attention_mask"] = v["attention_mask"][:per_protein_max]

        # Use the actual maximum clipped pair length (no separator; <=2*per_protein_max).
        # __getitem__.pad_data only pads; it no longer truncates.
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
            # PPI bank: concatenate A/B embeddings without a separator, then pad to max_length.
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
