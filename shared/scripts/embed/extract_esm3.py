"""
ESM3 protein-chain embedding extractor.

This script loads the locally cached ESM3 small-open weights and extracts a
per-residue / mean-pooled embedding from a protein sequence and (optionally)
structure information.  Structure can be passed as a PDB file path or as an
atom37 coordinate array.

Usage:
    # Run from the project root (or anywhere; the script resolves its own paths).
    python shared/scripts/embed/extract_esm3.py

Environment:
    Requires torch, numpy, and a biotite version compatible with `shared/esm`.
    The script includes a small compatibility shim so it also works with
    biotite versions that no longer ship `biotite.structure.io.npz`.

The standalone smoke test will pick a random PDB from
`data/POOR_benchmark/func_prediction/pdbs` and print the resulting embedding
shapes.
"""

from __future__ import annotations

import os
import random
import sys
import warnings
from pathlib import Path
from typing import Sequence

import numpy as np
import torch

# Make `shared/esm` importable as the top-level package `esm` regardless of
# the current working directory.
_PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(_PROJECT_ROOT / "shared"))

# Compatibility shim: some `biotite` versions do not provide
# `biotite.structure.io.npz`. The class is only used for optional NPZ export,
# so a no-op stub is sufficient for embedding extraction.
import types as _types  # noqa: E402

try:
    from biotite.structure.io.npz import NpzFile as _  # noqa: F401
except Exception:
    _npz_stub = _types.ModuleType("biotite.structure.io.npz")

    class _NpzFile:
        def set_structure(self, *args, **kwargs):
            pass

        def write(self, *args, **kwargs):
            pass

    _npz_stub.NpzFile = _NpzFile
    sys.modules["biotite.structure.io.npz"] = _npz_stub

from shared.esm.models.esm3 import ESM3
from shared.esm.models.function_decoder import FunctionTokenDecoder, FunctionTokenDecoderConfig
from shared.esm.models.vqvae import StructureTokenDecoder, StructureTokenEncoder
from shared.esm.sdk.api import ESMProtein, LogitsConfig
from shared.esm.tokenization import get_esm3_model_tokenizers
from shared.esm.utils import encoding
import shared.esm.utils.residue_constants as RC
from shared.esm.utils.sampling import _BatchedESMProteinTensor

from Bio import BiopythonWarning
from Bio.Data import PDBData
from Bio.PDB import MMCIFParser, PDBParser

warnings.filterwarnings("ignore", category=BiopythonWarning)

DEFAULT_MODEL_DIR = Path("/root/autodl-tmp/model_cache/esm3")


def _make_structure_encoder_fn(model_dir: Path):
    """Return a loader for the ESM3 structure encoder that reads from ``model_dir``."""

    def _load(device: torch.device | str):
        with torch.device(device):
            model = StructureTokenEncoder(
                d_model=1024, n_heads=1, v_heads=128, n_layers=2, d_out=128, n_codes=4096
            ).eval()
        state_dict = torch.load(
            model_dir / "esm3_structure_encoder_v0.pth",
            map_location=device,
            weights_only=True,
        )
        model.load_state_dict(state_dict)
        return model

    return _load


def _make_structure_decoder_fn(model_dir: Path):
    """Return a loader for the ESM3 structure decoder that reads from ``model_dir``."""

    def _load(device: torch.device | str):
        with torch.device(device):
            model = StructureTokenDecoder(d_model=1280, n_heads=20, n_layers=30).eval()
        state_dict = torch.load(
            model_dir / "esm3_structure_decoder_v0.pth",
            map_location=device,
            weights_only=True,
        )
        model.load_state_dict(state_dict)
        return model

    return _load


def _make_function_decoder_fn(model_dir: Path):
    """Return a loader for the ESM3 function decoder that reads from ``model_dir``."""

    def _load(device: torch.device | str):
        with torch.device(device):
            config = FunctionTokenDecoderConfig(
                interpro_entry_list=str(model_dir / "data" / "entry_list_safety_29026.list"),
                keyword_vocabulary_path=str(
                    model_dir / "data" / "keyword_vocabulary_safety_filtered_58641.txt"
                ),
            )
            model = FunctionTokenDecoder(config).eval()
        state_dict = torch.load(
            model_dir / "esm3_function_decoder_v0.pth",
            map_location=device,
            weights_only=True,
        )
        model.load_state_dict(state_dict)
        return model

    return _load


def load_esm3_sm_open_v1(
    model_dir: str | os.PathLike = DEFAULT_MODEL_DIR,
    device: torch.device | str | None = None,
) -> ESM3:
    """
    Load the ESM3 small-open model from a local weight directory.

    Args:
        model_dir: Directory containing the ``esm3_*.pth`` weight files.
        device: Torch device. Defaults to CUDA if available.

    Returns:
        An ESM3 model in evaluation mode.
    """
    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model_dir = Path(model_dir)
    if not model_dir.is_dir():
        raise FileNotFoundError(f"ESM3 model directory not found: {model_dir}")

    tokenizers = get_esm3_model_tokenizers()
    with torch.device(device):
        model = ESM3(
            d_model=1536,
            n_heads=24,
            v_heads=256,
            n_layers=48,
            structure_encoder_fn=_make_structure_encoder_fn(model_dir),
            structure_decoder_fn=_make_structure_decoder_fn(model_dir),
            function_decoder_fn=_make_function_decoder_fn(model_dir),
            tokenizers=tokenizers,
        ).eval()

    main_ckpt = model_dir / "esm3_sm_open_v1.pth"
    state_dict = torch.load(main_ckpt, map_location=device, weights_only=True)
    model.load_state_dict(state_dict)

    # Match the dtype behaviour of ESM3.from_pretrained on CUDA.
    if device.type != "cpu":
        model = model.to(torch.bfloat16)

    return model


def _residue_to_one_letter(
    residue,
    restype_3to1: dict[str, str],
) -> str | None:
    """Convert a Biopython residue to a one-letter code, or None if not an amino acid."""
    resname = residue.resname.strip().upper()
    return restype_3to1.get(resname)


def load_atom37_from_structure_file(
    path: str | os.PathLike,
    chain_id_hint: str | None = None,
) -> tuple[np.ndarray, str]:
    """
    Load a PDB/mmCIF file with Biopython and return atom37 coordinates
    (shape ``[L, 37, 3]``) together with the observed amino-acid sequence.

    Missing atoms are filled with NaN. HETATM / water / non-amino-acid residues
    are skipped.
    """
    path = Path(path)
    suffix = path.suffix.lower()

    if suffix in (".cif", ".mmcif"):
        parser = MMCIFParser(QUIET=True)
    elif suffix == ".pdb":
        parser = PDBParser(QUIET=True)
    else:
        raise ValueError(f"Unsupported structure file format: {suffix}")

    structure = parser.get_structure(path.stem, str(path))
    model = next(structure.get_models())

    chains = list(model.get_chains())
    if not chains:
        raise ValueError("No chains found in structure file")

    # Choose chain.
    chain = None
    if chain_id_hint and chain_id_hint in [c.id for c in chains]:
        chain = next(c for c in chains if c.id == chain_id_hint)
    else:
        # Pick the chain with the most amino-acid residues.
        def _chain_length(c):
            return sum(
                1 for r in c.get_residues() if not r.id[0].strip()
            )

        chain = max(chains, key=_chain_length)

    restype_3to1 = {v: k for k, v in RC.restype_1to3.items()}
    restype_3to1.update(PDBData.protein_letters_3to1)

    observed_residues = []
    observed_seq_chars = []
    for r in chain.get_residues():
        if r.id[0].strip():
            # Skip hetero/water residues.
            continue
        one_letter = _residue_to_one_letter(r, restype_3to1)
        if one_letter is None:
            continue
        observed_residues.append(r)
        observed_seq_chars.append(one_letter)

    if not observed_residues:
        raise ValueError(f"No amino-acid residues found in chain {chain.id!r}")

    L = len(observed_residues)
    coords = np.full((L, RC.atom_type_num, 3), np.nan, dtype=np.float32)

    for pos, residue in enumerate(observed_residues):
        for atom in residue.get_atoms():
            atom_name = atom.get_name().strip().upper()
            if atom_name in RC.atom_order:
                coords[pos, RC.atom_order[atom_name]] = atom.get_coord().astype(np.float32)

    observed_seq = "".join(observed_seq_chars)
    return coords, observed_seq


def _chain_id_hint(struct_file: str) -> str | None:
    """Extract chain id hint from filenames like ``7PLO_S.pdb`` or ``A.cif``."""
    base = Path(struct_file).stem
    if "_" in base:
        return base.rsplit("_", 1)[1]
    return None


class ESM3Embedder:
    """
    Convenience wrapper around ESM3 for extracting chain-level embeddings.

    The model is loaded once in ``__init__`` and can be reused for many calls.
    """

    def __init__(
        self,
        model_dir: str | os.PathLike = DEFAULT_MODEL_DIR,
        device: torch.device | str | None = None,
    ):
        self.device = (
            torch.device(device)
            if device is not None
            else torch.device("cuda" if torch.cuda.is_available() else "cpu")
        )
        self.model = load_esm3_sm_open_v1(model_dir, self.device)
        self.model.eval()

    @torch.no_grad()
    def extract(
        self,
        sequence: str | None = None,
        structure: str | os.PathLike | np.ndarray | torch.Tensor | None = None,
        apply_model_norm: bool = True,
    ) -> dict[str, torch.Tensor]:
        """
        Extract an embedding for one protein chain.

        Args:
            sequence: Amino-acid sequence (one-letter code). If ``structure`` is a
                PDB path and ``sequence`` is not given, the sequence is read from
                the PDB file.
            structure: One of the following:
                - ``None``: sequence-only embedding.
                - ``str`` / ``Path``: path to a PDB file.
                - ``np.ndarray`` / ``torch.Tensor``: atom37 coordinates of shape
                  ``(L, 37, 3)``.
            apply_model_norm: If ``True`` (default), apply the model's trained
                final ``LayerNorm`` (``model.transformer.norm``) to the embeddings
                before removing special tokens. This matches the post-norm
                representation used by the model's output heads. If ``False``,
                return the raw pre-norm embeddings directly from the last
                transformer block.

        Returns:
            Dictionary with keys:
                - ``residue_embeddings``: ``[L, D]`` tensor (special tokens removed).
                - ``mean_embedding``: ``[D]`` mean-pooled embedding.
                - ``attention_mask``: ``[L]`` boolean mask of retained residues.
        """
        if sequence is None and structure is None:
            raise ValueError("At least one of `sequence` or `structure` must be provided.")

        # ------------------------------------------------------------------
        # Build an ESMProteinTensor from the inputs.
        # ------------------------------------------------------------------
        if structure is None:
            # Sequence-only mode.
            if sequence is None:
                raise ValueError("`sequence` is required when `structure` is None.")
            protein = ESMProtein(sequence=sequence)
            protein_tensor = self.model.encode(protein)

        elif isinstance(structure, (str, os.PathLike)):
            # Try the fast ESM loader first; fall back to a robust Biopython-based
            # parser on failure.
            try:
                protein = ESMProtein.from_pdb(structure)
                if sequence is not None:
                    protein = ESMProtein(
                        sequence=sequence,
                        coordinates=protein.coordinates,
                        plddt=protein.plddt,
                    )
                protein_tensor = self.model.encode(protein)
            except Exception as e:
                # Use robust Biopython loader; it returns coordinates and the
                # observed sequence directly from the structure file.
                coords, observed_seq = load_atom37_from_structure_file(
                    structure,
                    _chain_id_hint(str(structure)),
                )
                if sequence is None:
                    sequence = observed_seq
                structure = coords
                # Fall through to the ndarray branch.

        if isinstance(structure, (np.ndarray, torch.Tensor)):
            coords = torch.as_tensor(structure, dtype=torch.float32, device=self.device)
            if coords.ndim != 3 or coords.shape[1:] != (37, 3):
                raise ValueError(
                    f"Coordinate array must have shape (L, 37, 3), got {tuple(coords.shape)}"
                )
            if sequence is None:
                raise ValueError(
                    "`sequence` is required when `structure` is a coordinate array."
                )
            seq_tok = encoding.tokenize_sequence(
                sequence, self.model.tokenizers.sequence, add_special_tokens=True
            )
            coords_tok, _plddt, struct_tok = encoding.tokenize_structure(
                coords,
                self.model.get_structure_encoder(),
                self.model.tokenizers.structure,
                reference_sequence=sequence,
                add_special_tokens=True,
            )
            protein_tensor = _BatchedESMProteinTensor(
                sequence=seq_tok.unsqueeze(0),
                structure=struct_tok.unsqueeze(0),
                coordinates=coords_tok.unsqueeze(0),
            ).to(self.device)

        else:
            raise TypeError(
                f"`structure` must be a PDB path, ndarray/tensor, or None; got {type(structure)}"
            )

        # Add batch dimension if it is missing.
        if protein_tensor.sequence is not None and protein_tensor.sequence.ndim == 1:
            protein_tensor = _BatchedESMProteinTensor.from_protein_tensor(protein_tensor).to(
                self.device
            )

        # ------------------------------------------------------------------
        # Forward pass and extract embeddings.
        # ------------------------------------------------------------------
        logits_output = self.model.logits(
            protein_tensor,
            LogitsConfig(
                sequence=True,
                structure=True,
                secondary_structure=True,
                sasa=True,
                function=True,
                residue_annotations=True,
                return_embeddings=True,
            ),
        )

        embeddings = logits_output.embeddings  # [B, L, D]  (pre-norm)
        if embeddings is None:
            raise RuntimeError("ESM3 did not return embeddings.")

        # Apply the model's trained final LayerNorm to get the post-norm
        # representation used by the output heads. The raw embeddings bypass
        # this norm (see TransformerStack.forward), so we apply it here when
        # requested to match the model's canonical representation.
        if apply_model_norm:
            embeddings = self.model.transformer.norm(embeddings)

        seq_tokens = protein_tensor.sequence[0]  # [L]
        tokenizer = self.model.tokenizers.sequence
        special_ids = {tokenizer.pad_token_id, tokenizer.cls_token_id, tokenizer.eos_token_id}
        # Some tokenizers expose BOS instead of CLS; cover both names safely.
        if getattr(tokenizer, "bos_token_id", None) is not None:
            special_ids.add(tokenizer.bos_token_id)

        valid_mask = torch.tensor(
            [tid.item() not in special_ids for tid in seq_tokens],
            dtype=torch.bool,
            device=self.device,
        )

        residue_embeddings = embeddings[0][valid_mask]
        mean_embedding = residue_embeddings.mean(dim=0)

        return {
            "residue_embeddings": residue_embeddings,
            "mean_embedding": mean_embedding,
            "attention_mask": valid_mask,
        }

    def extract_batch(
        self,
        items: Sequence[tuple[str | None, str | os.PathLike | np.ndarray | torch.Tensor | None]],
        apply_model_norm: bool = True,
    ) -> list[dict[str, torch.Tensor]]:
        """
        Extract embeddings for multiple chains (one chain per call).

        This is a simple sequential wrapper; for true batching, use the lower-level
        ``model.logits`` API with a manually built ``_BatchedESMProteinTensor``.

        Args:
            items: Sequence of ``(sequence, structure)`` tuples.
            apply_model_norm: Passed through to :meth:`extract`.
        """
        return [
            self.extract(sequence=seq, structure=struct, apply_model_norm=apply_model_norm)
            for seq, struct in items
        ]


def _smoke_test():
    pdb_root = _PROJECT_ROOT / "data" / "POOR_benchmark" / "func_prediction" / "pdbs"
    if not pdb_root.is_dir():
        raise FileNotFoundError(f"PDB directory not found: {pdb_root}")

    pdb_files = list(pdb_root.glob("*.pdb"))
    if not pdb_files:
        raise FileNotFoundError(f"No PDB files found in {pdb_root}")

    pdb_path = random.choice(pdb_files)
    print(f"Selected PDB for smoke test: {pdb_path}")

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Using device: {device}")

    embedder = ESM3Embedder(device=device)
    result = embedder.extract(structure=pdb_path)

    print(f"Residue embeddings shape: {result['residue_embeddings'].shape}")
    print(f"Mean embedding shape: {result['mean_embedding'].shape}")
    print(f"Attention mask sum (retained residues): {result['attention_mask'].sum().item()}")

    # Also test sequence-only mode.
    protein = ESMProtein.from_pdb(pdb_path)
    seq_only = embedder.extract(sequence=protein.sequence)
    print(f"Sequence-only residue embeddings shape: {seq_only['residue_embeddings'].shape}")


if __name__ == "__main__":
    _smoke_test()
