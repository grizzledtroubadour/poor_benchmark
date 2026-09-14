"""
ESM2 protein sequence embedding extractor.

This script loads the locally cached ESM2 weights using transformers library
and extracts per-residue and mean-pooled embeddings from protein sequences.

Usage:
    # Run from the project root (or anywhere; the script resolves its own paths).
    python shared/scripts/embed/extract_esm2.py

Environment:
    Requires torch and transformers.
    The script uses local model weights and is optimized for memory efficiency
    with dynamic batching support.

The standalone smoke test will generate embeddings for a few sample sequences
and print the resulting shapes.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Sequence

import torch

# Make the project root importable regardless of the current working directory.
_PROJECT_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(_PROJECT_ROOT))

from transformers import AutoTokenizer, AutoModelForMaskedLM

# Default model configurations
DEFAULT_MODEL_PATH = _PROJECT_ROOT / "model_cache" / "esm2_650m"


def load_esm2_local(
    model_path: str | os.PathLike = DEFAULT_MODEL_PATH,
    device: torch.device | str | None = None,
) -> tuple:
    """
    Load an ESM2 model from a local weight directory using transformers.

    Args:
        model_path: Path to the ESM2 model directory containing config.json and weights.
        device: Torch device. Defaults to CUDA if available.

    Returns:
        A tuple of (model, tokenizer):
            - model: ESM2 model in evaluation mode.
            - tokenizer: ESM2 tokenizer for sequence encoding.
    """
    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    
    model_path = Path(model_path)
    if not model_path.is_dir():
        raise FileNotFoundError(f"ESM2 model directory not found: {model_path}")

    model = AutoModelForMaskedLM.from_pretrained(str(model_path)).to(device)
    tokenizer = AutoTokenizer.from_pretrained(str(model_path))
    
    model.eval()
    
    return model, tokenizer


class ESM2Embedder:
    """
    Convenience wrapper around ESM2 for extracting sequence-level embeddings.

    The model is loaded once in ``__init__`` and can be reused for many calls.
    Supports both single sequence and batch processing with automatic memory
    management.
    """

    def __init__(
        self,
        model_path: str | os.PathLike = DEFAULT_MODEL_PATH,
        device: torch.device | str | None = None,
        max_length: int = 1022,
    ):
        """
        Initialize the ESM2 embedder.

        Args:
            model_path: Path to the ESM2 model directory.
            device: Torch device. Defaults to CUDA if available.
            max_length: Maximum sequence length (excluding special tokens).
        """
        self.device = (
            torch.device(device)
            if device is not None
            else torch.device("cuda" if torch.cuda.is_available() else "cpu")
        )
        self.max_length = max_length
        
        self.model, self.tokenizer = load_esm2_local(model_path, self.device)
        self.model.eval()

    @torch.no_grad()
    def extract(
        self,
        sequence: str,
        name: str = "protein",
    ) -> dict[str, torch.Tensor]:
        """
        Extract an embedding for one protein sequence.

        Args:
            sequence: Amino-acid sequence (one-letter code).
            name: Optional name/identifier for the sequence.

        Returns:
            Dictionary with keys:
                - ``residue_embeddings``: ``[L, D]`` tensor (per-residue embeddings).
                - ``mean_embedding``: ``[D]`` mean-pooled embedding.
                - ``sequence_length``: int, the actual sequence length.
        """
        if not sequence:
            raise ValueError("Sequence cannot be empty.")

        # Tokenize sequence (truncate if too long)
        seq_token = torch.tensor(self.tokenizer.encode(sequence))[:self.max_length]
        
        # Create batch with attention mask
        seq_token = seq_token.unsqueeze(0).to(self.device)  # [1, L]
        attention_mask = torch.ones_like(seq_token).to(self.device)

        # Forward pass with automatic mixed precision if on CUDA
        if self.device.type == "cuda":
            with torch.amp.autocast(device_type="cuda"):
                outputs = self.model.esm(
                    seq_token,
                    attention_mask=attention_mask,
                    return_dict=True,
                )
        else:
            outputs = self.model.esm(
                seq_token,
                attention_mask=attention_mask,
                return_dict=True,
            )

        # Extract embeddings: [batch, seq_len, hidden_dim]
        # ESM2 includes <cls> and <eos> tokens
        embeddings = outputs.last_hidden_state  # [1, L, D]

        # Remove special tokens: <cls> at position 0, <eos> at the end
        # The actual sequence is from index 1 to -1
        seq_len = attention_mask[0].sum().item() - 2  # Subtract <cls> and <eos>
        residue_embeddings = embeddings[0, 1:1+seq_len].cpu()
        mean_embedding = residue_embeddings.mean(dim=0)

        return {
            "residue_embeddings": residue_embeddings,
            "mean_embedding": mean_embedding,
            "sequence_length": seq_len,
        }

    @torch.no_grad()
    def extract_batch(
        self,
        sequences: Sequence[tuple[str, str]],
        max_tokens: int | None = None,
    ) -> list[dict[str, torch.Tensor]]:
        """
        Extract embeddings for multiple sequences.

        Args:
            sequences: List of (name, sequence) tuples.
            max_tokens: Optional maximum tokens per batch for memory management.
                If None, processes all sequences in one batch.
                Recommended: 12000 for 24GB GPU, 25000 for 40GB GPU.

        Returns:
            List of embedding dictionaries, one per input sequence.
        """
        if not sequences:
            return []

        if max_tokens is None:
            # Process all at once
            return self._extract_batch_internal(sequences)
        else:
            # Dynamic batching based on token limit
            results = []
            current_batch = []
            current_max_len = 0

            for name, seq in sequences:
                seq_len = min(len(seq), self.max_length) + 2  # +2 for special tokens
                
                if not current_batch:
                    current_max_len = seq_len
                
                # Estimate tokens: batch_size * max_length_in_batch
                estimated_tokens = (len(current_batch) + 1) * max(current_max_len, seq_len)
                
                if estimated_tokens > max_tokens and current_batch:
                    # Process current batch
                    results.extend(self._extract_batch_internal(current_batch))
                    current_batch = [(name, seq)]
                    current_max_len = seq_len
                else:
                    current_batch.append((name, seq))
                    current_max_len = max(current_max_len, seq_len)
            
            # Process remaining batch
            if current_batch:
                results.extend(self._extract_batch_internal(current_batch))
            
            return results

    def _pad_data(self, tensor: torch.Tensor, dim: int, max_length: int, pad_value: int = 0) -> torch.Tensor:
        """Pad a tensor to max_length along the specified dimension."""
        pad_size = list(tensor.shape)
        pad_size[dim] = max_length - tensor.shape[dim]
        if pad_size[dim] > 0:
            pad_tensor = torch.full(pad_size, pad_value, dtype=tensor.dtype, device=tensor.device)
            return torch.cat([tensor, pad_tensor], dim=dim)
        return tensor

    def _extract_batch_internal(
        self,
        sequences: Sequence[tuple[str, str]],
    ) -> list[dict[str, torch.Tensor]]:
        """Internal batch processing without dynamic splitting."""
        # Tokenize all sequences
        tokenized = []
        actual_lengths = []
        
        for name, seq in sequences:
            seq_token = torch.tensor(self.tokenizer.encode(seq))[:self.max_length]
            tokenized.append(seq_token)
            # Actual sequence length (excluding special tokens)
            actual_lengths.append(len(seq_token) - 2)
        
        # Find max length in this batch and pad
        max_length_batch = max(len(t) for t in tokenized)
        
        seq_batch = []
        attention_masks = []
        
        for seq_token in tokenized:
            attention_mask = torch.zeros(max_length_batch)
            attention_mask[:len(seq_token)] = 1
            
            seq_token_padded = self._pad_data(seq_token, dim=0, max_length=max_length_batch)
            
            seq_batch.append(seq_token_padded)
            attention_masks.append(attention_mask)
        
        seq_batch = torch.stack(seq_batch).to(self.device)
        attention_masks = torch.stack(attention_masks).to(self.device)

        # Forward pass
        if self.device.type == "cuda":
            with torch.amp.autocast(device_type="cuda"):
                outputs = self.model.esm(
                    seq_batch,
                    attention_mask=attention_masks,
                    return_dict=True,
                )
        else:
            outputs = self.model.esm(
                seq_batch,
                attention_mask=attention_masks,
                return_dict=True,
            )

        embeddings = outputs.last_hidden_state  # [B, L, D]

        # Extract embeddings for each sequence
        results = []
        for i, (name, seq) in enumerate(sequences):
            seq_len = actual_lengths[i]
            # Remove <cls> (index 0) and extract actual sequence
            residue_embeddings = embeddings[i, 1:1+seq_len].cpu()
            mean_embedding = residue_embeddings.mean(dim=0)

            results.append({
                "residue_embeddings": residue_embeddings,
                "mean_embedding": mean_embedding,
                "sequence_length": seq_len,
            })

        return results


def _smoke_test():
    """Run a simple smoke test with sample sequences."""
    print("Running ESM2Embedder smoke test...")
    
    # Sample protein sequences
    test_sequences = [
        ("protein1", "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL"),
        ("protein2", "MAFSAEDVLKEYDRRRRMEALLLSLYYPNDRKLLDYKEWSPPRVQVECPKAPVEWNNPPSEKGLIVGHFSGIKYKGEKAQASEVDVNKMCCWVSKFKDAMRRYQGIQTCKIPGKVLSDLDAKIKAYNLTVEGVEGFVRYSRVTKQHVAAFLKELRHSKQYENVNLIHYILTDKRVDIQHLEKDLVKDFKALVESAHRMRQGHMINVKYILYQLLKKHGHGPDGPDILTVKTGSKGVLYDDSFRKIYTDLGWKFTPL"),
        ("protein3", "MKKLVLSLSLVLAFSSATAAF"),
    ]

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Using device: {device}")

    # Test single extraction
    print("\n=== Test 1: Single sequence extraction ===")
    embedder = ESM2Embedder(device=device)
    result = embedder.extract(sequence=test_sequences[0][1], name=test_sequences[0][0])
    
    print(f"Sequence: {test_sequences[0][0]}")
    print(f"Length: {result['sequence_length']}")
    print(f"Residue embeddings shape: {result['residue_embeddings'].shape}")
    print(f"Mean embedding shape: {result['mean_embedding'].shape}")
    print(f"Embedding dtype: {result['residue_embeddings'].dtype}")

    # Test batch extraction
    print("\n=== Test 2: Batch extraction ===")
    batch_results = embedder.extract_batch(test_sequences)
    
    for i, (name, seq) in enumerate(test_sequences):
        result = batch_results[i]
        print(f"Sequence: {name}, Length: {result['sequence_length']}, "
              f"Residue shape: {result['residue_embeddings'].shape}")

    # Test dynamic batching
    print("\n=== Test 3: Dynamic batching (max_tokens=5000) ===")
    batch_results_dynamic = embedder.extract_batch(test_sequences, max_tokens=5000)
    print(f"Processed {len(batch_results_dynamic)} sequences with dynamic batching")
    
    for i, (name, seq) in enumerate(test_sequences):
        result = batch_results_dynamic[i]
        print(f"Sequence: {name}, Length: {result['sequence_length']}")

    print("\nSmoke test completed successfully!")


if __name__ == "__main__":
    _smoke_test()
