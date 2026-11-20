# shared/npu — NPU (Ascend) compatibility layer

All NPU adaptations are centralized here and **enabled by device type**. Apply workarounds on NPU; on CUDA/CPU, use no-ops or performant defaults to avoid GPU overhead from eager attention, disabled MHA fastpaths, or disabled nested tensors.

## Usage

Call once in the training entry point (`KPLM/tasks/main.py`) after importing torch:

```python
from shared.npu import apply_global_patches
apply_global_patches()  # Detect the device automatically; no-op outside NPU environments.
```

Select the attention implementation in model classes instead of hard-coding `"eager"`:

```python
from shared.npu import attn_implementation
config._attn_implementation = attn_implementation(self.device)  # NPU -> eager, CUDA -> sdpa
```

## Included patches (NPU only)

| patch | Reason |
|-------|------|
| `suppress_torch_npu_warnings` | torch_npu path-ownership warnings are harmless |
| `disable_nested_tensor` | NPU does not support `torch._nested_tensor_from_mask` |
| `disable_mha_fastpath` | MHA fastpath is very slow on long sequences during evaluation (111s → 0.28s/batch) |
| `patch_torch_save` | Move NPU tensors to CPU before serialization for compatibility with torch_npu checks |

Also provides `use_additive_attn_mask()`: TransformerAdapter in `finetune_model.py` uses an additive float mask instead of `src_key_padding_mask` on NPU.
