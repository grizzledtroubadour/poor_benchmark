"""NPU (Ascend) compatibility layer.

All NPU-specific workarounds live here and are applied conditionally by
device type, so CUDA runs keep their fast paths (sdpa attention, MHA
fastpath, nested tensors).
"""
from .patches import (
    apply_global_patches,
    attn_implementation,
    disable_mha_fastpath,
    disable_nested_tensor,
    is_npu_device,
    patch_torch_save,
    suppress_torch_npu_warnings,
    use_additive_attn_mask,
)

__all__ = [
    "apply_global_patches",
    "attn_implementation",
    "disable_mha_fastpath",
    "disable_nested_tensor",
    "is_npu_device",
    "patch_torch_save",
    "suppress_torch_npu_warnings",
    "use_additive_attn_mask",
]
