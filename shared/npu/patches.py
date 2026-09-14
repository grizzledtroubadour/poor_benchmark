"""NPU (Ascend) compatibility patches, applied only when running on NPU.

Centralizes every NPU workaround so training code can call a single entry
point instead of activating patches unconditionally (which slows down NVIDIA
GPU runs: forced eager attention, disabled MHA fastpath, etc.).

Patches included (NPU only):
  - suppress torch_npu environment-check warnings
  - disable nested-tensor fastpath (NPU lacks torch._nested_tensor_from_mask)
  - disable nn.MultiheadAttention fastpath (extremely slow on NPU for long
    sequences: 2047-len batch=64 went 111s -> 0.28s per batch)
  - wrap torch.save so NPU tensors are moved to CPU before serialization
  - set TORCH_SAVE_ZIPFILE_SERIALIZATION=1
"""
import os
import warnings

import torch


def is_npu_device(device_or_type=None) -> bool:
    """Return True when the target device is an NPU.

    Accepts a torch.device, a string ('npu', 'npu:0', ...), or None.
    None means auto-detect from the runtime environment.
    """
    if device_or_type is None:
        try:
            from shared.utils.device import detect_npu_available
        except ImportError:
            from ..utils.device import detect_npu_available
        return bool(detect_npu_available()) and not torch.cuda.is_available()
    if isinstance(device_or_type, torch.device):
        return device_or_type.type == "npu"
    return str(device_or_type).lower().startswith("npu")


def suppress_torch_npu_warnings():
    """torch_npu emits non-critical path-ownership warnings at import."""
    warnings.filterwarnings("ignore")
    warnings.filterwarnings("ignore", message=".*owner does not match the current owner.*")
    warnings.filterwarnings("ignore", category=UserWarning, module="torch_npu.*")


def disable_nested_tensor():
    """NPU backend doesn't support torch._nested_tensor_from_mask in
    PyTorch Transformers; fall back to the non-nested path."""
    try:
        torch.backends.cuda.enable_nested_tensor(False)
    except (AttributeError, RuntimeError):
        pass  # May not be available in all PyTorch versions


def disable_mha_fastpath():
    """Under eval/no_grad, the nn.TransformerEncoder MHA fastpath is extremely
    slow on NPU for long sequences (measured 2047-len batch=64: 111s -> 0.28s
    per batch after disabling)."""
    try:
        torch.backends.mha.set_fastpath_enabled(False)
    except (AttributeError, RuntimeError):
        pass  # May not be available in all PyTorch versions


def patch_torch_save():
    """Move NPU tensors to CPU before serialization (the torch_npu C-level
    hook may otherwise complain), and keep the new zipfile serialization for
    compatibility with torch_npu's check logic."""
    _original_torch_save = torch.save

    def _npu_safe_torch_save(obj, f, *args, **kwargs):
        if isinstance(obj, dict):
            obj = {k: v.cpu() if isinstance(v, torch.Tensor) and v.device.type == "npu" else v
                   for k, v in obj.items()}
        # May be a no-op in PyTorch 2.7+, but torch_npu still checks for it
        kwargs.setdefault("_use_new_zipfile_serialization", True)
        return _original_torch_save(obj, f, *args, **kwargs)

    torch.save = _npu_safe_torch_save
    os.environ["TORCH_SAVE_ZIPFILE_SERIALIZATION"] = "1"


def apply_global_patches(device_type=None):
    """Apply all NPU global patches. No-op on non-NPU environments.

    Call once, right after `import torch` and before building models.
    """
    if not is_npu_device(device_type):
        return
    suppress_torch_npu_warnings()
    disable_nested_tensor()
    disable_mha_fastpath()
    patch_torch_save()


def attn_implementation(device_or_type=None, default: str = "sdpa") -> str:
    """Attention implementation for HF model configs.

    NPU kernels lack efficient sdpa/flash paths for some of these models, so
    we force "eager" there; on CUDA we keep the fast default ("sdpa").
    """
    if is_npu_device(device_or_type):
        return "eager"
    return default


def use_additive_attn_mask(device_or_type=None) -> bool:
    """Whether TransformerAdapter must use an additive float mask instead of
    src_key_padding_mask (the latter triggers the unsupported nested-tensor
    path on NPU)."""
    return is_npu_device(device_or_type)
