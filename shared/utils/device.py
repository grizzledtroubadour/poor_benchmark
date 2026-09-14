"""
Device utilities: flexible selection of available compute device.
Supports CUDA if available, attempts to detect vendor NPUs when possible,
and falls back to CPU. Also provides helper to set random seeds.

All comments in this file are in English by default.
"""
import os
import random
import torch
import numpy as np


def detect_npu_available():
    """Try to detect if an NPU runtime is available.

    This tries a few known extension points (non-exhaustive) and returns True
    if an NPU-like backend is present. It's intentionally conservative.
    """
    # Some NPU builds expose torch.npu or provide vendor-specific packages.
    try:
        if hasattr(torch, "npu"):
            return True
    except Exception:
        pass

    # Try importing vendor-specific extensions (best-effort)
    try:
        import torch_npu  # type: ignore

        return True
    except Exception:
        pass

    try:
        import torch_mlu  # type: ignore

        return True
    except Exception:
        pass

    return False


def get_device(preferred: str = None):
    """Return a `torch.device` chosen from availability and preference.

    Args:
        preferred: Optional string preference, e.g. 'cuda', 'cuda:0', 'cpu', 'npu'

    Returns:
        (device, device_str) where device is `torch.device` and device_str is a
        normalized string describing the chosen device.
    """
    # If a preferred device is passed, try to honor it
    if preferred:
        pref = str(preferred).lower()
        if pref.startswith("cuda") and torch.cuda.is_available():
            return torch.device(pref), "cuda"
        if pref == "cpu":
            return torch.device("cpu"), "cpu"
        if pref.startswith("npu") and detect_npu_available():
            return torch.device(pref), "npu"

    # Auto-select: prefer CUDA, then NPU (best-effort), else CPU
    if torch.cuda.is_available():
        return torch.device("cuda"), "cuda"

    if detect_npu_available():
        return torch.device("npu"), "npu"

    return torch.device("cpu"), "cpu"


def setup_device(devices: str):
    device, device_type = get_device()
    
    # Parse devices argument
    if devices:
        device_list = [int(x.strip()) for x in devices.split(',')]
        device_count = len(device_list)
        # Set appropriate environment variable based on device type
        os.environ['CUDA_VISIBLE_DEVICES'] = devices
        # For Huawei NPU (torch_npu)
        os.environ['ASCEND_VISIBLE_DEVICES'] = devices
        os.environ['ASCEND_RT_VISIBLE_DEVICES'] = devices
    else:
        device_list = [0]  # Default to first device
        device_count = 1
    
    print(f"Device Type: {device_type}, Devices: {device_list}, Device Count: {device_count}")
    return device_type, device_count


def to_device(tensor, device: torch.device):
    """Move a tensor (or nested structure) to the target device.

    This helper keeps a small surface area so callers don't need to branch.
    """
    try:
        return tensor.to(device)
    except Exception:
        # If tensor is not a torch Tensor, raise to inform caller
        raise


def seed_everything(seed: int = 42, deterministic: bool = True):
    """Set random seeds for reproducibility across modules.

    Args:
        seed: integer seed
        deterministic: set cudnn/deterministic behavior when possible
    """
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    # Try to configure deterministic behavior if requested
    try:
        if deterministic:
            torch.backends.cudnn.deterministic = True
            torch.backends.cudnn.benchmark = False
    except Exception:
        pass


__all__ = ["get_device", "to_device", "seed_everything"]
