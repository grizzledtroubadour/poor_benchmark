"""
Shared utility helpers for the ProtHub monorepo.

Centralized imports for commonly used utilities across subprojects:
  - Device management (CUDA / NPU / CPU auto-detection)
  - Reproducibility (seed_everything)
  - Logging
  - Checkpointing
  - Registry system (models, datasets, runners, optimizers, callbacks)
"""

# Lazy imports to avoid hard dependencies on optional packages (e.g. torch_npu)
__all__ = [
    "seed_everything",
    "get_device",
    "setup_device",
    "to_device",
    "detect_npu_available",
    "setup_logger",
    "SetupCallback",
    "process_args",
    "AverageMeter",
    "save_checkpoint",
    "set_seed",
    "clear_device_cache",
    "load_object",
    "pmap_multi",
    "MODELS",
    "DATASETS",
    "OPTIMIZERS",
    "SCHEDULERS",
    "RUNNERS",
    "CALLBACKS",
    "Config",
    "load_config",
]


def __getattr__(name):
    """Lazy import to avoid circular dependencies and optional packages."""
    # device.py exports
    if name in ("seed_everything", "get_device", "setup_device", "to_device", "detect_npu_available"):
        from shared.utils.device import (
            seed_everything,
            get_device,
            setup_device,
            to_device,
            detect_npu_available,
        )
        return locals()[name]

    # logger.py
    if name == "setup_logger":
        from shared.utils.logger import setup_logger
        return setup_logger

    # logger.py
    if name == "SetupCallback":
        from shared.utils.logger import SetupCallback
        return SetupCallback

    # args.py
    if name == "process_args":
        from shared.utils.args import process_args
        return process_args

    # common.py
    if name in ("AverageMeter", "save_checkpoint", "set_seed", "clear_device_cache", "load_object", "pmap_multi"):
        from shared.utils.common import (
            AverageMeter,
            save_checkpoint,
            set_seed,
            clear_device_cache,
            load_object,
            pmap_multi,
        )
        return locals()[name]

    # registry.py
    if name in ("MODELS", "DATASETS", "OPTIMIZERS", "SCHEDULERS", "RUNNERS", "CALLBACKS"):
        from shared.utils.registry import (
            MODELS,
            DATASETS,
            OPTIMIZERS,
            SCHEDULERS,
            RUNNERS,
            CALLBACKS,
        )
        return locals()[name]

    # config.py
    if name in ("Config", "load_config"):
        from shared.utils.config import Config, load_config
        return locals()[name]

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
