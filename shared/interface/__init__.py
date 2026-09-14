"""
Shared interface base classes for PyTorch Lightning modules and data modules.

Provides MInterface_base and DInterface_base, extracted from KPLM so that
all subprojects (KPLM, PLM-SAE, model_fusion) can share the same optimizer,
scheduler, and config infrastructure.
"""
from shared.interface.model_interface import MInterface_base
from shared.interface.data_interface import DInterface_base

__all__ = ["MInterface_base", "DInterface_base"]
