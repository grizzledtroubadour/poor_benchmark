import logging
import os
import sys

import pytorch_lightning as pl


def setup_logger(name, log_dir):
    os.makedirs(log_dir, exist_ok=True)
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    logger.propagate = False # Prevent double logging if root logger is configured
    
    # Check if handlers already exist to avoid duplication
    if not logger.handlers:
        # File handler
        fh = logging.FileHandler(os.path.join(log_dir, "train.log"))
        fh.setLevel(logging.INFO)
        
        # Console handler
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(logging.INFO)
        
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        fh.setFormatter(formatter)
        ch.setFormatter(formatter)
        
        logger.addHandler(fh)
        logger.addHandler(ch)

    return logger


class SetupCallback(pl.Callback):
    """Save config YAML and argv to disk on fit start."""

    def __init__(self, now, logdir, ckptdir, cfgdir, config, argv_content=None):
        super().__init__()
        self.now = now
        self.logdir = logdir
        self.ckptdir = ckptdir
        self.cfgdir = cfgdir
        self.config = config
        self.argv_content = argv_content

    def on_fit_start(self, trainer, pl_module):
        os.makedirs(self.logdir, exist_ok=True)
        os.makedirs(self.ckptdir, exist_ok=True)
        os.makedirs(self.cfgdir, exist_ok=True)

        print("Project config")
        from omegaconf import OmegaConf
        print(OmegaConf.to_yaml(self.config))
        OmegaConf.save(self.config,
                       os.path.join(self.cfgdir, "{}-project.yaml".format(self.now)))

        with open(os.path.join(self.logdir, "argv_content.txt"), "w") as f:
            f.write(str(self.argv_content))
