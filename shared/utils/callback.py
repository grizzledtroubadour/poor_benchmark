import json
import os
import torch
import pytorch_lightning as pl

from shared.utils.logger import SetupCallback

class EpochMetricsCallback(pl.Callback):
    def __init__(self, save_path):
        self.save_path = save_path
        self.epochs = []

    def on_validation_end(self, trainer, pl_module):
        try:
            metrics = {}
            for k, v in trainer.callback_metrics.items():
                # keep numeric scalars only
                if isinstance(v, torch.Tensor):
                    metrics[k] = v.item()
                elif isinstance(v, (int, float)):
                    metrics[k] = float(v)
            metrics['epoch'] = int(trainer.current_epoch)
            self.epochs.append(metrics)
            os.makedirs(os.path.dirname(self.save_path), exist_ok=True)
            with open(self.save_path, 'w', encoding='utf-8') as f:
                json.dump(self.epochs, f, ensure_ascii=False, indent=2)
        except Exception:
            pass

