"""
Base LightningDataModule with DataLoader infrastructure.

Extracted from KPLM's src/interface/data_interface.py.
Subclasses override data_setup() and optionally data_process_fn().
"""
import inspect
import importlib
import pytorch_lightning as pl
from torch.utils.data import DataLoader


class DInterface_base(pl.LightningDataModule):
    def __init__(self, num_workers=8, dataset='', **kwargs):
        super().__init__()
        self.save_hyperparameters()
        self.num_workers = num_workers
        self.dataset = dataset
        self.kwargs = kwargs
        self.batch_size = kwargs.get('batch_size', 4)
        self.task_name = kwargs.get("task_name")
        self.finetune_type = kwargs.get("finetune_type")
        print("batch_size", self.batch_size)
        print("task_name", self.task_name)

    def train_dataloader(self):
        return DataLoader(self.train_set, batch_size=self.batch_size,
                         num_workers=self.num_workers, shuffle=True,
                         prefetch_factor=3,
                         collate_fn=getattr(self, 'data_process_fn', None))

    def val_dataloader(self):
        return DataLoader(self.val_set, batch_size=self.batch_size,
                         num_workers=self.num_workers, shuffle=False,
                         collate_fn=getattr(self, 'data_process_fn', None))

    def test_dataloader(self):
        return DataLoader(self.test_set, batch_size=self.batch_size,
                         num_workers=self.num_workers, shuffle=False,
                         collate_fn=getattr(self, 'data_process_fn', None))
