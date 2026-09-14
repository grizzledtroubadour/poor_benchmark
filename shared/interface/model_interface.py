"""
Base LightningModule with optimizer, scheduler, and config infrastructure.

Extracted from KPLM's src/interface/model_interface.py so that all
subprojects can share the same training infrastructure.
"""
import os
import torch
import pytorch_lightning as pl
import torch.nn as nn
import torch.optim.lr_scheduler as lrs
import inspect


class MInterface_base(pl.LightningModule):
    def __init__(self, model_name=None, loss=None, lr=None, **kargs):
        super().__init__()
        self.save_hyperparameters()
        self.load_model()
        self.configure_loss()
        os.makedirs(os.path.join(self.hparams.res_dir, self.hparams.ex_name), exist_ok=True)

    def forward(self, input):
        pass

    def training_step(self, batch, batch_idx, **kwargs):
        pass

    def validation_step(self, batch, batch_idx):
        pass

    def test_step(self, batch, batch_idx):
        return self.validation_step(batch, batch_idx)

    def on_validation_epoch_end(self):
        self.print('')

    def get_schedular(self, optimizer, lr_scheduler='onecycle'):
        if lr_scheduler == 'step':
            # StepLR 按 optimizer step 计数：把"每 N 个 epoch 衰减一次"换算成步数
            # lr_decay_steps 是总优化步数，steps_per_epoch = lr_decay_steps / epoch
            decay_epochs = getattr(self.hparams, 'lr_decay_epochs', 10)
            gamma = getattr(self.hparams, 'lr_decay_rate', 0.1)
            steps_per_epoch = max(1, int(self.hparams.lr_decay_steps / self.hparams.epoch))
            scheduler = lrs.StepLR(optimizer,
                                   step_size=max(1, int(steps_per_epoch * decay_epochs)),
                                   gamma=gamma)
        elif lr_scheduler == 'cosine':
            scheduler = lrs.CosineAnnealingLR(optimizer,
                                               T_max=self.hparams.lr_decay_steps)
        elif lr_scheduler == 'onecycle':
            scheduler = lrs.OneCycleLR(optimizer, max_lr=self.hparams.lr,
                                       total_steps=self.hparams.lr_decay_steps, three_phase=False)
        elif lr_scheduler == 'plateau':
            scheduler = lrs.ReduceLROnPlateau(optimizer, mode='min', factor=0.1,
                                              patience=3, verbose=True)
        else:
            raise ValueError('Invalid lr_scheduler type!')
        return scheduler

    def configure_optimizers(self):
        if hasattr(self.hparams, 'weight_decay'):
            weight_decay = self.hparams.weight_decay
        else:
            weight_decay = 0

        optimizer_g = torch.optim.AdamW(
            filter(lambda p: p.requires_grad, self.model.parameters()),
            lr=self.hparams.lr, weight_decay=weight_decay,
            betas=(0.9, 0.98), eps=1e-8)

        schecular_g = self.get_schedular(optimizer_g, self.hparams.lr_scheduler)
        if self.hparams.lr_scheduler == 'plateau':
            return [optimizer_g], [{"scheduler": schecular_g, "interval": "epoch", "monitor": "val_loss"}]
        else:
            return [optimizer_g], [{"scheduler": schecular_g, "interval": "step"}]

    def lr_scheduler_step(self, *args, **kwargs):
        scheduler = self.lr_schedulers()
        if self.hparams.lr_scheduler != 'plateau':
            scheduler.step()

    def configure_loss(self):
        self.loss_function = nn.CrossEntropyLoss(reduction='none')

    def load_model(self):
        self.model = None

    def instancialize(self, Model, **other_args):
        """Instancialize a model using the corresponding parameters
        from self.hparams dictionary. You can also input any args
        to overwrite the corresponding value in self.hparams."""
        class_args = inspect.getargspec(Model.__init__).args[1:]
        inkeys = self.hparams.keys()
        args1 = {}
        for arg in class_args:
            if arg in inkeys:
                args1[arg] = getattr(self.hparams, arg)
        args1.update(other_args)
        return Model(**args1)
