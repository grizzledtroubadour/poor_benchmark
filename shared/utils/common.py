import torch
import numpy as np
import random
import os
import importlib
from joblib import Parallel, delayed, parallel_backend, cpu_count
from tqdm import tqdm

def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False

def clear_device_cache(device: str):
    """
    Clear cache dynamically based on the device type.
    """
    if "cuda" in device:
        torch.cuda.empty_cache()
    elif "npu" in device:
        # Assuming torch_npu is imported elsewhere if using Ascend
        torch.npu.empty_cache()
    elif "mps" in device:
        torch.mps.empty_cache()

def save_checkpoint(state, save_dir, filename="checkpoint.pth"):
    os.makedirs(save_dir, exist_ok=True)
    torch.save(state, os.path.join(save_dir, filename))

def load_object(module_path, class_name, **kwargs):
    """
    Dynamically load a class from a module and instantiate it.
    """
    module = importlib.import_module(module_path)
    cls = getattr(module, class_name)
    return cls(**kwargs)

class AverageMeter(object):
    """Computes and stores the average and current value"""
    def __init__(self):
        self.reset()

    def reset(self):
        self.val = 0
        self.avg = 0
        self.sum = 0
        self.count = 0

    def update(self, val, n=1):
        self.val = val
        self.sum += val * n
        self.count += n
        self.avg = self.sum / self.count


def pmap_multi(pickleable_fn, data, n_jobs=None, verbose=1, desc=None, **kwargs):
    """并行 map：对 data 中每个元素 d 调用 pickleable_fn(*d, **kwargs)。

    使用 joblib loky 后端，fn 及参数需可 pickle；d 必须是 tuple，
    单参数时传成 (d,) 即可。
    """
    if n_jobs is None:
        n_jobs = cpu_count()

    # 定义一个真正可以 pickling 的函数，避免 lambda 引起问题
    def _wrapped(d):
        return pickleable_fn(*d, **kwargs)

    # tqdm 外部包裹，不要嵌入 generator 里
    data_iter = list(tqdm(data, desc=desc))

    with parallel_backend('loky'):  # 或 'multiprocessing'
        results = Parallel(n_jobs=n_jobs, verbose=verbose, timeout=None)(
            delayed(_wrapped)(d) for d in data_iter
        )
    return results
