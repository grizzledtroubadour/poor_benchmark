class Registry:
    def __init__(self, name):
        self._name = name
        self._module_dict = dict()

    def __repr__(self):
        return f"{self._name} Registry"

    def register(self, name=None):
        def _register(cls):
            module_name = name if name else cls.__name__
            if module_name in self._module_dict:
                raise KeyError(f"{module_name} is already registered in {self._name}")
            self._module_dict[module_name] = cls
            return cls
        return _register

    def get(self, name):
        if name not in self._module_dict:
            raise KeyError(f"{name} is not found in {self._name} Registry")
        return self._module_dict[name]

    def build(self, config, **kwargs):
        """
        Build an instance from config. 
        Config can be a dict or a DictConfig (Hydra/OmegaConf).
        It must contain a 'name' key (or 'type' key, but we standardize on 'name' or passed explicitly).
        """
        from omegaconf import OmegaConf, DictConfig
        
        # If config is None, we can't do anything unless kwargs has enough info? 
        # But usually we call build(cfg).
        if config is None:
            config = {}

        # Convert DictConfig to container (dict) for manipulation
        if isinstance(config, DictConfig):
            params = OmegaConf.to_container(config, resolve=True)
        elif isinstance(config, dict):
            params = config.copy()
        else:
            raise ValueError("Config must be a dict or DictConfig")

        # Determine the name of the class to instantiate
        # Priority: 1. kwargs['name'] (override) 2. config['name']
        name = kwargs.pop('name', params.pop('name', None))
        
        if name is None:
             raise ValueError(f"Config must contain 'name' key to build from {self._name} registry.")

        cls = self.get(name)
        
        # Merge kwargs into params (kwargs overwrite config)
        params.update(kwargs)
        
        return cls(**params)

MODELS = Registry("Models")
DATASETS = Registry("Datasets")
OPTIMIZERS = Registry("Optimizers")
SCHEDULERS = Registry("Schedulers")
RUNNERS = Registry("Runners")
CALLBACKS = Registry("Callbacks")

# Register standard PyTorch optimizers
import torch.optim
for k, v in torch.optim.__dict__.items():
    if isinstance(v, type) and issubclass(v, torch.optim.Optimizer) and v is not torch.optim.Optimizer:
        try:
            OPTIMIZERS.register(k)(v)
        except KeyError:
            pass # Already registered

# Register standard PyTorch Schedulers
import torch.optim.lr_scheduler
for k, v in torch.optim.lr_scheduler.__dict__.items():
     if isinstance(v, type) and k[0].isupper(): # heuristic for classes
        try:
            SCHEDULERS.register(k)(v)
        except KeyError:
            pass
