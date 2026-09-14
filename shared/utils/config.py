import yaml
import argparse
from types import SimpleNamespace
import os

class Config(dict):
    """
    A dictionary that allows access via dot notation.
    """
    def __getattr__(self, key):
        try:
            val = self[key]
        except KeyError:
            return None
        
        if isinstance(val, dict):
            return Config(val)
        return val

    def __setattr__(self, key, value):
        self[key] = value

    def to_dict(self):
        """Recursively convert back to dict."""
        output = {}
        for k, v in self.items():
            if isinstance(v, Config):
                output[k] = v.to_dict()
            else:
                output[k] = v
        return output

def load_config(config_path):
    with open(config_path, 'r') as f:
        config_dict = yaml.safe_load(f)
    
    def dict_to_config(d):
        if isinstance(d, dict):
            return Config({k: dict_to_config(v) for k, v in d.items()})
        elif isinstance(d, list):
            return [dict_to_config(v) for v in d]
        else:
            return d

    return dict_to_config(config_dict)
