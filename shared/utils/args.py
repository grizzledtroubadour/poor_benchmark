"""
Config merging utilities: flatten_dict + process_args.

Extracted from KPLM's src/utils/utils.py so that all subprojects
can share the same 3-way merge: CLI > Hydra YAML > argparse defaults.
"""
import sys


def flatten_dict(d, parent_key='', sep='.', level=0):
    """
    Recursively flatten a nested dict into a flat dict.
    Level <= 1 strips the parent key (so top-level YAML namespaces
    like Output:, Data:, Training: become flat top-level keys).
    """
    items = {}
    for k, v in d.items():
        if level <= 1:
            new_key = k
        else:
            new_key = f"{parent_key}{sep}{k}" if parent_key else k

        if isinstance(v, dict):
            items.update(flatten_dict(v, new_key, sep=sep, level=level + 1))
        else:
            items[new_key] = v
    return items


def process_args(parser, config_path):
    """
    Merge argparse defaults, Hydra YAML config, and CLI overrides.

    Priority: CLI > config_file > parser_default.

    Args:
        parser: argparse.ArgumentParser with all defined arguments
        config_path: relative path to the Hydra config directory

    Returns:
        args Namespace with all keys merged
    """
    from hydra import initialize, compose
    from omegaconf import DictConfig, OmegaConf

    def eval_resolver(expr: str):
        return eval(expr, {}, {})

    OmegaConf.register_new_resolver("eval", eval_resolver, use_cache=False)

    # 1) parser defaults only (no CLI)
    defaults = parser.parse_args([])
    defaults_dict = vars(defaults)

    # 2) parse with CLI overrides
    args = parser.parse_args()
    args_dict = vars(args)

    # 3) load Hydra config and flatten
    with initialize(config_path=config_path):
        cfg: DictConfig = compose(config_name=args.config_name)
    config_dict = flatten_dict(OmegaConf.to_container(cfg, resolve=True))

    # 4) detect which keys were explicitly passed on CLI
    passed = set()
    for tok in sys.argv[1:]:
        if not tok.startswith('--'):
            continue
        key = tok.lstrip('-').split('=')[0].replace('-', '_')
        passed.add(key)

    # 5) merge: CLI > config_file > parser_default
    merged = {}
    processed_config = {}
    for k, v in config_dict.items():
        new_k = k.split('.')[-1] if '.' in k else k
        processed_config[new_k] = v

    all_keys = set(defaults_dict.keys()).union(set(processed_config.keys()))

    for key in all_keys:
        if key in passed:
            merged[key] = args_dict[key]
        elif key in processed_config:
            merged[key] = processed_config[key]
        else:
            merged[key] = defaults_dict[key]

    args.__dict__.update(merged)
    return args
