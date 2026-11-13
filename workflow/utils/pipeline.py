"""
Utils that are specific to the overall pipeline but not the modules
"""

from pathlib import Path
from pprint import pformat

from .config import _get_or_default_from_config, set_defaults
from .ModuleConfig import ModuleConfig
from .InputFiles import InputFiles


# deprecated
def update_module_configs(config, params):
    """
    Update config with parameters from modules TSV
    """
    # Process params per module
    for dataset in config['DATASETS'].keys():

        def _get(params, dataset, module):
            submodules = params.query('dataset == @dataset and module == @module')['submodules']
            if len(submodules) == 0 or submodules.isna().all():
                submodules = _get_or_default_from_config(config['DATASETS'], config['defaults'], dataset, module)
            else:
                submodules = submodules.to_list()[0]
            if len(submodules) == 0:
                raise ValueError(
                    f'No submodules defined for {dataset}, {module}.'
                    ' Check that they are either defined either in the config or modules.tsv'
                )
            return submodules

        config['DATASETS'][dataset]['integration'] = _get(params, dataset, 'integration')
        config['DATASETS'][dataset]['metrics'] = _get(params, dataset, 'metrics')

    config = set_defaults(config, ['integration', 'metrics'])
    config['dataset_meta'] = str(Path(config['dataset_meta']).resolve())

    return config


# deprecated
def config_for_module(config, module):
    if 'DATASETS' not in config:
        config['DATASETS'] = {}
    datasets = config['DATASETS']
    for dataset in datasets.keys():
        # for module in datasets[dataset]['input'].keys():
        input_slot = datasets[dataset]['input']
        if input_slot is None or module not in input_slot:
            continue
        input_file = input_slot[module]
        if input_file in config['output_map']:
            file_name = config['output_map'][input_file].format(dataset=dataset)
            config['DATASETS'][dataset]['input'][module] = file_name
    return config


def update_input_files_per_dataset(
    dataset: str,
    module_name: str,
    config: dict,
    first_module: str = None,
    config_class_map: dict[str: ModuleConfig] = None,
    config_kwargs: dict[str: dict] = None,
):
    original_input = config['DATASETS'][dataset]['input'][module_name]

    if module_name == first_module:
        raise ValueError(
            f"Circular dependency detected in input file resolution: "
            f"\n{pformat(original_input)}"
            f"\nThe module '{first_module}' is being used as an input module for itself "
            f"in dataset '{dataset}'. This indicates a cycle in the module dependencies, "
            f"likely due to a module as input for '{first_module}' that directly or indirectly depends on '{first_module}' itself."
            f"Please check your configuration for circular references in the input modules."
        )
    if first_module is None:
        first_module = module_name
    if config_class_map is None:
        config_class_map = {}
    if config_kwargs is None:
        config_kwargs = {}
    
    # Check if user provided explicit file name mapping (dict) vs simple module reference (string)
    user_provided_mapping = isinstance(original_input, dict)
    
    # Parse original input into file name to input specification mapping
    file_map = InputFiles.parse(original_input)

    # Parse input specification to file paths
    input_files = {}
    for file_name, input_specifier in file_map.items():
        if '.zarr' in input_specifier or '.h5ad' in input_specifier or ',' not in input_specifier:
            # single file
            input_specifiers = [input_specifier]
        else:
            input_specifiers = input_specifier.split(',')
        
        for input_specifier in input_specifiers:  # iterate file paths (e.g. when more than 1 input module specified)
            input_specifier = input_specifier.strip()
            if '/' in input_specifier:
                # assert Path(input_specifier).exists(), f'Missing input file "{input_specifier}"'
                input_files |= {file_name: input_specifier}
                continue

            # Get output files for input module via recursion
            input_module = input_specifier  # rename for easier readability
            config = update_input_files_per_dataset(
                dataset=dataset,
                module_name=input_module,
                config=config,
                first_module=first_module,
                config_class_map=config_class_map,
                config_kwargs=config_kwargs,
            )
            
            # Create module config for input module to get output files
            ModuleConfigClass = config_class_map.get(input_module, ModuleConfig)
            input_cfg = ModuleConfigClass(
                module_name=input_module,
                config=config,
                warn=False,
                **config_kwargs.get(input_module, {})
            )
            output_files = input_cfg.get_output_files(
                subset_dict={'dataset': dataset},
                as_dict=True
            )
            
            # Handle file naming based on whether user provided explicit mapping
            if user_provided_mapping:
                # User explicitly provided file_name mapping, use their naming
                if len(output_files) == 1:
                    # Single output: use user-provided file name
                    output_path = list(output_files.values())[0]
                    input_files |= {file_name: output_path}
                else:
                    # Multiple outputs: prepend user-provided file name
                    for output_id, output_path in output_files.items():
                        new_file_id = f"{file_name}:{output_id}"
                        input_files |= {new_file_id: output_path}
            else:
                # Use output file IDs from ModuleConfig's output_map as file names
                input_files |= output_files
    
    config['DATASETS'][dataset]['input'][module_name] = input_files
    return config


def update_file_for_module_param(
    dataset: str,
    module_name: str,
    config: dict,
    key: str,
    subset_dict: dict = None,
    config_class_map: dict[str: ModuleConfig] = None,
    config_kwargs: dict[str: dict] = None,
):
    """
    :param dataset: dataset name
    :param module_name: module name
    :param config: config dict
    :param key: key in module that should specify a file name
    """
    if config_class_map is None:
        config_class_map = {}
    if config_kwargs is None:
        config_kwargs = {}
    
    dataset_config = config['DATASETS'].get(dataset, {})
    file_pattern = dataset_config.get(module_name, {}).get(key)
    if file_pattern in dataset_config['input']:
        input_module = file_pattern
        ModuleConfigClass = config_class_map.get(input_module, ModuleConfig)
        input_cfg = ModuleConfigClass(
            module_name=input_module,
            config=config,
            **config_kwargs.get(input_module, {})
        )
        if subset_dict is None:
            subset_dict = {}
        subset_dict |= {'dataset': dataset}
        dataset_config[module_name][key] = input_cfg.get_output_files(subset_dict=subset_dict)[0]
    return config