from .misc import ifelse, get_use_gpu


def get_env(
    config: dict,
    env_name: str,
    mode: str = None,
    env_dir: str = '../envs',
    gpu_env: str = None,
    no_gpu: bool = False,
):
    """
    Get the conda environment for a rule's ``conda:`` directive

    :param config: complete Snakemake config
    :param env_name: name of the environment
    :param mode: ``'local'`` returns the environment name (environment must already exist),
        ``'from_yaml'`` returns the path to ``<env_dir>/<env_name>.yaml`` for Snakemake to create.
        Defaults to ``config['env_mode']``, or ``'local'``
    :param env_dir: directory containing the environment YAML files, relative to the Snakefile
    :param gpu_env: environment to use instead of ``env_name`` when ``use_gpu`` is set in the config
    :param no_gpu: always use ``env_name``, even when ``use_gpu`` is set
    :return: environment name or path to environment file
    """
    if mode is None:
        mode = config.get('env_mode', 'local')
    if no_gpu or gpu_env is None:
        gpu_env = env_name
    env_name = ifelse(get_use_gpu(config), _if=gpu_env, _else=env_name)

    if mode == 'from_yaml':
        return f'{env_dir}/{env_name}.yaml'
    elif mode == 'local':
        return env_name
    else:
        raise ValueError(f'Unknown environment mode: {mode}')

