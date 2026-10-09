.. _architecture:

🏗️ Architecture
================

This page describes how the toolbox turns a configuration file into Snakemake jobs.
The mechanisms described here are shared by all modules; the module pages under *Modules* describe what each module computes.

Overview
--------

The toolbox consists of

* a top-level Snakefile (``workflow/Snakefile``) that imports all modules,
* one directory per module (``workflow/<module>/``) containing a ``Snakefile``, ``rules/*.smk``, ``scripts/`` and optionally a ``params.tsv``,
* shared Python utilities (``workflow/utils/``),
* conda environment definitions (``envs/*.yaml``).

When Snakemake is called, the following happens:

.. code-block:: text

   1. load config (YAML/JSON) given by --configfile
   2. for each task in config['DATASETS']:
          for each module in task['input']:
              resolve the module's input files            # utils/pipeline.py::update_input_files_per_dataset
              (module references are replaced by the output files of that module)
   3. import every module Snakefile with the (updated) config
          each module creates a ModuleConfig object        # utils/ModuleConfig.py
              collect module parameters for each task      # utils/WildcardParameters.py
              build the parameter space (one row per job configuration)
              compute output file names
   4. Snakemake builds the job DAG from the requested targets (e.g. integration_all)
   5. each job runs a script in its conda environment     # utils/environments.py::get_env
          scripts read/write AnnData objects as zarr       # utils/io.py

Input resolution
----------------

Each task lists its modules under ``input``.
A module either receives file paths or the name of another module.
``update_input_files_per_dataset`` resolves module names recursively: if ``integration: preprocessing`` is configured, the output files of the ``preprocessing`` module for this task become the input files of ``integration``.
Each input file gets a *file id*, which is part of all downstream file names.
See :ref:`input-file-mapping` for all supported formats and how file ids are derived.

Parameter space and wildcards
-----------------------------

Every module defines a ``ModuleConfig`` in its Snakefile with

* ``config_params``: the configuration keys of the module that are read for each task (e.g. ``methods``, ``batch``, ``label`` for integration),
* ``wildcard_names``: the parameters that become Snakemake wildcards and hence part of the output path,
* ``explode_by``: list-valued parameters that are expanded to one job per entry (e.g. one job per integration method),
* ``parameters`` (optional): a module-specific ``params.tsv`` with static per-method information (e.g. conda environment, output types, CPU/GPU resources).

``WildcardParameters`` builds a table with one row per task and parameter combination:

.. code-block:: text

   for task in tasks:
       row = {dataset: task}
       for param in config_params:
           row[param] = task[module][param] if set else defaults[module][param]
   table = rows
   for column in explode_by:
       table = table.explode(column)              # one row per list entry
   table = table.merge(params.tsv, on=method)     # add static method information
   table = table.merge(input file ids)            # one row per input file

Rules query this table with ``mcfg.get_from_parameters(wildcards, key)`` to retrieve the parameters of a given job.
Default values for all tasks can be set under the ``defaults`` key of the configuration (see :ref:`advanced-configuration`).

Output files
------------

All module outputs are written to ``<output_dir>/<module>/`` and plots to ``<images>/<module>/``.
Within these, paths follow the wildcard pattern of the module, e.g.

.. code-block:: text

   data/out/integration/dataset~my_task/file_id~preprocessing:file_1/batch~batch/var_mask~None/method~scanorama--hyperparams~<hash>--label~None--output_type~embed.zarr

* ``dataset~<task>``: task name from the configuration
* ``file_id~<module>:<id>``: input file id, prefixed by the module that produced it
* other ``<wildcard>~<value>`` entries: the module parameters
* ``<hash>``: dictionaries (e.g. method hyperparameters) and long file ids are replaced by a short hash

When the output of one module is the input to the next, the downstream file id is composed of the upstream module name and its wildcards (e.g. ``integration:preprocessing:file_1--integration=<hash>``).
Long names are shortened by hashing groups of wildcards (``ModuleConfig.get_output_files``).
The mapping between file ids and files is written to ``<output_dir>/<module>/input_files.tsv`` and ``output_files.tsv`` during the dry run.

Data format and linking
-----------------------

Intermediate and output data are stored as AnnData objects in `zarr <https://zarr.readthedocs.io>`_ format.
Inputs can be ``.h5ad`` or ``.zarr``.
To avoid copying large matrices, scripts only write the slots they changed and symlink all unchanged slots (e.g. ``X``, ``layers``) to the input zarr store (``utils/io.py::write_zarr_linked``).
Scripts read only the slots they need (``utils/io.py::read_anndata``), optionally as `dask <https://www.dask.org>`_ arrays for out-of-core processing of large datasets.

Software environments
---------------------

Each rule runs in a conda environment defined under ``envs/``.
The environment is chosen by ``utils/environments.py::get_env``:

.. code-block:: text

   get_env(config, env_name, gpu_env=None):
       if config['use_gpu'] and gpu_env is set:
           env_name = gpu_env                      # e.g. rapids_singlecell instead of scanpy
       if config['env_mode'] == 'local':
           return env_name                         # use pre-installed named environment
       if config['env_mode'] == 'from_yaml':
           return f'envs/{env_name}.yaml'          # Snakemake creates the environment

For method-based modules (e.g. integration, metrics), the environment of each method is defined in the module's ``params.tsv`` (column ``env``, with ``cpu_env`` used instead when ``use_gpu: false``).
See :ref:`advanced-configuration` for how to choose the environment mode.

Computational resources
-----------------------

Rules request resources (memory, partition, GPU) from profiles defined under ``resources`` in the configuration (``ModuleConfig.get_resource``).
Each rule declares whether it uses the ``cpu`` or ``gpu`` profile; with ``use_gpu: false`` all rules fall back to the ``cpu`` profile.
Memory requests increase with each retry of a failed job, and GPU jobs fall back to CPU after a configurable number of attempts.
Example resource configurations are under ``configs/computational_resources/``.

Common rules
------------

Every module imports the ``common`` module (``workflow/common``), which provides rules shared across modules:
saving the resolved configuration, plotting the rule graph and job DAG, generic plotting rules (dot plots, embeddings, UMAPs, bar and swarm plots), adding ``obs`` annotations and converting zarr to h5ad.
