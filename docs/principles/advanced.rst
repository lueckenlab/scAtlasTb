.. _advanced-configuration:

⚙️ Advanced configuration
=========================


Set defaults
------------

You can set module-specific defaults that will be used for all tasks (under ``configs['DATASETS']``), if the parameters have not been specified for those tasks.
This can shorten the configuration file, make it more readable and help avoid misconfiguration if you want to reuse the same configurations for multiple tasks.

Under the ``defaults`` directive, you can set the defaults in the same way as the task-specific configuration.

.. dropdown:: Example defaults for modules
   :icon: beaker

   .. code-block:: yaml

      defaults:
        preprocessing:
          highly_variable_genes:
            n_top_genes: 2000
          pca:
            n_comps: 50
          assemble:
            - normalize
            - highly_variable_genes
            - pca
        integration:
          raw_counts: raw/X
          norm_counts: X
          batch: batch
          methods:
            unintegrated:
            scanorama:
              batch_size: 100
            scvi:
              max_epochs: 10
              early_stopping: true
        metrics:
          unintegrated: layers/norm_counts
          batch: batch
          label: bulk_labels
          methods:
            - nmi
            - graph_connectivity

Additionally to the module defaults, you can set which datasets you want to include in your workflow, without having to remove or comment out any entries in ``configs['DATASETS']``.

.. code-block:: yaml

   defaults:
   ...
     datasets:
     # list of dataset/task names that you want your workflow to be restricted to
       - test
       - test2

Environment management
----------------------

The toolbox has 2 modes of managing conda environments for Snakemake:

* **Option 1:** ``env_mode: local`` (default) — Manually pre-install conda environments locally. This requires you to manually update conda environments when the YAML specifications change, but gives you full control over your environments and does not keep outdated environment copies.
* **Option 2:** ``env_mode: from_yaml`` — For maximum reproducibility (but potentially more computational overhead), let Snakemake manage the environments.

Option 1: ``env_mode: local`` (default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This option is convenient when you have limited space or a slow file system, where installing conda environments is quite expensive.
By setting ``env_mode`` to ``local``, you can manage the conda environments yourself, which is particularly convenient when you need to debug or try out multiple different versions of a package that fits your system.
This is the recommended approach if you are developing new features for the toolbox.

Set the global parameter in your configuration file. This is the default, so it will be used even when ``env_mode`` is not configured:

.. code-block:: yaml

   env_mode: local

Installing environments
^^^^^^^^^^^^^^^^^^^^^^^

To keep environment management overhead minimal, consider creating only the environments you need for your specific workflow (which is recommended for small workflows).
Each module should have a section on which environments it needs.
You can install each environment directly with the following conda command:

.. code-block:: bash

   conda env create -f envs/<env_name>.yaml

Alternatively, you can use ``install_environment.sh``, which automatically creates the environment when it doesn't exist yet, or updates it when it does:

.. code-block:: bash

   bash envs/install_environment.sh -h  # help message
   bash envs/install_environment.sh -f envs/<env_name>.yaml

If you want to pre-install all environments, ``envs/install_all_environments.sh`` provides a convenient wrapper:

.. code-block:: bash

   bash envs/install_all_environments.sh -h  # help message
   bash envs/install_all_environments.sh -n  # dry run
   bash envs/install_all_environments.sh

.. note::

   1. The script will create new environments for each file in the ``envs`` directory if they don't yet exist, and update any pre-existing environments.
   2. The environment names correspond to their respective file names and are documented under the ``name:`` directive in the ``envs/<env_name>.yaml`` file.
   3. If an environment creation fails, the script will skip that environment and you might need to troubleshoot the installation manually.
   4. Some environments require the channel priority to be set to ``flexible``. If your installation command fails, try setting ``conda config --set channel_priority flexible`` before restarting the command.

Updating environments
^^^^^^^^^^^^^^^^^^^^^

In cases where the environments have been updated, you might want a clean install of the new environment instead of updating existing ones.

You can manually remove conda environments via:

.. code-block:: bash

   conda env remove -n <env_name>

If you want to remove all toolbox-related environments:

.. code-block:: bash

   bash envs/install_all_environments.sh -r -n  # dry run (recommended)
   bash envs/install_all_environments.sh -r

This will remove all environments defined under ``envs/*.yaml``. After removing all environments, recreate your environments as needed.

Option 2: ``env_mode: from_yaml``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Snakemake supports automatically creating conda environments for each rule, which is convenient because you don't need to worry about environment management — Snakemake ensures you are always using the most up-to-date environment.

.. warning::

   Any updates to the environment YAML files will trigger a new environment to be installed without removing the old ones, which must be done manually. Additionally, if you are setting up multiple projects in different working directories, each project will require its own set of environments, increasing computational overhead. For small workflows requiring only a few environments, these concerns are minor compared to the convenience of Snakemake handling environments for you.

Set the global parameter in your configuration file:

.. code-block:: yaml

   env_mode: from_yaml

.. note::

   In this mode, **do not** pre-install the environments locally (i.e. do not combine Option 1 with Option 2) if you want to avoid redundant copies of environments.

Creating environments
^^^^^^^^^^^^^^^^^^^^^

You can trigger Snakemake to install all environments required for your workflow in advance by adding the following parameters:

.. code-block:: bash

   <snakemake_cmd> --use-conda --conda-create-envs-only --cores 1

Environments will be saved under ``.snakemake/conda`` from wherever you call the Snakemake commands, so make sure that directory has sufficient space or create a symlink for ``.snakemake/conda`` to a different location (e.g. scratch).
Refer to the `Snakemake documentation <https://snakemake.readthedocs.io/en/v7.31.1/snakefiles/deployment.html#integrated-package-management>`_ for more on configuring the location of the resulting conda environments.

Updating environments
^^^^^^^^^^^^^^^^^^^^^

Any updates to the environment specification will trigger creating a new environment under ``.snakemake/conda/envs``, but the old environment persists. You may need to clean up old environments periodically using `--conda-cleanup-envs <https://snakemake.readthedocs.io/en/v7.31.1/executing/cli.html#conda>`_:

.. code-block:: bash

   snakemake <target_rule> --conda-cleanup-envs


.. _snakemake_profiles:

Snakemake profiles
------------------

Snakemake profiles help you manage the many flags and options of a snakemake command in a single file, which will simplify the Snakemake call considerably.
The toolbox provides some example Snakemake profiles under ``.profiles``, which you can copy and adapt to your needs.

To use a profile (e.g. the local profile), add ``--profile .profiles/<profile_name>`` to your Snakemake command.
You can read more about profiles in `Snakemake's documentation <https://snakemake.readthedocs.io/en/v7.31.1/executing/cli.html#profiles>`_.


.. _cluster_execution:

Cluster execution
-----------------

Snakemake supports scheduling rules as jobs on a cluster and ``scAtlasTb`` has been developed and tested to work with SLURM.
In order for Snakemakek  to deploy jobs through SLURM create a Snakemake profile under ``.profiles/<your_profile>/config.yaml``.

.. dropdown:: Example profile for SLURM
   :icon: beaker

   Adapted from https://github.com/jdblischak/smk-simple-slurm

   .. code-block:: yaml

      cluster:
        mkdir -p logs/{rule} &&
        sbatch
          --partition={resources.partition}
          --qos={resources.qos}
          --gres=gpu:{resources.gpu}
          --cpus-per-task={threads}
          --mem={resources.mem_mb}
          --job-name={rule}
          --output=logs/%j-{rule}.out
          --parsable
      default-resources:
        - partition=cpu
        - qos=normal
        - gpu=0
        - mem_mb=90000
        - disk_mb=20000
      restart-times: 0
      max-jobs-per-second: 10
      max-status-checks-per-second: 1
      local-cores: 1
      latency-wait: 30
      jobs: 20
      keep-going: True
      rerun-incomplete: True
      printshellcmds: True
      scheduler: ilp
      use-conda: True
      cluster-cancel: scancel
      rerun-triggers:
        - mtime
        - params
        - input
        - software-env
        - code
      show-failed-logs: True


You can find detailed information on cluster execution in the `Snakemake documentation <https://snakemake.readthedocs.io/en/v7.31.1/executing/cluster.html>`_.

Additionally, you need to configure ``cpu`` and ``gpu`` resource settings in your workflow config file (NOT the snakemake profile config), under ``resources``.
These settings will be used by scAtlasTb to determine how to schedule the different rules on the cluster, depending on whether they require GPU or not.
For each resource profile, you need to set the memory requirements (``mem_mb``), ``partition``, ``qos``, and ``gpu`` count.
If your system does define any ``qos``, you can set it to ``normal``, which is the default.

.. code-block:: yaml

  use_gpu: true  # assuming you have and want to use GPU nodes
  resources:
    cpu: # set CPU resource settings for rules that do not require GPU; must be called "cpu"
      partition: cpu
      qos: normal
      gpu: 0
      mem_mb: 100000
    gpu: # set GPU resource settings for rules that require GPU; must be called "gpu"
      partition: gpu
      qos: normal
      gpu: 1
      mem_mb: 100000

In order for jobs to make use of the GPU, make sure that your GPU environments are installed correctly (see `Working with GPUs <troubleshooting.html#working-with-gpus>`_).
If you don't have GPU nodes, you can configure the gpu resources to be the same as the cpu resources.

.. code-block:: yaml

  use_gpu: false
  resources:
    cpu:
      partition: cpu
      qos: normal
      gpu: 0
      mem_mb: 100000
    gpu:
      partition: cpu
      qos: normal
      gpu: 1
      mem_mb: 100000

Finally, include the snakemake profile in your Snakemake command:

.. code-block:: bash

   <snakemake_cmd> --profile .profiles/<your_profile>

with ``<snakemake_cmd>`` being either ``snakemake`` or  your runner script (recommended, e.g. ``run.sh``).