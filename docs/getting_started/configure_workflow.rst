.. _configure-your-workflow:

🔧 Configure Your Workflow
==========================

Configuring your workflow requires setting global options as well as subworkflows consisting of modules.
Global configuration allows you to set output locations, computational resources, and other settings used across all modules, while module settings affect the behaviour of a module for a given task.

.. note::
   The recommended way to manage your workflow configuration files is to save them outside of the toolbox directory in a directory dedicated to your project. That way you can guarantee the separation of the toolbox and your own configuration.

You can find example configuration files under ``configs/``.

.. dropdown:: TL;DR Complete configuration for quickstart (click to expand)
   :icon: code

   .. literalinclude:: ../../configs/quickstart.yaml
      :language: yaml
      :caption: configs/quickstart.yaml

1. Global configuration
-----------------------

These settings define your project's output structure and hardware utilization.
The directories will be created automatically if they do not already exist.

**Output & Computational Settings**

* **Output Paths:** Intermediate and large files are stored under ``output_dir``, while images and summary files are stored under ``images``.
* **Hardware Acceleration:** If you have a supported GPU, set ``use_gpu`` to ``true``. In the backend, this ensures Snakemake utilizes GPU-enabled conda environments for supported rules.

.. note::
   Relative paths must be relative to the project root (where you call the pipeline), not the directory of the configuration file itself.

.. code-block:: yaml

   # Output locations
   output_dir: data/out
   images: images

   # Hardware settings
   use_gpu: false

2. Input configuration
----------------------

You can select and combine modules to create a custom workflow by specifying the input and module configuration in a YAML file.
Each instance of a workflow needs a user-defined task name and it can take any number of inputs consisting of modules.

Under each task, the ``input`` section lists module names, and each module is mapped to either its input files or the output of a previous module.

.. code-block:: yaml

   DATASETS:

     my_task:
       input:
         preprocessing:
           file_1: data/pbmc68k.h5ad
           # file_2: ... # more files if required
         integration: preprocessing # all outputs of module will automatically be used as input
         metrics: integration

     another_dataset:
       ...

See :ref:`documentation <input-file-mapping>` for more details on how to specify input mappings and the different formats you can use, as well as how they are resolved to file ids for downstream modules.

.. warning::
   There can only be one instance of a module as a key in the input mapping (in the backend this is a dictionary). But you can reuse the same module output as input for multiple other modules. The order of the entries in the input mapping doesn't matter.

3. Module configuration
-----------------------

You can configure the behaviour of each module by specifying their parameters under the same dataset name.

.. code-block:: yaml

   DATASETS:
     my_task:
       input:
         ...

       # module configuration
       preprocessing:
         highly_variable_genes:
           n_top_genes: 2000
         pca:
           n_comps: 50
         assemble:
           - normalize
           - highly_variable_genes
           - pca

       # module configuration
       integration:
         raw_counts: layers/counts
         norm_counts: layers/normcounts
         batch: batch
         methods:
           unintegrated:
           harmonypy:
           scanorama:
             batch_size: 100

       # module configuration
       metrics:
         unintegrated: layers/normcounts
         batch: batch
         label: bulk_labels
         metrics:
           - nmi
           - graph_connectivity

Each module has a specific set of parameters that can be configured.
Read more about the specific parameters in the README of the module you want to use.
