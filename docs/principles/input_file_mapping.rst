.. _input-file-mapping:

📥 Input file mapping guidelines
================================

Input mappings define how each module receives its input files for a dataset.

This page describes the current behavior implemented in
``workflow/utils/pipeline.py::update_input_files_per_dataset``
and ``workflow/utils/InputFiles.py::InputFiles.parse``.

Each dataset has an ``input`` section with one key per module.
For one module key, the value can be:

- an explicit mapping (dictionary),
- a single string,
- or a list of strings.

.. _input-mapping-formats:

Input mapping formats
---------------------

After defining input mappings and running a Snakemake dry run, the pipeline writes TSV maps
under the module output directory:

- ``{output_dir}/{module}/input_files.tsv``
- ``{output_dir}/{module}/output_files.tsv``

These files help debug which input file ids were resolved for a module and which output file ids
were produced.

Below are examples of how different input mapping formats are resolved to file ids for downstream modules.


.. _explicit-file-id-mapping:

1. Explicit file id mapping (recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use a dictionary when you want stable, readable downstream file ids.

.. code-block:: yaml

   DATASETS:
     my_task:
       input:
         preprocessing:
           donor_a: data/a.h5ad
           donor_b: data/b.h5ad
         integration: preprocessing

**Behavior:**

- The explicit mapping under ``preprocessing`` defines stable file ids: ``donor_a`` and ``donor_b``.
- Each key maps directly to one input file path.
- This format is recommended when you want predictable, human-readable ids.

In this example, ``integration`` uses a plain module reference,
so the resolved integration input ids are:

- ``preprocessing:donor_a``
- ``preprocessing:donor_b``


.. _single-or-list-file-paths:

1. Single string or list input (implicit ids)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A module input can be a single string or a list of strings.
These forms are first normalized by ``InputFiles.parse``:

- ``str`` -> ``[str]`` -> ``{hash(str): str}``
- ``list[str]`` -> ``{hash(item): item for item in list}``

So initial file ids are short hashes of the input strings.

.. important::

   - For plain file paths, those hashed ids remain final ids.
   - For module references, those hashed ids are replaced by upstream output ids.
   - File-path detection is currently based on the presence of ``/``.
     Strings without ``/`` are treated as module references.

.. code-block:: yaml

   DATASETS:
     my_task:
       input:
         preprocessing: data/pbmc68k.h5ad  # single file

.. code-block:: yaml

   DATASETS:
     my_task:
       input:
         merge:
           - data/a.h5ad  # list of files
           - data/b.h5ad

.. code-block:: yaml

   DATASETS:
     my_task:
      input:
        preprocessing: data/pbmc68k.h5ad  # single file
        integration: preprocessing  # module reference as plain string

**Behavior:**

- Single file path string: one hashed id.
- List of file path strings: one hashed id per path.
- Plain module reference string: upstream output ids are used directly, for example ``preprocessing:{hash("data/pbmc68k.h5ad")}``.


.. _comma-separated-input:

3. Comma-separated multiple input specifiers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Comma handling depends on a heuristic in ``update_input_files_per_dataset``.
The string is split on commas only when all conditions are true:

- the string does not contain ``.zarr``
- the string does not contain ``.h5ad``
- the string contains ``,``

.. code-block:: yaml

   DATASETS:
     my_task:
       input:
         preprocessing:
           pbmc: data/pbmc68k.h5ad
         integration: preprocessing
         collect: preprocessing,integration

**Behavior:**

- Values are split by comma and trimmed.
- Each item can be either a file path or a module reference.
- Resolved entries are merged into one mapping.
- If two entries resolve to the same file id, the later one overwrites the earlier one.

.. caution::

  A string like ``preprocessing,integration,data/external.zarr`` is not split,
  because it contains ``.zarr``. It is treated as a single path-like string due to the current heuristic.

.. _mixed-input-styles:

4. Mixed input styles
~~~~~~~~~~~~~~~~~~~~~

Different input styles can be combined within the same module input definition.
This is most useful when one downstream module should consume a mix of explicit file ids,
resolved module outputs, and direct file paths.

**Case A: Single upstream module output**

.. code-block:: yaml

   DATASETS:
     my_task:
       input:
         reference_mapping:
           pbmc: data/pbmc68k.h5ad
         merge:
           query: reference_mapping
           reference: data/reference.h5ad

Here, ``query`` comes from an upstream module reference and ``reference`` from a direct path.
If ``reference_mapping`` has one output, resolved ids for ``merge`` are:

- ``query``
- ``reference``

If ``reference_mapping`` had multiple outputs, ``query`` would expand to
``query:{upstream_output_id}``.

**Case B: Multiple upstream module outputs**

.. code-block:: yaml

   DATASETS:
     my_task:
       input:
         preprocessing:
           donor_a: data/a.h5ad
           donor_b: data/b.h5ad
         collect:
           batch: preprocessing
           external: data/external.zarr

In this example, ``collect`` has two separate explicit keys:

- ``batch`` points to a module reference (``preprocessing``).
- ``external`` points to a direct file path.

If ``preprocessing`` resolves to ``preprocessing:donor_a`` and ``preprocessing:donor_b``,
resolved ids for ``collect`` are:

- ``batch:preprocessing:donor_a`` and ``batch:preprocessing:donor_b``
  (because ``batch`` is an explicit key over a multi-output module reference),
- ``external`` (because direct file paths in an explicit mapping keep their key).

If two resolved entries produce the same file id, the later entry overwrites the earlier one.

Output naming for multi-output rules
------------------------------------

Internally, output ids are created in ``ModuleConfig.get_output_files`` and may later be
propagated through ``update_input_files_per_dataset``.

If a rule expands to multiple outputs for one input file id, each output gets its own output id entry that contains wildcards.
This affects the following modules:

- ``split_data`` where data is split by ``key`` into ``values`` specified by the user
- ``integration`` where each combination of integration method, batch variable, highly variable gene mask (var_mask) results in a separate output file
- ``sample_representation`` where each combination of representation and method results in a separate output file

Below is an example for the ``split_data`` module:

.. code-block:: yaml

   DATASETS:
     my_task:
       input:
         split_data:
           pbmc1: data/pbmc68k.h5ad
           pbmc2: data/pbmc68k.h5ad # same input file for demonstration purposes only
         collect:
           lineage: split_data
    
    split_data:
      key: lineage
      values:
        - TCD4
        - TCD8_NK

The input of ``split_data`` consists of 2 files named ``pbmc1`` and ``pbmc2`` that will each be split by the ``lineage`` key into 2 output files.

The output files include wildcards, resulting in the following output file names:
- ``split_data:filter:relabel:pbmc1--split_data_key=lineage--split_data_value=TCD4``, 
- ``split_data:filter:relabel:pbmc1--split_data_key=lineage--split_data_value=TCD8_NK``
- ``split_data:filter:relabel:pbmc2--split_data_key=lineage--split_data_value=TCD4``
- ``split_data:filter:relabel:pbmc2--split_data_key=lineage--split_data_value=TCD8_NK``

In cases of very long output ids e.g. for ``integration``, where more wildcards are in the output file name, they are shortened by ``shorten_name`` in ``ModuleConfig.get_output_files``.
In that case, wildcards from the same module are collapsed into ``module={hash({module_parameter}={module_value})}``.
This keeps ids deterministic while preventing path/name explosion.

.. _practical-recommendations:

Practical recommendations
-------------------------

1. Prefer explicit dictionary mappings for stable file ids.
2. Use module references to implicitly compose workflows from module outputs.
3. Avoid depending on auto-hashed ids when file names matter for downstream analysis.
4. Be cautious when mixing pure file paths and module references in the same input definition, as it can lead to unexpected id resolution.