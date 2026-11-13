.. _api:

API reference
=============

Helper classes and functions from ``workflow/utils``.
Snakefiles and scripts import them via ``utils``, e.g. ``from utils.ModuleConfig import ModuleConfig``.
See :doc:`/principles/architecture` for how they fit together.

Building the workflow
---------------------

Used in Snakefiles to resolve input files, expand parameters into wildcards and define rules.
Every module creates one :class:`~utils.ModuleConfig.ModuleConfig`, which wraps
:class:`~utils.InputFiles.InputFiles` (input files per dataset) and
:class:`~utils.WildcardParameters.WildcardParameters` (one row per job and its parameters).

.. autosummary::
   :toctree: generated
   :nosignatures:

   utils.ModuleConfig.ModuleConfig
   utils.WildcardParameters.WildcardParameters
   utils.InputFiles.InputFiles
   utils.environments.get_env

Reading and writing data
------------------------

Used in scripts to read only the slots they need and to write results without copying unchanged data.

.. autosummary::
   :toctree: generated
   :nosignatures:

   utils.io.read_anndata
   utils.io.write_zarr_linked
   utils.io.write_zarr
   utils.io.link_zarr

Matrix helpers
--------------

Apply conversions to ``X``, ``raw`` and ``layers`` at once, for in-memory and dask arrays.

.. autosummary::
   :toctree: generated
   :nosignatures:

   utils.misc.dask_compute
   utils.misc.ensure_sparse
   utils.misc.ensure_dense
   utils.misc.apply_layers
   utils.accessors.select_layer
   utils.accessors.subset_hvg
