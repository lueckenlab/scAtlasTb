"""
Helpers for integration methods that run in R on Seurat objects (``seurat_rpca``, ``stacas``).

The method script reads the prepared AnnData with :func:`read_counts`, hands the data to its R
script with :func:`run_r_method` and adds the returned embedding to the output.
Data and parameters are exchanged through files in a temporary directory: an HDF5 file with the
counts (read in R by ``read_seurat_input()`` in ``seurat_utils.R``), a JSON parameter file, and an
HDF5 file with the embedding written by R (``write_embedding()``).

Normalization follows Seurat's ``LogNormalize`` (``NormalizeData`` with its default
``scale.factor = 1e4``) on the raw counts, with each cell's library size taken over **all** genes
before the counts are subset to the integration features, exactly as ``NormalizeData`` on the full
object followed by subsetting the features would compute it. Only the exchanged genes are passed
to R, so library sizes are computed here.
"""
import json
import logging
import os
import subprocess
import tempfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from dask import array as da
from scipy import sparse

from utils.io import read_anndata
from utils.accessors import subset_hvg
from utils.misc import dask_compute


def r_arguments(hyperparams: dict) -> dict:
    """
    Translate hyperparameter names to R argument names: underscores become dots
    (e.g. ``k_weight`` -> ``k.weight``); names without underscores (e.g. ``genesBlockList``)
    are kept as they are.
    """
    return {key.replace('_', '.'): value for key, value in hyperparams.items()}


def library_size(X) -> np.ndarray:
    """
    Total counts per cell as float64, for a cells x genes matrix (numpy, scipy sparse or dask)
    """
    if isinstance(X, da.Array):
        block_sums = X.map_blocks(
            lambda block: np.asarray(block.sum(axis=1, dtype='float64')).reshape(-1, 1),
            chunks=(X.chunks[0], (1,) * len(X.chunks[1])),
            dtype='float64',
        )
        return np.asarray(block_sums.sum(axis=1).compute()).ravel()
    return np.asarray(X.sum(axis=1, dtype='float64')).ravel()


def read_counts(input_file: str, var_column: str = 'integration_features'):
    """
    Read the raw counts of a prepared integration input and compute library sizes on all genes

    :param input_file: prepared zarr of the integration module (raw counts in ``layers/counts``)
    :param var_column: boolean column in ``.var`` with the genes to pass to R; None keeps all genes
    :return: AnnData with the counts of the selected genes in ``.X`` (in memory),
        library size per cell over all genes
    """
    logging.info(f'Read {input_file}...')
    adata = read_anndata(
        input_file,
        X='layers/counts',
        obs='obs',
        var='var',
        uns='uns',
        dask=True,
        backed=True,
    )
    logging.info(f'Compute library sizes on all {adata.n_vars} genes...')
    size = library_size(adata.X)
    if var_column is None:
        adata = dask_compute(adata, layers='X')
    else:
        if adata.var[var_column].all():
            logging.warning(
                f'All {adata.n_vars} genes in layers/counts are in "{var_column}": if the counts were '
                'already subset to the integration features, library sizes cover only those genes.'
            )
        adata, _ = subset_hvg(adata, var_column=var_column)
    return adata, size


def _string_array(values) -> np.ndarray:
    """Strings for HDF5; missing values become empty strings (NA in R)"""
    values = pd.Series(values, dtype='object')
    values = values.where(values.notna(), '').astype(str)
    return values.to_numpy(dtype=object)


def write_input(file: Path, adata, library_size: np.ndarray, obs: dict):
    """
    Write counts (cells x genes CSR = genes x cells CSC as read by R), library sizes,
    cell and gene names and the obs columns needed by the R script
    """
    X = sparse.csr_matrix(adata.X, dtype='float64')
    X.sum_duplicates()
    X.sort_indices()
    if X.nnz > np.iinfo(np.int32).max:
        raise ValueError(
            f'{X.nnz} non-zero counts exceed the 2^31 - 1 entries that a sparse R matrix '
            '(dgCMatrix) can hold. Pass fewer genes to R (e.g. a smaller var_mask).'
        )
    string_dtype = h5py.string_dtype()
    with h5py.File(file, 'w') as f:
        counts = f.create_group('counts')
        counts.create_dataset('data', data=X.data)
        counts.create_dataset('indices', data=X.indices.astype(np.int32))
        counts.create_dataset('indptr', data=X.indptr.astype(np.int32))
        counts.create_dataset('shape', data=np.array(X.shape, dtype=np.int32))
        f.create_dataset('library_size', data=np.asarray(library_size, dtype='float64'))
        f.create_dataset('obs_names', data=_string_array(adata.obs_names), dtype=string_dtype)
        f.create_dataset('var_names', data=_string_array(adata.var_names), dtype=string_dtype)
        group = f.create_group('obs')
        for key, values in obs.items():
            group.create_dataset(key, data=_string_array(values), dtype=string_dtype)


def read_embedding(file: Path, obs_names) -> np.ndarray:
    """Read the embedding written by ``write_embedding()`` in R, in the order of obs_names"""
    with h5py.File(file, 'r') as f:
        names = f['obs_names'].asstr()[:]
        n_cells, n_dims = f['shape'][:]
        embedding = f['X_emb'][:].reshape(n_cells, n_dims)
    if not np.array_equal(names, np.asarray(obs_names, dtype=str)):
        raise ValueError('Cells of the R embedding do not match the input cells')
    return embedding


def run_r_method(
    script: Path,
    adata,
    library_size: np.ndarray,
    obs: dict,
    params: dict,
    threads: int = 1,
    tmpdir: str = None,
) -> np.ndarray:
    """
    Run an R integration script and return its embedding (cells x dimensions)

    The script is called as ``Rscript --vanilla <script> <params.json>``. Besides ``params``,
    the JSON file holds the paths of the input file (``input``), of the output file (``output``)
    and of ``seurat_utils.R`` (``utils``).

    :param script: R script
    :param adata: AnnData with the counts to pass to R in ``.X``
    :param library_size: library size per cell (over all genes)
    :param obs: obs columns to pass to R, as {name in R: values}
    :param params: parameters of the R script (JSON serialisable)
    :param threads: number of BLAS/OpenMP threads for R
    :param tmpdir: directory for the exchange files (default: system temporary directory)
    """
    script = Path(script)
    with tempfile.TemporaryDirectory(prefix=f'{script.stem}_', dir=tmpdir) as workdir:
        workdir = Path(workdir)
        input_file = workdir / 'input.h5'
        logging.info(f'Write R input to {input_file}...')
        write_input(input_file, adata, library_size, obs)

        params = params | {
            'input': str(input_file),
            'output': str(workdir / 'output.h5'),
            'utils': str(script.with_name('seurat_utils.R')),
        }
        params_file = workdir / 'params.json'
        params_file.write_text(json.dumps(params, indent=2))

        env = os.environ.copy()
        for var in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']:
            env[var] = str(threads)
        command = ['Rscript', '--vanilla', str(script), str(params_file)]
        logging.info(f'Run {" ".join(command)} with parameters:\n{json.dumps(params, indent=2)}')
        subprocess.run(command, check=True, env=env)

        return read_embedding(workdir / 'output.h5', adata.obs_names)
