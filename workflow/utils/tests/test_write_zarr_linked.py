"""
Tests for write_zarr_linked with h5ad inputs, which cannot be symlinked.
Each scenario runs on an h5ad and a zarr copy of the same data and checks that both give the same output.

Run with pytest:
    python -m pytest workflow/utils/tests
or without pytest (e.g. in the scanpy conda environment):
    python workflow/utils/tests/test_write_zarr_linked.py
"""
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # workflow directory
from utils.io import read_anndata, write_zarr_linked


N_OBS, N_VARS = 40, 15


def make_adata(seed=0):
    rng = np.random.default_rng(seed)
    counts = sparse.random(N_OBS, N_VARS, density=0.3, format='csr', random_state=seed, dtype=np.float32)
    counts.data = np.ceil(counts.data * 10)
    obs = pd.DataFrame(
        {
            'batch': pd.Categorical(rng.choice(['a', 'b', 'c'], N_OBS)),
            'score': rng.normal(size=N_OBS),
        },
        index=[f'cell_{i}' for i in range(N_OBS)],
    )
    var = pd.DataFrame(
        {'gene_symbol': [f'SYM{i}' for i in range(N_VARS)]},
        index=[f'gene_{i}' for i in range(N_VARS)],
    )
    adata = ad.AnnData(
        X=counts.copy(),
        obs=obs,
        var=var,
        layers={
            'counts': counts.copy(),
            'normcounts': counts.multiply(1 / 10).tocsr(),
        },
        obsm={'X_pca': rng.normal(size=(N_OBS, 3))},
        varm={'PCs': rng.normal(size=(N_VARS, 3))},
        obsp={'distances': sparse.random(N_OBS, N_OBS, density=0.1, format='csr', random_state=seed)},
        uns={'dataset': 'test', 'params': {'a': 1}},
    )
    return adata


def make_inputs(tmp_path):
    """Write the same data as h5ad and zarr"""
    adata = make_adata()
    h5ad_file = tmp_path / 'input.h5ad'
    zarr_file = tmp_path / 'input.zarr'
    adata.write_h5ad(h5ad_file)
    adata.write_zarr(zarr_file)
    return {'h5ad': h5ad_file, 'zarr': zarr_file}


def to_dense(x):
    if hasattr(x, 'compute'):
        x = x.compute()
    if sparse.issparse(x):
        x = x.toarray()
    return np.asarray(x)


def read_output(file):
    adata = read_anndata(str(file), verbose=False)
    return adata


def assert_adata_equal(a, b):
    assert a.shape == b.shape, (a.shape, b.shape)
    pd.testing.assert_frame_equal(a.obs, b.obs, check_categorical=False)
    pd.testing.assert_frame_equal(a.var, b.var, check_categorical=False)
    if a.X is None or b.X is None:
        assert a.X is None and b.X is None
    else:
        np.testing.assert_allclose(to_dense(a.X), to_dense(b.X))
    for slot in ['layers', 'obsm', 'varm', 'obsp', 'varp']:
        slot_a, slot_b = getattr(a, slot), getattr(b, slot)
        assert sorted(slot_a.keys()) == sorted(slot_b.keys()), (slot, list(slot_a.keys()), list(slot_b.keys()))
        for key in slot_a.keys():
            np.testing.assert_allclose(to_dense(slot_a[key]), to_dense(slot_b[key]), err_msg=f'{slot}/{key}')
    assert sorted(a.uns.keys()) == sorted(b.uns.keys())
    assert (a.raw is None) == (b.raw is None)
    if a.raw is not None:
        np.testing.assert_allclose(to_dense(a.raw.X), to_dense(b.raw.X))
        pd.testing.assert_frame_equal(a.raw.var, b.raw.var, check_categorical=False)


def run_both(tmp_path, func, **kwargs):
    """Run func(input_file, output_file) on h5ad and zarr input, return outputs read back from disk"""
    inputs = make_inputs(tmp_path)
    outputs = {}
    for fmt, in_file in inputs.items():
        out_file = tmp_path / f'output_{fmt}.zarr'
        func(in_file, out_file, **kwargs)
        outputs[fmt] = read_output(out_file)
    assert_adata_equal(outputs['h5ad'], outputs['zarr'])
    return outputs['h5ad'], make_adata()


def test_partial_load_keep_uns(tmp_path):
    def func(in_file, out_file):
        adata = read_anndata(str(in_file), obs='obs', uns='uns', verbose=False)
        adata.uns['new'] = 'value'
        write_zarr_linked(adata, in_file, out_file, files_to_keep=['uns'], verbose=False)

    out, orig = run_both(tmp_path, func)
    orig.uns['new'] = 'value'
    assert out.uns['new'] == 'value'
    assert_adata_equal(out, orig)


def test_pseudobulk_keep_uns(tmp_path):
    """e.g. marker_genes empty result: object with different n_obs, only uns is kept"""
    def func(in_file, out_file):
        adata = ad.AnnData(X=np.zeros((3, N_VARS)), uns={'rank_genes_groups': {'x': 1}})
        write_zarr_linked(adata, in_file, out_file, files_to_keep=['uns'], verbose=False)

    out, orig = run_both(tmp_path, func)
    assert out.n_obs == N_OBS
    assert 'rank_genes_groups' in out.uns
    # top-level keep replaces the slot of the input file
    orig.uns = {'rank_genes_groups': {'x': 1}}
    assert_adata_equal(out, orig)


def test_nested_layer_keep(tmp_path):
    def func(in_file, out_file):
        adata = read_anndata(str(in_file), X='X', obs='obs', var='var', layers='layers', verbose=False)
        adata.layers['x'] = adata.X * 2
        adata.var['new_col'] = np.arange(adata.n_vars)
        write_zarr_linked(
            adata, in_file, out_file,
            files_to_keep=['layers/x', 'var/new_col'],
            verbose=False,
        )

    out, orig = run_both(tmp_path, func)
    assert sorted(out.layers.keys()) == ['counts', 'normcounts', 'x']
    np.testing.assert_allclose(to_dense(out.layers['x']), orig.X.toarray() * 2)
    assert 'new_col' in out.var.columns
    assert out.obsm.keys() == orig.obsm.keys()


def test_nested_layer_keep_partial_load_merges(tmp_path):
    """
    h5ad only: a nested keep merges into the layers of the input file, even when the
    in-memory object holds no other layers.
    """
    in_file = make_inputs(tmp_path)['h5ad']
    out_file = tmp_path / 'output.zarr'
    adata = read_anndata(str(in_file), X='X', obs='obs', var='var', verbose=False)
    adata.layers['x'] = adata.X * 2
    write_zarr_linked(adata, in_file, out_file, files_to_keep=['layers/x'], verbose=False)
    out = read_output(out_file)
    assert sorted(out.layers.keys()) == ['counts', 'normcounts', 'x']


def test_slot_map_X(tmp_path):
    def func(in_file, out_file):
        adata = read_anndata(str(in_file), obs='obs', uns='uns', verbose=False)
        write_zarr_linked(
            adata, in_file, out_file,
            files_to_keep=['uns'],
            slot_map={'X': 'layers/normcounts'},
            verbose=False,
        )

    out, orig = run_both(tmp_path, func)
    np.testing.assert_allclose(to_dense(out.X), orig.layers['normcounts'].toarray())
    assert sorted(out.layers.keys()) == ['counts', 'normcounts']


def test_subset_mask(tmp_path):
    rng = np.random.default_rng(1)
    obs_mask = rng.random(N_OBS) > 0.4
    var_mask = rng.random(N_VARS) > 0.3

    def func(in_file, out_file):
        adata = read_anndata(str(in_file), obs='obs', var='var', uns='uns', verbose=False)
        adata = adata[obs_mask, var_mask].copy()
        adata.obs['new_col'] = 'x'
        write_zarr_linked(
            adata, in_file, out_file,
            files_to_keep=['obs', 'uns'],
            slot_map={'X': 'layers/counts'},
            subset_mask=(obs_mask, var_mask),
            verbose=False,
        )

    out, orig = run_both(tmp_path, func)
    assert out.shape == (obs_mask.sum(), var_mask.sum())
    expected = orig[obs_mask, var_mask]
    np.testing.assert_allclose(to_dense(out.X), expected.layers['counts'].toarray())
    np.testing.assert_allclose(to_dense(out.layers['normcounts']), expected.layers['normcounts'].toarray())
    np.testing.assert_allclose(to_dense(out.obsm['X_pca']), expected.obsm['X_pca'])
    np.testing.assert_allclose(to_dense(out.obsp['distances']), expected.obsp['distances'].toarray())
    assert (out.obs['new_col'] == 'x').all()
    assert out.obs_names.tolist() == expected.obs_names.tolist()


def test_in_dir_map_output(tmp_path):
    """e.g. preprocessing normalize: X is linked to a layer that was written from memory"""
    def func(in_file, out_file):
        adata = read_anndata(str(in_file), X='X', obs='obs', var='var', uns='uns', verbose=False)
        adata.layers['normcounts'] = adata.X * 3
        write_zarr_linked(
            adata, in_file, out_file,
            files_to_keep=['uns', 'var', 'layers/normcounts'],
            slot_map={'X': 'layers/normcounts', 'layers/counts': 'X'},
            in_dir_map={'layers/normcounts': out_file},
            verbose=False,
        )

    out, orig = run_both(tmp_path, func)
    np.testing.assert_allclose(to_dense(out.X), orig.X.toarray() * 3)
    np.testing.assert_allclose(to_dense(out.layers['normcounts']), orig.X.toarray() * 3)
    np.testing.assert_allclose(to_dense(out.layers['counts']), orig.X.toarray())


def test_in_dir_map_other_file(tmp_path):
    other = make_adata(seed=5)
    other_file = tmp_path / 'other.zarr'
    other.write_zarr(other_file)

    def func(in_file, out_file):
        adata = read_anndata(str(in_file), obs='obs', uns='uns', verbose=False)
        write_zarr_linked(
            adata, in_file, out_file,
            files_to_keep=['uns'],
            slot_map={'layers/other': 'X'},
            in_dir_map={'X': other_file},
            verbose=False,
        )

    out, orig = run_both(tmp_path, func)
    np.testing.assert_allclose(to_dense(out.layers['other']), other.X.toarray())
    # in_dir_map is keyed by input slot, so the default link of X is redirected as well
    np.testing.assert_allclose(to_dense(out.X), other.X.toarray())
    np.testing.assert_allclose(to_dense(out.layers['counts']), orig.layers['counts'].toarray())


def test_shape_mismatch(tmp_path):
    """A kept slot from a pseudobulked object cannot be combined with the input file"""
    in_file = make_inputs(tmp_path)['h5ad']
    adata = ad.AnnData(
        X=np.zeros((3, N_VARS)),
        obs=pd.DataFrame(index=['a', 'b', 'c']),
        var=pd.DataFrame(index=[f'gene_{i}' for i in range(N_VARS)]),
    )
    try:
        write_zarr_linked(adata, in_file, tmp_path / 'output.zarr', files_to_keep=['obs', 'uns'], verbose=False)
    except ValueError as e:
        assert 'Slot "obs"' in str(e), str(e)
    else:
        raise AssertionError('Expected ValueError for obs shape mismatch')
    assert not (tmp_path / 'output.zarr').exists()


if __name__ == '__main__':
    tests = [
        (name, func) for name, func in sorted(globals().items())
        if name.startswith('test_') and callable(func)
    ]
    for name, func in tests:
        with tempfile.TemporaryDirectory() as tmp:
            print(f'{name}...', flush=True)
            func(Path(tmp))
    print(f'{len(tests)} tests passed')
