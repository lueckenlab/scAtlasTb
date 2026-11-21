from pathlib import Path
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata, write_zarr_linked

input_files = snakemake.input.dfs
input_zarr = snakemake.input.zarr
output_file = snakemake.output.zarr

# original data, read lazily so that only .uns is changed
adata = read_anndata(input_zarr, dask=True, backed=True)

for file in input_files:
    key = f'marker_genes_{Path(file).stem}'  # group=<group>.zarr -> marker_genes_group=<group>
    uns = read_anndata(file, uns='uns').uns
    if key in uns:
        adata.uns[key] = uns[key]
    else:
        # all cells were filtered out in this run: drop any previous result
        adata.uns.pop(key, None)

logging.info(f'Writing {output_file}...')
# zarr input: only .uns is written, all other slots are linked to the input
# h5ad input: cannot be linked, so the full data is written to zarr
write_zarr_linked(
    adata,
    in_dir=input_zarr,
    out_dir=output_file,
    files_to_keep=['uns'],
)
