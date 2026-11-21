import anndata as ad
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata, write_zarr_linked

input_files = snakemake.input.dfs
input_zarr = snakemake.input.zarr
output_file = snakemake.output.zarr


uns_list = [
    read_anndata(file, uns='uns').uns
    for file in [input_zarr]+input_files
]
uns_all = uns_list[0]
empty_marker_keys = set()
for uns in uns_list[1:]:
    empty_marker_keys.update(
        key for key, value in uns.items()
        if key.startswith('marker_genes_group=') and not value
    )
    uns_all |= uns
for key in empty_marker_keys:
    uns_all[key] = {}
adata = ad.AnnData(uns=uns_all)

logging.info(f'Writing {output_file}...')
write_zarr_linked(
    adata,
    in_dir=input_files[0],
    out_dir=output_file,
    files_to_keep=['uns'],
)