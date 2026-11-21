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
for uns in uns_list[1:]:
    uns_all |= uns
adata = ad.AnnData(uns=uns_all)

logging.info(f'Writing {output_file}...')
write_zarr_linked(
    adata,
    in_dir=input_files[0],
    out_dir=output_file,
    files_to_keep=['uns'],
)