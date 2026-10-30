import logging
logging.basicConfig(level=logging.INFO)
import scanpy as sc
import anndata as ad
import pandas as pd
from pprint import pformat
from tqdm import tqdm
from scipy.sparse import csr_matrix
import re

from utils.io import read_anndata, get_file_reader, write_zarr_linked, link_file
from utils.misc import dask_compute
from collect_utils import analyze_slot_columns, merge_df, align_file_indices, LARGE_SLOTS


dataset = snakemake.wildcards.dataset
files = snakemake.input
output_file = snakemake.output.zarr
sep = snakemake.params.get('sep', '--')
obs_index_col_map = snakemake.params.get('obs_index_col', {})
same_slots = snakemake.params.get('same_slots', [])
merge_slots = snakemake.params.get('merge_slots', [])
skip_slots = snakemake.params.get('skip_slots', [])
read_kwargs = {k: k for k in same_slots+merge_slots if k not in skip_slots}

# parse obs index
if isinstance(obs_index_col_map, str):
    obs_index_col_map = {file_id: obs_index_col_map for file_id in files.keys()}
assert isinstance(obs_index_col_map, dict), 'obs_index_col_map must be a dict'

# parse regex keys
obs_index_col_map = {
    file_id: next(
        (v for k, v in obs_index_col_map.items() if re.fullmatch(k, file_id)), 
        None  # Default if no match is found
    )
    for file_id in files.keys()
}

if len(files) == 1:
    logging.info('Single file, write unmodified...')
    adata = read_anndata(files[0], **read_kwargs, dask=True, backed=True)
    print(adata, flush=True)
    if files[0].endswith('.zarr'):
        write_zarr_linked(adata, in_dir=files[0], out_dir=output_file)
    else:
        adata.write_zarr(output_file)
    exit(0)


def read_file(file, **kwargs):
    if file.endswith('.zarr'):
        kwargs = {k: v for k, v in kwargs.items() if k not in LARGE_SLOTS+['X', 'raw']}
        
        adata = read_anndata(file, verbose=False, **kwargs)
        func, _ = get_file_reader(file)
        
        default_values = [
            csr_matrix(adata.shape), # for layers
            csr_matrix((adata.n_obs, 0)), # for obsm
            csr_matrix((adata.n_obs, adata.n_obs)), # for obsp
        ]
        for slot_name, default_value in zip(LARGE_SLOTS, default_values):
            slot_keys = func(file, 'r').get(slot_name, {}).keys()
            slot_values = {key: default_value for key in slot_keys}
            setattr(adata, slot_name, slot_values)
        
        return adata
    return read_anndata(file, verbose=False, **kwargs)


def get_slot(adata, slot_name):
    slot = getattr(adata, slot_name, None)
    if slot is not None:
        return slot
    return adata.__dict__.get(f'_{slot_name}')


adatas = {
    file_id: read_file(file, **read_kwargs, dask=True, backed=True)
    for file_id, file in tqdm(
        files.items(),
        desc='Read files',
        miniters=1,
        total=len(files),
    )
}
logging.info(f'File ids: \n{pformat(list(adatas.keys()))}')

for file_id, adata in adatas.items():
    assert adata.obsm is not None, f'adata.obsm is None for file_id={file_id}'

# Ensure all adatas have the same obs_names and var_names, and match their order to the first adata
write_copy_for = align_file_indices(adatas, files, merge_slots)
column_suffix_map = {}

if 'obs' in merge_slots:
    logging.info('Determine which columns need suffixes for obs...')
    column_suffix_map['obs'] = analyze_slot_columns(
        {file_id: adata.obs for file_id, adata in adatas.items()},
        n_threads=snakemake.threads
    )

if 'var' in merge_slots:
    logging.info('Determine which columns need suffixes for var...')
    column_suffix_map['var'] = analyze_slot_columns(
        {file_id: adata.var for file_id, adata in adatas.items()},
        n_threads=snakemake.threads
    )


# TODO: link slots with new slot names for merge_slots

slots = dict()
slot_link_map = dict()
in_dir_map = dict()
file_to_link = None

for file_id, _ad in adatas.items():

    file_name = files[file_id]
    # intialize file_to_link
    if not file_to_link and file_name.endswith('.zarr'):
        file_to_link = file_name

    for slot_name in merge_slots:
        logging.info(f'Merge slot "{slot_name}" for file_id={file_id}...')
        slot = get_slot(_ad, slot_name)
        update_slot_link_map = dict()
        to_link = file_name.endswith('.zarr') and file_id not in write_copy_for.get(slot_name, [])
        
        if isinstance(slot, pd.DataFrame):
            slots[slot_name] = merge_df(
                df_current=slot,
                file_id=file_id,
                df_previous=slots.get(slot_name),
                column_suffix_map=column_suffix_map.get(slot_name, {}),
                sep=sep,
            )
        elif slot_name == 'X':
            if to_link:
                update_slot_link_map[f'layers/{slot_name}{sep}{file_id}'] = f'{slot_name}'
            else:
                new_slot = {f'{slot_name}{sep}{file_id}': slot}
                slots['layers'] = slots.get('layers', {}) | new_slot
        elif hasattr(slot, 'items'):
            if to_link:
                update_slot_link_map = {
                    f'{slot_name}/{key}{sep}{file_id}': f'{slot_name}/{key}'
                    for key in slot.keys()
                }
            else:
                new_slot = {
                    f'{key}{sep}{file_id}': value
                    for key, value in slot.items()
                }
                slots[slot_name] = slots.get(slot_name, {}) | new_slot
        else:
            logging.info(f'file: {file_name}')
            logging.info(slot)
            raise NotImplementedError(f'Slot "{slot_name}" not supported')
        
        slot_link_map |= update_slot_link_map
        in_dir_map |= {
            key: file_name for key in update_slot_link_map.keys()
        }

# deal with same slots
if file_to_link:
    slot_link_map |= {slot: slot for slot in same_slots}
else:
    for slot_name in same_slots:
        slots[slot_name] = get_slot(_ad, slot_name)

logging.info('Create AnnData object...')
adata = ad.AnnData(**slots)
print(adata, flush=True)

logging.info('Compute matrix...')
compute_layers = [
    key
    for slot in LARGE_SLOTS
    for key in (getattr(adata, slot, {}) or {}).keys()
]
print(f'Compute layers: {compute_layers}', flush=True)
adata = dask_compute(adata, layers=compute_layers)

logging.info(f'Write to {output_file}...')
write_zarr_linked(
    adata=adata,
    in_dir=file_to_link,
    out_dir=output_file,
    files_to_keep=merge_slots+skip_slots,
    # slot_map=slot_link_map,
    # in_dir_map=in_dir_map,
)

# workaround hack - need to fix bug in write_zarr_linked
print('in_dir_map\n', pformat(in_dir_map), flush=True)
for slot, file in in_dir_map.items():
    link_file(in_file=f'{file}/{slot_link_map[slot]}', out_file=f'{output_file}/{slot}')
