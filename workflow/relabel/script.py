import logging
logging.basicConfig(level=logging.INFO)
from tqdm import tqdm
from tqdm.dask import TqdmCallback
from pathlib import Path
import pandas as pd
from dask import dataframe as dd
from pprint import pformat
from dask import config as da_config

n_threads = snakemake.threads
da_config.set(num_workers=n_threads)
logging.info(f"Dask using {da_config.get('num_workers')} workers")

from utils.io import read_anndata, write_zarr_linked


def make_unique_index(df, index_name):
    if df.index.has_duplicates:
        logging.warning(
            f'Found {df.index.duplicated().sum()} duplicate entries in index "{index_name}";'
             ' making unique by appending suffixes.'
        )
        counts = df.groupby(level=0).cumcount()
        df.index = [f"{label}-{i}" if i > 0 else str(label) for label, i in zip(df.index, counts)]
        df.index.name = index_name
    return df


def convert_dtypes(col):
    # check if can be cast to boolean
    if set(col.dropna().astype(str).str.lower()) <= {"true", "false"}:
        bool_map = {"true": True, "false": False}
        return col.astype(str).str.lower().map(bool_map).astype('boolean')
    try:
        return pd.to_numeric(col)
    except Exception:
        return col


def read_table(file, index_name, obs_names, usecols):
    kwargs = dict(
        usecols=usecols,
        comment='#',
        na_values=['nan', 'NaN', 'NA', 'na', '', 'None', None],
        keep_default_na=True,
        engine='c',
        dtype=str,
        low_memory=False,
        blocksize='200MB',
    )
    
    # define file reader
    if file.endswith('.csv'):
        read_func = dd.read_csv
    elif file.endswith('.tsv'):
        read_func = dd.read_table
    elif file.endswith('.parquet'):
        read_func = pd.read_parquet
        kwargs = dict(
            engine='pyarrow',
            columns=usecols,
            use_threads=True,
        )
    else:
        raise ValueError(f'Unsupported file type: {file}')
    
    df = read_func(file, **kwargs)
    
    # load to memory
    if isinstance(df, dd.DataFrame):
        if index_name in df.columns:
            # subset to entries in adata.obs
            overlapping_idx = df[index_name].isin(obs_names)
            df = df.loc[overlapping_idx, :].copy()
        with TqdmCallback(desc=f'Loading {file}', total=df.shape[0], miniters=1, leave=True):
            df = df.compute()
    
    if index_name in df.columns:
        df = df.set_index(index_name)
        df = make_unique_index(df, index_name)
    
    if df.index.name == index_name and index_name is not None:
        # make sure shape matches adata.obs
        index_old = df.index
        df = df.reindex(obs_names)
    
        # check that values aren't all NaN
        if df.empty or df.isna().all().all():
            raise ValueError(
                'No matching data found in mapping file.\n'
                f'Mapping file index:\n{index_old}\nobs_names:\n{obs_names}\n'
                f'Mapping file "{file}" might be mismatched with adata.obs'
            )
    
    return df.apply(convert_dtypes)


def _parse_merge_columns_config(merge_config, file_id):
    required_cols = ['file_id', 'column_name', 'columns']

    if merge_config is None:
        return pd.DataFrame(columns=required_cols)

    if isinstance(merge_config, list):
        rows = merge_config
    elif isinstance(merge_config, dict):
        if 'entries' in merge_config:
            rows = merge_config['entries']
        elif set(required_cols).issubset(merge_config.keys()):
            rows = [merge_config]
        else:
            rows = []
            for col_name, cfg in merge_config.items():
                if col_name in {'file', 'sep'}:
                    continue
                row = {'column_name': col_name}
                if isinstance(cfg, dict):
                    row |= cfg
                else:
                    row['columns'] = cfg
                rows.append(row)
    else:
        raise ValueError(
            'merge_columns must be one of: '
            'TSV file config, list[dict], dict with entries, or dict keyed by output column.'
        )

    normalized_rows = []
    for row in rows:
        if not isinstance(row, dict):
            raise ValueError(f'Each merge_columns entry must be a dictionary, got: {type(row)}')

        out_col = row.get('column_name')
        cols = row.get('columns')
        row_file_id = row.get('file_id', file_id)

        assert out_col is not None, f'"column_name" missing in merge_columns entry: {row}'
        assert cols is not None, f'"columns" missing in merge_columns entry: {row}'

        if isinstance(cols, str):
            cols = [x.strip() for x in cols.split(',') if x.strip()]
        elif isinstance(cols, (list, tuple)):
            cols = [str(x).strip() for x in cols if str(x).strip()]
        else:
            raise ValueError(f'"columns" must be a comma-separated string or list, got: {type(cols)}')

        normalized_rows.append(
            {
                'file_id': str(row_file_id),
                'column_name': str(out_col),
                'columns': ','.join(cols),
            }
        )

    return pd.DataFrame(normalized_rows, columns=required_cols)
    

input_file = snakemake.input.anndata
input_new_cols = snakemake.input.get('new_columns')
input_merge_cols = snakemake.input.get('merge_columns')
output_file = snakemake.output.zarr
file_id = snakemake.wildcards.file_id
rename_columns = snakemake.params.get('rename_columns')
merge_columns_config = snakemake.params.get('merge_columns_config')
selective_update = snakemake.params.get('selective_update')
rename_obsm_keys = snakemake.params.get('rename_obsm_keys')
constant_columns = snakemake.params.get('constant_columns')
dask = snakemake.params.get('dask', True)

files_to_keep = ['obs']
if rename_obsm_keys:
    files_to_keep.append('obsm')

logging.info('Read adata...')
kwargs = dict(dask=True, backed=True)
if input_file.endswith(('.zarr', '.zarr/')):
    kwargs |= {x: x for x in files_to_keep}
adata = read_anndata(input_file, **kwargs)
logging.info(adata.__str__())

# add constant columns
if constant_columns:
    assert isinstance(constant_columns, dict), 'constant_columns must be a dictionary'
    logging.info(f'Add constant columns:\n{pformat(constant_columns)}')
    for col_name, col_value in tqdm(constant_columns.items(), desc='Adding constant columns'):
        adata.obs[col_name] = col_value
        logging.info(f'Added column "{col_name}" with value: {col_value}')

# merge new columns
if rename_columns:
    assert isinstance(rename_columns, dict), 'rename_columns must be a dictionary'
    logging.info(f'Rename columns:\n{pformat(rename_columns)}')
    for old_col, new_col in tqdm(rename_columns.items(), desc='Renaming'):
        if old_col not in adata.obs.columns:
            logging.warning(f'Column "{old_col}" not found in adata.obs. Skipping.')
            logging.info(f'{new_col} already in adata.obs.columns: {new_col in adata.obs.columns}')
            continue
        adata.obs[new_col] = adata.obs[old_col]

if rename_obsm_keys:
    assert isinstance(rename_obsm_keys, dict), 'rename_obsm_keys must be a dictionary'
    logging.info(f'Rename obsm keys:\n{pformat(rename_obsm_keys)}')
    for old_key, new_key in tqdm(rename_obsm_keys.items(), desc='Renaming obsm keys'):
        if old_key not in adata.obsm.keys():
            logging.warning(f'Key "{old_key}" not found in adata.obsm. Skipping.')
            continue
        adata.obsm[new_key] = adata.obsm[old_key]
        del adata.obsm[old_key]


if input_new_cols:
    index_col = snakemake.params.get('index_col')
    mapping_order = snakemake.params.get('mapping_order')
    mapping_order = list(dict.fromkeys(mapping_order))
    obs_names = adata.obs_names # save copy to preserve index

    if index_col in adata.obs.columns:
        adata.obs[index_col] = adata.obs[index_col].astype(str)
        adata.obs = adata.obs.reset_index().set_index(index_col)
    elif index_col is not None:
        adata.obs_names.name = index_col
    
    logging.info(f'Read new columns from "{input_new_cols}" with index "{index_col}"...')
    label_mapping = read_table(
        file=input_new_cols,
        index_name=index_col,
        obs_names=adata.obs_names,
        usecols=mapping_order
    )

    map_first_col = mapping_order[0]
    use_index = index_col is not None and label_mapping.index.name == map_first_col

    # Critical for large data: ensure 1:1 or M:1 join by dropping duplicates in mapping
    if use_index:
        dup = label_mapping.index.duplicated()
        key = f"{input_new_cols} index"
    else:
        dup = label_mapping[map_first_col].duplicated()
        key = f"{input_new_cols} column '{map_first_col}'"

    if dup.any():
        logging.warning(f"Duplicates in {key}; keeping first.")
        label_mapping = label_mapping[~dup]
    
    # set up merging
    join_kwargs = dict(how='left')

    # handle index-based merging
    if use_index:
        # set join parameters
        join_kwargs['left_index'] = True
        join_kwargs['right_index'] = True

        if index_col in mapping_order:
            # no need to merge on index column if already set as index
            mapping_order.remove(index_col)

        # remove any duplicate columns from adata.obs
        dup_cols = [col for col in mapping_order if col in adata.obs.columns]
        logging.warning(f'Columns "{dup_cols}" already exist in adata.obs; dropping before merge.')
        adata.obs = adata.obs.drop(columns=dup_cols)

        # make unique index if needed to prevent merge issues
        adata.obs = make_unique_index(adata.obs, index_col)        

    else: # default merging on columns will be used on all present columns
        join_kwargs['on'] = adata.obs.columns.intersection(mapping_order).tolist()

    # merge
    logging.info(f'Merging new columns with adata.obs using join_kwargs:\n{pformat(join_kwargs)}')
    adata.obs = adata.obs.merge(label_mapping[mapping_order], **join_kwargs)
    
    # restore original obs_names
    adata.obs_names = obs_names

    logging.info(f"Successfully mapped columns:\n{pformat(mapping_order)}")
    if use_index:
        index_col = adata.obs_names.to_series()
    else:
        assert map_first_col in adata.obs.columns, f'Mapping column "{map_first_col}" not found in adata.obs after merge'
        index_col = adata.obs[map_first_col]
    logging.info(f"Mapping rate: {index_col.notna().mean():.1%}")


# merge existing columns
if input_merge_cols or merge_columns_config:
    sep = snakemake.params.get('merge_sep', '-')
    file_id = file_id.split(':')[-1]  # in case file_id contains ':'

    if input_merge_cols:
        merge_source = input_merge_cols
        merge_config_df = pd.read_table(input_merge_cols, comment='#')
    else:
        merge_source = 'inline merge_columns config'
        merge_config_df = _parse_merge_columns_config(merge_columns_config, file_id=file_id)
    
    merge_config_df = merge_config_df[merge_config_df['file_id'] == file_id]
    
    if merge_config_df.shape[0] == 0:
        logging.warning(f'No entries found in {merge_source} for {file_id}')
    
    else:
        logging.info(f'Merge columns:\n{pformat(merge_config_df)}')
        
        # check if required columns are present
        for col in ['file_id', 'column_name', 'columns']:
            assert col in merge_config_df.columns, \
                f'"{col}" not found in {merge_source}\n{merge_config_df}'
        assert not merge_config_df.duplicated().any(), \
            f'Duplicated rows in {merge_source} for {file_id}\n{merge_config_df[merge_config_df.duplicated()]}'
        
        for _, row in tqdm(
            merge_config_df.iterrows(),
            total=merge_config_df.shape[0],
            desc=f'Merge values with sep="{sep}"',
        ):
            col_name = row['column_name']
            cols = [x.strip() for x in str(row['columns']).split(',') if x.strip()]
            adata.obs[col_name] = adata.obs[cols].apply(
                lambda x: sep.join(x.values.astype(str)), axis=1
            )


if selective_update:
    base_column = selective_update.get('base_column')
    new_column = selective_update.get('new_column', base_column)
    update_map = selective_update.get('update_map', {})
    query = selective_update.get('query', 'True')

    assert isinstance(update_map, dict)
    logging.info(f'query: {query}')
    
    adata.obs[new_column] = adata.obs[base_column].astype(str)
    
    logging.info(f'Selective update of column "{base_column}" to "{new_column}" using update_map:\n{pformat(update_map)}')
    for col, _dict in update_map.items():
        assert col in adata.obs.columns, f'"{col}" not found in adata.obs'
        
        # parse to dictionary if not already
        if isinstance(_dict, str) and Path(_dict).exists():
            logging.info(f'Read mapping from file: {_dict}')
            df = pd.read_table(_dict, comment='#')
            _dict = df[[col, base_column]].drop_duplicates().set_index(col)[base_column].to_dict()
        
        # build query mask
        labels = list(_dict.keys())
        col_query = f'`{col}`.isin(@labels)'
        
        # combine with user query if provided
        combined_query = f'({query}) & ({col_query})'

        # log mapping before update
        columns = list(dict.fromkeys([col, base_column]))
        value_counts = adata.obs.query(combined_query)[columns].value_counts(dropna=False, sort=False)
        logging.info(f'Before mapping "{col}"\nquery={combined_query}\n{labels=}\nn={value_counts.sum()}:\n{pformat(value_counts)}')
        
        # determine new values
        dtype = adata.obs[new_column].dtype
        s = adata.obs[col].map(_dict).astype(dtype)
        
        # evaluate query
        mask = adata.obs.eval(combined_query) # & ~s.isna()
        
        # map values selectively
        adata.obs.loc[mask, new_column] = s.loc[mask]

        # log mapping results
        columns = list(dict.fromkeys([col, new_column]))
        value_counts = adata.obs.query(combined_query)[columns].value_counts(dropna=False, sort=False)
        logging.info(f'After mapping "{col}"\nquery={combined_query}\n{labels=}\nn={value_counts.sum()}:\n{pformat(value_counts)}')


logging.info(f'Write to {output_file}...')
write_zarr_linked(
    adata,
    in_dir=input_file,
    out_dir=output_file,
    files_to_keep=files_to_keep,
    compute=not dask,
)
