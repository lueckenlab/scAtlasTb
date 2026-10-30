# Relabel Module

This module is used to relabel and merge columns in an AnnData object based on a configuration file.
It is useful for reproducibly ingesting or mapping obs-level metadata based on existing columns or the index of `.obs`.

## Config file

**TL;DR:** Here is an example config showcasing some possible options of the `relabel` module.

```yaml
DATASETS:
  
  dataset_name:  # can replace with a representative name
    
    input:
      # input file mapping for relabel module
      relabel:
        file_1: test/input/pbmc68k.h5ad
        file_2: test/input/pbmc68k.h5ad
    
    # module configuration
    relabel:
      
      # add new columns with constant values
      constant_columns:
        dataset: "pbmc68k"
        processing_version: "1.0"
      
      # mapping setup for adding new columns
      new_columns:
        file:  test/input/mapping.tsv
        # index_col: # optional
        order:
          - cell_type
          - harmonized_label
          - lineage
      
      # mapping setup for merging existing columns
      merge_columns:
        file:  test/input/merge_test.tsv
        sep: '-'
      
      # rename existing obs columns
      rename_columns:
        bulk_labels: cell_type
        phase: cell_cycle_phase
      
      # rename existing obsm keys
      rename_obsm_keys:
        X_pca: X_emb
        old_embedding: new_embedding
      
      # mapping setup for selective update of existing columns
      selective_update:
        base_column: bulk_labels
        new_column: bulk_labels_new
        query: 'n_counts > 1000'  # optional filter
        update_map:
          louvain:
            '6': Dendritic
            '8': 'CD19+ B'
```

## Mapping modes

This module supports different modes of mapping metadata to the AnnData object.
You decide on any number of mapping modes in a single setup to map different metadata columns in a single run.
In the following, the different mapping modes are explained.

### Mapping new columns: `new_columns`

Adding new columns to the AnnData object requires a mapping file that specifies how to map existing columns to new ones.
In the configuration, you need to specify a TSV file with the mapping as well as the order of the columns to be matched.

```yaml
...
      new_columns:
        file:  test/input/mapping.tsv  # mapping file
        order:
          - cell_type  # first entry MUST be an existing column in the anndata object
          - harmonized_label  # metadata that hierarchically relates to first entry ('cell_type')
          - lineage  # metadata that hierarchically relates to first entry ('cell_type')
```

* `file`: the path to the mapping file in a tabular format (TSV, CSV, or Parquet).
* `order`: a list of columns in the mapping file that should be matched to the existing columns in the AnnData object.
  The first entry must be an existing column in the AnnData object, and the subsequent entries are the new columns to be created based on the mapping.

The new column mapping file must contain at least one column that exists in the anndata object.
The columns must include at least all columns that are specified in the `order` list.

Example for `test/input/mapping.tsv`:

```
bulk_labels	lineage
CD14+ Monocyte	myeloid
Dendritic	lymphoid
CD56+ NK	lymphoid
CD4+/CD25 T Reg	lymphoid
```

Alternatively, you can map values via index (`.obs_names`) of the AnnData object.
In that case, the first entry of the `order` list must be the name of the index column in the mapping file.

```yaml
...
      new_columns:
        index_col: index  # column in the mapping file that should be used as index
        file: test/input/mapping_index_test.tsv  # mapping file with metadata mapped to index
        order:
          - index  # first entry MUST be the index column in the mapping file
          - relabel_by_index  # new column to be added
          - another_column  # optional: additional columns can be added
```

* `index_col`: the name of the index column in the mapping file that should be used as index.
  If not specified, index mode is not used, so `index_col` must be set explicitly when mapping by index.
* `file`: the path to the mapping file in a tabular format (TSV, CSV, or Parquet).
* `order`: a list of columns in the mapping file. The first entry must be the index column name, followed by the new columns to be added.

Example for `test/input/mapping_index_test.tsv`:

```
index	relabel_by_index	cell_quality
AAACCTGAGAAACCAT-1	typeA	high
AAACCTGAGAAACCGC-1	typeB	medium
AAACCTGAGAAACCTA-1	typeA	high
```

**Important:** If the index values in the mapping file do not match the index in your AnnData object (`adata.obs_names`), the pipeline will raise an error.

#### Mapping file format

The file format needs to be tabular and can be either one of TSV, CSV, or Parquet.
Parquet files are fully [compatible with Pandas](https://pandas.pydata.org/pandas-docs/stable/user_guide/10min.html#parquet) and are more computationally efficient than text-based files, both in storage and processing speed.

```yaml
...
      new_columns:
        index_col: index
        file: test/input/mapping_index_test.parquet  # use parquet file for mapping
        order:
          - index
          - relabel_by_index
```


### Merge existing columns: `merge_columns`

The values of existing columns can be merged into a new column.
This is most useful when you want more granular categories in a single column, e.g. when a batch consists of a sample + pool.

```yaml
...
      merge_columns:  # mapping setup for merging existing columns
        file:  test/input/merge_test.tsv  # TSV file with cell label mapping
        sep: '-'  # separator for merging columns
      
```

* `file`: the path to the mapping file in a TSV format.
* `sep`: the separator used to merge the column values. The default is `'-'`, but you can specify any separator you like.

As an alternative to a TSV file, you can define merge rules directly in YAML.
In this case, each key is the output column and each value defines `columns` and optionally `file_id`:

```yaml
...
      merge_columns:
        sep: '-'
        cc_joint:
          file_id: file_1
          columns: [S_score, G2M_score]
        ct_joint:
          columns: bulk_labels,louvain
```

For inline YAML mode:

* `columns` can be either a comma-separated string or a list.
* `file_id` is optional; if omitted, the current pipeline `file_id` is used.

The merge column mapping file (`file`) must be in TSV format with the following columns:

* `file_id`: input file id that is used in the pipeline to match which file to add the new column to
* `column_name`: new column name
* `columns`: columns to merge, separated by `,` (no whitespaces)

Example for `test/input/merge_test.tsv`:

```
file_id	column_name	columns
file_1	cc_joint	S_score,G2M_score
file_1	ct_joint	bulk_labels,louvain
file_2	counts_joint	n_genes,n_counts
```

### Selective update of an existing column: `selective_update`

The selective update allows you to create a new column based on an existing column and a mapping of values.

```yaml
...
      selective_update:
        base_column: bulk_labels # column for which values should be updated
        new_column: bulk_labels_new # new column to be created
        query: 'n_counts > 1000' # optional pandas query to restrict update to specific cells
        update_map: # mapping of values to be updated based on column
          bulk_labels: # column to condition which cells get remapped
            'Dendritic':  # empty value means np.nan
            'CD19+ B': 'nan' # can map to specific values including 'nan'
          louvain:
            '6': Dendritic # reannotate cells where 'louvain == "6"' to "Dendritic"
            '8': 'CD19+ B'
          
```

* `base_column`: the starting column of which the values should be modified
* `new_column`: the name of the new column to be created. If it is not defined, it will be the same as `base_column`.
* `query`: optional pandas query string to restrict the update to specific cells. The query is combined with column-specific queries using AND logic: `(user_query) & (column.isin(mapped_values))`. Defaults to `'True'` if not specified.
* `update_map`: a mapping of values that should be updated. The keys are the columns in the AnnData object, and the values are dictionaries with the values to be updated. Alternatively, the values can be a file path to a TSV file with the mapping. An empty value (no value after `:`) means the mapped cells will be set to `np.nan`.

The selective update works by:
1. Starting with a copy of `base_column` as `new_column`
2. For each column in `update_map`, creating a query that combines the user `query` with `column.isin(mapped_keys)`
3. Only updating cells that match both the user query AND have values that are keys in the mapping dictionary
4. Applying the mapped values to matching cells

Below is an example where one update entry is defined by a TSV file instead of a dictionary:

```yaml
...
      selective_update:
        base_column: bulk_labels
        new_column: lineage
        update_map:
          lineage: test/input/mapping_test.tsv # mapping defined in a TSV file
          louvain:
            '6': myeloid
            '8': lymphoid
```

When using a mapping TSV file for `update_map`, the file must contain the `base_column` (in this example `bulk_labels`) as well as the `update_map` key column (in this example `lineage`).
The following example `test/input/mapping_test.tsv` shows a valid TSV definition of the `update_map` key `lineage`: 

```
bulk_labels	lineage
CD14+ Monocyte	myeloid
Dendritic	myeloid
CD56+ NK	lymphoid
CD4+/CD25 T Reg	lymphoid
CD19+ B	lymphoid
CD8+ Cytotoxic T	lymphoid
CD4+/CD45RO+ Memory	lymphoid
CD8+/CD45RA+ Naive Cytotoxic	lymphoid
CD4+/CD45RA+/CD25- Naive T	lymphoid
CD34+	HSCs
```

### Rename existing columns: `rename_columns`

You can rename existing columns in `adata.obs` by specifying a mapping of old to new column names.

```yaml
...
      rename_columns:
        old_column_name: new_column_name
        another_old: another_new
```

This creates new columns with the specified names containing the data from the old columns. If an old column doesn't exist in the data, a new column filled with `'nan'` will be created.

### Rename existing obsm keys: `rename_obsm_keys`

You can rename existing obsm keys in the AnnData object by specifying a mapping of old to new keys.

```yaml
...
      rename_obsm_keys:
        X_pca: X_emb  # rename 'X_pca' to 'X_emb'
        old_embedding: new_embedding
```

This will rename the obsm key `X_pca` to `X_emb` in the AnnData object, without modifying any data. If a key doesn't exist, a warning will be logged and that key will be skipped.

### Add constant columns: `constant_columns`

The constant columns feature allows you to add new columns with a single constant value across all cells.
This is useful for adding metadata such as dataset name, processing version, quality control status, or any other fixed value that applies to all observations.

```yaml
...
      constant_columns:
        column_name_1: value_1
        column_name_2: value_2
        column_name_3: value_3
```

* Each key is the name of the new column to be created
* The value is the constant value to be assigned to all cells in that column

The values can be of any type (strings, numbers, booleans, etc.) and will be broadcast to all rows in the AnnData object.

Example:

```yaml
relabel:
  constant_columns:
    dataset: "pbmc68k"
    processing_version: "1.0"
    quality_control: "passed"
```

This would add three new columns to the AnnData object:
- `dataset`: all cells would have the value "pbmc68k"
- `processing_version`: all cells would have the value "1.0"
- `quality_control`: all cells would have the value "passed"

## Other parameters

* `dask`: Whether to use Dask for dask-support when reading and writing files, only applies to inputs that are `h5ad`, because the complete input needs to be written to a new `zarr` store. Default is `True`.
* `threads`: Number of threads to use for reading and writing files. Default is `1`.

## Testing

Activate the snakemake environment and call `test/run_test.sh` with run specific Snakemake parameters.

```
conda activate snakemake
bash test/run_test.sh -n  # dry run
bash test/run_test.sh -c2  # actual run with max 2 cores
```
