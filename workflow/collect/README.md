# Collect files

Collect and merge slots from multiple files into a single AnnData file.
This module handles cases where data is split across multiple files and needs to be consolidated, with support for merging metadata and preserving file-specific information.

## Parameters

```yaml
DATASETS:
  SchulteSchrepping2020:
    input:
      collect:
        file_1: test/input/load_data/download/SchulteSchrepping2020.h5ad
        file_2: test/input/load_data/harmonize_metadata/SchulteSchrepping2020_local.zarr
    collect:
      same_slots:
        - X
        - var
      merge_slots:
        - obs
        - obsm
      obs_index_col:
        file_2: barcode
      sep: "--"
```

## Configuration Options

* **`same_slots`**: List of AnnData slots that are identical across all input files. These slots are copied directly from the first file without modification. Common slots include `X` (expression matrix) and `var` (gene metadata).

* **`merge_slots`**: List of AnnData slots that need to be merged from all input files:
  - For **DataFrame slots** (like `obs`): Columns are merged, with file-specific columns getting suffixed with the separator and file ID (e.g., `cell_type--file_1`)
  - For **dictionary slots** (like `obsm`, `obsp`): Keys are merged, with file-specific keys getting suffixed with the separator and file ID

* **`sep`**: Separator string used to distinguish file-specific columns/keys (default: `"--"`). For example, if `sep="--"` and a column exists in file_1, it becomes `column_name--file_1`.

* **`obs_index_col`**: Specifies which column to use as the obs index for each file:
  - Can be a string (same column for all files) 
  - Can be a dict mapping file IDs to column names (as shown above)
  - Can use regex patterns as keys to match multiple file IDs

* **`skip_slots`**: List of slots to exclude from processing, even if listed in `same_slots` or `merge_slots`.

## Behavior

### Single File
If only one input file is provided, it's written directly to the output without modification.

### Multiple Files
1. **Index Alignment**: All files are aligned to have matching `obs_names` and `var_names` in the same order as the first file
2. **Large Slot Handling**: For memory efficiency, large slots (`layers`, `obsm`, `obsp`) are reloaded from disk when reordering is needed
3. **Obs Index Setting**: If `obs_index_col` is specified, the designated column becomes the new obs index
4. **Column Merging**: For `obs` merging, columns that differ across files are renamed with file suffixes, while identical columns are preserved as-is
5. **Data Linking**: For `same_slots`, data is linked/copied from the first file to avoid duplication

## Output

A single zarr file containing the merged AnnData object with:
- Aligned indices across all input files
- File-specific metadata preserved with appropriate suffixes
- Shared data slots copied from the reference file
- Merged slots combining information from all input files