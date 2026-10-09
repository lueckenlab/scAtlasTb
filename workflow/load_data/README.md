# Workflow

Given a TSV file and a schema mapping, the pipeline does the following:

1. [Load data from CELLxGENE, DCP, URL or file](#load-data)
2. [Aggregate metadata](#aggregate-metadata)
3. [Merge datasets per study, organ or custom subset](#merge-data)
4. [Filter cells](#filter)

Examples of configuration files are under `test/configs/` (e.g. `test/configs/cellxgene.yaml` and `test/configs/dataset_info/cellxgene.tsv`).
By default, dataset configurations are available under `configs` at the top-level pipeline (git root of this repository), but they can be modified or replaced by custom files.

## Load data

Files can either be read from a specified input file or downloaded from an URL, CELLxGENE or DCP directly.
Which way a dataset is loaded, depends on the [dataset mapping](#dataset-file-mapping).

## Aggregate Metadata

The `harmonize_metadata` rule adds additional dataset-level information that is included from the input TSV file.
This steps expects the data to follow
the [CELLxGENE schema 3.0.0](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/3.0.0/schema.md)
and is extended as described below.
Other steps include:

+ adding external annotations if annotation file and columns are available in the TSV
+ saving donor IDs under `.obs['donor']`
+ building a technical sample ID from the columns listed in `tech_id` and saving it under `.obs['tech_id']` and `.obs['sample']`

The input `AnnData` objects must contain:

+ `.X` raw counts, sparse format
+ `.obs` columns from the schema defined in the schema 
+ `.obs.index` cell barcode

The output `AnnData` objects will contain:
+ `.X` raw counts, sparse format
+ `.uns['meta']` metadata from TSV file
+ `.obs` columns from CELLxGENE schema 3.0.0 and a subset of information in `.uns['meta']` from `EXTRA_COLUMNS`
  from `scripts/load_data_utils.py`
+ `.obs['dataset']`, `.uns['dataset']` name of task/dataset
+ `.obs['organ']`, `.uns['organ']` organ
+ `.obs['donor']` donor ID
+ `.obs['sample']` sample ID (same as `.obs['tech_id']`)
+ `.obs['barcode']` cell barcodes as declared in index
+ `.obs['author_annotation']` author annotation under the `author_annotation` column of the input TSV
+ `.var` gene information as specified in CELLxGENE schema 3.0.0
+ `.obs.index` unique cell identifiers `<barcode>-<tech_id>`

The `AnnData` is saved as a zarr file for a better speed to compression tradeoff compared ot gzipped h5ad files.

## Merge Data

This operation is applied by the following rules

+ `merge_study`: merge all datasets of a study if multiple datasets are available, else create a symlink
+ `merge_organ`: merge all studies that belong to an organ
+ `merge_organ_filter`: merge all cells removed by filtering per organ
+ `merge_subset`: merge datasets by subset defined in the input TSV under `subset` (overlapping subsets allowed)

The `AnnData` must contain all the slots described in [Aggregate Metadata](#aggregate-metadata) apart from the `.uns` slot.

## Filter

Filter cells per study depending on the `config.yaml` specification.
Two keys are available for controlling the filtering behaviour, `filter_per_organ` specifies global filter paramters
for all datasets per organ and `filter_per_study` allows for study specific filter options.
An example of an organ-level filter specification is shown below:

```yaml
filter_per_organ:
  blood:
    cells_per_sample:
      min: 50
      max: 10000
    mito_pct: 30
    remove_by_column:
      dataset:
        - Lee2020_2
```

All organ-level filtering decisions are applied per study.
The `remove_by_column` key can include any columns that are available in the anndata objects per study.
Currently only `remove_by_column` is applied (other keys such as `cells_per_sample` or `mito_pct` are ignored), and cells are only flagged in `.obs['filtered']`, not removed.

An example of per study filters shows that only the studies that require further filtering need to be overwritten.
The filter options are the same as for the organ-level filters.

```yaml
filter_per_study:
  SchulteSchrepping2020:
    remove_by_column:
      sample:
        - Schulte-Schrepping_C2P01H_d0
        - Schulte-Schrepping_C2P05F_d0
        - Schulte-Schrepping_C2P07H_d0
        - Schulte-Schrepping_C2P10H_d0
        - Schulte-Schrepping_C2P13F_d0
        - Schulte-Schrepping_C2P15H_d0
        - Schulte-Schrepping_C2P16H_d0
        - Schulte-Schrepping_C2P19H_d0
      donor:
        - C19-CB-0008
      disease:
        - influenza
```

Both keys can be empty or missing from the config file.
In that case, no filtering is applied.

# Preparing the input data

In order to use the data loader module, you need to define the following files with the file locations and dataset-level metadata.

* Dataset definiton file (`configs/datasets.tsv`)
* Schema mapping (e.g. `configs/schema_mapping.tsv`)
* Configuration file (e.g. `configs/imported/config.yaml`)
* DCP metadata (optional)

Additionally, you need to prepare your input `AnnData` files to contain the metadata that is defined in your dataset definition file.

## Dataset definition file

The dataset definition file should specify which datasets you want to include for your analyses together with any additional dataset-level metadata.

| Column            | Description                                                                                                                                                                                                |
| ----------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| dataset           | Name of the dataset, multiple datasets make up a study                                                                                                                                                     |
| study             | Name of the study, used for aggregating datasets to study level                                                                                                                                            |
| organ             | Name of the organ, can be tissue or any other name for aggregating the atlas                                                                                                                               |
| donor_column      | Column with donor IDs                                      a donor is an individual who provided the sample                                                                                                |
| tech_id           | Column(s) that identify a technical sample/library, multiple columns joined by `+` (e.g. `sample+lane`); stored as `.obs['tech_id']` and `.obs['sample']` |
| keep_covariates   | Optional. Comma-separated additional `.obs` columns to keep in addition to the schema columns |
| author_annotation | Column with author annotations (needed for annotation quality assessment and label harmonisation)                                                                                                          |
| cell_type         | Optional. Column with cell ontology labels (needed to for different versions of the CELLxGENE schema, TODO: deprecate). If cell_type is missing from the file, it will be generated from author_annotation |
| schema            | Name of schema to be mapped to `cellxgene`. Naming must match the columns in schema mapping.                                                                                                               |
| url               | URL or path to the `h5ad` file. For data from CELLxGENE data: `cellxgene` when `collection_id` and `dataset_id` are defined. For HCA DCP data `dcp` when `project_uuid` is defined.                        |
| collection_id     | Only for `url=cellxgene`, if data should be downloaded directly from CELLxGENE.                                                                                                                            |
| dataset_id        | Only for `url=cellxgene`, if data should be downloaded directly from CELLxGENE.                                                                                                                            |
| project_uuid      | Only for `url=dcp`, if data should be downloaded directly from HCA DCP data portal                                                                                                                         |
| annotation_file   | Optional. Any additional annotations that are not in the `AnnData` object. Needs to have a matching barcode column.                                                                                        |
| barcode_column    | Optional. Column in `AnnData.obs` for merging external annotations                                                                                                                                         |

All other columns are optional and will be added to `AnnData.uns['meta']`.

## Schema Mapping

The data loader ensures that the data adheres to the [CELLxGENE schema 3.0.0](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/3.0.0/schema.md) specifications.
For datasets that do not adhere to that schema, the schema mapping file allows to provide a mapping of custom `AnnData.obs` columns to the ones defined in CELLxGENE.
Below is an example of a schema mapping for the schemas `custom` and `dcp` to `cellxgene`.
The column names are used to map the dataset to its corresponding schema.

```
cellxgene                custom                       dcp                                                            
study                    study_PI                     project.contributors.name                                      
sample                   sample_ID                    specimen_from_organism.biomaterial_core.biomaterial_id         
donor_id                 subject_ID                   donor_organism.biomaterial_core.biomaterial_id                 
development_stage        subject_developmental_state  donor_organism.development_stage.text                          
sex                      sex                          donor_organism.sex                                             
self_reported_ethnicity  ethnicity_free_text          donor_organism.human_specific.ethnicity.text                   
suspension_type          biological_unit              library_preparation_protocol.nucleic_acid_source               
assay                    library_platform             library_preparation_protocol.library_construction_method.text  
organism                 species                      donor_organism.genus_species.text                              
cell_type                cell_type                                                                                   
...                      ...                          ...                                                            
```


## DCP metadata (optional)

This file is optional and used for datasets for which users want to map additional HCA DCP metadata annotations.
The mapping should contain a `study` and a `filename` column, where `filename` is a TSV file that follows the [DCP metadata schema](https://data.humancellatlas.org/metadata).
The mapping does not have to include all the studies that you want to include in your analysis.

TODO: extend to other metadata input.

## Configuration file

The `config.yaml` file is the main configuration file of the pipeline and is included in the top-level workflow by default.
For data loading, you just need to define the location of the files defined above.
By default, `configs/imported/config.yaml` should already contain the correct paths to the previously described input files.
If you configured your datasets with different files than listed below, you need to update them accordingly.

```yaml
dataset_meta: configs/datasets.tsv
schema_file: data/input/schema_mapping.tsv
dcp_metadata: configs/dcp_metadata.tsv

filter_per_organ:
...

filter_per_study:
...
```

For more information on Snakemake configuration files, please refer to the [documentation](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html).

# Testing

The test configuration and command are under `test/`.
All paths in the following are relative to the module root directory.

## Prepare test data
Before running the test pipeline for the first time, you must download a test dataset.
The following script downloads the "SchulteSchrepping" dataset and then copies it to the location that is defined in `dataset.tsv`.

```
bash test/download_test_data.sh -c1
```

This needs to be done only once.

## Run pipeline on test configuration

Activate the snakemake environment and call `test/run_test_*.sh` with run specific Snakemake parameters.

```
conda activate snakemake
bash test/run_test_cellxgene.sh -n
bash test/run_test_cellxgene.sh -c
```
