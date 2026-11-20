# Exploration

Exploratory summaries of one or more input files: dataset statistics, cells per sample/donor, marker gene dot plots and barcode overlap between samples.

## Configuration

```yaml
DATASETS:
  test_exploration: # name the exploration task
    input:
      exploration:
        # file name to file path mapping for the files you want to explore
        Lee2020: test/input/load_data/harmonize_metadata/Lee2020.zarr
        SchulteSchrepping2020: test/input/load_data/harmonize_metadata/SchulteSchrepping2020.zarr
    exploration:
      marker_genes: marker_gene_key  # specify which marker gene set/mapping you want to use
      sample: sample # sample ID
      donor: donor # donor ID
      summary_columns:
        # define colors for the summary plot
        - disease

MARKER_GENES:
  marker_gene_key: # provide marker gene mapping definition
    T:
      - CD3D
      - CD4
      - CD8A
    Naive T:
      - CCR7
      - LEF1
      - CD27
```