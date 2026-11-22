# Batch Integration

This module provides different scRNA integration methods for batch correction.

Available methods (as defined in `params.tsv`): `unintegrated`, `bbknn`, `combat`, `harmony_pytorch`, `harmonypy`, `scanorama`, `scvi`, `scanvi`, `scvi_vitkl`, `drvi`, `sysvi`, `scpoli` and `scgen`.

## Configuration
Here is an example configuration with all available parameters.
Parameters that are optional can be ommitted in the config and the pipeline will use default values.

```yaml
DATASETS:

  dataset_name:

    # input file configuration
    input:
      integration:
        # file_id to file path mapping
        pbmc: test/input/pbmc68k.h5ad

    # module configuration
    integration:
      raw_counts: layers/counts  # optional
      norm_counts: layers/normcounts  # optional

      label:
        - bulk_labels
        - louvain

      batch:
        - phase
        - batch_2

      output_types:
        - full
        - embed
        - knn

      var_mask:  # optional
        - highly_variable
        - highly_variable_2

      neighbors:  # optional
        n_neighbors: 30

      # random seed used by each integration method, default: 0
      seed: 0

      # number of CPU threads to use for preparation and running integration methods
      threads: 3

      methods: # mandatory, needs at least 1 value

        unintegrated: # no hyperparameters defined, will use defaults

        bbknn:
          # define hyperparamters if desired
          neighbors_within_batch: 3

        combat:
          covariates:
            - bulk_labels

        scanorama:
          batch_size: 100

        # Seurat v5 RPCA, runs in R: raw_counts are log-normalized on all genes (Seurat's LogNormalize),
        # then subset to var_mask; see scripts/methods/seurat_rpca.py for all hyperparameters
        seurat_rpca:
          dims: 30
          k_weight: 100  # IntegrateLayers() argument k.weight (underscores stand for dots)
          reference:  # optional: batches the other batches are integrated onto (a list inside a list)
            -
            - [G1, S]

        scvi:
          max_epochs: 100
          early_stopping: true
          categorical_covariate_keys:
            - []
            - [batch_2, phase]
          continuous_covariate_keys:
            - percent_mito

        scgen:
          n_epochs: 10

        scpoli:
          cell_type_keys:
            - 
            - bulk_labels
          embedding_dims: 5
          recon_loss: nb
          n_epochs: 50
          early_stopping_kwargs:
            early_stopping_metric: val_prototype_loss
            mode: min
            threshold: 0
            patience: 20
            reduce_lr: true
            lr_patience: 13
            lr_factor: 0.1

        harmonypy:
          sigma: 0.1
          key: batch_2
          n_comps:
            - 10
            - 30

        harmony_pytorch:
          sigma: 0.1
          batch_key: batch_2
          n_comps: 30

        drvi:
          max_epochs: 100
          n_latent: 32

        sysvi:
          system_key: phase  # mandatory
          max_epochs: 100
```

### Module specification
Mandatory parameters that must be defined are the input file definition, `label`, `batch` and `methods`. All other paramters are optional. 
The AnnData file in h5ad or zarr requires the following input:

* `raw_counts`: matrix under `.X` or `.layers` containing counts that have not been normalized or log-transformed, needed by e.g. scVI, scANVI, scPoli. Must be set explicitly (use `X` if the counts are stored in `.X`).
* `norm_counts`: matrix under `.X` or `.layers` containing counts that have been normalized or log-transformed, needed by e.g. Scanorama, scGEN, unintegrated. Must be set explicitly (use `X` if the normalised counts are stored in `.X`).
* `output_types`: filters which methods and method modes to run. There are 3 different abstraction layers that an integration method can work on. Integration could either correct counts (`full`), a low-dimensional representation of the data (`embed`) or the kNN graph (`knn`). Some methods (e.g. Scanorama, unintegrated) provide multiple outputs, which can be utilised independently.
* `batch`, `label`: columns in `.obs` that are used to correct the batch effect or inform the integration method with cell type information (`label`, only for semi-supervised methods like scANVI, scPoli)
* `neighbors`: arguments that get passed to the [`sc.pp.neighbors`](https://scanpy.readthedocs.io/en/stable/generated/generated/scanpy.pp.neighbors.html/) function. The kNN is computed for `full` and `embed` outputs after integration
* `methods`: methods configuration to determine which methods should be used, as well as which hyperparameters those methods should use
  * hyperparamters correspond to the parameters that the integration method supports. For e.g. scvi-tools methods, you can define all parameters for the module setup and training functions at the same level
  * The module also supports parameter exploration. If you pass a list to a hyperparamter, all combinations of that hyperparameter with the other parameters will be computed as a separate hyperparameter computation
* `var_mask`: boolean column(s) in `.var` (e.g. `highly_variable`) defining the features used for integration. Each value becomes a separate run. The column must exist in the input file. If not set, all genes are used.
* `save_subset`: if `true`, the prepared file stores the count matrices physically subset to the integration features instead of linking the full matrices (default: `false`)
* `seed`, `threads`: random seed passed to the methods (default: 0) and number of CPU threads (default: 1)
* `umap_colors`, `plots` (`colors`, `plot_centroids`, `plot_gene_chunk_size`): additional `.obs` columns or genes to colour the UMAP plots by (label and batch columns are always plotted)

## Output

The output will be saved under `<output_dir>/integration/`.
There are mappings of inputs (`input_files.tsv`) and outputs (`output_files.tsv`) as well as the mapping of hyperparameters to hex code used in the output files.

### Metadata
+ `.uns['wildcards']` wildcards of the run, prefixed with `integration_` (e.g. `integration_method`, `integration_batch`)
+ `.uns['integration']` integration specific entries
    + `.uns['integration']['method']` intgration method name
    + `.uns['integration']['label_key']` label used for integration
    + `.uns['integration']['batch_key']` batch used for integration
    + `.uns['integration']['output_type']` output types of the method (list of 'knn', 'embed' and/or 'full')
    + `.uns['integration']['hyperparams']` hyperparameters used
    + `.uns['integration']['model_history']` training history (deep learning methods only)
+ `.uns['output_type']` output type of the processed file (processed output only)


### Integrated output
The direct method output is under `{out_dir}/integration/run_method/dataset~{dataset}/file_id~{file_id}/batch~{batch}/var_mask~{var_mask}/method~{method}--hyperparams~{hyperparams}--label~{label}/adata.zarr`
Different integration methods provide different types of output.
Depending on the output type, the object must contain the following slots:

1. corrected graph (`knn`)
   + `.obsp['connectivities']`, `.obsp['distances']` integrated kNN graph returned by integration method
2. corrected embedding (`embed`)
   + `.obsm['X_emb']` integrated embedding returned by integration method
3. corrected features (`full`)
   + `.X` corrected feature counts returned by integration method
   + `.obsm['X_pca']` PCA on corrected feature counts (computed during postprocessing)

### Processed integrated output
The integrated output files require further processing to be used in downstream analysis.
Files with computed embeddings (for full feature outputs) and kNN graphs (for full feature and embedding outputs) are stored under `{out_dir}/integration/dataset~{dataset}/file_id~{file_id}/batch~{batch}/var_mask~{var_mask}/method~{method}--hyperparams~{hyperparams}--label~{label}--output_type~{output_type}.zarr` (including the UMAP in `.obsm['X_umap']`).

## Contributing to the module

Adding new integration methods to the module simply requires adding a new script under `scripts/methods/` and adding meta information of the method in the `params.tsv`.
The name of the script must match the method name in the `params.tsv` and the method must be defined in the `params.tsv`, otherwise it won't be recognised by the pipeline.
When developing the script, you can make use of the global variables and functions from the `integration_utils.py` script.
This is particlarly useful if you want to manage different parameter assignments to e.g. train and model setup functions (see `scripts/methods/scvi.py` for an example).
The best way to get started, is to look at one of the scripts (ideally one that has a similar model to what you want to implement).

Additionally, you might want to add a new conda environment YAML file, if the new integration requires different dependencies than what is already provided in `envs`.
The new environment must then be specified in `params.tsv` and must be installed if you're using the default running mode of the pipeline.

Once you have implemented the new integration method, you can test if the the method gets recognised via a dry run. See below about testing.

## Testing
Activate the snakemake environment and call `test/run_test.sh` with run specific Snakemake parameters.

```
conda activate snakemake
bash test/run_test.sh -n  # dry run
bash test/run_test.sh -c2  # actual run with max 2 cores
```

The script will call the pipeline and run a test script.
If the input files don't yet exist, the test script might throw a `FileNotFoundError`.
You can ignore it, if you haven't yet executed the pipeline completely and therefore don't yet have all outputs.

### Add new test scenarios

You can configure a new test case in the `test/config.yaml`, which will be recognised in the bash command presented above.

### Test model loading
Some integration methods store pytorch model.
Below are tests if these models can be loaded.

```
conda run --live-stream -n scarches python test/run_scarches_model_loading.py 
conda run --live-stream -n scvi-tools python test/run_scvi-tools_model_loading.py 
```
This requires the `scarches` and `scvi-tools` conda environments to be installed, which are defined under `envs/`. 
