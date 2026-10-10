# Majority voting

This module allows majority-voting per cell for multiple cell-assignments of labels.
The main computation is the mode per cell across multiple columns in the AnnData object.

## Environments

The following environment is needed:

- [`scanpy`](https://github.com/HCA-integration/scAtlasTb/blob/main/envs/scanpy.yaml)


## Configuration
```yaml
DATASETS:
  test:
    input:
      majority_voting:
        file_1: test/input/pbmc68k.h5ad
        file_2: test/input/pbmc68k.h5ad
    majority_voting:
      columns:
        - bulk_labels
        - na_.*
      threshold: 0.6
```

* `columns`: columns or patterns (parsable by [`re`](https://docs.python.org/3/library/re.html)) of columns in `.obs` to consider for majority voting.
* `threshold`: threshold to determine whether a label has low consensus agreement. If the majority agreement (#major_label/#columns) is lower than or equal to `threshold`, it is considered to have low consensus agreement. Must be set explicitly.

> Note: The main assumption of the `columns` is that the labels across columns are from the same (or mostly same) set of unique labels/


## Output

Additionally to the predicted labels, the module evaluates how strongly cells agree on their consensus label, and produces summary files and plots describing agreement levels.

* `<out_dir>/majority_voting/dataset~<datasets>/file_id~<file_id>.zarr`: Updated AnnData object with consensus label (`majority_consensus`) based on the most common value across selected columns, alongside a measure of how many columns agree (`majority_consensus_agreement`), and a Boolean flag for low agreement (`majority_consensus_low_agreement`). This table is written to a TSV file named majority_consensus_agreement.tsv in the output plots folder.

* `<image_dir>/majority_voting/dataset~<datasets>/file_id~<file_id>/majority_consensus_agreement.tsv`: Summary statistic of majority voting containing fractions of cells for every consensus label and agreement level.

* `<image_dir>/majority_voting/dataset~<datasets>/file_id~<file_id>/majority_consensus_frac.png`: Agreement bar plot showing the fraction of cells that have low agreement for each consensus label.
