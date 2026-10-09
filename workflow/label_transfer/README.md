# Label transfer

```yaml
DATASETS:
  test:
    input:
      label_transfer:
        file_1: test/input/pbmc68k.h5ad
        file_2: test/input/pbmc68k.h5ad
    label_transfer:
      majority_reference:
        reference_key: bulk_labels
        query_key: louvain
```

* `majority_reference`: Majority voting to assign reference labels to query clusters.
  * `reference_key`: The key in the reference dataset that contains the labels to be transferred.
  * `query_key`: The key in the query dataset that contains the clusters that will be assigned.
  * `crosstab_kwargs` (optional): keyword arguments passed to `pandas.crosstab` when counting reference labels per query group (e.g. `dropna: false`).
