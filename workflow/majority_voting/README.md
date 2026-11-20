# Majority voting

This module allows majority-voting per cell for multiple cell-assignments of labels.
The main computation is the mode per cell across multiple columns in the AnnData object.

```mermaid
flowchart TD
  in(["Input AnnData .obs"])
  sel["Select label columns<br/>regex match, at least 2 columns"]
  cats["Unify label categories across columns"]
  vote["Per cell: most frequent label across columns<br/>missing labels count as votes"]
  agr["Agreement = votes for winning label / number of columns<br/>all votes different: consensus NaN, agreement 0"]
  thr["Flag low agreement<br/>agreement <= threshold"]
  out(["Consensus label, agreement + low-agreement flag in .obs"])
  stats(["Summary table + plot<br/>fraction of low-agreement cells per label"])
  in --> sel
  sel --> cats
  cats --> vote
  vote --> agr
  agr --> thr
  thr --> out
  out --> stats
```

*Conceptual overview of the main steps of the module. See the [functional description](#functional-description) below for details.*

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

> Note: The main assumption of the `columns` is that the labels across columns are from the same (or mostly same) set of unique labels.

## Functional description

### Inputs

* **File formats:** `.h5ad` or `.zarr` (AnnData), configured under `input: majority_voting:` (see [architecture](https://scatlastb.readthedocs.io/en/latest/principles/architecture.html)). Only `.obs` is read.
* **`.obs` label columns:** selected by `columns`, a list of regular expressions matched with `re.fullmatch` against all `.obs` column names (order of `.obs`, duplicates removed). At least two matching columns are required, otherwise the job fails.
* **`threshold`:** agreement threshold for flagging low consensus. The rule passes the configured value without a default, so it must be set in the config (the script's fallback of `0.5` only applies when the parameter is absent).

### Processing steps

1. **`majority_voting_majority_voting`** (script `scripts/majority_voting.py`, environment `scanpy`).
   ```text
   columns = [c for c in obs.columns if any(re.fullmatch(p, c) for p in column_patterns)]
   assert len(columns) >= 2
   convert each column to categorical; warn if category sets differ between columns
   categories = union of all columns' categories
   codes = obs[columns] recoded to the shared categories (NaN -> -1)
   for each cell: mode = most frequent code (ties -> smallest code, i.e. NaN, then
                  alphabetically first category); count = frequency of the mode
   majority_consensus           = category of mode (NaN if mode is NaN)
   majority_consensus_agreement = count / n_columns
   if agreement <= 1 / n_columns:           # every vote different
       majority_consensus = NaN; agreement = 0
   majority_consensus_low_agreement = agreement <= threshold
   ```
   Missing labels count as votes for "NaN", so they lower the agreement of the other labels and can win the vote.
2. Summary statistics: counts of cells per (`majority_consensus`, `majority_consensus_agreement`) combination, with the total and fraction of cells per consensus label, written as TSV.
3. Plot: `pandas.crosstab(majority_consensus, majority_consensus_low_agreement, normalize='index')` drawn as a horizontal stacked bar chart (sorted by fraction of low-agreement cells).

### Outputs

Additionally to the consensus labels, the module evaluates how strongly the label columns agree on the consensus label, and produces summary files and plots describing agreement levels.

* `<output_dir>/majority_voting/dataset~<dataset>/file_id~<file_id>.zarr` — `.obs` with added `majority_consensus` (categorical), `majority_consensus_agreement` (float in [0, 1]) and `majority_consensus_low_agreement` (bool); input label columns are converted to categorical. All other slots are linked to the input.
* `<images>/majority_voting/dataset~<dataset>/file_id~<file_id>/majority_consensus_agreement.tsv` — columns `majority_consensus`, `majority_consensus_agreement`, `count`, `total_cells_majority_consensus_label`, `frac_cells_in_majority_consensus_label`.
* `<images>/majority_voting/dataset~<dataset>/file_id~<file_id>/majority_consensus_frac.png` — fraction of low-agreement cells per consensus label.

### Environments

* [`scanpy`](https://github.com/HCA-integration/scAtlasTb/blob/main/envs/scanpy.yaml) — the only environment used.
