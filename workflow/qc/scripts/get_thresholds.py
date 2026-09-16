import sctk
import numpy as np
import pandas as pd
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from utils.io import read_anndata, write_zarr_linked
from qc_utils import get_thresholds, apply_thresholds, parse_autoqc, update_thresholds


def thresholds_to_df(df, wildcards, qc_type=None, **kwargs):
    if qc_type is None: 
        if 'user_thresholds' in kwargs:
            qc_type = 'user'
            if 'autoqc_thresholds' in kwargs:
                qc_type = 'updated'
        elif 'autoqc_thresholds' in kwargs:
            qc_type = 'sctk_autoqc'
        else:
            qc_type = str(qc_type)
    
    # get thresholds
    thresholds = get_thresholds(**kwargs, transform=False)
    qc_thresholds = {
        key: (thresholds[f'{key}_min'], thresholds[f'{key}_max'])
        for key in kwargs['threshold_keys']
    }
    has_defined_thresholds = any(
        not (pd.isna(lower) and pd.isna(upper))
        for lower, upper in qc_thresholds.values()
    )
    
    # sort threshold columns
    threshold_keys = sorted(list(thresholds.keys()), reverse=True)
    thresholds = {key: thresholds[key] for key in threshold_keys}
    
    if df.shape[0] == 0:
        thresholds = {key: np.nan for key in thresholds}
        thresholds = dict(
            dataset=wildcards.dataset,
            file_id=wildcards.file_id,
            threshold_type=qc_type,
            passed_frac=np.nan,
            removed_frac=np.nan,
            n_passed=np.nan,
            n_removed=np.nan,
            n_total=np.nan,
        ) | thresholds
        return pd.DataFrame(thresholds, index=[0])

    if not has_defined_thresholds:
        thresholds = dict(
            dataset=wildcards.dataset,
            file_id=wildcards.file_id,
            threshold_type=qc_type,
            passed_frac=np.nan,
            removed_frac=np.nan,
            n_passed=np.nan,
            n_removed=np.nan,
            n_total=np.nan,
        ) | thresholds
        return pd.DataFrame(thresholds, index=[0])
    
    passed_qc = apply_thresholds(
        df,
        thresholds=qc_thresholds,
        threshold_keys=kwargs['threshold_keys'],
        column_name='passed_qc',
    )
    passed_qc = pd.Categorical(df['passed_qc'], categories=[True, False]).value_counts()
    
    thresholds = dict(
        dataset=wildcards.dataset,
        file_id=wildcards.file_id,
        threshold_type=qc_type,
        passed_frac=passed_qc[True] / passed_qc.sum(),
        removed_frac=passed_qc[False] / passed_qc.sum(),
        n_passed=passed_qc[True],
        n_removed=passed_qc[False],
        n_total=passed_qc.sum(),
    ) | thresholds
    
    return pd.DataFrame(thresholds, index=[0])


input_file = snakemake.input[0]
output_file = snakemake.output[0]
output_tsv = snakemake.output.tsv
output_qc_stats = snakemake.output.qc_stats
metrics_params_file = snakemake.input.get('metrics_params')

adata = read_anndata(input_file, obs='obs', uns='uns', verbose=False)

threshold_keys = snakemake.params.get('scautoqc_metrics')
user_thresholds = snakemake.params.get('thresholds')
alternative_thresholds = snakemake.params.get('alternative_thresholds')
has_alternative_thresholds = bool(alternative_thresholds)
autoqc_thresholds = adata.uns.get('scautoqc_ranges')
if autoqc_thresholds is None:
    autoqc_thresholds = pd.DataFrame()

# Add metric_params to autoqc_thresholds
metrics_params = sctk.default_metric_params_df
if metrics_params_file:
    user_params = pd.read_table(metrics_params_file, index_col=0)
    metrics_params.update(user_params)
print(f"autoqc_thresholds (before)=\n{autoqc_thresholds}")
autoqc_thresholds = autoqc_thresholds.merge(
    metrics_params,
    left_index=True,
    right_index=True,
    how='left',
)
autoqc_thresholds = autoqc_thresholds
threshold_keys = [key for key in threshold_keys if key in adata.obs.columns]
print(f"{threshold_keys=}")
print(f"autoqc_thresholds (after)=\n{autoqc_thresholds}")
if has_alternative_thresholds:
    print(f"{alternative_thresholds=}")
adata.uns['scautoqc_ranges'] = autoqc_thresholds

# Calculate threshold stats
threshold_frames = [
    # autoqc thresholds
    thresholds_to_df(
        df=adata.obs,
        wildcards=snakemake.wildcards,
        threshold_keys=threshold_keys,
        autoqc_thresholds=autoqc_thresholds,
        init_nan=True,
    ),
    # user thresholds
    thresholds_to_df(
        df=adata.obs,
        wildcards=snakemake.wildcards,
        threshold_keys=threshold_keys,
        user_thresholds=user_thresholds,
        init_nan=True,
    ),
    # updated thresholds
    thresholds_to_df(
        df=adata.obs,
        wildcards=snakemake.wildcards,
        threshold_keys=threshold_keys,
        autoqc_thresholds=autoqc_thresholds,
        user_thresholds=user_thresholds,
        init_nan=True,
    ),
]

if has_alternative_thresholds:
    threshold_frames.append(
        thresholds_to_df(
            df=adata.obs,
            wildcards=snakemake.wildcards,
            qc_type='alternative',
            threshold_keys=threshold_keys,
            user_thresholds=alternative_thresholds,
            init_nan=True,
        )
    )

df = pd.concat(threshold_frames)
adata.uns['qc'] = df

df.to_csv(output_tsv, sep='\t', index=False)

apply_thresholds(
    adata,
    thresholds=get_thresholds(
        threshold_keys,
        # parse_autoqc to make sure it follows autoqc's "side"
        user_thresholds=parse_autoqc(autoqc_thresholds),
    ),
    threshold_keys=threshold_keys,
    column_name='sctk_qc_status',
)

# add QC status column to .obs with 'passed', 'failed', and 'ambiguous'
apply_thresholds(
    adata,
    thresholds=get_thresholds(
        threshold_keys,
        user_thresholds=user_thresholds,
        autoqc_thresholds=autoqc_thresholds,
    ),
    threshold_keys=threshold_keys,
    column_name='user_qc_status',
)

user_status = adata.obs['user_qc_status']

if has_alternative_thresholds:
    # set defaults for alternative thresholds
    updated_thresholds = update_thresholds(user_thresholds, autoqc_thresholds)
    alternative_thresholds = updated_thresholds | alternative_thresholds
    print(f"{alternative_thresholds=}")
    apply_thresholds(
        adata,
        thresholds=get_thresholds(
            threshold_keys,
            user_thresholds=alternative_thresholds,
            # autoqc_thresholds=autoqc_thresholds,
        ),
        threshold_keys=threshold_keys,
        column_name='alternative_qc_status',
    )

    # get 'passed' if user_qc_status and alternative_qc_status are both True
    alt_status = adata.obs['alternative_qc_status']
    adata.obs['qc_status'] = 'ambiguous'
    adata.obs.loc[user_status & alt_status, 'qc_status'] = 'passed'
    adata.obs.loc[~(user_status | alt_status), 'qc_status'] = 'failed'
else:
    adata.obs['qc_status'] = np.where(user_status, 'passed', 'failed')

# convert QC status to ordered categorical
adata.obs['qc_status'] = pd.Categorical(
    adata.obs['qc_status'],
    categories=['passed', 'failed', 'ambiguous'],
    ordered=True
)

print("calculate QC stats")
qc_status_counts = adata.obs['qc_status'].value_counts()
qc_status_counts = pd.DataFrame(qc_status_counts).T
qc_status_counts['file_id'] = snakemake.wildcards.file_id
qc_status_counts = qc_status_counts[['file_id', 'passed', 'failed', 'ambiguous']]

qc_status_counts.to_csv(output_qc_stats, sep='\t', index=False)

write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=['obs', 'uns'],
    verbose=False,
)
