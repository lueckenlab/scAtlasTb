import ast
import logging

import pandas as pd


def read_threshold_file(file: str):
    """
    Read user threshold file in TSV format
    """
    df = pd.read_table(file, comment='#')
    assert 'file_id' in df.columns
    # keep a threshold column for every QC metric (e.g. percent_hb_max, percent_ribo_min),
    # not only the default metrics, so that no user threshold is silently dropped
    threshold_columns = [col for col in df.columns if col.endswith(('_min', '_max'))]
    if len(threshold_columns) == 0:
        logging.warning(f'WARNING: No QC threshold columns (<metric>_min or <metric>_max) found in {file}.')
    if 'threshold_type' in df.columns:
        df = df[df['threshold_type'].isin(['user', 'alternative'])]
    else:
        df['threshold_type'] = 'user'
    columns = ['file_id', 'threshold_type'] + threshold_columns
    return df[columns].drop_duplicates()


def unpack_thresholds(row: pd.Series):
    thresholds = row.thresholds
    if isinstance(thresholds, str):
        thresholds = ast.literal_eval(thresholds)
    if isinstance(thresholds, dict):
         thresholds = thresholds.get(row.file_id, thresholds)
    return thresholds
