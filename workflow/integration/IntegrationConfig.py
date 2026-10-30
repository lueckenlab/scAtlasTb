from pathlib import Path
import numpy as np
import pandas as pd
from snakemake.rules import Rule

from utils.ModuleConfig import ModuleConfig, get_use_gpu
from utils.config import _get_or_default_from_config
from utils.misc import expand_dict_and_serialize, unique_dataframe


class IntegrationConfig(ModuleConfig):

    def __init__(
        self,
        **kwargs
    ):
        """
        :param kwargs: parameters for ModuleConfig
        """
        # parse output type info
        parameters = kwargs.get('parameters')
        if isinstance(parameters, str):
            parameters = pd.read_table(parameters)
        if isinstance(parameters, pd.DataFrame):
            profile = 'gpu' if get_use_gpu(kwargs['config']) else 'cpu'
            if profile == 'cpu':
                mask = ~parameters['cpu_env'].isna()
                parameters.loc[mask, 'env'] = parameters.loc[mask, 'cpu_env']
                parameters['resources'] = 'cpu'

            # set output types
            parameters['output_type'] = parameters['output_type'].str.split(',')
            parameters['output_types'] = parameters['output_type']
            parameters = parameters.explode('output_type')

            kwargs['parameters'] = parameters
        
        kwargs['dont_inherit'] = ['methods']

        super().__init__(**kwargs, write_output_files=False)
        
        # set hyperparameters
        self.set_hyperparams()

        # remove redundant label wildcards
        wildcards_df = self.parameters.wildcards_df
        if 'use_cell_type' in wildcards_df.columns:
            wildcards_df['label'] = np.where(
                wildcards_df['use_cell_type'],
                wildcards_df['label'],
                'None'
            )
        
        # subset to user defined output types
        remove_indices = []
        for dataset in wildcards_df['dataset'].unique():
            output_type_col = wildcards_df.query('dataset == @dataset')['output_type']
            output_types = self.get_for_dataset(
                dataset=dataset,
                query=[self.module_name, 'output_types'],
                default=output_type_col.unique().tolist()
            )
            remove_indices.extend(output_type_col[~output_type_col.isin(output_types)].index)
        wildcards_df = wildcards_df.drop(remove_indices)
        
        # set default paramspace arguments
        if kwargs.get('paramspace_kwargs') is None:
            paramspace_kwargs = dict(
                filename_params=['method', 'hyperparams', 'label', 'output_type'],
                filename_sep='--',
            )
        else:
            paramspace_kwargs = {}
        self.update_parameters(
            wildcards_df=unique_dataframe(wildcards_df),
            wildcard_names=self.parameters.wildcard_names + ['hyperparams'],
            **paramspace_kwargs,
        )

        # write output files
        self.write_output_files()


    def set_hyperparams(self):
        """
        Set hyperparameters specific to each method of a module for all datasets

        :param config: config containing dataset specific information
        :param module: name of module, key must be present for each dataset entry
        """
        defaults = self.parameters.default_config.get(self.module_name)
        records = []
        for dataset, dataset_dict in self.parameters.dataset_config.items():
            methods_config = _get_or_default_from_config(
                config=dataset_dict,
                defaults=defaults,
                key=self.module_name,
                value='methods',
            )
            if methods_config is None:
                continue
            for method, hyperparams_dict in methods_config.items():
                if isinstance(hyperparams_dict, dict):
                    # make hyperparameters invariant to order of configuration
                    hyperparams_dict = dict(sorted(hyperparams_dict.items()))
                    
                    df = self.parameters.wildcards_df.query('method == @method and dataset == @dataset')
                    if df.empty:
                        continue
                    do_not_expand = df.iloc[0]['do_not_expand']
                    do_not_expand = do_not_expand.split(',') if isinstance(do_not_expand, str) else [do_not_expand]
                    records.extend(
                        (dataset, method, *rec)
                        for rec in expand_dict_and_serialize(hyperparams_dict, do_not_expand=do_not_expand)
                    )
                else:
                    records.append((dataset, method, str(hyperparams_dict), hyperparams_dict))
        self.hyperparams_df = pd.DataFrame(
            records,
            columns=['dataset', 'method', 'hyperparams', 'hyperparams_dict']
        )
        
        # write to file 
        if self.out_dir.exists():
            unique_dataframe(
                self.hyperparams_df[['method', 'hyperparams', 'hyperparams_dict']]
            ).to_csv(self.out_dir / 'hyperparams.tsv', sep='\t', index=False)
        
        # merge to wildcards
        self.parameters.wildcards_df = self.parameters.wildcards_df.merge(
            self.hyperparams_df,
            on=['dataset', 'method'],
            how='left'
        )


    def get_hyperparams(self):
        return self.hyperparams_df
