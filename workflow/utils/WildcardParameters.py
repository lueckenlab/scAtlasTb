import logging
from collections.abc import Iterable
from collections import defaultdict
from pprint import pformat, pprint
import numpy as np
import pandas as pd
from snakemake.utils import Paramspace
from snakemake.io import Wildcards

from .config import _get_or_default_from_config
from .misc import unique_dataframe

logger = logging.getLogger('WildcardParameters')
logger.setLevel(logging.WARNING)


class WildcardParameters:
    
    def __init__(
        self,
        module_name: str,
        parameters: pd.DataFrame,
        input_file_wildcards: dict,
        dataset_config: dict,
        default_config: dict,
        default_datasets: list,
        dont_inherit: list,
        wildcard_names: list,
        mandatory_wildcards: list = None,
        config_params: list = None,
        rename_config_params: dict = None,
        explode_by: [str, list] = None,
        paramspace_kwargs: dict = None,
        dtypes: dict = None,
    ):
        """
        :param module_name: name of module
        :param parameters: dataframe with parameters for each dataset
        :param input_file_wildcards: dictionary with input file wildcards for each dataset
        :param dataset_config: dictionary with dataset configuration for module
        :param default_config: dictionary with default configuration for module
        :param default_datasets: list of default datasets from config
        :param wildcard_names: list of wildcard names for expanding rules
        :param config_params: list of parameters that a module should consider as wildcards, order and length must match wildcard_names, by default will take wildcard_names
        :param explode_by: column(s) to explode wildcard_names extracted from config by
        :param paramspace_kwargs: additional arguments for snakemake.utils.Paramspace
        :param dtypes: dictionary with dtypes for wildcard columns
        """
        self.module_name = module_name
        
        # determine wildcards
        self.dataset_config = dataset_config
        self.default_config = default_config
        self.default_datasets = default_datasets
        self.set_wildcards(
            config_params=config_params,
            wildcard_names=wildcard_names,
            explode_by=explode_by,
            rename_config_params=rename_config_params,
            dont_inherit=dont_inherit,
            dtypes=dtypes,
        )
        if mandatory_wildcards is None:
            mandatory_wildcards = ['dataset', 'file_id']
        for wildcard in mandatory_wildcards:
            while wildcard in wildcard_names:
                wildcard_names.remove(wildcard)
        self.wildcard_names = mandatory_wildcards + wildcard_names
        
        if parameters is not None:
            self.parameters_df = parameters
            if parameters.shape[1] > 0:
                self.wildcards_df = self.wildcards_df.merge(self.parameters_df, how='left')
        else:
            self.parameters_df = pd.DataFrame()
        
        # add input file wildcards
        self.wildcards_df = self.wildcards_df.merge(
            pd.DataFrame(input_file_wildcards, dtype=str),
            on='dataset',
            how='left',
        )
        
        # set wildcards to string
        self.wildcards_df.loc[:, self.wildcard_names] = self.wildcards_df.loc[:, self.wildcard_names].astype(str)

        # set paramspace
        if paramspace_kwargs is None:
            paramspace_kwargs = {}
        self.paramspace_kwargs = paramspace_kwargs
        # self.paramspace = Paramspace(
        #     self.wildcards_df[self.wildcard_names],
        #     **paramspace_kwargs
        # )


    def subset_by_query(
        self,
        query_dict: [dict, Wildcards],
        columns: list = None,
        verbose: bool = False,
    ) -> pd.DataFrame:
        """
        Subset ``wildcards_df`` using scalar, iterable, or operator filters.
        
        Scalar values use equality; lists, tuples, and sets select rows whose
        values are members of the collection. A one-entry operator dictionary
        can specify ``eq``, ``ne``, ``in``, or ``notin`` (also ``==``, ``!=``,
        ``isin``, ``not_in``, ``nin``, and ``!in`` respectively). Equality
        and inequality follow pandas `Series.eq
        <https://pandas.pydata.org/docs/reference/api/pandas.Series.eq.html>`_
        and `Series.ne
        <https://pandas.pydata.org/docs/reference/api/pandas.Series.ne.html>`_
        semantics. Membership follows `Series.isin
        <https://pandas.pydata.org/docs/reference/api/pandas.Series.isin.html>`_
        semantics.

        Examples::

            # Keep rows where the dataset is "study_a" (scalar equality).
            params.subset_by_query({"dataset": "study_a"})

            # Select any of several values (iterable membership).
            params.subset_by_query({"batch": ["batch_1", "batch_2"]})

            # Select with explicit equality, inequality, or membership.
            params.subset_by_query({"n_vars": {"eq": 2000}})
            params.subset_by_query({"batch": {"ne": "batch_1"}})
            params.subset_by_query({"dataset": {"in": ["study_a", "study_b"]}})
            params.subset_by_query({"batch": {"notin": ["batch_1", "batch_3"]}})

        :param query_dict: Mapping of column names to scalar values, iterable
            membership values, or one-entry operator dictionaries. Supported
            operators are ``eq``, ``ne``, ``in``, and ``notin`` and their
            aliases listed above.
        :param columns: columns to subset wildcards_df by
        :return: subset of wildcards_df
        """
        if not query_dict and not columns:
            return self.wildcards_df
        
        if columns is None:
            columns = self.wildcards_df.columns
        
        df = self.wildcards_df.copy()

        def _as_iterable(x):
            if isinstance(x, (list, tuple, set, pd.Index, np.ndarray)):
                return x
            return [x]

        def _apply_operator_filter(_df, _key, _op, _val):
            if _op in ('eq', '=='):
                return _df[_df[_key] == _val]
            if _op in ('ne', '!='):
                return _df[_df[_key] != _val]
            if _op in ('in', 'isin'):
                return _df[_df[_key].isin(_as_iterable(_val))]
            if _op in ('notin', 'not_in', 'nin', '!in'):
                return _df[~_df[_key].isin(_as_iterable(_val))]
            raise ValueError(
                f'Unsupported operator "{_op}" for query key "{_key}". '
                'Use one of: eq, ne, in, notin.'
            )

        # subset by query
        for key, value in query_dict.items():
            if key not in df.columns:
                continue
            if isinstance(value, dict):
                if len(value) != 1:
                    raise ValueError(
                        f'Operator query for key "{key}" must have exactly one operator. '
                        f'Got: {list(value.keys())}'
                    )
                op, op_value = next(iter(value.items()))
                df = _apply_operator_filter(df, key, op, op_value)
            elif isinstance(value, (list, tuple, set)):
                df = df[df[key].isin(value)]
            else:
                df = df[df[key] == value]
            
            if verbose:
                print(f'subset by {key} ==  {value}')
                print(df.shape)
        
        if verbose:
            print(df.transpose())
            print(df.dtypes)
        
        # subset by columns
        return unique_dataframe(df[columns]).reset_index(drop=True)


    def set_wildcards(
        self,
        config_params: list = None,
        wildcard_names: list = None,
        explode_by: list = None,
        config_entries: list = None,
        rename_config_params: dict = None,
        dtypes: dict = None,
        dont_inherit: list = None,
        warn: bool = False,
    ):
        """
        Collect wildcards and parameters from an configuration instance (e.g. dataset) for a given module.
        This function assumes that the keys of the given config keys are the different instances that contain specific parameters for different modules.

        :param config_params: List of parameters for each config entry of a module
            e.g. ['integration', 'label', 'batch']
        :param wildcard_names: names of wildcards to be extracted.
            Must map to config keys, and prepended by a wildcard name for the config entries
            e.g. ['dataset', 'method', 'label', 'batch']
        :param explode_by: column to explode by, expecting list entry for that column
        :param config_entries: list of entries to subset the config by, otherwise use all keys
        """
        if dtypes is None:
            dtypes = {}
        
        if not wildcard_names:
            wildcard_names = []
        
        if not config_params:  # len(self.dataset_config) == 0 or 
            config_params = wildcard_names
        
        # if len(config_params) != len(wildcard_names):
        #     raise ValueError('config_params and wildcard_names must be of same length')
            
        if config_entries is None:
            config_entries = self.dataset_config.keys()

        # collect entries for dataframe
        defaults = self.default_config.get(self.module_name, {})
        records = [
            (
                key,
                *[
                    _get_or_default_from_config(
                        config=self.dataset_config.get(key, {}),
                        defaults=defaults,
                        key=self.module_name,
                        value=param,
                        dont_inherit=dont_inherit,
                        warn=warn,
                    )
                    for param in config_params
                ]
            )
            for key in config_entries
        ]
        # create dataframe
        columns = [*['dataset']+config_params]
        df = pd.DataFrame.from_records(records, columns=columns)
        
        # rename columns
        if rename_config_params is None:
            rename_config_params = {}
        rename_config_params = {k: v for k, v in rename_config_params.items() if k in df.columns}
        df = df.rename(columns=rename_config_params)
        
        # deal with duplicate columns after renaming columns
        duplicated_columns = df.columns[df.columns.duplicated(keep=False)].unique().tolist()
        if len(duplicated_columns) > 0:
            reverse_map = defaultdict(list)
            for k, v in rename_config_params.items():
                if v not in duplicated_columns:
                    continue
                reverse_map[v].append(k)
            logger.warn(f'WARNING: Duplicated columns: {dict(reverse_map)}')
        for col in duplicated_columns:
            merged_col = df[duplicated_columns].bfill(axis=1).iloc[:, 0]
            del df[col]
            df[col] = merged_col
        
        # explode columns
        if explode_by is not None:
            explode_by = [explode_by] if isinstance(explode_by, str) else explode_by
            for column in explode_by:
                df = df.explode(column)
        
        # infer and set dtypes
        dtypes = df.convert_dtypes().dtypes.to_dict() | dtypes
        
        def get_default_value(x, dtype):
            if pd.api.types.is_bool_dtype(x) or dtype in (bool, np.bool_):
                return False
            if pd.api.types.is_integer(x) or dtype in (int, np.int32, np.int64):
                return 0
            if pd.api.types.is_float(x) or dtype in (float, np.float32, np.float64):
                return 0.0
            return None

        na_map = {
            col: get_default_value(df[col], dtype)
            for col, dtype in dtypes.items()
        }
        na_map = {k: v for k, v in na_map.items() if v is not None}
        df = df.fillna(value=na_map).astype(dtypes)
        
        if df.empty:
            self.wildcards_df = df
            return
        
        # set dtypes
        df = df.replace({np.nan: None})
        # for i, v in df.items():
        #     if isinstance(v[0], (list, dict)):
        #         continue
        #     else:
        #         df[i] = df[i].astype(str)
        self.wildcards_df = df.reset_index(drop=True)


    def update(
        self,
        wildcards_df: pd.DataFrame = None,
        parameters_df: pd.DataFrame = None,
        wildcard_names: list = None,
        **kwargs
    ):
        """
        :param wildcards_df: dataframe with updated wildcards
        :param parameters_df: dataframe with updated parameters
        :param wildcard_names: list of wildcard names to subset the paramspace by
        :param kwargs: additional arguments for snakemake.utils.Paramspace
        """
        if wildcards_df is not None:
            self.wildcards_df = unique_dataframe(wildcards_df)
        if parameters_df is not None:
            self.parameters_df = unique_dataframe(parameters_df)
        if wildcard_names is not None:
            self.wildcard_names = wildcard_names
        self.paramspace_kwargs = kwargs
        # self.paramspace = Paramspace(
        #     self.wildcards_df[self.wildcard_names],
        #     **kwargs
        # )


    def get_wildcards(
        self,
        subset_dict: [dict, Wildcards] = None,
        exclude: [list, str] = None,
        wildcard_names: list = None,
        all_params: bool = False,
        as_df: bool = False,
        default_datasets: bool = True,
        verbose: bool = False,
    ) -> [dict, pd.DataFrame]:
        """
        Retrieve wildcard instances as dictionary
        
        :param exclude: list of wildcard names to exclude
        :param subset_dict: dictionary with column (must be present in parameters_df) to value mapping
        :param wildcard_names: list of wildcard names to subset the wildcards by
        :param all_params: whether to include all parameters. If False (default), used defined wilcard names
        :param as_df: whether to return a dataframe instead of a dictionary
        :param default_datasets: whether to subset to default datasets (default: True)
        :return: dictionary of wildcards that can be applied directly for expanding target files
        """
        if exclude is None:
            exclude = []
        elif isinstance(exclude, str):
            exclude = [exclude]
        if subset_dict is None:
            subset_dict = {}
        if wildcard_names is None:
            wildcard_names = self.wildcards_df.columns if all_params else self.wildcard_names
            if verbose:
                print(f'wildcard_names: {wildcard_names}')
        
        if default_datasets:
            if self.default_datasets is not None:
                query_dict = {'dataset': self.default_datasets}
            else:
                query_dict = {}
        else:
            query_dict = {}
        query_dict |= subset_dict
        
        if verbose:
            print(f'query_dict: {query_dict}')
        
        df = self.subset_by_query(
            query_dict=query_dict,
            columns=[w for w in wildcard_names if w not in exclude]
        )
        
        if verbose:
            print(f'wildcards_df:\n{df}')
        
        return df if as_df else df.to_dict('list')


    def get_parameters(self) -> pd.DataFrame:
        return self.parameters_df


    def get_paramspace(
        self,
        wildcard_names: list = None,
        exclude: list = None,
        **kwargs
    ) -> Paramspace:
        """
        :param default: whether to return the default paramspace (all wildcards)
        :param wildcard_names: list of wildcard names to subset the paramspace by
        :param kwargs: additional arguments for snakemake.utils.Paramspace
        :return: snakemake.utils.Paramspace object
        """
        if exclude is None:
            exclude = []
        if wildcard_names is None:
            wildcard_names = [w for w in self.wildcard_names if w not in exclude]
        if not kwargs:
            kwargs = self.paramspace_kwargs
        try:
            paramspace = Paramspace(
                self.wildcards_df[wildcard_names],
                **kwargs
            )
        except KeyError as e:
            raise ValueError(f'{e}\nwildcard names: {wildcard_names}, kwargs: {kwargs}') from e
        return paramspace


    def get_from_parameters(
        self,
        query_dict: [dict, Wildcards],
        parameter_key: str,
        wildcards_sub: [list, None] = None,
        exclude: [list, str] = None,
        check_query_keys: bool = False,
        check_null: bool = False,
        default: [str, None] = None,
        single_value: bool = True,
        verbose: bool = False,
        as_type: type = None,
    ):
        """
        Get parameter value(s) from wildcards dataframe by querying with wildcards.
        
        Subsets the wildcards_df by query_dict conditions and retrieves the specified parameter.
        Supports single and multiple value returns, type conversion, and null checking.

        :param query_dict: Dictionary mapping wildcard column names to values for filtering
        :param parameter_key: Name of the column to retrieve values from
        :param wildcards_sub: List of wildcard columns to use for subsetting. Defaults to all columns.
        :param exclude: Wildcard name(s) to exclude from wildcards_sub
        :param check_query_keys: If True, verify all query_dict keys are in wildcards_sub
        :param check_null: If True, raise error if parameter is None/NaN/'None'
        :param default: Default value to return if parameter is None/NaN/'None'
        :param single_value: If True, return single value; if False, return list
        :param verbose: If True, print debug information
        :param as_type: Optional type to cast parameter to (e.g., int, str, list)
        :return: Parameter value(s) matching the query, or default if null
        """
        if wildcards_sub is None:
            wildcards_sub = self.wildcards_df.columns.tolist()
        if exclude is not None:
            exclude = [exclude] if isinstance(exclude, str) else exclude
            wildcards_sub = [w for w in wildcards_sub if w not in exclude]
        
        assert parameter_key in self.wildcards_df.columns, f'"{parameter_key}" not in wildcards_df.columns'
        if check_query_keys:  # redundant?
            for key in query_dict:
                assert key in wildcards_sub, f'Query key "{key}" is not in wildcards_sub'

        wildcards_sub = list(set(wildcards_sub + [parameter_key]))
        query_dict = {k: v for k, v in query_dict.items() if k in wildcards_sub}
        
        try:
            params_sub = self.subset_by_query(
                query_dict=query_dict,
                columns=[parameter_key],
                verbose=verbose,
            )
            assert not params_sub.empty, 'No wildcard combination found'
            if single_value:
                assert params_sub.shape[0] == 1, f'More than 1 row after subsetting'
        
        except AssertionError as e:
            columns = list(query_dict.keys())+[parameter_key]
            
            if params_sub.empty:
                params_sub = self.wildcards_df[columns]
            else:
                # params_sub = self.subset_by_query(query_dict, columns, verbose=False)
                params_sub = self.subset_by_query(query_dict, verbose=False)
                # remove columns that are the same for all rows
                n_unique_columns = params_sub.astype(str).nunique() > 1
                non_unique_columns = n_unique_columns[n_unique_columns].index.tolist()
                # subset to columns
                non_unique_columns = [col for col in non_unique_columns if col not in columns]
                params_sub = params_sub[columns + non_unique_columns]

            raise AssertionError(
                f'{e} for:\n\tparameter_key="{parameter_key}"\n\twildcards_sub={wildcards_sub}'
                f'\nquery_dict:\n{pformat(query_dict)}'
                f'\n{params_sub}'
                f'\n{params_sub.dtypes}'
                f'\nall columns: {self.wildcards_df.columns.tolist()}'
            ) from e
        
        parameter = params_sub[parameter_key].tolist()
        if single_value:
            parameter = parameter[0]
        
        # check if NULL
        if isinstance(parameter, Iterable):
            is_null = parameter == 'None'
        else:
            is_null = pd.isna(parameter) or pd.isnull(parameter) or parameter is None
        if is_null:
            if check_null:
                df = self.subset_by_query(
                    query_dict=query_dict,
                    columns=[parameter_key],
                    verbose=True,
                )
                print('query_dict:')
                pprint(query_dict)
                print('parameter_key:', parameter_key)
                raise ValueError(f'Parameter should not be null. parameter: {parameter_key}, value:{parameter}')
            else:
                parameter = default
        if as_type is not None:
            parameter = as_type(parameter)
        if verbose:
            print(f'parameter: {parameter}, {type(parameter)}')
        return parameter
