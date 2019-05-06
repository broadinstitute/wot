# -*- coding: utf-8 -*-

import wot

from .ot_model import *


def initialize_ot_model(matrix, **kwargs):
    """
    Initializes an OTModel from a list of files.

    Parameters
    ----------
    matrix : str
        Path to a gene expression matrix file.
    **kwargs : dict
        Other keywords arguments, will be passed to OT configuration.

    Returns
    -------
    ot_model : wot.OTModel
        The OTModel with all data from the input files available.

    Example
    -------
    >>> # Basic initialization
    >>> initialize_ot_model('matrix.txt', 'days.txt')

    >>> # Tweaking unbalanced parameters
    >>> initialize_ot_model('matrix.txt', 'days.txt', lambda1=50, lambda2=80, epsilon=.01)
    """
    ds = wot.io.read_dataset(matrix)
    if kwargs.pop('transpose', False):
        ds = ds.T

    wot.io.add_row_metadata_to_dataset(dataset=ds, days=kwargs.pop('cell_days', None),
                                       growth_rates=kwargs.pop('cell_growth_rates', None),
                                       covariate=kwargs.pop('covariate', None))
    return OTModel(ds, **kwargs)


def parse_configuration(config):
    """
    Build a valid configuration from the argument.

    Parameters
    ----------
    config : None, str, pandas.DataFrame, or dict
        None to always use default, or a path to configuration, or parseable configuration

    Returns
    -------
    parse_config : None, or dict of (float, float): dict
        The configuration for each timepair as a dictionnary.

    Raises
    ------
    ValueError
        If the input cannot be converted to a valid configuration.

    Notes
    -----
    When passing a DataFrame, it must have a column 't', or at least columns 't0' and 't1'
    When passing a path to a file, the same constraints apply.
    When passing a dict, all keys must be castable to float or all castable to (float, float).
    """
    if config is None:
        return None
    elif isinstance(config, str):
        return parse_configuration(wot.io.read_day_pairs(config))
    elif isinstance(config, pd.DataFrame):
        if 't' in config.columns:
            return parse_per_timepoint_configuration(config)
        elif 't0' in config.columns and 't1' in config.columns:
            return parse_per_timepair_configuration(config)
        else:
            raise ValueError("Configuration must have at least a column 't' or two columns 't0' and 't1'")
    elif isinstance(config, dict):
        raise ValueError("Not implemented")
    else:
        raise ValueError("Unrecognized argument type for configuration. Use DataFrame, dict, str or None")


def parse_per_timepoint_configuration(config):
    """
    Parses per-timepoint configurations to per-timepair OTModel-compatible.

    Parameters
    ----------
    config : pandas.DataFrame or dict of str: dict
        The configuration to parse.

    Returns
    -------
    parsed_config : dict of (float, float): dict
        The configuration for each timepair as a dictionnary.

    Raises
    ------
    ValueError
        If the input cannot be coverted to a valid configuration

    Notes
    -----
    When passing a DataFrame, the column t must be present to indicate timepoints.
    When passing a dictionnary, all keys must be castable to float, and values must be dictionnaries.
    """
    if isinstance(config, pd.DataFrame):
        if 't' not in config.columns:
            raise ValueError("Invalid per-timepoint configuration : must have column t")
        types = {'t': float, 'epsilon': float, 'lambda1': float, 'lambda2': float}
        config = config.sort_values(by='t').astype({x: types[x] for x in types if x in config.columns})
        fields = [x for x in config.columns if x != 't' and x in types]
        day_pairs = {}
        for i in range(len(config) - 1):
            t0c = config.iloc[i]
            t1c = config.iloc[i + 1]
            day_pairs[(t0c['t'], t1c['t'])] = {x: (t0c[x] + t1c[x]) / 2 for x in fields}
        return day_pairs
    elif isinstance(config, dict):
        raise ValueError("Not implemented")
    else:
        raise ValueError("Unrecognized argument type for per-timepoint configuration. Use DataFrame, str or None")


def parse_parameter_file(path):
    df = pd.read_csv(path, engine='python', sep=None, header=None)
    #  two column file containing parameter and value
    result = {}
    for i in range(len(df)):
        result[df.iloc[i, 0]] = df.iloc[i, 1]
    return result


def parse_per_timepair_configuration(config):
    """
    Build a valid per-timepair configuration from the argument

    Parameters
    ----------
    config : pandas.DataFrame or dict of (float, float): dict
        The configuration to be checked or parsed

    Returns
    -------
    parsed_config : dict of (float, float): dict
        The configuration for each timepair as a dictionnary

    Raises
    ------
    ValueError
        If the input cannot be converted to a valid configuration.

    Notes
    -----
    When passing a DataFrame, the columns t0 and t1 must be present to indicate timepoints
    When passing a dict, all keys must be castable to (float, float), all values must be dictionnaries.
    """
    if isinstance(config, pd.DataFrame):
        if 't0' not in config.columns or 't1' not in config.columns:
            raise ValueError("Invalid per-timepair configuration : must have columns t0 and t1")
        result = {}
        for i in range(len(config)):
            h = config.loc[i].to_dict()
            t0, t1 = h.pop('t0'), h.pop('t1')
            result[(t0, t1)] = h
        return parse_per_timepair_configuration(result)
    elif isinstance(config, dict):
        try:
            if not all(isinstance(z, (int, float)) for x, y in config for z in [x, y]):
                raise ValueError("Dictionary keys for config must be pairs of floats")
        except:
            # `x, y in config` failed, so a key is not a pair (wrong unpack count)
            raise ValueError("Dictionnary keys for config must be pairs")
        # At this point, we know all keys are pairs of float-castable scalars
        valid_fields = ['epsilon', 'lambda1', 'lambda2']
        for key in config:
            if not isinstance(config[key], dict):
                raise ValueError("Dictionnary values for config must be dictionnaries")
            selected = config[key]
            # Drop all keys that are not valid configuration options for day pairs
            config[key] = {x: selected[x] for x in valid_fields if x in selected}
        return config
    else:
        raise ValueError("Unrecognized argument type for config. Use DataFrame or dict")
