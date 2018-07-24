# -*- coding: utf-8 -*-

import wot.io
import pandas as pd
from .ot_model import *

def initialize_ot_model(matrix, days, tmap_dir = '.', tmap_prefix = None, **kwargs):
    """
    Initializes an OTModel from a list of files.

    Parameters
    ----------
    matrix : str
        Path to a gene expression matrix file.
    days : str
        Path to a days file for the matrix.
    tmap_dir : str, optional
        Path to the transport maps directory for the OTModel.
    tmap_prefix : str, optional
        Prefix for transport maps cached by the OTModel.
    **kwargs : dict
        Other keywords arguments, will be passed to OT configuration.

    Returns
    -------
    ot_model : wot.OTModel
        The OTModel with all data from the input files available.

    Example
    -------
    # Basic initialization
    initialize_ot_model('matrix.txt', 'days.txt')

    # Tweaking unbalanced parameters
    intialize_ot_model('matrix.txt', 'days.txt', lambda1=50, lambda2=80, epsilon=.01)
    """
    ds = wot.io.read_dataset(matrix)
    wot.io.incorporate_days_information_in_dataset(ds, days)
    ot_model = OTModel(ds, tmap_dir, tmap_prefix, **kwargs)
    return ot_model

def parse_day_pairs(day_pairs):
    """
    Build a valid day pairs configuration from the argument, or return None

    Parameters
    ----------
    day_pairs : pandas.DataFrame or dict or str or None
        Either a path to a file, a DataFrame, or None (to always use defaults)

    Returns
    -------
    parsed_day_pairs : dict of (float, float): dict or None
        The configuration for each day pair as a dictionnary, or None

    Raises
    ------
    ValueError
        If the input cannot be converted to a valid day pairs configuration.

    Notes
    -----
    When passing a DataFrame, the columns t0 and t1 must be present to indicate timepoints
    When passing a path to a file, the same columns must appear in the file.
    When passing a dict, all keys must be castable to (float, float), all values must be dictionnaries.
    """
    if day_pairs is None:
        return None
    elif isinstance(day_pairs, str):
        # The argument is a path. Parse it to a DataFrame and continue
        return parse_day_pairs(wot.io.read_day_pairs(day_pairs))
    elif isinstance(day_pairs, pd.DataFrame):
        if 't0' not in day_pairs.columns or 't1' not in day_pairs.columns:
            raise ValueError("Invalid day pairs : must have columns t0 and t1")
        result = {}
        for i in range(len(day_pairs)):
            h = day_pairs.loc[i].to_dict()
            t0, t1 = h.pop('t0'), h.pop('t1')
            result[(t0, t1)] = h
        return parse_day_pairs(result)
    elif isinstance(day_pairs, dict):
        try:
            if not all(isinstance(z, (int, float)) for x,y in day_pairs for z in [x,y]):
                raise ValueError("Dictionnary keys for day_pairs must be pairs of floats")
        except:
            # `x, y in day_pairs` failed, so a key is not a pair (wrong unpack count)
            raise ValueError("Dictionnary keys for day_pairs must be pairs")
        # At this point, we know all keys are pairs of float-castable scalars
        valid_fields = [ 'epsilon', 'lambda1', 'lambda2' ]
        for key in day_pairs:
            if not isinstance(day_pairs[key], dict):
                raise ValueError("Dictionnary values for day_pairs must be dictionnaries")
            config = day_pairs[key]
            # Drop all keys that are not valid configuration options for day pairs
            day_pairs[key] = { x: config[x] for x in valid_fields if x in config }
        return day_pairs
    else:
        raise ValueError("Unrecognized argument type for day_pairs. Use DataFrame, str or None")
