# -*- coding: utf-8 -*-

import wot.io
from .core import *

def initialize_core(matrix, days, tmap_dir = '.', tmap_prefix = None, **kwargs):
    """
    Initializes a core from a list of files.

    Parameters
    ----------
    matrix : str
        Path to a gene expression matrix file.
    days : str
        Path to a days file for the matrix.
    tmap_dir : str, optional
        Path to the transport maps directory for the Core.
    tmap_prefix : str, optional
        Prefix for transport maps cached by the Core.
    **kwargs : dict
        Other keywords arguments, will be passed to OT configuration.

    Returns
    -------
    core : wot.Core
        The Core with all data from the input files available.

    Example
    -------
    # Basic initialization
    initialize_core('matrix.txt', 'days.txt')

    # Tweaking unbalanced parameters
    intialize_core('matrix.txt', 'days.txt', lambda1=50, lambda2=80, epsilon=.01)
    """
    ds = wot.io.read_dataset(matrix)
    wot.io.incorporate_days_information_in_dataset(ds, days)
    core = Core(ds, tmap_dir, tmap_prefix)
    core.set_ot_config(**kwargs)
    return core
