# -*- coding: utf-8 -*-

import wot.io
from .core import *

def initialize_core(matrix, days, tmap_dir = '.', tmap_prefix = None):
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

    Returns
    -------
    core : wot.Core
        The Core with all data from the input files available.
    """
    ds = wot.io.read_dataset(matrix)
    wot.io.incorporate_days_information_in_dataset(ds, days)
    return Core(ds, tmap_dir, tmap_prefix)
