# -*- coding: utf-8 -*-

import wot.io
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
