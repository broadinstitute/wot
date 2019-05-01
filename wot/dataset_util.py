# -*- coding: utf-8 -*-

import math

import anndata
import numpy as np
import pandas as pd
import scipy.sparse
import scipy.sparse


def split_anndata(dataset, metadata):
    """
    Split AnnData into sub-datasets according to a metadata

    Parameters
    ----------
    metadata : str
        The metadata to use for the split

    Returns
    -------
    splits : dict of t: anndata.AnnData
        Dictionnary of datasets. t is the type of the 'metadata' column.
        Each cell in splits[k] has its 'metadata' column constant to k.

    Raises
    ------
    ValueError
        If the metadata is not present
    """

    if metadata not in dataset.obs.columns:
        raise ValueError("Cannot split on '{}' : column not present".format(metadata))

    def extract(group):
        indices = dataset.obs.index.get_indexer_for(group.index)
        return dataset[indices]

    return {name: extract(group) for name, group in dataset.obs.groupby(metadata)}


def mean_and_variance(x):
    if scipy.sparse.issparse(x):
        mean = x.mean(axis=0)
        mean_sq = x.multiply(x).mean(axis=0)
        mean = np.asarray(mean)
        mean_sq = np.asarray(mean_sq)
        var = (mean_sq - mean ** 2)
        return mean, var
    else:
        return x.mean(axis=0), x.var(axis=0)


def dataset_from_x(x, rows=None, columns=None,
                   row_prefix="cell_", column_prefix="gene_"):
    x = np.asarray(x, dtype=np.float64)
    if x.ndim == 1:
        x = np.asarray([x]).T
    if rows is None:
        row_count_len = math.floor(math.log10(x.shape[0])) + 1
        rows = ["{}{:0{}}".format(row_prefix, i, row_count_len) for i in range(x.shape[0])]
    if columns is None:
        col_count_len = math.floor(math.log10(x.shape[1])) + 1
        columns = ["{}{:0{}}".format(column_prefix, i, col_count_len) for i in range(x.shape[1])]
    return anndata.AnnData(x,
                           pd.DataFrame([], index=rows, columns=[]),
                           pd.DataFrame([], index=columns, columns=[])
                           )
