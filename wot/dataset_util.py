# -*- coding: utf-8 -*-

import anndata
import math
import numpy as np
import pandas as pd
import scipy.sparse
import scipy.sparse


def list_of_days_in_dataset(dataset):
    if 'day' not in dataset.obs.columns:
        raise ValueError("No day information available for this dataset")
    return sorted(list(set(dataset.obs['day'].values)))


def cell_indices_by_day(dataset):
    """Returns a dictionary mapping each day with the list of indices of cells from that day"""
    day_to_indices = {}
    unique_days = list_of_days_in_dataset(dataset)
    for day in unique_days:
        day_query = dataset.obs['day'] == day
        indices = np.where(day_query)[0]
        day_to_indices[day] = indices
    return day_to_indices


def get_cells_in_gene_sets(gene_sets, dataset, quantile=.99):
    cell_sets = {}
    for gene_set_index in range(gene_sets.X.shape[1]):
        gene_indices = list(np.where(gene_sets.X[:, gene_set_index] == 1)[0])
        extracted = dataset.X[:, gene_indices]
        thresholds = np.percentile(extracted, axis=0, q=quantile * 100)
        selected = []
        for i in range(extracted.shape[0]):
            if all(extracted[i] > thresholds):
                selected.append(i)
        cell_sets[gene_sets.var.index[gene_set_index]] = \
            dataset.obs.index[selected]
    return cell_sets


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
        return anndata.AnnData(dataset.X[indices], dataset.obs.iloc[indices].copy(), dataset.var.copy())

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


def extract_cells_at_indices(ds, indices):
    return anndata.AnnData(ds.X[indices], ds.obs.iloc[indices].copy(), ds.var.copy())


def add_cell_metadata(dataset, name, data):
    dataset.obs[name] = data


def set_cell_metadata(dataset, name, data, indices=None):
    if indices is None:
        dataset.obs[name] = data
    else:
        if isinstance(indices, set) or isinstance(indices[0], str):
            dataset.obs.loc[indices, name] = data
        else:
            dataset.obs.loc[dataset.obs.index[indices], name] = data


def merge_datasets(*args):
    datasets = list(args)
    merged_x = np.concatenate([d.X for d in datasets])
    row_columns = set(datasets[0].obs.columns)
    if not all([set(d.obs.columns) == row_columns for d in datasets]):
        raise ValueError("Unable to merge: incompatible metadata between datasets")
    merged_row_meta = pd.concat([d.obs for d in datasets], sort=True)
    if merged_row_meta.index.duplicated().any():
        raise ValueError("Unable to merge: duplicate rows between datasets, cannot lose information")
    col_index = datasets[0].var.index
    if not all([d.var.index.equals(col_index) for d in datasets]):
        raise ValueError("Unable to merge: incompatible genes between datasets")
    merged_col_meta = datasets[0].var
    return anndata.AnnData(merged_x, merged_row_meta, merged_col_meta)


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
