# -*- coding: utf-8 -*-

import math

import numpy
import pandas
import scipy.sparse
import scipy.stats

import wot


# Compute the z-score for each gene in the set. Truncate these z-scores at 5
# or −5, and define the signature of the cell to be the mean z-score over
# all genes in the gene set.
# ds and gs must be in the same order
def score_gene_sets(dataset_to_score, gs, background_ds=None, method='mean_z_score', use_dask=False):
    if use_dask:
        import dask.array as np
    else:
        import numpy as np
    # gene sets has genes on rows, sets on columns
    # ds has cells on rows, genes on columns
    gs_x = gs.x
    if background_ds is None:
        background_ds = dataset_to_score

    if method == 'mean_z_score':
        background_x = background_ds.x.toarray() if scipy.sparse.isspmatrix(background_ds.x) else background_ds.x
        gene_indices = (gs_x.sum(axis=1) > 0) & (
                background_x.std(axis=0) > 0)  # keep genes that are in gene sets and have standard deviation > 0
        gs_x = gs_x[gene_indices]
        background_x = background_x[:, gene_indices]
        dataset_to_score.x = dataset_to_score.x[:, gene_indices]
        std = np.std(background_x, axis=0)
        mean = np.mean(background_x, axis=0)

        dataset_to_score.x = (dataset_to_score.x - mean) / std
        dataset_to_score.x[dataset_to_score.x < -5] = -5
        dataset_to_score.x[dataset_to_score.x > 5] = 5
        dataset_to_score.x[dataset_to_score.x == np.nan] = 0

        scores = dataset_to_score.x.dot(gs_x)
        ngenes_in_set = gs_x.sum(axis=0)
        ngenes_in_set[ngenes_in_set == 0] = 1  # avoid divide by zero
        scores = scores / ngenes_in_set  # scores contains cells on rows, gene sets on columns
        return wot.Dataset(x=scores, row_meta=dataset_to_score.row_meta, col_meta=gs.col_meta)
    elif method == 'mean_rank':
        ranks = np.zeros(dataset_to_score.x.shape)
        x = dataset_to_score.x.toarray() if scipy.sparse.isspmatrix(dataset_to_score.x) else dataset_to_score.x
        for i in range(dataset_to_score.x.shape[0]):
            ranks[i] = scipy.stats.rankdata(x[i, :])
        scores = ranks.dot(gs_x)
        ngenes_in_set = gs_x.sum(axis=0)
        ngenes_in_set[ngenes_in_set == 0] = 1  # avoid divide by zero
        scores /= ngenes_in_set  # scores contains cells on rows, gene sets on columns
        scores /= dataset_to_score.x.shape[1]
        return wot.Dataset(x=scores, row_meta=dataset_to_score.row_meta, col_meta=gs.col_meta)


def list_of_days_in_dataset(dataset):
    if 'day' not in dataset.row_meta.columns:
        raise ValueError("No day information available for this dataset")
    return sorted(list(set(dataset.row_meta['day'].values)))


def cell_indices_by_day(dataset):
    """Returns a dictionary mapping each day with the list of indices of cells from that day"""
    day_to_indices = {}
    unique_days = list_of_days_in_dataset(dataset)
    for day in unique_days:
        day_query = dataset.row_meta['day'] == day
        indices = numpy.where(day_query)[0]
        day_to_indices[day] = indices
    return day_to_indices


def extract_cells_at_indices(ds, indices):
    return wot.Dataset(ds.x[indices], ds.row_meta.iloc[indices].copy(), ds.col_meta.copy())


def add_cell_metadata(dataset, name, data):
    dataset.row_meta[name] = data


def set_cell_metadata(dataset, name, data, indices=None):
    if indices is None:
        dataset.row_meta[name] = data
    else:
        if isinstance(indices, set):
            dataset.row_meta.loc[indices, name] = data
        else:
            dataset.row_meta.loc[dataset.row_meta.index[indices], name] = data


def merge_datasets(*args):
    datasets = list(args)
    merged_x = numpy.concatenate([d.x for d in datasets])
    row_columns = set(datasets[0].row_meta.columns)
    if not all([set(d.row_meta.columns) == row_columns for d in datasets]):
        raise ValueError("Unable to merge: incompatible metadata between datasets")
    merged_row_meta = pandas.concat([d.row_meta for d in datasets], sort=True)
    if merged_row_meta.index.duplicated().any():
        raise ValueError("Unable to merge: duplicate rows between datasets, cannot lose information")
    col_index = datasets[0].col_meta.index
    if not all([d.col_meta.index.equals(col_index) for d in datasets]):
        raise ValueError("Unable to merge: incompatible genes between datasets")
    merged_col_meta = datasets[0].col_meta
    return wot.Dataset(merged_x, merged_row_meta, merged_col_meta)


def dataset_from_x(x, rows=None, columns=None,
                   row_prefix="cell_", column_prefix="gene_"):
    if rows is None:
        row_count_len = math.floor(math.log10(x.shape[0])) + 1
        rows = ["{}{:0{}}".format(row_prefix, i, row_count_len) for i in range(x.shape[0])]
    if columns is None:
        col_count_len = math.floor(math.log10(x.shape[1])) + 1
        columns = ["{}{:0{}}".format(column_prefix, i, col_count_len) for i in range(x.shape[1])]
    return wot.Dataset(x,
                       pandas.DataFrame([], index=rows, columns=[]),
                       pandas.DataFrame([], index=columns, columns=[])
                       )
