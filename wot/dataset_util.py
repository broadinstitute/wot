# -*- coding: utf-8 -*-
import numpy as np
import wot
import scipy.sparse


def z_score(x, z_min=-5, z_max=5):
    x = x.todense() if scipy.sparse.isspmatrix(x) else x
    mean = np.mean(x, axis=0)  # 1 by ngenes
    variance = np.var(x, axis=0)
    stdev = np.sqrt(variance)
    x = (x - mean) / stdev

    if z_min is not None:
        x[x < z_min] = z_min
    if z_max is not None:
        x[x > z_max] = z_max
    return x


# Compute the z-score for each gene in the set. Truncate these z-scores at 5
# or −5, and define the signature of the cell to be the mean z-score over
# all genes in the gene set.
# ds and gs must be in the same order
def score_gene_sets(ds, gs, z_score_ds=True):
    # gene sets has genes on rows, sets on columns
    # ds has cells on rows, genes on columns
    gs_x = gs.x
    ds_x = ds.x

    indices = gs_x.sum(axis=1) > 0 # keep genes that are in gene sets
    gs_x = gs_x[indices, :]
    ds_x = ds_x[:, indices]
    if z_score_ds:
        ds_x = ds_x.todense() if scipy.sparse.isspmatrix(ds_x) else ds_x
        mean = np.mean(ds_x, axis=0)
        variance = np.var(ds_x, axis=0)
        indices = variance > 0
        gs_x = gs_x[indices, :]
        ds_x = ds_x[:, indices]

        ds_x = (ds_x - mean) / np.sqrt(variance)
        ds_x[ds_x < -5] = -5
        ds_x[ds_x > 5] = 5

    scores = ds_x.dot(gs_x)
    ngenes_in_set = gs_x.sum(axis=0)
    scores = scores / ngenes_in_set  # scores contains cells on rows, gene sets on columns
    return wot.Dataset(x=scores, row_meta=ds.row_meta, col_meta=gs.col_meta)
