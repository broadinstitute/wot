# -*- coding: utf-8 -*-

import wot
import scipy.sparse


# Compute the z-score for each gene in the set. Truncate these z-scores at 5
# or −5, and define the signature of the cell to be the mean z-score over
# all genes in the gene set.
# ds and gs must be in the same order
def score_gene_sets(ds, gs, z_score_ds=True, use_dask=False):
    if use_dask:
        import dask.array as np
    else:
        import numpy as np
    # gene sets has genes on rows, sets on columns
    # ds has cells on rows, genes on columns
    gs_x = gs.x
    ds_x = ds.x
    gene_indices = (gs_x.sum(axis=1) > 0)  # keep genes that are in gene sets
    gs_x = gs_x[gene_indices]
    ds_x = ds_x[:, gene_indices]
    if z_score_ds:
        ds_x = ds_x.toarray() if scipy.sparse.isspmatrix(ds_x) else ds_x
        std = np.std(ds_x, axis=0)
        if not use_dask:
            gene_indices = (std > 0)  # keep genes that have std > 0
            gs_x = gs_x[gene_indices]
            ds_x = ds_x[:, gene_indices]
            std = std[gene_indices]
        mean = np.mean(ds_x, axis=0)
        ds_x = (ds_x - mean) / std
        ds_x[ds_x < -5] = -5
        ds_x[ds_x > 5] = 5
        ds_x[ds_x == np.nan] = 0

    scores = ds_x.dot(gs_x)
    ngenes_in_set = gs_x.sum(axis=0)
    ngenes_in_set[ngenes_in_set == 0] = 1  # avoid divide by zero
    scores = scores / ngenes_in_set  # scores contains cells on rows, gene sets on columns
    return wot.Dataset(x=scores, row_meta=ds.row_meta, col_meta=gs.col_meta)
