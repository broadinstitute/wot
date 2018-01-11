# -*- coding: utf-8 -*-
import dask.array as da
import wot


def z_score(x, z_min=-5, z_max=5):
    mean = da.mean(x, axis=0)  # 1 by ngenes
    variance = da.var(x, axis=0)
    stdev = da.sqrt(variance)
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

    indices = gs_x.sum(axis=1) > 0
    gs_x = gs_x[indices, :]
    ds_x = ds_x[:, indices]
    if z_score_ds:
        ds_x = z_score(ds_x)

    scores = ds_x.dot(gs_x)
    ngenes_in_set = gs_x.sum(axis=0)
    scores = scores / ngenes_in_set
    return wot.Dataset(x=scores, row_meta=ds.row_meta, col_meta=gs.col_meta)
