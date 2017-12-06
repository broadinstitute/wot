import numpy as np
import wot
import scipy


def score_gene_sets(ds, gs, z_score=True):
    intersect = np.intersect1d(ds.col_meta.index, gs.row_meta.index,
                               assume_unique=True)
    # subset dataset and gene sets to include only genes found in both and
    # genes in dataset that have at least one non-zero entry
    gs_indices = np.zeros(len(intersect), dtype=np.uint32)
    ds_indices = np.zeros(len(intersect), dtype=np.uint32)
    for i in range(len(intersect)):
        gs_indices[i] = gs.row_meta.index.get_loc(intersect[i])
        ds_indices[i] = ds.col_meta.index.get_loc(intersect[i])

    # gene sets has genes on rows, sets on columns
    gs_x = gs.x[gs_indices, :]
    # ds has cells on rows, genes on columns
    ds_x = ds.x[:, ds_indices]

    # exclude genes that are all zeros
    sums = ds_x.sum(axis=0)
    not_zero = np.where(sums > 0)
    gs_x = gs_x[not_zero[1], :]
    ds_x = ds_x[:, not_zero[1]]

    # Compute the z-score for each gene in the set. Truncate these z-scores at 5
    # or −5, and define the signature of the cell to be the mean z-score over
    # all genes in the gene set.
    if z_score:
        x = ds_x.todense() if scipy.sparse.issparse(ds_x) else ds_x
        mean = np.mean(x, axis=0)  # 1 by ngenes
        variance = np.var(x, axis=0)
        stdev = np.sqrt(variance)
        x = (x - mean) / stdev
        x[np.where(x > 5)] = 5
        x[np.where(x < -5)] = -5
    else:
        x = ds_x

    scores = x.dot(gs_x.todense() if not scipy.sparse.issparse(x) else gs_x)
    ngenes_in_set = gs_x.sum(axis=0)
    scores = scores / ngenes_in_set
    return wot.Dataset(x=scores, row_meta=ds.row_meta, col_meta=gs.col_meta)
