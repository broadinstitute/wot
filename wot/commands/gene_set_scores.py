#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os

import wot.io
import numpy as np


def compute(matrix, gene_sets, out, format, cell_filter=None, background_cell_set=None,
            permutations=None, method='mean_z_score', nbins=None):
    use_dask = False
    ds = wot.io.read_dataset(matrix, use_dask=use_dask, chunks=(10000, None))
    if cell_filter is not None:
        cell_filter = wot.io.read_gene_sets(cell_filter)
        cells_ids = cell_filter.row_meta.index.values[np.where(cell_filter.x[:, 0] > 0)[0]]
        cell_filter = ds.row_meta.index.isin(cells_ids)
        ds = wot.Dataset(ds.x[cell_filter], ds.row_meta.iloc[cell_filter], ds.col_meta)

    background_ds = None
    if background_cell_set is not None:
        background_cells_ds = wot.io.read_gene_sets(background_cell_set)
        background_cells_ids = background_cells_ds.row_meta.index.values[np.where(background_cells_ds.x[:, 0] > 0)[0]]
        cell_filter = ds.row_meta.index.isin(background_cells_ids)
        background_ds = wot.Dataset(ds.x[cell_filter], ds.row_meta.iloc[cell_filter], ds.col_meta)

    gs = wot.io.read_gene_sets(gene_sets, ds.col_meta.index.values)
    if gs.x.shape[1] is 0:
        raise ValueError('No overlap of genes in gene sets and dataset')

    # scores contains cells on rows, gene sets on columns
    scores = []
    p_vals = []
    fdrs = []
    k = []
    n = []
    for j in range(gs.x.shape[1]):
        result = wot.score_gene_sets(dataset_to_score=ds, background_ds=background_ds,
                                     gs=wot.Dataset(gs.x[:, [j]], gs.row_meta, gs.col_meta.iloc[[j]]),
                                     permutations=permutations, method=method, nbins=nbins)
        scores.append(result['score'])
        if permutations is not None and permutations > 0:
            p_vals.append(result['p_value'])
            fdrs.append(result['fdr'])
            k.append(result['k'])
            n.append(result['n'])

    # import dask.array as da
    # da.to_npy_stack('/Users/jgould/git/wot/bin/data/', result.x, axis=0)
    scores = wot.Dataset(x=np.hstack(scores), row_meta=ds.row_meta, col_meta=gs.col_meta)
    wot.io.write_dataset(ds=scores, path=out + '_score', output_format=format, txt_full=False)
    if permutations is not None and permutations > 0:
        wot.io.write_dataset(ds=wot.Dataset(x=np.vstack(p_vals).T, row_meta=ds.row_meta, col_meta=gs.col_meta),
                             path=out + '_p_value', output_format=format, txt_full=False)
        wot.io.write_dataset(ds=wot.Dataset(x=np.vstack(fdrs).T, row_meta=ds.row_meta, col_meta=gs.col_meta),
                             path=out + '_fdr', output_format=format, txt_full=False)
        wot.io.write_dataset(ds=wot.Dataset(x=np.vstack(k).T, row_meta=ds.row_meta, col_meta=gs.col_meta),
                             path=out + '_k', output_format=format, txt_full=False)
        wot.io.write_dataset(ds=wot.Dataset(x=np.vstack(n).T, row_meta=ds.row_meta, col_meta=gs.col_meta),
                             path=out + '_n', output_format=format, txt_full=False)


def main(argv):
    parser = argparse.ArgumentParser(description='Compute cell gene set scores')
    parser.add_argument('--matrix', help=wot.commands.MATRIX_HELP, required=True)
    parser.add_argument('--gene_sets',
                        help='Gene sets in gmx or gmt format. If not specified gene sets for apoptosis and cell cycle are used')
    parser.add_argument('--cell_filter', help='Cells to include')
    parser.add_argument('--out', help='Output file name prefix')
    # parser.add_argument('--dask', help='Dask scheduler URL')
    parser.add_argument('--format', help=wot.commands.FORMAT_HELP, default='loom', choices=wot.commands.FORMAT_CHOICES)
    parser.add_argument('--nperm', help='Number of permutations', default=10000, type=int)
    parser.add_argument('--method', help='Method to compute gene set scores',
                        choices=['mean_z_score', 'mean', 'mean_rank'])
    parser.add_argument('--nbins', help='Number of bins for sampling', default=25, type=int)

    args = parser.parse_args(argv)
    if args.out is None:
        args.out = wot.io.get_filename_and_extension(os.path.basename(args.matrix))[0] + '_gene_set_scores'

    use_dask = False
    # use_dask = args.dask is not None
    # if use_dask:
    #     from dask.distributed import Client
    #
    #     args.format = 'loom'
    #     file_path = args.out + '.loom'
    #     if os.path.exists(file_path):
    #         os.remove(file_path)
    #     client = Client(args.dask)

    gene_sets = args.gene_sets
    if gene_sets is None:
        import sys
        gene_sets = os.path.join(os.path.dirname(sys.argv[0]), 'resources', 'growth_scores_gene_sets.gmx')
    compute(matrix=args.matrix, cell_filter=args.cell_filter, gene_sets=gene_sets, out=args.out,
            format=args.format, permutations=args.nperm, method=args.method, nbins=args.nbins)
