#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os

import numpy as np
import pandas as pd
import wot.io


def compute(matrix, gene_sets, out, format, cell_filter=None, background_cell_set=None,
            permutations=None, method='mean_z_score', nbins=None, drop_frequency=None, drop_p_value_threshold=None,
            gene_set_filter=None):
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
    if gene_set_filter is not None:
        if os.path.exists(gene_set_filter):
            set_names = pd.read_table(gene_set_filter, header=None, index_col=0, engine='python', sep=None).index.values
        else:
            set_names = gene_set_filter.split(',')
        gs_filter = gs.col_meta.index.isin(set_names)
        gs = wot.Dataset(gs.x[:, gs_filter], gs.row_meta, gs.col_meta.iloc[gs_filter])

    # scores contains cells on rows, gene sets on columns
    for j in range(gs.x.shape[1]):
        result = wot.score_gene_sets(dataset_to_score=ds, background_ds=background_ds,
                                     gs=wot.Dataset(gs.x[:, [j]], gs.row_meta, gs.col_meta.iloc[[j]]),
                                     permutations=permutations, method=method, nbins=nbins,
                                     drop_frequency=drop_frequency, drop_p_value_threshold=drop_p_value_threshold)
        column_names = [str(gs.col_meta.index.values[j])]
        if permutations is not None and permutations > 0:
            column_names.append('p_value')
            column_names.append('FRD(BH)')
            column_names.append('k')
            column_names.append('n')
            if drop_frequency > 0:
                column_names.append('p_value_ci')
                x = np.hstack((result['score'], np.vstack(
                    (result['p_value'], result['fdr'], result['k'], result['n'], result['p_value_ci'])).T))
            else:
                x = np.hstack((result['score'], np.vstack(
                    (result['p_value'], result['fdr'], result['k'], result['n'])).T))

        else:
            x = result['score']
        wot.io.write_dataset(ds=wot.Dataset(x=x, row_meta=ds.row_meta, col_meta=pd.DataFrame(index=column_names)),
                             path=out + '_' + column_names[0], output_format=format, txt_full=False)
    # import dask.array as da
    # da.to_npy_stack('/Users/jgould/git/wot/bin/data/', result.x, axis=0)


def main(argv):
    parser = argparse.ArgumentParser(description='Compute cell gene set scores')
    parser.add_argument('--matrix', help=wot.commands.MATRIX_HELP, required=True)
    parser.add_argument('--gene_sets',
                        help='Gene sets in gmx or gmt format. If not specified gene sets for apoptosis and cell cycle are used')
    parser.add_argument('--cell_filter', help='Cells to include')
    parser.add_argument('--gene_set_filter', help='Gene sets to include')
    parser.add_argument('--out', help='Output file name prefix')
    # parser.add_argument('--dask', help='Dask scheduler URL')
    parser.add_argument('--format', help=wot.commands.FORMAT_HELP, default='loom', choices=wot.commands.FORMAT_CHOICES)
    parser.add_argument('--nperm', help='Number of permutations', default=10000, type=int)
    parser.add_argument('--method', help='Method to compute gene set scores',
                        choices=['mean_z_score', 'mean', 'mean_rank'])
    parser.add_argument('--nbins', help='Number of bins for sampling', default=25, type=int)
    parser.add_argument('--drop_p_value_threshold',
                        help='Exclude cells from further permutations when the estimated lower bound of the nominal p-value is >= drop_p_value_threshold',
                        default=0.05, type=float)
    parser.add_argument('--drop_frequency',
                        help='Check the estimated lower bound of the nominal p-value every drop_frequency permutations',
                        default=1000, type=int)

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
            format=args.format, permutations=args.nperm, method=args.method, nbins=args.nbins,
            drop_frequency=args.drop_frequency, drop_p_value_threshold=args.drop_p_value_threshold,
            gene_set_filter=args.gene_set_filter)
