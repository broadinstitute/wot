#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os

import anndata
import numpy as np
import pandas as pd

import wot.io


def main(argv):
    parser = argparse.ArgumentParser(description='Compute cell gene set scores')
    parser.add_argument('--matrix', help=wot.commands.MATRIX_HELP, required=True)
    parser.add_argument('--transpose', help='Transpose the matrix', action='store_true')
    parser.add_argument('--gene_sets',
                        help='Gene sets in gmx, gmt, or grp format', required=True)
    parser.add_argument('--cell_filter',
                        help='File with one cell id per line to include')
    parser.add_argument('--day_filter',
                        help='File with one cell id per line to include')
    parser.add_argument('--gene_set_filter', help='Gene sets to include')
    parser.add_argument('--nperm', help='Number of permutations to perform', type=int)
    parser.add_argument('--out', help='Output file name prefix', default='')
    parser.add_argument('--format', help=wot.commands.FORMAT_HELP, default='txt', choices=wot.commands.FORMAT_CHOICES)
    parser.add_argument('--method', help='Method to compute gene set scores',
                        choices=['mean_z_score', 'mean', 'mean_rank'], required=True)
    parser.add_argument('--verbose', action='store_true', help='Print verbose information')

    args = parser.parse_args(argv)
    if args.out is None or args.out == '':
        args.out = wot.io.get_filename_and_extension(os.path.basename(args.matrix))[0] + '_gene_set_scores'

    gene_sets = args.gene_sets

    ds = wot.io.read_dataset(args.matrix)
    if args.transpose:
        ds = ds.T
    if args.verbose:
        print('Read ' + args.matrix)
    ds = wot.io.filter_adata(ds, obs_filter=args.cell_filter)

    # background_ds = None
    # if background_cell_set is not None:
    #     background_cells_ds = wot.io.read_sets(background_cell_set)
    #     background_cells_ids = background_cells_ds.obs.index.values[np.where(background_cells_ds.X[:, 0] > 0)[0]]
    #     cell_filter = ds.obs.index.isin(background_cells_ids)
    #     background_ds = anndata.AnnData(ds.X[cell_filter], ds.obs.iloc[cell_filter], ds.var)

    gs = wot.io.read_sets(gene_sets, ds.var.index.values)
    if args.verbose:
        print('Read ' + gene_sets)

    if gs.shape[1] == 0:
        raise ValueError('No overlap of genes in gene sets and dataset')
    if args.gene_set_filter is not None:
        if os.path.exists(args.gene_set_filter):
            set_names = pd.read_csv(args.gene_set_filter, header=None, index_col=0, engine='python',
                                    sep='\n').index.values
        else:
            set_names = args.gene_set_filter.split(',')
        gs_filter = gs.var.index.isin(set_names)
        gs = gs[:, gs_filter]
    if gs.shape[1] == 0:
        raise ValueError('No gene sets')

    output_prefix = args.out + '_'

    # scores contains cells on rows, gene sets on columns
    permutations = args.nperm
    for j in range(gs.shape[1]):
        if args.verbose and gs.shape[1] > 1:
            print(gs.var.index.values[j])
        result = wot.score_gene_sets(ds=ds,
                                     gs=gs[:, [j]],
                                     permutations=permutations, method=args.method,
                                     verbose=args.verbose)
        column_names = [str(gs.var.index.values[j]) + '_score']
        if permutations is not None and permutations > 0:
            column_names.append('p_value')
            column_names.append('FDR_BH')
            column_names.append('k')
            x = np.hstack((np.array([result['score']]).T, np.array([result['p_value']]).T,
                           np.array([result['fdr']]).T, np.array([result['k']]).T))

        else:
            x = np.array([result['score']]).T

        # separate file for each gene set
        name = output_prefix + str(gs.var.index.values[j])
        name = name.replace(' ', '_')

        wot.io.write_dataset(ds=anndata.AnnData(X=x, obs=ds.obs, var=pd.DataFrame(index=column_names)),
                             path=name, output_format=args.format)
