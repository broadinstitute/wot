#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import logging
import os

import anndata
import pandas as pd

import wot

logger = logging.getLogger('wot')


def create_parser():
    parser = argparse.ArgumentParser(
        description='Score each cell according to its expression of input gene signatures')
    parser.add_argument('--matrix', help='A matrix with cells on rows and genes on columns', required=True)
    parser.add_argument('--gene_sets',
                        help='Gene sets in gmx, gmt, or grp format', required=True)
    parser.add_argument('--method', help='Method to compute gene set scores',
                        choices=['mean_z_score', 'mean', 'mean_rank'], default='mean_z_score')
    parser.add_argument('--cell_filter',
                        help='File with one cell id per line to include')
    parser.add_argument('--gene_set_filter', help='Gene sets to include')
    parser.add_argument('--max_z_score', help='Threshold z-scores at specified value', type=float, default=5)
    parser.add_argument('--nperm', help='Number of permutations to perform', type=int)
    parser.add_argument('--out', help='Output file name prefix', default='')
    parser.add_argument('--transpose', help='Transpose the matrix', action='store_true')
    parser.add_argument('--format', help=wot.commands.FORMAT_HELP, default='txt', choices=wot.commands.FORMAT_CHOICES)
    parser.add_argument('--verbose', action='store_true', help='Print verbose information')
    return parser


def main(args):
    if args.verbose:
        logger.setLevel(logging.DEBUG)
        logger.addHandler(logging.StreamHandler())
    if args.out is None or args.out == '':
        args.out = wot.io.get_filename_and_extension(os.path.basename(args.matrix))[0] + '_gene_set_scores'

    gene_sets = args.gene_sets
    max_z_score = abs(args.max_z_score)
    ds = wot.io.read_dataset(args.matrix)
    if args.transpose:
        ds = ds.T
    logger.info('Read ' + args.matrix)
    ds = wot.io.filter_adata(ds, obs_filter=args.cell_filter)

    # background_ds = None
    # if background_cell_set is not None:
    #     background_cells_ds = wot.io.read_sets(background_cell_set)
    #     background_cells_ids = background_cells_ds.obs.index.values[np.where(background_cells_ds.X[:, 0] > 0)[0]]
    #     cell_filter = ds.obs.index.isin(background_cells_ids)
    #     background_ds = anndata.AnnData(ds.X[cell_filter], ds.obs.iloc[cell_filter], ds.var)

    gs = wot.io.read_sets(gene_sets, ds.var.index.values)
    logger.info('Read ' + gene_sets)

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

    # scores contains cells on rows, gene sets on columns
    permutations = args.nperm
    result_df = pd.DataFrame(index=ds.obs)
    for j in range(gs.shape[1]):
        if gs.shape[1] > 1:
            logger.info(gs.var.index.values[j])
        result = wot.score_gene_sets(ds=ds,
                                     gs=gs[:, [j]],
                                     permutations=permutations, method=args.method, max_z_score=max_z_score)
        gene_set_name = str(gs.var.index.values[j])
        result_df[gene_set_name + '_score'] = result['score']
        if permutations is not None and permutations > 0:
            result_df[gene_set_name + '_p_value'] = result['p_value']
            result_df[gene_set_name + '_fdr_bh'] = result['fdr']
            result_df[gene_set_name + '_k'] = result['k']

        # separate file for each gene set
        # name = output_prefix + str(gs.var.index.values[j])
        # name = name.replace(' ', '_')

        # wot.io.write_dataset(ds=anndata.AnnData(X=x, obs=ds.obs, var=pd.DataFrame(index=column_names)),
        #                      path=name, output_format=args.format)
    wot.io.write_dataset(ds=anndata.AnnData(X=result_df.values, obs=ds.obs, var=pd.DataFrame(index=result_df.columns)),
                         path=args.out, output_format=args.format)
