#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os

import wot.io


def compute(matrix, gene_sets, no_zscore, out, format, use_dask=False):
    ds = wot.io.read_dataset(matrix, use_dask=use_dask, chunks=(10000, None))
    gs = wot.io.read_gene_sets(gene_sets, ds.col_meta.index.values)
    if gs.x.shape[1] is 0:
        raise ValueError('No overlap of genes in gene sets and dataset')

    # scores contains cells on rows, gene sets on columns
    result = wot.score_gene_sets(ds=ds, gs=gs, z_score_ds=not no_zscore, use_dask=use_dask)
    # import dask.array as da
    # da.to_npy_stack('/Users/jgould/git/wot/bin/data/', result.x, axis=0)
    wot.io.write_dataset(ds=result, path=out, output_format=format)


def main(argv):
    parser = argparse.ArgumentParser(description='Compute gene set scores for each cell')
    parser.add_argument('--matrix', help=wot.commands.MATRIX_HELP, required=True)
    parser.add_argument('--gene_sets',
                        help='Gene sets in gmx or gmt format. If not specified gene sets for apoptosis and cell cycle are used')
    parser.add_argument('--out', help='Output file name prefix')
    parser.add_argument('--no_zscore', action='store_true',
                        help='Do not z-score genes')
    parser.add_argument('--verbose', action='store_true',
                        help='Print progress information')
    # parser.add_argument('--dask', help='Dask scheduler URL')
    parser.add_argument('--format', help=wot.commands.FORMAT_HELP, default='loom', choices=wot.commands.FORMAT_CHOICES)

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
    compute(matrix=args.matrix, gene_sets=gene_sets, no_zscore=args.no_zscore, out=args.out,
            format=args.format)
