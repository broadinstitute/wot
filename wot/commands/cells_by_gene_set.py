#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys

import numpy as np

import wot.io


def create_parser():
    parser = argparse.ArgumentParser(description='Generate cell sets from gene set scores')
    parser.add_argument('--score',
                        help='Gene sets scores generated from the gene_set_scores command', required=True,
                        action='append')
    parser.add_argument('--filter', help='Comma separated list of column names to include')
    parser.add_argument('--quantile', default=99,
                        help='Quantile for cells to be considered a member of a cell set', type=float)
    parser.add_argument('--out', help='Output file name prefix', default="wot")
    return parser


def main(args):
    quantile = args.quantile
    column_filter = None
    if args.filter is not None:
        column_filter = args.filter.split(',')
    cell_set_name_to_ids = {}
    for path in args.score:
        adata = wot.io.read_dataset(path)
        if column_filter is not None:
            adata = adata[:, column_filter]
        percentiles = np.percentile(adata.X, axis=0, q=quantile)
        if len(percentiles.shape) == 0:
            percentiles = [percentiles]

        for j in range(adata.shape[1]):  # each gene set
            x = adata[:, j].X
            selected = x >= percentiles[j]
            cell_ids = adata[selected].obs.index
            if len(cell_ids) > 0:
                cell_set_name_to_ids[adata.var.index[j]] = cell_ids

    if len(cell_set_name_to_ids) == 0:
        sys.exit('No cells found')

    wot.io.write_sets(cell_set_name_to_ids, args.out)
