#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import wot.io


def main(argv):
    parser = argparse.ArgumentParser(description= \
                                         'Compute the list of cells in each gene set')
    parser.add_argument('--matrix',
                        help=wot.commands.MATRIX_HELP, required=True)
    parser.add_argument('--gene_sets',
                        help='Gene sets in gmx of gmt format', required=True)
    parser.add_argument('--out', help='Output file name prefix')
    parser.add_argument('--format', help=wot.commands.FORMAT_HELP,
                        default='gmt', choices=['gmt', 'gmx', 'txt'])
    parser.add_argument('--quantile', default=.99,
                        help='Proportion of cells considered to have high expression')

    args = parser.parse_args(argv)

    if args.out is None:
        args.out = 'stdout'

    dataset = wot.io.read_dataset(args.matrix)
    gene_sets = wot.io.read_sets(args.gene_sets, dataset.var.index.values)

    result = wot.get_cells_in_gene_sets(gene_sets, dataset, quantile=float(args.quantile))
    wot.io.write_gene_sets(result, args.out, args.format)
