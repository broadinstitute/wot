#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import numpy
import wot.io


def get_cells_in_gene_sets(gene_sets, dataset, quantile=.99):
    cell_sets = {}
    for gene_set_index in range(gene_sets.x.shape[1]):
        gene_indices = list(numpy.where(gene_sets.x[:,gene_set_index] == 1)[0])
        extracted = dataset.x[:, gene_indices]
        thresholds = numpy.percentile(extracted, axis=0, q=quantile * 100)
        selected = []
        for i in range(extracted.shape[0]):
            if all(extracted[i] > thresholds):
                selected.append(i)
        cell_sets[gene_sets.col_meta.index[gene_set_index]] = \
                dataset.row_meta.index[selected]
    return cell_sets

def main(argv):
    parser = argparse.ArgumentParser(description=\
            'Compute the list of cells in each gene set')
    parser.add_argument('--matrix',
            help=wot.commands.MATRIX_HELP, required=True)
    parser.add_argument('--gene_sets',
            help='Gene sets in gmx of gmt format', required=True)
    parser.add_argument('--out', help='Output file name prefix')
    parser.add_argument('--format', help=wot.commands.FORMAT_HELP,
            default='gmt', choices=['gmt', 'gmx', 'txt'])
    parser.add_argument('--quantile', default=.01,
            help='Proportion of cells considered to have high expression')

    args = parser.parse_args(argv)

    if args.out is None:
        args.out = 'stdout'

    dataset = wot.io.read_dataset(args.matrix)
    gene_sets = wot.io.read_gene_sets(args.gene_sets, dataset.col_meta.index.values)

    result = get_cells_in_gene_sets(gene_sets, dataset, quantile=float(args.quantile))
    wot.io.write_gene_sets(result, args.out, args.format)
