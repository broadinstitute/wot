#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import wot.io


def main(argv):
    parser = argparse.ArgumentParser(description='Compute neighborhood graph')
    parser.add_argument('--matrix', help=wot.commands.MATRIX_HELP, required=True)
    parser.add_argument('--gene_filter',
                        help='File with one gene id per line to include from the matrix')
    parser.add_argument('--cell_filter',
                        help='File with one cell id per line to include from the matrix')
    parser.add_argument('--transpose', help='Transpose the matrix', action='store_true')
    parser.add_argument('--pca_comps',
                        help='Number of PCA components.',
                        type=int, default=50)
    parser.add_argument('--diff_comps',
                        help='Number of diffusion components.',
                        type=int, default=15)
    parser.add_argument('--neighbors', help='Number of nearest neighbors',
                        type=int, default=15)
    parser.add_argument('--space', help='Space to compute the neighborhood graph in', choices=['dmap', 'pca', 'input'])
    parser.add_argument('--out',
                        help='Output file name. The file is saved in gexf format (https://gephi.org/gexf/format/)')

    args = parser.parse_args(argv)
    if args.out is None:
        args.out = 'wot-neighborhood-graph'
    space = args.space
    neighbors = args.neighbors
    pca_comps = args.pca_comps
    diff_comps = args.diff_comps
    adata = wot.io.read_dataset(args.matrix)
    if args.transpose:
        adata = adata.T
    adata = wot.io.filter_adata(adata, obs_filter=args.cell_filter, var_filter=args.gene_filter)
    wot.neighborhood_graph.compute_neighborhood_graph(adata, space=space, neighbors=neighbors, pca_comps=pca_comps,
                                                      diff_comps=diff_comps)
    wot.neighborhood_graph.graph_to_gexf(adata, args.out + '.gexf')
