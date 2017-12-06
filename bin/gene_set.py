#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import wot
import numpy as np

parser = argparse.ArgumentParser(
    description='Compute the gene set scores for the given gene sets')
parser.add_argument('--matrix',
                    help='Gene expression tab delimited file with cells on '
                         'rows and features on columns', required=True)
parser.add_argument('--gene_sets',
                    help='Gene sets in gmx format.',
                    required=True)
parser.add_argument('--prefix',
                    help='Output file name prefix',
                    required=True)
parser.add_argument('--transpose', action='store_true',
                    help='Transpose the matrix')
args = parser.parse_args()

ds = wot.read_sparse(args.matrix)

gs = wot.read_gmx(args.gene_sets)
if args.transpose:
    ds.transpose()
result = wot.score_gene_sets(ds=ds, gs=gs, z_score=True)
result.write(path=args.prefix, output_format='txt')
