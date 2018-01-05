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
parser.add_argument('--verbose', action='store_true',
                    help='Print progress information')
args = parser.parse_args()
gs = wot.read_gmx(args.gene_sets)
if args.verbose:
    print('Read gene sets')
ds = wot.read_dataset(args.matrix)
if args.verbose:
    print('Read dataset')

if args.transpose:
    ds.transpose()
result = wot.score_gene_sets(ds=ds, gs=gs)
if args.verbose:
    print('Computed gene set scores')
output_format = 'txt'
wot.write_dataset(ds=result, path=wot.io.check_file_extension(args.prefix,
                                                              format=output_format),
                  output_format=output_format)
