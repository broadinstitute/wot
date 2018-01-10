#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import wot.io
import dask

parser = argparse.ArgumentParser(
    description='Compute the gene set scores for the given gene sets')
parser.add_argument('--matrix',
                    help='Gene expression tab delimited file with cells on '
                         'rows and features on columns', required=True)
parser.add_argument('--gene_sets',
                    help='Gene sets in gmx or gmt format.',
                    required=True)
parser.add_argument('--prefix',
                    help='Output file name prefix',
                    required=True)
parser.add_argument('--no_zscore', action='store_true',
                    help='Do not z-score genes')
parser.add_argument('--verbose', action='store_true',
                    help='Print progress information')

args = parser.parse_args()
ds = wot.io.read_dataset(args.matrix)
if args.verbose:
    print('Read dataset')
gs = wot.io.read_gene_sets(args.gene_sets, ds.col_meta.index.values.compute())
if args.verbose:
    print('Read gene sets')

result = wot.score_gene_sets(ds=ds, gs=gs, z_score_ds=not args.no_zscore)

result.x, result.row_meta, result.col_meta = dask.compute(result.x,
                                                          result.row_meta,
                                                          result.col_meta)
if args.verbose:
    print('Computed gene set scores')
output_format = 'txt'
wot.io.write_dataset(ds=result, path=wot.io.check_file_extension(args.prefix,
                                                              format=output_format),
                  output_format=output_format)
