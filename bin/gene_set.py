#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import wot.io

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
parser.add_argument('--to_upper',
                    help='Convert matrix ids to uppercase', action='store_true')
parser.add_argument('--no_zscore', action='store_true',
                    help='Do not z-score genes')
parser.add_argument('--verbose', action='store_true',
                    help='Print progress information')

args = parser.parse_args()
use_dask = False

ds = wot.io.read_dataset(args.matrix, use_dask=use_dask)
if args.to_upper:
    ds.col_meta.index = ds.col_meta.index.str.upper()
if args.verbose:
    print('Read dataset')
gene_ids = ds.col_meta.index.values

if use_dask:
    gene_ids = gene_ids.compute()
gs = wot.io.read_gene_sets(args.gene_sets, gene_ids)

if args.verbose:
    print('Read gene sets')

result = wot.score_gene_sets(ds=ds, gs=gs, z_score_ds=not args.no_zscore)

if use_dask:
    import dask

    result.x, result.row_meta, result.col_meta = dask.compute(result.x, result.row_meta, result.col_meta)

output_format = 'txt'
path = wot.io.check_file_extension(args.prefix,
                                   format=output_format)
if args.verbose:
    print('Computed gene set scores, writing output to ' + path)

wot.io.write_dataset(ds=result, path=path,
                     output_format=output_format)
