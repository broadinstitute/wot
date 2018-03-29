#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import wot.io
import os.path


def compute():
    ds = wot.io.read_dataset(args.matrix, use_dask=use_dask, chunks=(10000, None))
    if args.to_upper:
        ds.col_meta.index = ds.col_meta.index.str.upper()
    if args.verbose and not use_dask:
        print('Read dataset')
    gene_ids = ds.col_meta.index.values
    gs = wot.io.read_gene_sets(args.gene_sets, gene_ids)

    if args.verbose and not use_dask:
        print('Read gene sets')
    # scores contains cells on rows, gene sets on columns
    result = wot.score_gene_sets(ds=ds, gs=gs, z_score_ds=not args.no_zscore, use_dask=args.dask)
    if args.verbose:
        print('Saving results...')
    # import dask.array as da
    # da.to_npy_stack('/Users/jgould/git/wot/bin/data/', result.x, axis=0)
    wot.io.write_dataset(ds=result, path=args.prefix, output_format=args.format,
                         progress=lambda progress: print(progress))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compute gene set scores')
    parser.add_argument('--matrix',
                        help='Gene expression tab delimited file with cells on rows and features on columns',
                        required=True)
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
    parser.add_argument('--dask', help='Dask scheduler URL')
    parser.add_argument('--format', help='Output file format. When using dask, output file format is set to loom',
                        default='loom')

    args = parser.parse_args()
    use_dask = args.dask is not None
    if use_dask:
        from dask.distributed import Client

        args.format = 'loom'
        file_path = args.prefix + '.loom'
        if os.path.exists(file_path):
            os.remove(file_path)
        client = Client(args.dask)

    compute()
