#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import wot.io
import pandas as pd


def main(argv):
    parser = argparse.ArgumentParser(description='Convert matrix data formats')

    parser.add_argument('--format', help=wot.commands.FORMAT_HELP, default='loom',
                        choices=wot.commands.FORMAT_CHOICES + ['json'])
    parser.add_argument('matrix', help='File(s) to convert', nargs='+')
    parser.add_argument('--row_meta', help='Row metadata to join with ids in matrix', action='append')
    parser.add_argument('--col_meta', help='Column metadata to join with ids in matrix', action='append')
    args = parser.parse_args(argv)
    files = args.matrix
    row_meta = []
    if args.row_meta is not None:
        for path in args.row_meta:
            row_meta.append(pd.read_table(path, index_col='id', engine='python', sep=None))
    col_meta = []
    if args.col_meta is not None:
        for path in args.col_meta:
            col_meta.append(pd.read_table(path, index_col='id', engine='python', sep=None))
    for f in files:
        name = wot.io.get_filename_and_extension(f)[0]
        ds = wot.io.read_dataset(f)
        for df in row_meta:
            ds.row_meta = ds.row_meta.join(df)
        for df in col_meta:
            ds.col_meta = ds.col_meta.join(df)
        wot.io.write_dataset(ds, name, output_format=args.format)
