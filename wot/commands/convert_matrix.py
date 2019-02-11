#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import pandas as pd

import wot.io


def main(argv):
    parser = argparse.ArgumentParser(description='Convert matrix data formats')

    parser.add_argument('--format', help=wot.commands.FORMAT_HELP, choices=wot.commands.FORMAT_CHOICES)
    parser.add_argument('matrix', help='File(s) to convert', nargs='+')
    parser.add_argument('--obs', help='Row metadata to join with ids in matrix', action='append')
    parser.add_argument('--var', help='Column metadata to join with ids in matrix', action='append')
    parser.add_argument('--transpose', help='Transpose the matrix before saving', action='store_true')
    args = parser.parse_args(argv)
    files = args.matrix
    obs = []
    if args.obs is not None:
        for path in args.obs:
            obs.append(pd.read_csv(path, index_col='id', engine='python', sep=None))
    var = []
    if args.var is not None:
        for path in args.var:
            var.append(pd.read_csv(path, index_col='id', engine='python', sep=None))
    for f in files:
        name = wot.io.get_filename_and_extension(f)[0]
        ds = wot.io.read_dataset(f)
        for df in obs:
            ds.obs = ds.obs.join(df)
        for df in var:
            ds.var = ds.var.join(df)
        wot.io.write_dataset(ds.T if args.transpose else ds, name, output_format=args.format)
