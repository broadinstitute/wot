#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import pandas as pd

import wot.io


def create_parser():
    parser = argparse.ArgumentParser(description='Convert matrix data format')
    parser.add_argument('--matrix', help='File to convert', required=True)
    parser.add_argument('--format', help=wot.commands.FORMAT_HELP, choices=wot.commands.FORMAT_CHOICES, required=True)
    parser.add_argument('--out', help='Output file name', required=True)
    parser.add_argument('--obs', help='Row metadata to join with ids in matrix', action='append')
    parser.add_argument('--var', help='Column metadata to join with ids in matrix', action='append')
    parser.add_argument('--transpose', help='Transpose the matrix before saving', action='store_true')
    return parser


def main(args):
    obs = []
    if args.obs is not None:
        for path in args.obs:
            obs.append(pd.read_csv(path, index_col='id', engine='python', sep=None))
    var = []
    if args.var is not None:
        for path in args.var:
            var.append(pd.read_csv(path, index_col='id', engine='python', sep=None))

    ds = wot.io.read_dataset(args.matrix)
    for df in obs:
        ds.obs = ds.obs.join(df)
    for df in var:
        ds.var = ds.var.join(df)
    wot.io.write_dataset(ds.T if args.transpose else ds, args.out, output_format=args.format)
