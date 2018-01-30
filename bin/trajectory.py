#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import pandas
import os.path
import wot
import wot.io
import wot.ot
import csv
import glob

parser = argparse.ArgumentParser(
    description='Compute the descendants given cell ids, a time t, '
                'and transport maps')
parser.add_argument('--dir',
                    help='Directory of transport maps as produced by ot',
                    required=True)
parser.add_argument('--time',
                    help='The time t',
                    required=True, type=float)
parser.add_argument('--prefix',
                    help='Prefix for ouput file names.',
                    required=True)
parser.add_argument('--id',
                    help='A file with one cell id per line to '
                         'compute '
                         'the '
                         'trajectory for. If unspecified, compute the trajectory for all cells.')
parser.add_argument('--normalize', action='store_true',
                    help='Normalize the total mass at each timepoint to one.',
                    required=False)
parser.add_argument('--compress', action='store_true',
                    help='gzip output files', required=False)

args = parser.parse_args()
input_dir = args.dir
time = args.time
prefix = args.prefix
ids = pandas.read_table(args.id, quoting=csv.QUOTE_NONE).iloc[:, 0] if args.id is not None else None

# transport map file names end with start_end.csv
# csv has prior time on rows, later time on columns


result = wot.ot.trajectory(wot.io.read_transport_maps(input_dir), time, ids=ids, normalize=args.normalize)


def output(f, is_descendants):
    f.to_csv(prefix + ('_descendants.txt' if is_descendants else
                       '_ancestors.txt') + (
                 '.gz' if
                 args.compress
                 else ''),
             index_label='id', sep='\t',
             doublequote=False, quoting=csv.QUOTE_NONE,
             compression=('gzip' if args.compress else None))


output(result['descendants_summary'], True)
output(result['ancestors_summary'], False)
