#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import pandas
import os.path
import wot
import wot.io
import wot.ot
import csv

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
                         'trajectory for.',
                    required=True)
parser.add_argument('--normalize', action='store_true',
                    help='Normalize the total mass at each timepoint to one.',
                    required=False)
parser.add_argument('--compress', action='store_true',
                    help='gzip output files', required=False)

args = parser.parse_args()
input_dir = args.dir
time = args.time
prefix = args.prefix
ids = pandas.read_table(args.id, quoting=csv.QUOTE_NONE).iloc[:, 0]

# transport map file names end with start_end.csv
# csv has prior time on rows, later time on columns
transport_maps_inputs = []  # file, start, end
for f in os.listdir(input_dir):
    if os.path.isfile(os.path.join(input_dir, f)):
        file_info = wot.io.get_file_basename_and_extension(f)
        basename = file_info[0]
        tokens = basename.split('_')
        path = os.path.join(input_dir, f)
        t = tokens[len(tokens) - 1]
        t_minus_1 = tokens[len(tokens) - 2]
        try:
            t_minus_1_f = float(t_minus_1)

        except ValueError:
            continue
        transport_map = pandas.read_table(path, index_col=0,
                                          quoting=csv.QUOTE_NONE)
        if str(t_minus_1_f) == str(time):
            # subset rows
            transport_map = transport_map[transport_map.index.isin(ids)]
        if str(float(t)) == str(time):
            # subset columns
            transport_map = transport_map[ids]
        transport_maps_inputs.append(
            {'transport_map': transport_map,
             't_minus_1': t_minus_1, 't': t, 't_minus_1_f': t_minus_1_f})
transport_maps_inputs.sort(
    key=lambda x: x['t_minus_1_f'])  # sort by t_minus_1 (start time)

result = wot.ot.trajectory(ids, transport_maps_inputs, time, args.normalize)


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
