#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import pandas
import os.path
import wot

parser = argparse.ArgumentParser(
    description='Compute ancestors and descendants given cell ids, a time t, '
                'and transport maps')
parser.add_argument('--dir',
                    help='Directory of transport maps as produced by ot',
                    required=True)
parser.add_argument('--time',
                    help='The time t',
                    required=True)
parser.add_argument('--prefix',
                    help='Prefix for ouput file names. Command produces '
                         'prefix_ancestors.txt and prefix_descendants.txt',
                    required=True)
parser.add_argument('--id',
                    help='A list of cell ids at time t to compute the '
                         'trajectory for.',
                    required=True)
parser.add_argument('--compress', action='store_true',
                    help='gzip output files')

args = parser.parse_args()
input_dir = args.dir
time = args.time
prefix = args.prefix
ids = pandas.read_table(args.id).iloc[:, 0]

# transport map file names end with start_end.csv
# csv has prior time on rows, later time on columns
transport_maps_inputs = []  # file, start, end
for f in os.listdir(input_dir):
    if os.path.isfile(os.path.join(input_dir, f)):
        file_info = wot.get_file_basename_and_extension(f)
        basename = file_info['basename']
        tokens = basename.split('_')
        path = os.path.join(input_dir, f)
        t = tokens[len(tokens) - 1]
        t_minus_1 = tokens[len(tokens) - 2]
        transport_map = pandas.read_table(path, index_col=0)
        if t_minus_1 == time:
            transport_map = transport_map[transport_map.index.isin(ids)]
        if t == time:
            transport_map = transport_map[ids]
        transport_maps_inputs.append(
            {'transport_map': transport_map,
             't_minus_1': t_minus_1, 't': t, 't_minus_1_f': float(
                t_minus_1)})
transport_maps_inputs.sort(
    key=lambda x: x['t_minus_1_f'])  # sort by t_minus_1 (start time)

result = wot.trajectory(ids, transport_maps_inputs, time)
result['ancestors'].to_csv(prefix + '_ancestors.txt' + ('.gz' if
                                                        args.compress else ''),
                           index_label='id', sep='\t',
                           compression='gzip' if args.compress else None)
result['descendants'].to_csv(prefix + '_descendants.txt' + ('.gz' if
                                                            args.compress
                                                            else ''),
                             index_label='id', sep='\t',
                             compression='gzip' if args.compress else None)
