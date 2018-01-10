#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import wot
import pandas as pd
import os

parser = argparse.ArgumentParser(
    description='Computes the distance between a reference transport map and '
                'one or more transport maps. The columns of the transport '
                'maps must be in the same order.')

parser.add_argument('--ref',
                    help='The reference cluster transport map',
                    required=True)
parser.add_argument('--dir',
                    help='Directory of reference transport maps as produced by '
                         'ot. Used to compute weights for each column.',
                    required=True)
parser.add_argument('--cluster',
                    help='Two column tab delimited file without header with '
                         'cell id and '
                         'cluster id',
                    required=True)
parser.add_argument('file', nargs='+', help='One or more transport '
                                            'maps to compare against '
                                            'the reference')
parser.add_argument('--prefix',
                    help='Prefix for ouput file names',
                    required=True)
parser.add_argument('--compress', action='store_true',
                    help='gzip output files')
parser.add_argument('--verbose', action='store_true',
                    help='Print progress information')

args = parser.parse_args()
ref_transport_map = pd.read_table(args.ref, index_col=0)
clusters = pd.read_table(args.cluster, index_col=0, header=None,
                             names=['cluster'])
grouped_by_cluster = clusters.groupby(clusters.columns[0], axis=0)
cluster_ids = list(grouped_by_cluster.groups.keys())
all_cell_ids = set()
for f in os.listdir(args.dir):
    if os.path.isfile(os.path.join(args.dir, f)):
        transport_map = pd.read_table(os.path.join(args.dir, f),
                                          index_col=0)
        all_cell_ids.update(transport_map.columns)
        all_cell_ids.update(transport_map.index)

total_cluster_size = wot.get_column_weights(all_cell_ids, grouped_by_cluster,
                                            cluster_ids)
names = list()
distances = list()

for f in args.file:
    compare_to = pd.read_table(f, index_col=0)
    if args.verbose:
        print(f)
    d = wot.transport_map_distance(transport_map_1=ref_transport_map,
                                   transport_map_2=compare_to,
                                   column_weights=total_cluster_size)
    names.append(f)
    distances.append(d)

pd.DataFrame({"names": names, "distances": distances}).to_csv(
    args.prefix +
    '.txt' + ('.gz' if
    args.compress
    else
    ''),
    sep='\t', header=False,
    compression='gzip' if args.compress
    else None)
