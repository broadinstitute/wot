#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import wot
import pandas as pd
import os
import csv

parser = argparse.ArgumentParser(
    description='Convert a series of cell transport maps to a single '
                'summarized cluster '
                'by cluster transport map')

parser.add_argument('--dir',
                    help='Directory of transport maps as produced by ot',
                    required=True)
parser.add_argument('--clusters',
                    help='Two column tab delimited file without header with '
                         'cell id and '
                         'cluster id',
                    required=True)
parser.add_argument('--prefix',
                    help='Prefix for ouput file names',
                    required=True)
parser.add_argument('--details', action='store_true',
                    help='Save cluster by cluster transport maps for each '
                         'input transport map')
parser.add_argument('--compress', action='store_true',
                    help='gzip output files')

args = parser.parse_args()
input_dir = args.dir
cluster_transport_maps = []
clusters = pd.read_table(args.clusters, index_col=0, header=None,
                             names=['cluster'], quoting=csv.QUOTE_NONE)
grouped_by_cluster = clusters.groupby(clusters.columns[0], axis=0)
cluster_ids = list(grouped_by_cluster.groups.keys())
column_cell_ids_by_time = []
all_cell_ids = set()

for f in os.listdir(input_dir):
    if os.path.isfile(os.path.join(input_dir, f)):
        file_info = wot.get_file_basename_and_extension(f)
        basename = file_info['basename']
        path = os.path.join(input_dir, f)
        transport_map = pd.read_table(path, index_col=0,
                                          quoting=csv.QUOTE_NONE)
        all_cell_ids.update(transport_map.columns)
        all_cell_ids.update(transport_map.index)
        column_cell_ids_by_time.append(transport_map.columns)
        cluster_transport_map = wot.transport_map_by_cluster(
            transport_map, grouped_by_cluster, cluster_ids)
        if args.details:
            cluster_transport_map.to_csv(args.prefix + basename + '.txt' + (
                '.gz' if
                args.compress
                else ''),
                                         index_label="id",
                                         sep='\t',
                                         compression='gzip' if args.compress
                                         else None)
        cluster_transport_maps.append(cluster_transport_map)

weights = wot.get_weights(all_cell_ids, column_cell_ids_by_time,
                          grouped_by_cluster, cluster_ids)
cluster_weights_by_time = weights['cluster_weights_by_time']
combined_cluster_map = wot.transport_maps_by_time(cluster_transport_maps,
                                                  cluster_weights_by_time)
combined_cluster_map.to_csv(args.prefix + '.txt' + ('.gz' if
                                                    args.compress else ''),
                            index_label="id",
                            sep='\t',
                            compression='gzip' if args.compress
                            else None)
