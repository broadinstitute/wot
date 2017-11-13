#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import wot
import pandas
import os

parser = argparse.ArgumentParser(
    description='Convert a series of cell transport maps to a single '
                'summarized cluster '
                'by cluster transport map')

parser.add_argument('--dir',
                    help='Directory of transport maps as produced by ot',
                    required=True)
parser.add_argument('--clusters',
                    help='Two column file without header with cell id and '
                         'cluster id',
                    required=True)
parser.add_argument('--prefix',
                    help='Prefix for ouput file names',
                    required=True)
parser.add_argument('--save', action='store_true',
                    help='Save cluster by cluster transport maps for each '
                         'input transport map')
parser.add_argument('--compress', action='store_true',
                    help='gzip output files')

args = parser.parse_args()
input_dir = args.dir
cluster_transport_maps = []
clusters = pandas.read_table(args.clusters, index_col=0, header=None,
                             names=['cluster'])
grouped_by_cluster = clusters.groupby(clusters.columns[0], axis=0)
cluster_ids = list(grouped_by_cluster.groups.keys())
data_frames_with_column_ids_only = []
for f in os.listdir(input_dir):
    if os.path.isfile(os.path.join(input_dir, f)):
        file_info = wot.get_file_basename_and_extension(f)
        if file_info['ext'] == '.txt':
            basename = file_info['basename']
            path = os.path.join(input_dir, f)
            transport_map = pandas.read_table(path, index_col=0)
            cluster_transport_map = wot.transport_map_by_cluster(
                transport_map, grouped_by_cluster, cluster_ids)
            if args.save:
                cluster_transport_map.to_csv(args.prefix + basename + '.txt' + (
                    '.gz' if
                    args.compress
                    else ''),
                                             index_label="id",
                                             sep='\t',
                                             compression='gzip' if args.compress
                                             else None)
            cluster_transport_maps.append(cluster_transport_map)
            data_frames_with_column_ids_only.append(pandas.DataFrame(index=[],
                                                                     columns=transport_map.columns))

weights = wot.get_weights(data_frames_with_column_ids_only, grouped_by_cluster,
                          cluster_ids)
cluster_weights_by_time = weights['cluster_weights_by_time']
combined_cluster_map = wot.transport_maps_by_time(cluster_transport_maps,
                                                  cluster_weights_by_time)
combined_cluster_map.to_csv(args.prefix + '.txt' + ('.gz' if
                                                    args.compress else ''),
                            index_label="id",
                            sep='\t',
                            compression='gzip' if args.compress
                            else None)
