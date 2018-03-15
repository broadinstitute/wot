#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import wot.ot
import wot.io
import csv
import pandas as pd

parser = argparse.ArgumentParser(
    description='Compute the trajectory')
parser.add_argument('--dir',
                    help='Directory of transport maps as produced by ot',
                    required=True)
parser.add_argument('--start_time',
                    help='The start time',
                    required=True, type=float)
parser.add_argument('--end_time',
                    help='The end time',
                    required=True, type=float)
parser.add_argument('--prefix',
                    help='Prefix for ouput file names.',
                    required=True)

parser.add_argument('--clusters',
                    help='Two column tab delimited file without header with '
                         'cell id and '
                         'cluster id',
                    required=True)

parser.add_argument('--cell_sets',
                    help='Grouping of clusters into cell sets',
                    required=True)

args = parser.parse_args()
start_time = args.start_time
end_time = args.end_time
prefix = args.prefix
clusters = pd.read_table(args.clusters, index_col=0, header=None, names=['cluster'], engine='python', sep=None)
cell_sets = wot.io.read_gene_sets(args.cell_sets)
transport_maps = wot.io.list_transport_maps(args.dir)
start_time_index = None
end_time_index = None
for i in range(len(transport_maps)):
    if transport_maps[i]['t1'] == start_time:
        start_time_index = i
    if transport_maps[i]['t2'] == end_time:
        end_time_index = i

if start_time_index is None:
    raise RuntimeError(
        'Transport transport_map for time ' + str(start_time) + ' not found.')
if end_time_index is None:
    raise RuntimeError(
        'Transport transport_map for time ' + str(end_time) + ' not found.')
tmap = None
for i in range(start_time_index, end_time_index + 1):
    print(transport_maps[i]['path'])
    ds = wot.io.read_dataset(transport_maps[i]['path'])
    tmap_i = pd.DataFrame(index=ds.row_meta.index, columns=ds.col_meta.index, data=ds.x)
    if tmap is None:
        tmap = tmap_i
    else:
        tmap = tmap.dot(tmap_i / np.sqrt(np.sum(tmap_i.values)))
nsets = cell_sets.x.shape[1]
cell_set_ids = cell_sets.col_meta.index.values
cell_set_id_to_row_indices = {}
cell_set_id_to_column_indices = {}
for i in range(len(cell_set_ids)):
    cell_set_id = cell_set_ids[i]
    cluster_ids = cell_sets.row_meta.index[cell_sets.x[:, i] > 0]
    cell_ids = clusters.index.values[clusters['cluster'].isin(cluster_ids)]
    cell_set_id_to_row_indices[cell_set_id] = np.where(tmap.index.isin(cell_ids))
    cell_set_id_to_column_indices[cell_set_id] = np.where(np.isin(tmap.columns, cell_ids))
summary = pd.DataFrame(index=cell_set_ids, columns=cell_set_ids, data=0.0)
for i in range(nsets):
    row_indices = cell_set_id_to_row_indices[cell_set_ids[i]]
    tmap_r = tmap.values[row_indices]
    for j in range(nsets):
        column_indices = cell_set_id_to_row_indices[cell_set_ids[j]]
        summary.iloc[i, j] += tmap_r[:, column_indices].sum()

# tmap.to_csv(prefix + '_' + str(start_time) + '_' + str(end_time) + '_transition.txt', index_label='id', sep='\t',
#             doublequote=False, quoting=csv.QUOTE_NONE)
summary.to_csv(prefix + '_' + str(start_time) + '_' + str(end_time) + '_transition_summary.txt', index_label='id',
               sep='\t', doublequote=False, quoting=csv.QUOTE_NONE)
