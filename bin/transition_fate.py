#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import wot.ot
import wot.io
import csv
import pandas as pd

parser = argparse.ArgumentParser(
    description='Compute cell trajectories')
parser.add_argument('--dir',
                    help='Directory of transport maps as produced by ot',
                    required=True)
parser.add_argument('--start_time',
                    help='The start time',
                    required=True, type=float, action='append')
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

parser.add_argument('--format', help='Output file format', default='loom')

args = parser.parse_args()
start_times = args.start_time
end_time = args.end_time
prefix = args.prefix
clusters = pd.read_table(args.clusters, index_col=0, header=None, names=['cluster'], engine='python', sep=None)
cell_sets = wot.io.read_gene_sets(args.cell_sets)
transport_maps = wot.io.list_transport_maps(args.dir)

for start_time in start_times:
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
    start_time_ncells = None
    start_time_g = None
    end_time_ncells = None

    for i in range(start_time_index, end_time_index + 1):
        ds = wot.io.read_dataset(transport_maps[i]['path'])
        if i == start_time_index:
            start_time_ncells = ds.x.shape[0]
            if ds.row_meta.get('g') is not None:
                start_time_g = ds.row_meta['g'].values
        elif i == end_time_index:
            end_time_ncells = ds.x.shape[0]
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
        row_indices = np.where(tmap.index.isin(cell_ids))[0]
        cell_set_id_to_row_indices[cell_set_id] = row_indices
        column_indices = np.where(np.isin(tmap.columns, cell_ids))[0]
        cell_set_id_to_column_indices[cell_set_id] = column_indices
    summary = np.zeros(shape=(nsets, nsets))

    for i in range(nsets):
        row_indices = cell_set_id_to_row_indices[cell_set_ids[i]]
        if row_indices is None or len(row_indices) == 0:
            continue
        tmap_r = tmap.values[row_indices]
        for j in range(nsets):
            column_indices = cell_set_id_to_column_indices[cell_set_ids[j]]
            if column_indices is None or len(column_indices) == 0:
                continue
            summary[i, j] += tmap_r[:, column_indices].sum()

    row_meta = pd.DataFrame(index=cell_set_ids)
    col_meta = pd.DataFrame(index=cell_set_ids)
    cells_start = np.zeros(nsets)
    g = np.zeros(nsets)
    for i in range(nsets):
        row_indices = cell_set_id_to_row_indices[cell_set_ids[i]]
        if row_indices is None or len(row_indices) == 0:
            continue
        cells_start[i] = len(row_indices)
        if start_time_g is not None:
            g[i] = start_time_g[row_indices].sum()
    cells_start /= start_time_ncells
    row_meta['cells_start'] = cells_start
    cells_end = np.zeros(nsets)
    for j in range(nsets):
        column_indices = cell_set_id_to_column_indices[cell_set_ids[j]]
        if column_indices is None or len(column_indices) == 0:
            continue
        cells_end[j] = len(column_indices)

    cells_end /= end_time_ncells
    col_meta['cells_end'] = cells_end
    if start_time_g is not None:
        row_meta['g'] = g

    tmap_sum = tmap.values.sum()
    row_sums = summary.sum(axis=1)
    row_sums /= tmap_sum
    row_meta['sum'] = row_sums
    # row_filter = row_sums > 0

    column_sums = summary.sum(axis=0)
    column_sums /= tmap_sum
    col_meta['sum'] = column_sums
    # column_filter = column_sums > 0

    # summary = summary[row_filter]
    # summary = summary[:, column_filter]
    # row_meta = row_meta.iloc[row_filter]
    # col_meta = col_meta.iloc[column_filter]

    summary /= tmap_sum

    # tmap.to_csv(prefix + '_' + str(start_time) + '_' + str(end_time) + '_transition.txt', index_label='id', sep='\t',
    #             doublequote=False, quoting=csv.QUOTE_NONE)
    # summary.to_csv(prefix + '_' + str(start_time) + '_' + str(end_time) + '_transition_summary.txt', index_label='id',
    #                sep='\t', doublequote=False, quoting=csv.QUOTE_NONE)

    wot.io.write_dataset(
        wot.Dataset(summary, row_meta, col_meta),
        args.prefix + '_' + str(start_time) + '_' + str(args.end_time) + '_transition_summary',
        output_format=args.format)
