#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import numpy as np
import pandas as pd
import wot.io


def get_set_id_to_indices(cell_sets, tmap, is_columns):
    cell_set_id_to_indices = {}
    for i in range(len(cell_sets)):
        cell_set_name = cell_sets[i]['name']
        cell_ids = cell_sets[i]['set']
        if is_columns:
            indices = np.where(np.isin(tmap.columns, cell_ids))[0]
        else:
            indices = np.where(tmap.index.isin(cell_ids))[0]
        if len(indices) is 0:
            print(cell_set_name + ' has zero members in dataset')
        else:
            cell_set_id_to_indices[cell_set_name] = indices
    return cell_set_id_to_indices


def summarize_transport_map(transport_maps, start_cell_sets, end_cell_sets, start_time, end_time):
    transport_results = multiply_tmaps(transport_maps, start_time, end_time)
    tmap = transport_results['tmap']
    start_time_ncells = transport_results['start_time_ncells']
    start_time_g = transport_results['start_time_g']
    end_time_ncells = transport_results['end_time_ncells']
    cell_set_id_to_row_indices = get_set_id_to_indices(start_cell_sets, tmap, False)
    cell_set_id_to_column_indices = get_set_id_to_indices(end_cell_sets, tmap, True)
    summary = np.zeros(shape=(len(start_cell_sets), len(end_cell_sets)))
    nrows = len(start_cell_sets)
    ncols = len(end_cell_sets)
    rids = []
    for s in start_cell_sets:
        rids.append(s['name'])
    cids = []
    for s in end_cell_sets:
        cids.append(s['name'])

    for i in range(nrows):
        row_indices = cell_set_id_to_row_indices[rids[i]]
        if len(row_indices) == 0:
            continue
        tmap_r = tmap.values[row_indices]
        for j in range(ncols):
            column_indices = cell_set_id_to_column_indices[cids[j]]
            if len(column_indices) == 0:
                continue
            summary[i, j] += tmap_r[:, column_indices].sum()

    row_meta = pd.DataFrame(index=rids)
    col_meta = pd.DataFrame(index=cids)
    cells_start = np.zeros(nrows)
    g = np.zeros(nrows)
    for i in range(nrows):
        row_indices = cell_set_id_to_row_indices[rids[i]]
        if len(row_indices) == 0:
            continue
        cells_start[i] = len(row_indices)
        if start_time_g is not None:
            g[i] = start_time_g[row_indices].sum()
    g /= start_time_g.sum()
    cells_start /= start_time_ncells
    row_meta['cells_start'] = cells_start
    cells_end = np.zeros(ncols)
    for j in range(ncols):
        column_indices = cell_set_id_to_column_indices[cids[j]]
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

    return wot.Dataset(summary, row_meta=row_meta, col_meta=col_meta)


def multiply_tmaps(start_time, end_time, transport_maps, store=False):
    start_time_index = None
    end_time_index = None
    for i in range(len(transport_maps)):
        if transport_maps[i]['t1'] == start_time:
            start_time_index = i
        if transport_maps[i]['t2'] == end_time:
            end_time_index = i

    if start_time_index is None and end_time_index is None:
        raise RuntimeError(
            'Transport map for time ' + str(start_time) + ' and time ' + str(end_time) + ' not found.')
    elif start_time_index is None:
        raise RuntimeError(
            'Transport map for time ' + str(start_time) + ' not found.')
    elif end_time_index is None:
        raise RuntimeError(
            'Transport map for time ' + str(end_time) + ' not found.')
    tmap = None
    start_time_ncells = None
    start_time_g = None
    end_time_ncells = None
    tmaps = []
    for i in range(start_time_index, end_time_index + 1):
        ds = wot.io.read_dataset(transport_maps[i]['path'])
        if i == start_time_index:
            start_time_ncells = ds.x.shape[0]
            if ds.row_meta.get('g') is not None:
                start_time_g = ds.row_meta['g'].values
        elif i == end_time_index:
            end_time_ncells = ds.x.shape[1]
        tmap_i = pd.DataFrame(index=ds.row_meta.index, columns=ds.col_meta.index, data=ds.x)
        if tmap is None:
            tmap = tmap_i
        else:
            tmap = tmap.dot(tmap_i / np.sqrt(np.sum(tmap_i.values)))
        if store:
            tmaps.append(tmap)

    return {'start_time_ncells': start_time_ncells,
            'start_time_g': start_time_g,
            'end_time_ncells': end_time_ncells,
            'tmap': tmap,
            'tmaps': tmaps}


def main(argv):
    parser = argparse.ArgumentParser(
        description='Generate a transition table from one cell set to another cell set')
    parser.add_argument('--tmap', help=wot.commands.TMAP_HELP, required=True)
    parser.add_argument('--cell_set', help=wot.commands.CELL_SET_HELP, required=True, action='append')
    parser.add_argument('--cell_days', help=wot.commands.CELL_DAYS_HELP, required=True)
    parser.add_argument('--start_time',
                        help='The start time for the cell sets to compute the transitions to cell sets at end_time',
                        required=True, type=float)
    parser.add_argument('--end_time', help='The end time', required=True, type=float)
    parser.add_argument('--out', help='Prefix for ouput file.')
    parser.add_argument('--format', help=wot.commands.FORMAT_HELP, default='loom', choices=wot.commands.FORMAT_CHOICES)

    args = parser.parse_args(argv)

    time_to_cell_sets = wot.io.group_cell_sets(args.cell_set,
                                               pd.read_table(args.cell_days, index_col='id',
                                                             engine='python', sep=None,
                                                             dtype={'day': np.float64}))
    nsets = 0
    for t in time_to_cell_sets:
        nsets += len(time_to_cell_sets[t])

    transport_maps = wot.io.list_transport_maps(args.tmap)
    if len(transport_maps) == 0:
        print('No transport maps found in ' + args.tmap)
        exit(1)
    ds = summarize_transport_map(transport_maps=transport_maps, start_cell_sets=time_to_cell_sets[args.start_time],
                                 end_cell_sets=time_to_cell_sets[args.end_time],
                                 start_time=args.start_time, end_time=args.end_time)
    wot.io.write_dataset(ds, args.out + '_' + str(args.start_time) + '_' + str(args.end_time) + '_transition_table',
                         output_format=args.format)
