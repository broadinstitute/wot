#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import wot.io
import numpy as np
import pandas as pd


def main(argv):
    parser = argparse.ArgumentParser(
        description='Generate ancestors and descendants given a starting cell cet and transport maps')
    parser.add_argument('--tmap',
                        help='Directory of transport maps', required=True)
    parser.add_argument('--cell_set',
                        help='One or more gmt or gmx files containing cell sets.',
                        required=True, action='append')
    parser.add_argument('--cell_days',
                        help='File with headers "id" and "day" corresponding to cell id and days',
                        required=True)
    parser.add_argument('--out', help='Output file name', default='wot_trajectory')

    args = parser.parse_args(argv)

    time_to_cell_sets = wot.ot.TrajectorySampler.group_cell_sets(args.cell_set,
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
    transport_map_times = set()
    for tmap in transport_maps:
        transport_map_times.add(tmap['t1'])
        transport_map_times.add(tmap['t2'])

    sampled_results = wot.ot.TrajectorySampler.sample_all_timepoints(transport_maps=transport_maps,
                                                                     time_to_cell_sets=time_to_cell_sets)
    sampled_results = sampled_results['results']
    # create a matrix with cell set on columns, cell ids on rows

    time_to_results = {}
    cell_set_name_to_column_index = {}
    for result_dict in sampled_results:  # cell set
        for p in result_dict['pvecs']:  # time
            t = p['t']
            results = time_to_results.get(t)
            if results is None:
                results = []
                time_to_results[t] = results
            results.append(p)
            cindex = cell_set_name_to_column_index.get(p['cell_set'])
            if cindex is None:
                cindex = len(cell_set_name_to_column_index)
                cell_set_name_to_column_index[p['cell_set']] = cindex

    data_to_stack = []
    ids_to_stack = []
    for t in time_to_results:
        results = time_to_results[t]
        arrays = [None] * len(cell_set_name_to_column_index)
        cell_ids_t = None
        for r in results:
            cell_set = r['cell_set']
            v = r['v']
            arrays[cell_set_name_to_column_index[cell_set]] = v
            cell_ids = r['cell_ids']
            if cell_ids_t is None:
                cell_ids_t = cell_ids
            else:
                print((cell_ids_t != cell_ids).sum())
        array_for_t = np.array(arrays).T  # create a 2d array. Each array is one cell set
        data_to_stack.append(array_for_t)
        ids_to_stack.append(cell_ids_t)
    for x in data_to_stack:
        print(x.shape)
    rids = np.concatenate(ids_to_stack)
    data = np.concatenate(data_to_stack, axis=0)

    cids = [None] * len(cell_set_name_to_column_index)
    for cid in cell_set_name_to_column_index:
        cids[cell_set_name_to_column_index[cid]] = cid
    wot.io.write_dataset(wot.Dataset(data, row_meta=pd.DataFrame(index=rids), col_meta=pd.DataFrame(index=cids)),
                         args.out)
