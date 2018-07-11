#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import numpy as np
import pandas as pd
import wot.io


def main(argv):
    parser = argparse.ArgumentParser(
        description='Generate ancestors and descendants given a starting cell cet and transport maps')
    parser.add_argument('--tmap', help=wot.commands.TMAP_HELP, required=True)
    parser.add_argument('--cell_set', help=wot.commands.CELL_SET_HELP, required=True, action='append')
    parser.add_argument('--cell_days', help=wot.commands.CELL_DAYS_HELP, required=True)
    parser.add_argument('--out', help='Output file name', default='wot_trajectory')
    parser.add_argument('--progress', action='store_true', help='Print a progress bar while computing')

    args = parser.parse_args(argv)

    time_to_cell_sets = wot.io.group_cell_sets(args.cell_set, wot.io.read_days_data_frame(args.cell_days))
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

    if args.progress:
        wot.io.output_progress(0)
    trajectories = \
            wot.ot.Trajectory.trajectory_for_cell_sets(
                    transport_maps=transport_maps,
                    time_to_cell_sets=time_to_cell_sets,
                    cache_transport_maps=False,
                    print_progress=args.progress)

    # create a matrix with cell set on columns, cell ids on rows

    time_to_results = {}
    cell_set_name_to_column_index = {}
    for trajectory in trajectories:  # cell set
        t = trajectory['t']
        results = time_to_results.get(t)
        if results is None:
            results = []
            time_to_results[t] = results
        results.append(trajectory)
        cindex = cell_set_name_to_column_index.get(trajectory['cell_set'])
        if cindex is None:
            cindex = len(cell_set_name_to_column_index)
            cell_set_name_to_column_index[trajectory['cell_set']] = cindex

    data_to_stack = []
    ids_to_stack = []
    for t in time_to_results:
        results = time_to_results[t]
        arrays = [None] * len(cell_set_name_to_column_index)
        cell_ids_t = None
        for r in results:
            cell_set = r['cell_set']
            p = r['p']
            arrays[cell_set_name_to_column_index[cell_set]] = p
            cell_ids = r['cell_ids']
            if cell_ids_t is None:
                cell_ids_t = cell_ids
            else:
                print((cell_ids_t != cell_ids).sum())
        array_for_t = np.array(arrays).T  # create a 2d array. Each array is one cell set
        data_to_stack.append(array_for_t)
        ids_to_stack.append(cell_ids_t)
    rids = np.concatenate(ids_to_stack)
    data = np.concatenate(data_to_stack, axis=0)

    cids = [None] * len(cell_set_name_to_column_index)
    for cid in cell_set_name_to_column_index:
        cids[cell_set_name_to_column_index[cid]] = cid
    wot.io.write_dataset(wot.Dataset(data, row_meta=pd.DataFrame(index=rids), col_meta=pd.DataFrame(index=cids)),
                         args.out)
