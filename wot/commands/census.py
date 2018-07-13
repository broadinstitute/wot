#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import numpy
import pandas
import wot.io
import wot.ot

def main(argv):
    parser = argparse.ArgumentParser(
            description='Generate ancestor census for each time point given an initial cell set')
    parser.add_argument('--tmap', help=wot.commands.TMAP_HELP, required=True)
    parser.add_argument('--cell_days', help=wot.commands.CELL_DAYS_HELP, required=True)
    parser.add_argument('--cell_set', help=wot.commands.CELL_SET_HELP, required=True)
    parser.add_argument('--matrix', help=wot.commands.MATRIX_HELP, required=True)
    parser.add_argument('--out', help='Output files prefix', default='census')
    parser.add_argument('--progress', action='store_true', help='Print a progress bar while computing')

    args = parser.parse_args(argv)
    time_to_cell_sets = wot.io.group_cell_sets(args.cell_set, wot.io.read_days_data_frame(args.cell_days))
    if len(time_to_cell_sets) == 0:
        print("No cell sets found")
        exit(1)

    nsets = 0
    for t in time_to_cell_sets:
        nsets += len(time_to_cell_sets[t])
    transport_maps = wot.io.list_transport_maps(args.tmap)
    if len(transport_maps) == 0:
        print("No transport maps found in " + args.tmap)
        exit(1)

    transport_map_times = set()
    for tmap in transport_maps:
        transport_map_times.add(tmap['t1'])
        transport_map_times.add(tmap['t2'])

    ds = wot.io.read_dataset(args.matrix)
    if args.progress:
        wot.io.init_progress()
    trajectories = wot.ot.Trajectory.trajectory_for_cell_sets(transport_maps, time_to_cell_sets,
            cache_transport_maps=False, print_progress=args.progress)
    if args.progress:
        wot.io.finalize_progress()
    cell_set_name_to_trajectories = wot.ot.Trajectory.group_trajectories_by_cell_set(trajectories)

    cell_set_ds = wot.io.read_gene_sets(args.cell_set)
    for timestamped_cellset in cell_set_name_to_trajectories.keys() :
        census_by_time = {}
        for distribution in cell_set_name_to_trajectories[timestamped_cellset]:
            distribution_df = pandas.DataFrame(distribution['p'],
                    index=distribution['cell_ids'], columns=['p'])
            cell_ids = distribution_df.index.intersection(cell_set_ds.row_meta.index)
            ds_order = cell_set_ds.row_meta.index.get_indexer_for(cell_ids)
            ds_order = ds_order[ds_order != -1]
            if len(cell_ids) == 0:
                census = [ 0 ] * cell_set_ds.x.shape[1]
            else:
                census = numpy.dot(distribution_df.loc[cell_ids,:]['p'], cell_set_ds.x[ds_order])
            census_by_time[distribution['t']] = census
        census_by_time = sorted(census_by_time.items())
        row_meta = pandas.DataFrame([], index=[ k for k, v in census_by_time ], columns=[])
        col_meta = cell_set_ds.col_meta
        result = wot.Dataset(numpy.asarray([ v for k, v in census_by_time ],  dtype=numpy.float64),
                row_meta, col_meta)
        wot.io.write_dataset(result, args.out + '_' + timestamped_cellset, output_format='txt', txt_full=False)
