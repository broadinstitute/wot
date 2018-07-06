#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os

import h5py
import numpy as np
import wot.io
import scipy


def main(argv):
    parser = argparse.ArgumentParser(
        description='Generate mean expression profiles given a starting cell cet and transport maps.')
    parser.add_argument('--tmap', help=wot.commands.TMAP_HELP, required=True)
    parser.add_argument('--cell_days', help=wot.commands.CELL_DAYS_HELP, required=True)
    parser.add_argument('--cell_set', help=wot.commands.CELL_SET_HELP, required=True, action='append')
    parser.add_argument('--matrix', help=wot.commands.MATRIX_HELP, required=True)
    parser.add_argument('--matrix_transform', help='Transform matrix values', action='append',
                        choices=['expm1', 'log1p', 'rank'])

    args = parser.parse_args(argv)
    time_to_cell_sets = wot.io.group_cell_sets(args.cell_set, wot.io.read_days_data_frame(args.cell_days))
    if len(time_to_cell_sets) == 0:
        print('No cell sets found')
        exit(1)

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
    datasets = []
    dataset_names = []
    dataset_names.append(wot.io.get_filename_and_extension(os.path.basename(args.matrix))[0])
    ds = wot.io.read_dataset(args.matrix)
    nfeatures = ds.x.shape[1]
    datasets.append(ds)
    value_transform_functions = []
    if args.matrix_transform is not None and len(args.matrix_transform) > 0:
        for matrix_transform in args.matrix_transform:
            if matrix_transform == 'expm1':
                value_transform_functions.append(np.expm1)
            elif matrix_transform == 'log1p':
                value_transform_functions.append(np.log1p)
            elif matrix_transform == 'rank':
                value_transform_functions.append(scipy.stats.rankdata)

        def value_transform(values):
            for op in value_transform_functions:
                values = op(values)
            return values

    trajectories = wot.ot.Trajectory.trajectory_for_cell_sets(transport_maps=transport_maps,
                                                              time_to_cell_sets=time_to_cell_sets)

    dataset_name_to_traces = wot.ot.TrajectoryTrends.compute(trajectories, datasets, dataset_names,
                                                             value_transform=value_transform if len(
                                                                 value_transform_functions) > 0 else None)

    transport_map_times = list(transport_map_times)
    transport_map_times.sort()
    ds_shape = (len(transport_map_times), nfeatures)
    for ds_name in dataset_name_to_traces:
        all_traces = dataset_name_to_traces[ds_name]  # each trace is a cell set/feature combo
        cell_set_name_to_traces = {}
        for trace in all_traces:
            cell_set_name = trace['cell_set_name']
            traces = cell_set_name_to_traces.get(cell_set_name)
            if traces is None:
                traces = []
                cell_set_name_to_traces[cell_set_name] = traces
            traces.append(trace)

        for cell_set_name in cell_set_name_to_traces:
            traces = cell_set_name_to_traces[cell_set_name]
            traces.sort(key=lambda x: x['feature'])
            # for each dataset, output a matrix with time on rows and features/cell sets on columns. Values in matrix are mean expression
            f = h5py.File(cell_set_name + '_' + 'trajectory_trends.loom', 'w')
            f.create_group('/layers')
            f.create_group('/row_graphs')
            f.create_group('/col_graphs')
            f.create_dataset('/row_attrs/id', data=transport_map_times)
            # cids = []
            features = []
            # cell_set_names = []

            dset = f.create_dataset('/matrix', shape=ds_shape,
                                    chunks=(1000, 1000) if ds_shape[0] >= 1000 and ds_shape[1] >= 1000 else None,
                                    maxshape=(None, ds_shape[1]),
                                    compression='gzip', compression_opts=9,
                                    data=None)
            vdset = f.create_dataset('/layers/variance', shape=ds_shape,
                                     chunks=(1000, 1000) if ds_shape[0] >= 1000 and ds_shape[1] >= 1000 else None,
                                     maxshape=(None, ds_shape[1]),
                                     compression='gzip', compression_opts=9,
                                     data=None)

            for i in range(len(traces)):
                trace = traces[i]
                # cell_set_name = np.string_(trace['cell_set_name'])
                feature = np.string_(trace['feature'])
                # cell_set_names.append(cell_set_name)
                features.append(feature)
                dset[:, i] = trace['y']
                vdset[:, i] = trace['variance']
                # cids.append(np.string_(trace['name']))
                # cids.append()

            f.create_dataset('/col_attrs/id', data=features)
            # f.create_dataset('/col_attrs/feature', data=features)
            # f.create_dataset('/col_attrs/cell_set', data=cell_set_names)
            f.close()
