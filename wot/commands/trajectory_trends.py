#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import wot.io
import h5py
import numpy as np
import os


def main(argv):
    parser = argparse.ArgumentParser(
        description='Generate mean expression profiles given a starting cell cet and transport maps.')
    parser.add_argument('--dir',
                        help='Directory of transport maps', required=True)
    parser.add_argument('--cell_sets',
                        help='One or more gmt or gmx files containing cell sets. Each set id should end with _time (e.g. my_set_9)',
                        required=True, action='append')
    parser.add_argument('--matrix',
                        help='A matrix with cells on rows and features, such as genes or pathways on columns',
                        required=True)

    args = parser.parse_args(argv)

    cell_set_info = wot.ot.TrajectorySampler.create_time_to_cell_sets(args.cell_sets)
    time_to_cell_sets = cell_set_info['time_to_cell_sets']
    nsets = 0
    for t in time_to_cell_sets:
        nsets += len(time_to_cell_sets[t])

    transport_maps = wot.io.list_transport_maps(args.dir)
    if len(transport_maps) == 0:
        print('No transport maps found in ' + args.dir)
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

    sampled_results = wot.ot.TrajectorySampler.sample_all_timepoints(transport_maps=transport_maps,
                                                                     time_to_cell_sets=time_to_cell_sets,
                                                                     datasets=datasets,
                                                                     dataset_names=dataset_names)

    dataset_name_to_traces = sampled_results['dataset_name_to_traces']
    transport_map_times = list(transport_map_times)
    transport_map_times.sort()
    ds_shape = (len(transport_map_times), nsets * nfeatures)
    for ds_name in dataset_name_to_traces:
        # for each cell set, output a matrix with time on rows and features on columns. Values in matrix are mean expression
        f = h5py.File(ds_name + '_trajectory_trends.loom', 'w')
        f.create_group('/layers')
        f.create_group('/row_graphs')
        f.create_group('/col_graphs')
        f.create_dataset('/row_attrs/id', data=transport_map_times)
        cids = []
        features = []
        cell_set_names = []

        dset = f.create_dataset('/matrix', shape=ds_shape,
                                chunks=(1000, 1000) if ds_shape[0] >= 1000 and ds_shape[1] >= 1000 else None,
                                maxshape=(None, ds_shape[1]),
                                compression='gzip', compression_opts=9,
                                data=None)
        traces = dataset_name_to_traces[ds_name]  # each trace is a cell set/feature combo

        for i in range(len(traces)):
            trace = traces[i]
            y = trace['y']  # across time
            trace_name = trace['name']  # cell_set_name_feature
            # cell_set_name = np.string_(trace['cell_set_name'])
            sep = trace_name.rfind('_')
            feature = trace_name[sep + 1:]
            cell_set_name = trace_name[0:sep]
            cell_set_names.append(cell_set_name)
            features.append(feature)
            dset[:, i] = y
            cids.append(np.string_(trace_name))

        f.create_dataset('/col_attrs/id', data=cids)
        f.create_dataset('/col_attrs/feature', data=features)
        f.create_dataset('/col_attrs/cell_set', data=cell_set_names)
        f.close()
