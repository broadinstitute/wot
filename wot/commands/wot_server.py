#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import wot.io
import wot.ot
import numpy as np
import pandas as pd
import flask
import sys


def main(argsv):
    parser = argparse.ArgumentParser(description='Visualize cell set trajectories')
    parser.add_argument('--dir',
                        help='Directory of transport maps as produced by ot',
                        required=True)
    parser.add_argument('--cell_days',
                        help='Two column tab delimited file without header with cell ids and days', required=True)
    parser.add_argument('--cell_sets',
                        help='One or more gmt or gmx files containing cell sets. Each set id should end with _time (e.g. my_set_9)',
                        required=True, action='append')
    parser.add_argument('--cell_filter',
                        help='File with one cell id per line to include or or a python regular expression of cell ids to include')
    parser.add_argument('--coords',
                        help='Three column tab delimited file with header fields id, x, and y that contains 2-d cell coordinates',
                        required=True)

    parser.add_argument('--matrix',
                        help='One or more matrices with cells on rows and features, such as genes or pathways on columns',
                        action='append')

    args = parser.parse_args(argsv)
    cell_set_info = wot.ot.TrajectorySampler.create_time_to_cell_sets(args.cell_sets)
    cell_set_group_to_names = cell_set_info['cell_set_group_to_names']
    time_to_cell_sets = cell_set_info['time_to_cell_sets']

    days_data_frame = pd.read_table(args.cell_days, index_col=0, header=None, names=['day'], engine='python', sep=None,
                                    dtype={'day': np.float64})

    coords = pd.read_csv(args.coords, index_col='id', engine='python', sep=None)
    nx = 400
    ny = 400
    xmin = np.min(coords['x'])
    xmax = np.max(coords['x'])
    ymin = np.min(coords['y'])
    ymax = np.max(coords['y'])
    coords['px'] = np.floor(np.interp(coords['x'].values, [xmin, xmax], [0, nx])).astype(int)
    coords['py'] = np.floor(np.interp(coords['y'].values, [ymin, ymax], [0, ny])).astype(int)
    coords = coords.drop(['x', 'y'], axis=1)
    coords = coords.join(days_data_frame)
    coords[np.isnan(coords['day'].values)] = -1
    transport_maps = wot.io.list_transport_maps(args.dir)
    if len(transport_maps) == 0:
        raise ValueError('No transport maps found')
    transport_map_times = set()
    for tmap in transport_maps:
        transport_map_times.add(tmap['t1'])
        transport_map_times.add(tmap['t2'])
    datasets = []
    dataset_names = []
    if args.matrix is not None:
        for path in args.matrix:
            dataset_names.append(wot.io.get_filename_and_extension(os.path.basename(path))[0])
            ds = wot.io.read_dataset(path)

            if args.cell_filter is not None:
                if not os.path.isfile(args.cell_filter):
                    import re
                    expr = re.compile(args.cell_filter)
                    cell_ids = [elem for elem in ds.row_meta.index.values if expr.match(elem)]
                else:
                    cell_ids = pd.read_table(args.cell_filter, index_col=0, header=None).index.values

                # row_indices = np.isin(ds.row_meta.index.values, cell_ids, assume_unique=True)
                row_indices = ds.row_meta.index.isin(cell_ids)

                ds = wot.Dataset(ds.x[row_indices], ds.row_meta.iloc[row_indices], ds.col_meta)
            # ds_order = ds.row_meta.index.get_indexer_for(days_data_frame.index.values)
            # ds = wot.Dataset(ds.x[ds_order], ds.row_meta.iloc[ds_order], ds.col_meta)
            # np.sum(ds_order == -1)
            datasets.append(ds)

    static_folder = os.path.join(os.path.dirname(sys.argv[0]), 'web')
    app = flask.Flask(__name__, static_folder=static_folder)
    app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0

    @app.route("/cell_info/", methods=['GET'])
    def get_coords():
        return flask.jsonify(
            {'id': coords.index.values.tolist(), 'x': coords['px'].values.tolist(),
             'y': coords['py'].values.tolist(), 't': coords['day'].values.tolist()})

    @app.route("/list_times/", methods=['GET'])
    def list_times():
        return flask.jsonify(transport_map_times)

    @app.route("/list_features/", methods=['GET'])
    def list_features():
        features = set()
        for dataset_index in range(len(datasets)):
            features.update(datasets[dataset_index].col_meta.index.values.astype(
                str).tolist())
        return flask.jsonify(list(features))

    @app.route("/feature_value/", methods=['GET'])
    def feature_value():
        feature_ids = [flask.request.args.get('feature')]
        feature_ids = set(map(lambda x: x.upper(), feature_ids))
        for dataset_index in range(len(datasets)):
            column_indices = np.where(ds.col_meta.index.str.upper().isin(feature_ids))[0]
            if len(column_indices) == 1:
                return flask.jsonify(
                    {'ids': ds.row_meta.index.values.astype(str).tolist(),
                     'values': ds.x[:, column_indices[0]].tolist()})

        raise ValueError('Feature not found')

    @app.route("/list_cell_sets/", methods=['GET'])
    def list_cell_sets():
        return flask.jsonify(cell_set_group_to_names)

    @app.route("/cell_set_members/", methods=['GET'])
    def cell_set_members():
        cell_set_ids = set(flask.request.args.getlist('cell_set[]'))  # list of ids
        print(cell_set_ids)
        filtered_cell_sets = []
        for t in time_to_cell_sets:
            cell_sets = time_to_cell_sets[t]
            for cell_set in cell_sets:
                if cell_set['name'] in cell_set_ids:
                    json = {}
                    json['name'] = cell_set['name']
                    json['ids'] = list(cell_set['set'])
                    filtered_cell_sets.append(json)
        return flask.jsonify(filtered_cell_sets)

    @app.route("/trajectory/", methods=['POST'])
    def trajectory():
        feature_ids = flask.request.form.getlist('feature[]')  # list of ids
        cell_set_ids = set(flask.request.form.getlist('cell_set[]'))  # list of ids
        ncustom_cell_sets = int(flask.request.form.get('ncustom_cell_sets', '0'))
        print(flask.request.form)
        filtered_time_to_cell_sets = {}
        if ncustom_cell_sets > 0:
            for i in range(ncustom_cell_sets):
                cell_set_name = flask.request.form.get('cell_set_name' + str(i))
                cell_ids = flask.request.form.getlist('cell_set_ids' + str(i) + '[]')
                if len(cell_ids) == 0:
                    raise ValueError('No cell ids specified for custom cell set ' + cell_set_name)

                tokens = cell_set_name.split('_')
                t = float(tokens[len(tokens) - 1])
                cell_sets = filtered_time_to_cell_sets.get(t)
                if cell_sets is None:
                    cell_sets = []
                    filtered_time_to_cell_sets[t] = cell_sets
                cell_sets.append({'set': set(cell_ids), 'name': cell_set_name})

        for t in time_to_cell_sets:
            cell_sets = time_to_cell_sets[t]
            filtered_cell_sets = []
            for cell_set in cell_sets:
                if cell_set['name'] in cell_set_ids:
                    filtered_cell_sets.append(cell_set)
            if len(filtered_cell_sets) > 0:
                filtered_time_to_cell_sets[t] = filtered_cell_sets

        if len(filtered_time_to_cell_sets) == 0:
            raise ValueError('No cell sets specified')
        filtered_datasets = []
        filtered_dataset_names = []
        if feature_ids is not None and len(feature_ids) > 0:
            feature_ids = set(map(lambda x: x.upper(), feature_ids))
            for i in range(len(datasets)):
                ds = datasets[i]
                column_indices = np.where(ds.col_meta.index.str.upper().isin(feature_ids))[0]
                if len(column_indices) > 0:
                    filtered_datasets.append(
                        wot.Dataset(ds.x[:, column_indices], ds.row_meta, ds.col_meta.iloc[column_indices]))
                    filtered_dataset_names.append(dataset_names[i])

        data = wot.ot.TrajectorySampler.trajectory_plot(transport_maps=transport_maps, coords=coords,
                                                        time_to_cell_sets=filtered_time_to_cell_sets,
                                                        datasets=filtered_datasets,
                                                        dataset_names=filtered_dataset_names, cache_transport_maps=True,
                                                        smooth=True)
        dataset_name_to_traces = data['dataset_name_to_traces']
        for name in dataset_name_to_traces:
            traces = dataset_name_to_traces[name]
            for trace in traces:
                trace['x'] = trace['x'].tolist()
                trace['y'] = trace['y'].tolist()
        return flask.jsonify(data)

    app.run()


if __name__ == '__main__':
    main()
