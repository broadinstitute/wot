#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import wot.io
import wot.ot
import numpy as np
import pandas as pd
import sys
import flask
import scipy
import gunicorn.app.base
from gunicorn.six import iteritems


class StandaloneApplication(gunicorn.app.base.BaseApplication):

    def __init__(self, app, options=None):
        self.options = options or {}
        self.application = app
        super(StandaloneApplication, self).__init__()

    def load_config(self):
        config = dict([(key, value) for key, value in iteritems(self.options)
                       if key in self.cfg.settings and value is not None])
        for key, value in iteritems(config):
            self.cfg.set(key.lower(), value)

    def load(self):
        return self.application


def main(argsv):
    parser = argparse.ArgumentParser(description='Run wot server')
    parser.add_argument('--dir',
                        help='Directory of transport maps as produced by ot')
    parser.add_argument('--cell_sets',
                        help='One or more gmt or gmx files containing cell sets. Each set id should end with _time (e.g. my_cell_set_9)',
                        action='append')
    parser.add_argument('--cell_filter',
                        help='File with one cell id per line to include or or a python regular expression of cell ids to include')
    parser.add_argument('--cell_meta',
                        help='Should have headers "id", "x", "y", "day", and any additional metadata fields',
                        action='append', required=True)
    parser.add_argument('--matrix',
                        help='One or more matrices with cells on rows and features, such as genes or pathways on columns',
                        action='append')
    parser.add_argument('--workers', help='Number of worker processes', type=int, default=2)
    parser.add_argument('--port', help='Web server port', type=int, default=8080)
    parser.add_argument('--host', help='Web server host', default='127.0.0.1')

    args = parser.parse_args(argsv)
    if args.cell_sets is not None:
        cell_set_info = wot.ot.TrajectorySampler.create_time_to_cell_sets(args.cell_sets)
        cell_set_group_to_names = cell_set_info['cell_set_group_to_names']
        time_to_cell_sets = cell_set_info['time_to_cell_sets']
    else:
        cell_set_group_to_names = {}
        time_to_cell_sets = {}
    cell_metadata = None
    if args.cell_meta is not None:
        for f in args.cell_meta:
            df = pd.read_table(f, engine='python', sep=None, index_col='id')
            if cell_metadata is None:
                cell_metadata = df
            else:
                cell_metadata = cell_metadata.join(df)

    nx = 800
    ny = 800
    xmin = np.min(cell_metadata['x'])
    xmax = np.max(cell_metadata['x'])
    ymin = np.min(cell_metadata['y'])
    ymax = np.max(cell_metadata['y'])
    cell_metadata['x'] = np.floor(np.interp(cell_metadata['x'].values, [xmin, xmax], [0, nx])).astype(int)
    cell_metadata['y'] = np.floor(np.interp(cell_metadata['y'].values, [ymin, ymax], [0, ny])).astype(int)
    # coords = coords.drop(['x', 'y'], axis=1)
    # cell_metadata[np.isnan(cell_metadata['day'].values)] = -1

    if args.dir is not None:
        transport_maps = wot.io.list_transport_maps(args.dir)
        if len(transport_maps) == 0:
            raise ValueError('No transport maps found')
    else:
        transport_maps = []
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

    @app.route("/info/", methods=['GET'])
    def info():
        json = {'id': cell_metadata.index.values.tolist()}
        for column_name in cell_metadata:
            json[column_name] = cell_metadata[column_name].values.tolist()

        features = set()
        for dataset_index in range(len(datasets)):
            features.update(datasets[dataset_index].col_meta.index.values.astype(
                str).tolist())
        info_json = {}
        info_json['features'] = list(features)
        info_json['transport_map_times'] = list(transport_map_times)
        info_json['cell'] = json
        return flask.jsonify(info_json)

    @app.route("/feature_value/", methods=['GET'])
    def feature_value():
        feature_ids = [flask.request.args.get('feature')]
        feature_ids = set(map(lambda x: x.upper(), feature_ids))
        for dataset_index in range(len(datasets)):
            column_indices = np.where(ds.col_meta.index.str.upper().isin(feature_ids))[0]
            if len(column_indices) == 1:
                values = ds.x[:, column_indices[0]]
                if scipy.sparse.isspmatrix(values):
                    values = values.toarray()
                return flask.jsonify({'ids': ds.row_meta.index.values.astype(str).tolist(), 'values': values.tolist()})

        raise ValueError('Feature not found')

    @app.route("/list_cell_sets/", methods=['GET'])
    def list_cell_sets():
        return flask.jsonify(cell_set_group_to_names)

    @app.route("/cell_set_members/", methods=['GET'])
    def cell_set_members():
        cell_set_ids = set(flask.request.args.getlist('cell_set[]'))  # list of ids
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

        data = wot.ot.TrajectorySampler.trajectory_plot(transport_maps=transport_maps, coords=cell_metadata[['x', 'y']],
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

    options = {
        'bind': '%s:%s' % (args.host, args.port),
        'workers': args.workers,
    }
    print('Please go to http://127.0.0.1:' + str(args.port) + '/web/index.html')
    StandaloneApplication(app, options).run()
