#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import sys

import flask
import gunicorn.app.base
import numpy as np
import pandas as pd
import scipy
import simplejson as json
import wot.io
import wot.ot
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
    parser.add_argument('--tmap', help=wot.commands.TMAP_HELP, action='append')
    parser.add_argument('--cell_set', help=wot.commands.CELL_SET_HELP, action='append')
    parser.add_argument('--cell_filter',
                        help='File with one cell id per line to include or or a python regular expression of cell ids to include')
    parser.add_argument('--cell_meta',
                        help='Every file needs to have the header "id". One file should have "x" and "y", the x and y cell coordinates',
                        action='append')
    parser.add_argument('--matrix',
                        help=wot.commands.MATRIX_HELP + '. The 1st matrix needs to be the one used to compute transport maps',
                        required=True, action='append')
    parser.add_argument('--workers', help='Number of worker processes', type=int, default=2)
    parser.add_argument('--port', help='Web server port', type=int, default=8080)
    parser.add_argument('--host', help='Web server host', default='127.0.0.1')
    parser.add_argument('--timeout', help='Worker timeout', default=600, type=int)

    args = parser.parse_args(argsv)

    cell_metadata = None
    datasets = []
    dataset_names = []
    if args.cell_meta is not None:
        for f in args.cell_meta:
            df = pd.read_table(f, engine='python', sep=None, index_col='id')
            if cell_metadata is None:
                cell_metadata = df
            else:
                cell_metadata = cell_metadata.join(df, how='outer')
    else:
        cell_metadata = pd.DataFrame()

    if args.matrix is not None:

        for path in args.matrix:
            name_and_ext = wot.io.get_filename_and_extension(os.path.basename(path))
            dataset_names.append(name_and_ext[0])
            ds = wot.io.read_dataset(path, backed=True)

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
            ds.row_meta = ds.row_meta.join(cell_metadata)
            datasets.append(ds)
        # FIXME align datasets with first dataset
        # ref_dataset = datasets[0]
        # for i in range(1, len(datasets)):
        #     ds = datasets[i]
        #     ds_order = ds.row_meta.index.get_indexer_for(ref_dataset.row_meta.index.values)
        #     ds_order = ds_order[ds_order != -1]
        #     ds = wot.Dataset(ds.x[ds_order], ds.row_meta.iloc[ds_order], ds.col_meta)
        #     datasets[i] = ds

    if args.cell_set is not None:
        time_to_cell_sets = wot.io.group_cell_sets(args.cell_set, cell_metadata)
    else:
        time_to_cell_sets = {}
    name_to_transport_maps = {}
    if args.tmap is not None:
        for d in args.tmap:
            name_to_transport_maps[os.path.basename(os.path.abspath(d))] = d

    static_folder = os.path.join(os.path.dirname(sys.argv[0]), 'web')
    app = flask.Flask(__name__, static_folder=static_folder)
    app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0

    # @app.route("/coords/", methods=['GET'])
    # def coords():
    #     nx = flask.request.args.get('nx', 800, type=int)
    #     ny = flask.request.args.get('ny', 800, type=int)
    #
    #     xmin = np.min(cell_metadata['x'])
    #     xmax = np.max(cell_metadata['x'])
    #     ymin = np.min(cell_metadata['y'])
    #     ymax = np.max(cell_metadata['y'])
    #
    #     x = np.floor(np.interp(cell_metadata['x'].values, [xmin, xmax], [0, nx])).astype(int)
    #     y = np.floor(np.interp(cell_metadata['y'].values, [ymin, ymax], [0, ny])).astype(int)
    #     return flask.Response(json.dumps(info_json, ignore_nan=True), mimetype='application/json')

    @app.route("/info/", methods=['GET'])
    def info():

        cell_json = {'id': datasets[0].row_meta.index.values.tolist()}
        for column_name in datasets[0].row_meta:
            cell_json[column_name] = datasets[0].row_meta[column_name].values.tolist()

        features = set()
        for dataset_index in range(len(datasets)):
            ds = datasets[dataset_index]
            features.update(ds.col_meta.index.values.astype(str).tolist())
        info_json = {}
        info_json['transport_maps'] = list(name_to_transport_maps.keys())
        info_json['features'] = list(features)
        info_json['cell'] = cell_json
        return flask.Response(json.dumps(info_json, ignore_nan=True), mimetype='application/json')

    @app.route("/feature_value/", methods=['GET'])
    def feature_value():
        feature_ids = flask.request.args.getlist('feature[]')
        feature_ids = set(map(lambda x: x.upper(), feature_ids))
        result = {'v': []}
        for dataset_index in range(len(datasets)):
            ds = datasets[dataset_index]
            feature_indices = np.where(ds.col_meta.index.str.upper().isin(feature_ids))[0]
            if len(feature_indices) > 0:
                for j in feature_indices:
                    values = ds.x[:, j]
                    if scipy.sparse.isspmatrix(values):
                        values = values.toarray().flatten()
                    result['v'].append({'values': values.tolist(), 'name': ds.col_meta.index.values[j]})
                result['ids'] = ds.row_meta.index.values.astype(str).tolist()
                return flask.jsonify(result)  # only return values for one dataset as dataset might not be aligned
        # not found
        return flask.jsonify(result)

    @app.route("/list_cell_sets/", methods=['GET'])
    def list_cell_sets():
        names = []
        for name in time_to_cell_sets:
            cell_sets = time_to_cell_sets[name]
            for cell_set in cell_sets:
                names.append(cell_set['name'])
        return flask.jsonify(names)

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
    def compute_trajectory():
        feature_ids = flask.request.form.getlist('feature[]')  # list of ids
        cell_set_ids = set(flask.request.form.getlist('cell_set[]'))  # list of ids
        ncustom_cell_sets = int(flask.request.form.get('ncustom_cell_sets', '0'))
        transport_map_name = flask.request.form.get('transport_map', '')
        if transport_map_name is '' and len(name_to_transport_maps) == 1:
            transport_map_name = list(name_to_transport_maps.keys())[0]

        filtered_time_to_cell_sets = {}  # maps to to dict of sets where key is set name and values are list of set ids
        if ncustom_cell_sets > 0:
            for cell_set_idx in range(ncustom_cell_sets):
                cell_set_name = flask.request.form.get('cell_set_name' + str(cell_set_idx))
                cell_ids_in_set = flask.request.form.getlist('cell_set_ids' + str(cell_set_idx) + '[]')
                subset = cell_metadata[cell_metadata.index.isin(cell_ids_in_set)]
                group_by_day = subset.groupby('day')

                for group_name, group in group_by_day:
                    cell_sets_dict = filtered_time_to_cell_sets.get(group_name)
                    if cell_sets_dict is None:
                        cell_sets_dict = {}
                        filtered_time_to_cell_sets[group_name] = cell_sets_dict
                    full_name = cell_set_name + '_' + str(group_name)
                    cell_sets_dict[full_name] = list(group.index.values)

        # filter predefined sets
        for t in time_to_cell_sets:
            cell_sets = time_to_cell_sets[t]
            filtered_cell_sets = []
            for cell_set in cell_sets:
                if cell_set['name'] in cell_set_ids:
                    filtered_cell_sets.append(cell_set)
            if len(filtered_cell_sets) > 0:
                cell_sets_dict = filtered_time_to_cell_sets.get(t)
                if cell_sets_dict is None:
                    cell_sets_dict = {}
                    filtered_time_to_cell_sets[t] = cell_sets_dict = cell_sets_dict
                for s in filtered_cell_sets:
                    cell_sets_dict[s['name']] = list(s['set'])

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

        tmap_model = wot.model.TransportMapModel.from_directory(name_to_transport_maps[transport_map_name])
        trajectory_ds = None
        for t in filtered_time_to_cell_sets:  # compute trajectories for all unique starting times
            populations = tmap_model.population_from_cell_sets(filtered_time_to_cell_sets[t], at_time=t)
            trajectory_ds_t = wot.ot.compute_trajectories(tmap_model, populations)
            if trajectory_ds is None:
                trajectory_ds = trajectory_ds_t
            else:
                # concatentate columns on trajectory datasets
                trajectory_ds.x = np.concatenate((trajectory_ds.x, trajectory_ds_t.x), axis=1)
                trajectory_ds.col_meta = pd.concat((trajectory_ds.col_meta, trajectory_ds_t.col_meta), axis=1,
                                                   sort=False)

        trajectory_similarities = wot.ot.trajectory_similarities(trajectory_ds) if trajectory_ds.x.shape[
                                                                                       1] > 1 else None
        trajectory_trends_json = {}
        for i in range(len(filtered_datasets)):
            trends = wot.ot.compute_trajectory_trends_from_trajectory(trajectory_ds, filtered_datasets[i])
            traces = []
            trajectory_trends_json[filtered_dataset_names[i]] = traces
            for j in range(len(trends)):
                mean, variance = trends[j]
                for k in range(mean.shape[1]):
                    trace = {'x': mean.row_meta.index.values.tolist(), 'y': mean.x[:, k].tolist(),
                             'name': str(mean.col_meta.index.values[k]) + '_' + str(
                                 trajectory_ds.col_meta.index.values[k]), 'mode': 'lines+markers'}
                    traces.append(trace)

        trajectory_similarities_json = []
        if trajectory_similarities is not None:
            for key in trajectory_similarities:
                t = trajectory_similarities[key]
                trajectory_similarities_json.append(
                    {'name': str(key), 'x': t['time'].tolist(), 'y': t['similarity'].tolist(), 'mode': 'lines+markers'})

        trajectory_json = []
        for j in range(trajectory_ds.x.shape[1]):
            name = str(trajectory_ds.col_meta.index.values[j])
            x = trajectory_ds.x[:, j].tolist()
            trajectory_json.append({'name': name, 'p': x})
        data = {
            'trajectory_similarities': trajectory_similarities_json,
            'trajectory_trends': trajectory_trends_json,
            'trajectory': {'id': trajectory_ds.row_meta.index.values.tolist(), 'data': trajectory_json}
        }
        return flask.Response(json.dumps(data, ignore_nan=True), mimetype='application/json')

    options = {
        'bind': '%s:%s' % (args.host, args.port),
        'workers': args.workers,
        'timeout': args.timeout
    }
    print('WOT running at http://' + args.host + ':' + str(args.port) + '/web/index.html')
    StandaloneApplication(app, options).run()
