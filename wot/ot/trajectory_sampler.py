#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import wot.ot
import wot.io
import pandas as pd
import scipy.stats
import os


class TrajectorySampler:

    @staticmethod
    def create_time_to_cell_sets(cell_set_paths):
        time_to_cell_sets = {}
        cell_set_group_to_names = {}
        for path in cell_set_paths:
            cell_set_ds = wot.io.read_gene_sets(path)
            cell_set_names = []
            cell_set_group_to_names[wot.io.get_filename_and_extension(os.path.basename(path))[0]] = cell_set_names
            for i in range(cell_set_ds.x.shape[1]):
                cell_set_name = cell_set_ds.col_meta.index.values[i]
                tokens = cell_set_name.split('_')
                try:
                    t = float(tokens[len(tokens) - 1])
                except ValueError:
                    raise ValueError('Cell set name ' + cell_set_name + ' must end with _time')

                cell_ids_in_set = cell_set_ds.row_meta.index[cell_set_ds.x[:, i] > 0]
                cell_sets = time_to_cell_sets.get(t)
                if cell_sets is None:
                    cell_sets = []
                    time_to_cell_sets[t] = cell_sets
                cell_set_names.append(cell_set_name)
                cell_sets.append({'set': set(cell_ids_in_set), 'name': cell_set_name})

        return {'time_to_cell_sets': time_to_cell_sets, 'cell_set_group_to_names': cell_set_group_to_names}

    @staticmethod
    def trajectory_plot(transport_maps, time_to_cell_sets, coords=None,
                        datasets=None, dataset_names=None, cache_transport_maps=False, smooth=False):

        cell_set_name_to_force_layout_traces = {}
        sampled_results = TrajectorySampler.sample_all_timepoints(transport_maps=transport_maps,
                                                                  time_to_cell_sets=time_to_cell_sets,
                                                                  datasets=datasets,
                                                                  dataset_names=dataset_names,
                                                                  cache_transport_maps=cache_transport_maps,
                                                                  smooth=smooth)
        cell_set_names = sampled_results['cell_set_names']
        results = sampled_results['results']

        for result_dict in results:  # each cell set
            highs = []
            lows = []
            for p in result_dict['pvecs']:
                cell_set = p['cell_set']
                traces = cell_set_name_to_force_layout_traces.get(cell_set)
                if traces is None:
                    traces = []
                    cell_set_name_to_force_layout_traces[cell_set] = traces

                if coords is not None:
                    t = p['t']
                    v = p['v']
                    cell_ids = p['cell_ids']
                    joined = coords.join(pd.DataFrame(index=cell_ids, data={'v': v}), how='right')
                    df_sum = joined.groupby(['px', 'py']).sum()
                    p = np.percentile(df_sum['v'].values, [5, 95])
                    lows.append(p[0])
                    highs.append(p[1])
                    traces.append({'v': v, 't': t, 'x': df_sum.index.get_level_values(0).tolist(),
                                   'y': df_sum.index.get_level_values(1).tolist(),
                                   'marker': {'color': df_sum['v'].values.tolist()}})
            cmax = np.percentile(highs, 50)
            cmin = np.percentile(lows, 50)
            for cell_set_name in cell_set_name_to_force_layout_traces:
                traces = cell_set_name_to_force_layout_traces[cell_set_name]
                for trace in traces:
                    trace['marker']['cmin'] = cmin
                    trace['marker']['cmax'] = cmax

        for name in cell_set_name_to_force_layout_traces:
            traces = cell_set_name_to_force_layout_traces[name]
            traces.sort(key=lambda x: x['t'])

        def ancestry_divergence(ancestor_dist1, ancestor_dist2):
            return 1.0 - 0.5 * np.sum(np.abs(ancestor_dist1 - ancestor_dist2))

        ancestry_divergence_traces = []

        for i in range(1, len(cell_set_names)):
            cell_set_name_i = cell_set_names[i]
            traces1 = cell_set_name_to_force_layout_traces[cell_set_name_i]

            for j in range(i):
                cell_set_name_j = cell_set_names[j]
                traces2 = cell_set_name_to_force_layout_traces[cell_set_name_j]
                x = []
                y = []
                for k in range(len(traces1)):
                    d = ancestry_divergence(traces1[k]['v'], traces2[k]['v'])
                    x.append(traces1[k]['t'])
                    y.append(d)

                ancestry_divergence_traces.append(
                    {'x': x, 'y': y, 'name': cell_set_name_i + ' vs. ' + cell_set_name_j, 'mode': 'lines',
                     'type': 'scatter'})
        for key in cell_set_name_to_force_layout_traces:
            traces = cell_set_name_to_force_layout_traces[key]
            for t in traces:
                del t['v']

        dataset_name_to_traces = sampled_results['dataset_name_to_traces']
        return {'ancestry_divergence_traces': ancestry_divergence_traces,
                'dataset_name_to_traces': dataset_name_to_traces,
                'force': cell_set_name_to_force_layout_traces}

    @staticmethod
    def interpolate(x, xi, yi, sigma):
        diff = (x - xi)
        diff *= -diff
        sigma2 = 2 * np.power(sigma, 2)
        wi = np.exp(diff / sigma2)
        fx = np.sum(yi * wi) / np.sum(wi)
        return fx

    @staticmethod
    def kernel_smooth(xi, yi, stop, start=0, steps=1000, sigma=0.7):
        xlist = np.linspace(start, stop, steps)
        fhat = np.zeros(len(xlist))
        for i in range(len(xlist)):
            fhat[i] = TrajectorySampler.interpolate(xlist[i], xi, yi, sigma)

        return xlist, fhat

    @staticmethod
    def weighted_average(result, t, weights, cell_set_name, datasets=None, dataset_names=None):
        if datasets is not None:
            for ds_index in range(len(datasets)):
                ds = datasets[ds_index]
                values = ds.x
                values = values.toarray() if scipy.sparse.isspmatrix(values) else values
                for feature_index in range(ds.x.shape[1]):
                    array = values[:, feature_index]
                    # count = 100 * (np.count_nonzero(array) / len(array))

                    trace = {
                        "dataset": dataset_names[ds_index],
                        "cell_set_name": cell_set_name,
                        "name": cell_set_name + '_' + str(ds.col_meta.index.values[feature_index]),
                        "y": np.average(array, weights=weights),
                        # "size": count,
                        # "text": "{:.1f}".format(count) + ' %',
                        "x": t
                    }
                    result.append(trace)

    @staticmethod
    def do_sampling(result, t, sampled_indices, cell_set_name, datasets=None, dataset_names=None, summaries=None):
        if datasets is not None:
            for ds_index in range(len(datasets)):
                ds = datasets[ds_index]
                values = ds.x[sampled_indices] if sampled_indices is not None else ds.x
                values = values.toarray() if scipy.sparse.isspmatrix(values) else values
                for feature_index in range(ds.x.shape[1]):
                    array = values[:, feature_index]

                    if summaries[ds_index]:
                        # count = 100 * (np.count_nonzero(array) / len(array))
                        trace = {
                            "dataset": dataset_names[ds_index],
                            "cell_set_name": cell_set_name,
                            "name": cell_set_name + '_' + str(ds.col_meta.index.values[feature_index]),
                            "y": np.mean(array),
                            # "size": count,
                            # "text": "{:.1f}".format(count) + ' %',
                            "x": t
                        }
                    else:
                        trace = {
                            "dataset": dataset_names[ds_index],
                            "set_name": cell_set_name + '' + str(ds.col_meta.index.values[feature_index]),
                            "name": t,
                            "type": 'violin',
                            "median": np.median(array),
                            "boxpoints": False,
                            "line": {
                                "color": 'black'
                            },
                            "opacity": 0.7,
                            "y": array.tolist(),
                            "box": {
                                "visible": False
                            },
                            "meanline": {
                                "visible": True
                            },
                            "points": False
                        }
                    result.append(trace)

    @staticmethod
    def sample_all_timepoints(transport_maps, time_to_cell_sets, datasets=None, dataset_names=None,
                              smooth=False, cache_transport_maps=False):
        results = []
        cell_set_names = []
        for t in time_to_cell_sets:
            cell_sets = time_to_cell_sets[t]
            for cell_set_index in range(len(cell_sets)):
                cell_set_names.append(cell_sets[cell_set_index]['name'])

            result = wot.ot.TrajectorySampler.sample_for_one_timepoint(
                cell_sets=cell_sets,
                transport_maps=transport_maps,
                time=t, unaligned_datasets=datasets, dataset_names=dataset_names,
                cache_transport_maps=cache_transport_maps)
            results.append(result)

        trace_fields = ['x', 'y']
        line_traces = []
        for result_dict in results:  # each group of cell sets
            for trace in result_dict['traces']:
                line_traces.append(trace)
        dataset_name_to_all_traces = {}
        for trace in line_traces:
            dataset_traces = dataset_name_to_all_traces.get(trace['dataset'])
            if dataset_traces is None:
                dataset_traces = []
                dataset_name_to_all_traces[trace['dataset']] = dataset_traces
            dataset_traces.append(trace)
        dataset_name_to_line_traces = {}
        for dataset_name in dataset_name_to_all_traces:
            all_traces = dataset_name_to_all_traces[dataset_name]
            trace_name_to_line_trace = {}
            dataset_name_to_line_traces[dataset_name] = trace_name_to_line_trace
            for trace in all_traces:
                line_trace = trace_name_to_line_trace.get(trace['name'])
                if line_trace is None:
                    line_trace = {'name': trace['name'],
                                  "mode": "lines+markers",
                                  "showlegend": True,
                                  "type": 'scatter'}
                    for field in trace_fields:
                        line_trace[field] = np.array([trace[field]])
                    trace_name_to_line_trace[line_trace['name']] = line_trace
                else:
                    for field in trace_fields:
                        line_trace[field] = np.concatenate((line_trace[field], [trace[field]]))
        dataset_name_to_traces = {}
        for dataset_name in dataset_name_to_line_traces:
            trace_name_to_line_trace = dataset_name_to_line_traces[dataset_name]
            traces = []
            dataset_name_to_traces[dataset_name] = traces
            for trace_name in trace_name_to_line_trace:
                trace = trace_name_to_line_trace[trace_name]
                traces.append(trace)
                sort_order = np.argsort(trace['x'])
                # max_size = max(max_size, np.max(trace['size']))
                for field in trace_fields:
                    trace[field] = trace[field][sort_order]

                if smooth:
                    x = trace['x']
                    xsmooth, ysmooth = wot.ot.TrajectorySampler.kernel_smooth(x, trace['y'], stop=x[len(x) - 1])
                    trace['x'] = xsmooth
                    trace['y'] = ysmooth

                trace['mode'] = 'lines'

                # trace['size'] = trace['size'].tolist()
                # trace['text'] = trace['text'].tolist()
                trace['showlegend'] = True
                # trace['sizemode'] = 'area'
                # trace['sizemin'] = 4
                # trace['marker'] = {'size': trace['size'], 'sizeref': (2 * 100) / (4 * 4), 'size_min': 4}

        return {'results': results, 'dataset_name_to_traces': dataset_name_to_traces,
                'cell_set_names': cell_set_names}

    @staticmethod
    def sample_for_one_timepoint(cell_sets, transport_maps, time, unaligned_datasets=[], dataset_names=[],
                                 start=None, end=None, cache_transport_maps=False):
        """

        Args:
            cell_set_ds (list): A list of dicts containing "name" and "set"
            time (float): The time at which the cell sets were defined

        """

        if len(transport_maps) == 0:
            raise ValueError('No transport maps specified.')
        t2_index = None
        t1_index = None
        for transport_index in range(len(transport_maps)):
            if transport_maps[transport_index]['t1'] == time:
                t1_index = transport_index
            if transport_maps[transport_index]['t2'] == time:
                t2_index = transport_index

        traces = []
        pvecs = []
        n_cell_sets = len(cell_sets)

        ranges = []
        t0_loaded = False
        if t2_index is not None:
            if start is None:
                start = -1
            start = max(-1, start)
            ranges.append({'backward': True, 'range': range(t2_index, - 1, start)})

        if t1_index is not None:
            if end is None:
                end = len(transport_maps)
            end = min(len(transport_maps), end)
            ranges.append({'backward': False, 'range': range(t1_index, end)})

        if len(ranges) == 0:
            raise ValueError('No transport map starting or ending at ' + str(time) + ' found.')
        for r in ranges:
            back = r['backward']
            init = True
            for transport_index in r['range']:
                tmap_dict = transport_maps[transport_index]
                t = tmap_dict['t1'] if back else tmap_dict['t2']
                t_index = t2_index if back else t1_index
                tmap = tmap_dict.get('ds')
                if tmap is None:
                    tmap = wot.io.read_dataset(tmap_dict['path'])
                    if cache_transport_maps:
                        tmap_dict['ds'] = tmap

                # align dataset and tmap
                datasets = []
                if unaligned_datasets is not None:
                    for ds_index in range(len(unaligned_datasets)):
                        unaligned_ds = unaligned_datasets[ds_index]
                        if back:
                            ds_order = unaligned_ds.row_meta.index.get_indexer_for(tmap.row_meta.index.values)
                        else:
                            ds_order = unaligned_ds.row_meta.index.get_indexer_for(tmap.col_meta.index.values)

                        ds_order = ds_order[ds_order != -1]
                        ds = wot.Dataset(unaligned_ds.x[ds_order], unaligned_ds.row_meta.iloc[ds_order],
                                         unaligned_ds.col_meta)

                        datasets.append(ds)

                if init and transport_index == t_index:

                    init = False
                    pvec_array = []
                    filtered_cell_sets = []
                    for cell_set_index in range(n_cell_sets):

                        cell_ids_in_set = cell_sets[cell_set_index]['set']
                        membership = tmap.col_meta.index.isin(
                            cell_ids_in_set) if back else tmap.row_meta.index.isin(cell_ids_in_set)

                        if np.sum(membership) > 0:
                            membership = membership.astype(np.float)
                            membership /= membership.sum()
                            pvec_array.append(membership)
                            filtered_cell_sets.append(cell_sets[cell_set_index])
                            if back:
                                entropy = np.exp(scipy.stats.entropy(membership))
                                pvecs.append(
                                    {'cell_set': cell_sets[cell_set_index]['name'],
                                     'v': membership,
                                     'entropy': entropy,
                                     'normalized_entropy': entropy / len(membership), 't': time,
                                     'cell_ids': tmap.col_meta.index.values if back else tmap.row_meta.index.values})

                            if not t0_loaded:
                                if unaligned_datasets is not None:
                                    datasets0 = []

                                    for ds_index in range(len(unaligned_datasets)):
                                        unaligned_ds = unaligned_datasets[ds_index]
                                        ds0_order = unaligned_ds.row_meta.index.get_indexer_for(cell_ids_in_set)
                                        ds0_order = ds0_order[ds0_order != -1]
                                        ds0 = wot.Dataset(unaligned_ds.x[ds0_order],
                                                          unaligned_ds.row_meta.iloc[ds0_order],
                                                          unaligned_ds.col_meta)
                                        datasets0.append(ds0)

                                    TrajectorySampler.weighted_average(result=traces, t=time, weights=None,
                                                                       datasets=datasets0,
                                                                       cell_set_name=cell_sets[cell_set_index]['name'],
                                                                       dataset_names=dataset_names)
                        else:
                            print('No overlap found for ' + str(cell_sets[cell_set_index]['name']))
                    t0_loaded = True
                    cell_sets = filtered_cell_sets
                    n_cell_sets = len(cell_sets)
                    if n_cell_sets == 0:
                        raise ValueError('No cell sets found at time ' + str(time))
                new_pvec_array = []
                for cell_set_index in range(n_cell_sets):

                    v = pvec_array[cell_set_index]
                    if back:
                        v = tmap.x.dot(v)
                    else:
                        v = v.dot(tmap.x)

                    v /= v.sum()

                    entropy = np.exp(scipy.stats.entropy(v))
                    pvecs.append(
                        {'cell_set': cell_sets[cell_set_index]['name'], 'v': v,
                         'entropy': entropy,
                         'normalized_entropy': entropy / len(v), 't': t,
                         'cell_ids': tmap.row_meta.index.values if back else tmap.col_meta.index.values})
                    # n_choose = int(np.ceil(entropy))
                    # n_choose = min(ncells, n_choose)
                    # n_choose = ncells
                    # sampled_indices = np.random.choice(len(v), n_choose, p=v, replace=True)

                    TrajectorySampler.weighted_average(result=traces, t=t, weights=v,
                                                       datasets=datasets,
                                                       cell_set_name=cell_sets[cell_set_index]['name'],
                                                       dataset_names=dataset_names)
                    new_pvec_array.append(v)
                pvec_array = new_pvec_array

        return {'traces': traces, 'pvecs': pvecs}
