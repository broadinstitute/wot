#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import wot.ot
import wot.io
import csv
import pandas as pd
import scipy.stats


class TrajectorySampler:
    @staticmethod
    def plot(df):
        import seaborn as sns
        sns.set_style("whitegrid")
        sns.set_style("ticks", {"xtick.major.size": 2})
        # ax = sns.violinplot(x="t", y="values", data=df)
        g = sns.factorplot(x="t", y="value", row="cell_set", col='name', data=df, kind='violin')
        g.set_xticklabels(rotation=45)
        return g

    @staticmethod
    def do_sampling(result, t, sampled_indices, cell_set_name, datasets=None, summaries=None, color=None):
        if datasets is not None:
            for ds_index in range(len(datasets)):
                ds = datasets[ds_index]
                values = ds.x[sampled_indices] if sampled_indices is not None else ds.x
                values = values.toarray() if scipy.sparse.isspmatrix(values) else values
                for column_index in range(ds.x.shape[1]):
                    # #key_data = result.get(key)
                    # if key_data is None:
                    #     key_data = {'x': np.array([]), 'y': np.array([])}
                    #     result[key] = key_data
                    array = values[:, column_index]
                    if summaries[ds_index]:
                        count = 100 * (np.count_nonzero(array) / len(array))
                        trace = {
                            "name": cell_set_name + ' ' + str(ds.col_meta.index.values[column_index]),
                            "mode": "lines+markers",
                            "showlegend": True,
                            "type": 'scatter',
                            "y": np.mean(array),
                            "size": count,
                            "text": "{:.1f}".format(count) + ' %',
                            "x": t
                        }
                    else:
                        trace = {
                            "set_name": cell_set_name + ' ' + str(ds.col_meta.index.values[column_index]),
                            "name": t,
                            "type": 'violin',
                            "median": np.median(array),
                            "boxpoints": False,
                            "line": {
                                "color": 'black'
                            },
                            "fillcolor": color,
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
    def compute(cell_set_ds, transport_maps, time, unaligned_datasets=[], summaries=[], verbose=False,
                sampling_loader=None, ncells=1000, start=None, end=None):

        t2_index = None
        t1_index = None
        for transport_index in range(len(transport_maps)):
            if transport_maps[transport_index]['t1'] == time:
                t1_index = transport_index
            if transport_maps[transport_index]['t2'] == time:
                t2_index = transport_index

        traces = []
        pvecs = []
        n_cell_sets = cell_set_ds.x.shape[1]
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
        for r in ranges:
            back = r['backward']
            color = '#2c7bb6' if back else '#d7191c'
            for transport_index in r['range']:
                tmap_dict = transport_maps[transport_index]
                t = tmap_dict['t1'] if back else tmap_dict['t2']
                t_index = t2_index if back else t1_index
                if sampling_loader is not None:
                    import h5py
                    f = h5py.File(tmap_dict['path'], 'r')
                    ids = f['/row_attrs/id'][()]
                    if ids.dtype.kind == 'S':
                        ids = ids.astype(str)
                    tmap = wot.Dataset(None, pd.DataFrame(index=ids), None)
                    f.close()
                else:
                    path = tmap_dict['path']
                    tmap = tmap_dict.get('ds')

                    # if cache_getter is not None:
                    #     cached = cache_getter(path)
                    #     if cached is not None:
                    #         tmap = cached
                    if tmap is None:
                        if verbose:
                            print('Reading transport map ' + path)
                        tmap = wot.io.read_dataset(tmap_dict['path'])
                        tmap_dict['ds'] = tmap
                        # if cache_setter is not None:
                        #     cache_setter(path, tmap)

                # align ds and tmap
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

                if transport_index == t_index:
                    if sampling_loader is None:
                        pvec_array = []
                        cell_sets_to_keep = []
                        for cell_set_index in range(n_cell_sets):
                            cell_ids_in_set = cell_set_ds.row_meta.index[cell_set_ds.x[:, cell_set_index] > 0]
                            membership = tmap.col_meta.index.isin(
                                cell_ids_in_set) if back else tmap.row_meta.index.isin(cell_ids_in_set)

                            if np.sum(membership) > 0:
                                pvec_array.append(membership)
                                cell_sets_to_keep.append(cell_set_index)

                                if not t0_loaded:
                                    t0_loaded = True

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
                                        TrajectorySampler.do_sampling(result=traces, t=time, sampled_indices=None,
                                                                      datasets=datasets0, summaries=summaries,
                                                                      cell_set_name=cell_set_ds.col_meta.index.values[
                                                                          cell_set_index], color='#ffffbf')

                        cell_set_ds = wot.Dataset(cell_set_ds.x[:, cell_sets_to_keep], cell_set_ds.row_meta,
                                                  cell_set_ds.col_meta.iloc[cell_sets_to_keep])
                    n_cell_sets = cell_set_ds.x.shape[1]

                    if verbose:
                        print('Initializing ' + str(n_cell_sets) + ' cell sets')
                    if n_cell_sets is 0:
                        raise Exception('No cell sets')
                new_pvec_array = []
                for cell_set_index in range(n_cell_sets):
                    sampled_indices = None
                    if sampling_loader is None:
                        v = pvec_array[cell_set_index]
                        if back:
                            v = tmap.x.dot(v)
                        else:
                            v = v.dot(tmap.x)
                        v /= v.sum()

                        entropy = np.exp(scipy.stats.entropy(v))
                        pvecs.append(
                            {'cell_set': cell_set_ds.col_meta.index.values[cell_set_index], 'v': v, 'entropy': entropy,
                             'normalized_entropy': entropy / len(v), 't': t,
                             'cell_ids': tmap.row_meta.index.values if back else tmap.col_meta.index.values})
                        # n_choose = int(np.ceil(entropy))
                        # n_choose = min(ncells, n_choose)
                        n_choose = ncells
                        if verbose:
                            print('Sampling ' + str(n_choose) + ' cells')
                        sampled_indices = np.random.choice(len(v), n_choose, p=v, replace=True)

                    TrajectorySampler.do_sampling(result=traces, t=t, sampled_indices=sampled_indices,
                                                  datasets=datasets,
                                                  summaries=summaries,
                                                  cell_set_name=cell_set_ds.col_meta.index.values[cell_set_index],
                                                  color=color)
                    new_pvec_array.append(v)
                pvec_array = new_pvec_array

        return {'traces': traces, 'pvecs': pvecs}
