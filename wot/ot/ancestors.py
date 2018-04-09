#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import wot.ot
import wot.io
import csv
import pandas as pd
import scipy.stats


class Ancestors:
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
    def from_cmd_line(cmd_line=None):
        parser = argparse.ArgumentParser(
            description='Compute cell ancestors/descendants')
        parser.add_argument('--dir',
                            help='Directory of transport maps as produced by ot',
                            required=True)
        parser.add_argument('--time',
                            help='The time',
                            required=True, type=float)
        parser.add_argument('--prefix',
                            help='Prefix for ouput file names.',
                            required=True)
        parser.add_argument('--matrix', help='Gene expression matrix')
        parser.add_argument('--verbose', action='store_true',
                            help='Print progress information')
        parser.add_argument('--gene', help='List of genes', action='append')
        parser.add_argument('--gene_sets', help='Gene sets')
        parser.add_argument('--gene_set_filter',
                            help='A file with one gene set per line to include or a python regular expression of gene set ids to include')

        parser.add_argument('--cell_sets',
                            help='Grouping of cells into cell sets in gmt or gmx format',
                            required=True)
        parser.add_argument('--cell_set_filter',
                            help='A file with one cell set per line to include or a python regular expression of cell set ids to include')

        args = parser.parse_args(cmd_line)

        time = args.time

        cell_set_ds = wot.io.read_gene_sets(args.cell_sets)
        if args.cell_set_filter is not None:
            cell_set_ds = wot.ot.filter_sets(cell_set_ds, args.cell_set_filter)
        transport_maps = wot.io.list_transport_maps(args.dir)
        full_ds = None
        if args.gene is not None:
            genes = set(np.char.lower(np.array(args.gene)))

            full_ds = wot.io.read_dataset(args.matrix, col_filter={'id': lambda x: x.lower() in genes})

        gene_set_scores = None
        if args.gene_sets is not None:
            gene_set_scores = pd.read_table(args.gene_sets, index_col=0, quoting=csv.QUOTE_NONE, engine='python',
                                            sep=None)
            if args.gene_set_filter is not None:
                import os
                if not os.path.isfile(args.gene_set_filter):
                    import re
                    expr = re.compile(args.gene_set_filter)
                    set_ids = [elem for elem in gene_set_scores.columns if expr.match(elem)]
                else:
                    set_ids = pd.read_table(args.gene_set_filter, index_col=0, header=None).index.values

                gene_set_scores = gene_set_scores[set_ids]

        for t in range(len(transport_maps)):
            if transport_maps[t]['t2'] == time:
                time_index = t
                break

        if time_index is None:
            raise RuntimeError(
                'Transport transport_map for time ' + str(time) + ' not found.')

        return {'cell_set_ds': cell_set_ds, 'transport_maps': transport_maps, 'time': time,
                'full_ds': full_ds,
                'gene_set_scores': gene_set_scores, 'verbose': args.verbose, 'args': args}

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
                        trace = {
                            "name": cell_set_name + ' ' + str(ds.col_meta.index.values[column_index]),
                            "mode": "lines",
                            "type": 'scatter',
                            "y": float(np.mean(array)),
                            "std": float(np.std(array)),
                            "x": t
                        }
                    else:
                        trace = {
                            "group": cell_set_name + ' ' + str(ds.col_meta.index.values[column_index]),
                            "name": str(t),
                            "type": 'violin',
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
                sampling_loader=None, cache=False, ncells=1000):

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
            ranges.append({'backward': True, 'range': range(t2_index, - 1, -1)})
        if t1_index is not None:
            ranges.append({'backward': False, 'range': range(t1_index, len(transport_maps))})
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
                    tmap = tmap_dict.get('ds')
                    if tmap is None:
                        if verbose:
                            print('Reading transport map ' + tmap_dict['path'])
                        tmap = wot.io.read_dataset(tmap_dict['path'])
                        if cache:
                            tmap_dict['ds'] = tmap

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
                                        Ancestors.do_sampling(result=traces, t=time, sampled_indices=None,
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
                        pvecs.append({'v': v, 'entropy': entropy, 't': t,
                                      'cell_ids': tmap.row_meta.index.values if back else tmap.col_meta.index.values})
                        n_choose = ncells if ncells is not None else int(np.ceil(entropy))
                        if verbose:
                            print('Sampling ' + str(n_choose) + ' cells')
                        sampled_indices = np.random.choice(len(v), n_choose, p=v, replace=True)
                    # else:
                    #     sampled_indices = sampling_loader(t=t1, cell_set_name=cell_set_name)
                    # if save_sampling is not None:
                    #     save_sampling(t=t1, cell_set_name=cell_set_name, sampled_indices=sampled_indices)
                    Ancestors.do_sampling(result=traces, t=t, sampled_indices=sampled_indices, datasets=datasets,
                                          summaries=summaries,
                                          cell_set_name=cell_set_ds.col_meta.index.values[cell_set_index], color=color)
                    new_pvec_array.append(v)
                pvec_array = new_pvec_array

        return {'traces': traces, 'pvecs': pvecs}
