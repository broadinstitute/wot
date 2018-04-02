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
    def from_cmd_line(cmd_line=None, save_image=True):
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
        n_cell_sets = cell_set_ds.x.shape[1]
        full_ds = None
        if args.gene is not None:
            full_ds = wot.io.read_dataset(args.matrix)

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

        start_time_index = None
        time_index = None
        end_time_index = None
        min_time = np.inf
        max_time = -np.inf
        start_time = None
        end_time = None
        for t in range(len(transport_maps)):
            if start_time is None and transport_maps[t]['t1'] < min_time:
                min_time = transport_maps[t]['t1']
                start_time_index = t
            elif transport_maps[t]['t1'] == start_time:
                start_time_index = t

            if end_time is None and transport_maps[t]['t2'] > max_time:
                max_time = transport_maps[t]['t2']
                end_time_index = t
            elif transport_maps[t]['t2'] == end_time:
                end_time_index = t

            if transport_maps[t]['t2'] == time:
                time_index = t

        if time_index is None or end_time_index is None or start_time_index is None:
            raise RuntimeError(
                'Transport transport_map for time ' + str(time) + ' not found.')

        df = Ancestors.compute(cell_set_ds=cell_set_ds, transport_maps=transport_maps,
                               t2_index=start_time_index, time_index=time_index, full_ds=full_ds,
                               gene_set_scores=gene_set_scores,
                               genes=args.gene, verbose=args.verbose)

        g = Ancestors.plot(df)
        if save_image:
            g.savefig(args.prefix + '.png')
        return g

    @staticmethod
    def compute(cell_set_ds, transport_maps, time_index, t2_index, full_ds=None,
                gene_set_scores=None,
                genes=None, verbose=False, save_sampling=False):
        list_of_gene_indices = []

        if genes is not None:
            for gene in genes:
                order = full_ds.col_meta.index.get_indexer_for([gene])
                if order[0] != -1:
                    list_of_gene_indices.append((gene, order[0]))
            if len(list_of_gene_indices) == 0:
                print('No genes found')

        df_names = np.array([])
        df_cell_set_names = np.array([])
        df_times = np.array([])
        df_vals = np.array([])
        n_cell_sets = cell_set_ds.x.shape[1] if cell_set_ds is not None else 0
        load_sampling = False
        if save_sampling:
            import h5py
            saved_sampling_file = h5py.File('saved_sampling.h5', 'w')
        elif load_sampling:
            import h5py
            saved_sampling_file = h5py.File('saved_sampling.h5', 'r')
        for t in range(time_index, t2_index - 1, -1) if time_index > t2_index else range(time_index, t2_index + 1):
            if verbose:
                print('Reading transport map ' + transport_maps[t]['path'])
            t1 = transport_maps[t]['t1']
            if load_sampling:
                f = h5py.File(transport_maps[t]['path'], 'r')
                ids = f['/row_attrs/id'][()]
                if ids.dtype.kind == 'S':
                    ids = ids.astype(str)
                tmap = wot.Dataset(None, pd.DataFrame(index=ids), None)
                f.close()
            else:
                tmap = wot.io.read_dataset(transport_maps[t]['path'])
            # align ds and tmap
            if full_ds is not None:
                ds_order = tmap.row_meta.index.get_indexer_for(full_ds.row_meta.index)
                ds = wot.Dataset(full_ds.x[ds_order], full_ds.row_meta.iloc[ds_order], full_ds.col_meta)
            if gene_set_scores is not None:
                # put gene set scores in same order as tmap
                _gene_set_scores = tmap.row_meta.align(gene_set_scores, join='left', axis=0)[1]
            if t == time_index:
                pvec_array = []
                cell_sets_to_keep = []
                if verbose:
                    print('Initializing cell sets')
                for cell_set_index in range(n_cell_sets):
                    cell_ids = cell_set_ds.row_meta.index[cell_set_ds.x[:, cell_set_index] > 0]
                    membership = tmap.col_meta.index.isin(cell_ids)
                    if np.sum(membership) > 0:
                        pvec_array.append(membership)
                        cell_sets_to_keep.append(cell_set_index)
                cell_set_ds = wot.Dataset(cell_set_ds.x[:, cell_sets_to_keep], cell_set_ds.row_meta,
                                          cell_set_ds.col_meta.iloc[cell_sets_to_keep])
                n_cell_sets = cell_set_ds.x.shape[1]
            new_pvec_array = []
            for cell_set_index in range(n_cell_sets):
                cell_set_name = cell_set_ds.col_meta.index.values[cell_set_index]
                if not load_sampling:
                    v = pvec_array[cell_set_index]
                    v = tmap.x.dot(v)
                    v /= v.sum()
                    entropy = np.exp(scipy.stats.entropy(v))
                    cell_ids = tmap.row_meta.index.values
                    n = int(np.ceil(entropy))
                    if verbose:
                        print('Sampling ' + str(n) + ' cells')
                    sampled_indices = np.random.choice(len(cell_ids), n, p=v, replace=True)
                else:
                    sampled_indices = saved_sampling_file[str(t1) + '_' + str(cell_set_name)][()]
                if save_sampling:
                    saved_sampling_file[str(t1) + '_' + str(cell_set_name)] = sampled_indices
                if full_ds is not None:
                    values = ds.x[sampled_indices]
                    for gene_index in range(len(list_of_gene_indices)):
                        gene = list_of_gene_indices[gene_index]
                        df_vals = np.concatenate((df_vals, values[:, gene[1]]))
                        df_names = np.concatenate((df_names, np.repeat(gene[0], n)))
                        df_cell_set_names = np.concatenate((df_cell_set_names, np.repeat(cell_set_name, n)))
                        df_times = np.concatenate((df_times, np.repeat(t1, n)))

                if gene_set_scores is not None:
                    tmp_scores = _gene_set_scores.iloc[sampled_indices]
                    for gene_set_index in range(gene_set_scores.shape[1]):
                        vals = tmp_scores.iloc[:, gene_set_index].values
                        gene_set_name = gene_set_scores.columns[gene_set_index]
                        df_vals = np.concatenate((df_vals, vals))
                        df_names = np.concatenate((df_names, np.repeat(gene_set_name, n)))
                        df_cell_set_names = np.concatenate((df_cell_set_names, np.repeat(cell_set_name, n)))
                        df_times = np.concatenate((df_times, np.repeat(t1, n)))

                new_pvec_array.append(v)
            pvec_array = new_pvec_array

            if save_sampling:
                saved_sampling_file.close()

        return pd.DataFrame(data={'name': df_names,
                                  'cell_set': df_cell_set_names,
                                  'value': df_vals,
                                  't': df_times})


if __name__ == 'main':
    Ancestors.from_cmd_line()
