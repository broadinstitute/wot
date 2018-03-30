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
        # ax = sns.violinplot(x="t", y="values", data=df)
        return sns.factorplot(x="t", y="value", row="cell_set", col='name', data=df, kind="violin")

    @staticmethod
    def from_cmd_line(cmd_line=None):
        parser = argparse.ArgumentParser(
            description='Compute cell ancestors')
        parser.add_argument('--dir',
                            help='Directory of transport maps as produced by ot',
                            required=True)
        parser.add_argument('--start_time',
                            help='The start time',
                            required=True, type=float)
        parser.add_argument('--end_time',
                            help='The end time',
                            required=True, type=float)
        parser.add_argument('--prefix',
                            help='Prefix for ouput file names.',
                            required=True)
        parser.add_argument('--matrix', help='Gene expression matrix')

        parser.add_argument('--gene', help='List of genes', action='append')
        parser.add_argument('--gene_sets', help='Gene sets')

        parser.add_argument('--cell_sets',
                            help='Grouping of cells into cell sets in gmt or gmx format',
                            required=True)
        parser.add_argument('--cell_set_filter',
                            help='A file with one cell set per line to include or a python regular expression of cell set ids to include')

        args = parser.parse_args(cmd_line)

        start_time = args.start_time
        end_time = args.end_time

        cell_set_ds = wot.io.read_gene_sets(args.cell_sets)
        if args.cell_set_filter is not None:
            cell_set_ds = wot.ot.filter_sets(cell_set_ds, args.cell_set_filter)
        transport_maps = wot.io.list_transport_maps(args.dir)

        full_ds = None
        if args.gene is not None:
            full_ds = wot.io.read_dataset(args.matrix)

        gene_set_scores = None
        if args.gene_sets is not None:
            gene_set_scores = pd.read_table(args.gene_sets, index_col=0, quoting=csv.QUOTE_NONE, engine='python',
                                            sep=None)

        start_time_index = None
        end_time_index = None
        for t in range(len(transport_maps)):
            if transport_maps[t]['t1'] == start_time:
                start_time_index = t
            if transport_maps[t]['t2'] == end_time:
                end_time_index = t

        if start_time_index is None:
            raise RuntimeError(
                'Transport transport_map for time ' + str(start_time) + ' not found.')
        if end_time_index is None:
            raise RuntimeError(
                'Transport transport_map for time ' + str(end_time) + ' not found.')
        df = Ancestors.compute(cell_set_ds=cell_set_ds, transport_maps=transport_maps,
                               start_time_index=start_time_index,
                               end_time_index=end_time_index, full_ds=full_ds, gene_set_scores=gene_set_scores,
                               genes=args.gene)

        g = Ancestors.plot(df)
        g.savefig(args.prefix + '.png')

    @staticmethod
    def compute(cell_set_ds, transport_maps, start_time_index, end_time_index, full_ds=None, gene_set_scores=None,
                genes=None):
        list_of_gene_indices = []

        if genes is not None:
            for gene in genes:
                order = full_ds.col_meta.index.get_indexer_for([gene])
                if order[0] != -1:
                    list_of_gene_indices.append((gene, order[0]))
            if len(list_of_gene_indices) == 0:
                print('No genes found')

        n_cell_sets = cell_set_ds.x.shape[1]
        pvec_array = None
        data_table = []
        columns = ['name', 'value', 'cell_set', 't']
        for t in range(end_time_index, start_time_index - 1, -1):
            tmap = wot.io.read_dataset(transport_maps[t]['path'])
            # align ds and tmap
            if full_ds is not None:
                order = full_ds.row_meta.index.get_indexer_for(full_ds.row_meta.index)
                ds = wot.Dataset(full_ds.x[order], full_ds.row_meta.iloc[order], full_ds.col_meta)
            if gene_set_scores is not None:
                # put gene set scores in same order as tmap
                _gene_set_scores = tmap.row_meta.align(gene_set_scores, join='left', axis=0)[1]
            if t == end_time_index:
                pvec_array = []
                for cell_set_index in range(n_cell_sets):
                    cell_ids = cell_set_ds.row_meta.index[cell_set_ds.x[:, cell_set_index] > 0]
                    pvec_array.append(tmap.col_meta.index.isin(cell_ids))
            new_pvec_array = []
            for cell_set_index in range(n_cell_sets):
                v = pvec_array[cell_set_index]
                v = tmap.x.dot(v)
                v /= v.sum()
                entropy = np.exp(scipy.stats.entropy(v))
                cell_ids = tmap.row_meta.index.values
                n = int(entropy)
                sampled_indices = np.random.choice(len(cell_ids), n, p=v, replace=True)
                cell_set_name = cell_set_ds.col_meta.index.values[cell_set_index]
                if full_ds is not None:
                    values = ds.x[sampled_indices]
                    for gene_index in range(len(list_of_gene_indices)):
                        gene = list_of_gene_indices[gene_index]
                        data_table.append([gene[0], values[:, gene[1]], cell_set_name, t])
                if gene_set_scores is not None:
                    tmp_scores = _gene_set_scores.iloc[sampled_indices]
                    for gene_set_index in range(gene_set_scores.shape[1]):
                        vals = tmp_scores.iloc[:, gene_set_index].values
                        for i in range(len(vals)):  # FIXME
                            data_table.append(
                                [gene_set_scores.columns[gene_set_index], vals[i],
                                 cell_set_name,
                                 t])
                new_pvec_array.append(v)
            pvec_array = new_pvec_array

        return pd.DataFrame(data_table, columns=columns)


if __name__ == 'main':
    Ancestors.from_cmd_line()
