#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import wot.ot
import wot.io
import pandas as pd


parser = argparse.ArgumentParser(
    description='Compute cell trajectories')
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

parser.add_argument('--start_cell_set_filter',
                    help='A file with one cell set per line to include or a python regular expression of cell set ids to include')
parser.add_argument('--start_cell_sets',
                    help='Grouping of cells into cell sets in gmt or gmx format',
                    required=True)
parser.add_argument('--end_cell_set_filter',
                    help='A file with one cell set per line to include or a python regular expression of cell set ids to include')
parser.add_argument('--end_cell_sets', help='Grouping of cells into cell sets in gmt or gmx format', required=True)

parser.add_argument('--format', help='Output file format', default='loom')


class TransitionFate:

    @staticmethod
    def read_cell_sets(cell_set_filter, path):
        cell_sets = wot.io.read_gene_sets(path)
        if cell_set_filter is not None:
            return wot.ot.filter_sets(cell_sets, cell_set_filter)
        else:
            return cell_sets

    @staticmethod
    def compute_transport_map(args, store=False):
        start_time = args.start_time
        end_time = args.end_time
        transport_maps = wot.io.list_transport_maps(args.dir)
        start_time_index = None
        end_time_index = None
        for i in range(len(transport_maps)):
            if transport_maps[i]['t1'] == start_time:
                start_time_index = i
            if transport_maps[i]['t2'] == end_time:
                end_time_index = i

        if start_time_index is None and end_time_index is None:
            raise RuntimeError(
                'Transport map for time ' + str(start_time) + ' and time ' + str(end_time) + ' not found.')
        elif start_time_index is None:
            raise RuntimeError(
                'Transport map for time ' + str(start_time) + ' not found.')
        elif end_time_index is None:
            raise RuntimeError(
                'Transport map for time ' + str(end_time) + ' not found.')
        tmap = None
        start_time_ncells = None
        start_time_g = None
        end_time_ncells = None
        tmaps = []
        for i in range(start_time_index, end_time_index + 1):
            ds = wot.io.read_dataset(transport_maps[i]['path'])
            if i == start_time_index:
                start_time_ncells = ds.x.shape[0]
                if ds.row_meta.get('g') is not None:
                    start_time_g = ds.row_meta['g'].values
            elif i == end_time_index:
                end_time_ncells = ds.x.shape[1]
            tmap_i = pd.DataFrame(index=ds.row_meta.index, columns=ds.col_meta.index, data=ds.x)
            if tmap is None:
                tmap = tmap_i
            else:
                tmap = tmap.dot(tmap_i / np.sqrt(np.sum(tmap_i.values)))
            if store:
                tmaps.append(tmap)

        return {'start_time_ncells': start_time_ncells,
                'start_time_g': start_time_g,
                'end_time_ncells': end_time_ncells,
                'tmap': tmap,
                'tmaps': tmaps}

    @staticmethod
    def get_set_id_to_indices(cell_set_ds, tmap, is_columns):
        cell_set_ids = cell_set_ds.col_meta.index.values
        cell_set_id_to_indices = {}

        for i in range(len(cell_set_ids)):
            cell_set_id = cell_set_ids[i]
            cell_ids = cell_set_ds.row_meta.index[cell_set_ds.x[:, i] > 0]
            if is_columns:
                indices = np.where(np.isin(tmap.columns, cell_ids))[0]
            else:
                indices = np.where(tmap.index.isin(cell_ids))[0]
            if len(indices) is 0:
                raise Exception(cell_set_id + ' has zero members in dataset')
            cell_set_id_to_indices[cell_set_id] = indices
        return cell_set_id_to_indices

    @staticmethod
    def summarize_transport_map(args):
        transport_results = TransitionFate.compute_transport_map(args)
        tmap = transport_results['tmap']
        start_time_ncells = transport_results['start_time_ncells']
        start_time_g = transport_results['start_time_g']
        end_time_ncells = transport_results['end_time_ncells']

        start_cell_sets = TransitionFate.read_cell_sets(args.start_cell_set_filter, args.start_cell_sets)
        print(start_cell_sets)
        end_cell_sets = TransitionFate.read_cell_sets(args.end_cell_set_filter, args.end_cell_sets)
        if start_cell_sets.x.shape[1] == 0 and end_cell_sets.x.shape[1] == 0:
            print('No start or end cell sets')
            exit(1)
        elif start_cell_sets.x.shape[1] == 0:
            print('No start cell sets')
            exit(1)
        elif end_cell_sets.x.shape[1] == 0:
            print('No end cell sets')
            exit(1)
        cell_set_id_to_row_indices = TransitionFate.get_set_id_to_indices(start_cell_sets, tmap, False)
        cell_set_id_to_column_indices = TransitionFate.get_set_id_to_indices(end_cell_sets, tmap, True)
        summary = np.zeros(shape=(start_cell_sets.x.shape[1], end_cell_sets.x.shape[1]))

        rids = start_cell_sets.col_meta.index.values
        cids = end_cell_sets.col_meta.index.values
        nrows = start_cell_sets.x.shape[1]
        ncols = end_cell_sets.x.shape[1]

        for i in range(nrows):
            row_indices = cell_set_id_to_row_indices[rids[i]]
            if len(row_indices) == 0:
                continue
            tmap_r = tmap.values[row_indices]
            for j in range(ncols):
                column_indices = cell_set_id_to_column_indices[cids[j]]
                if len(column_indices) == 0:
                    continue
                summary[i, j] += tmap_r[:, column_indices].sum()

        row_meta = pd.DataFrame(index=rids)
        col_meta = pd.DataFrame(index=cids)
        cells_start = np.zeros(nrows)
        g = np.zeros(nrows)
        for i in range(nrows):
            row_indices = cell_set_id_to_row_indices[rids[i]]
            if len(row_indices) == 0:
                continue
            cells_start[i] = len(row_indices)
            if start_time_g is not None:
                g[i] = start_time_g[row_indices].sum()
        g /= start_time_g.sum()
        cells_start /= start_time_ncells
        row_meta['cells_start'] = cells_start
        cells_end = np.zeros(ncols)
        for j in range(ncols):
            column_indices = cell_set_id_to_column_indices[cids[j]]
            if column_indices is None or len(column_indices) == 0:
                continue
            cells_end[j] = len(column_indices)

        cells_end /= end_time_ncells
        col_meta['cells_end'] = cells_end
        if start_time_g is not None:
            row_meta['g'] = g

        tmap_sum = tmap.values.sum()
        row_sums = summary.sum(axis=1)
        row_sums /= tmap_sum
        row_meta['sum'] = row_sums
        # row_filter = row_sums > 0

        column_sums = summary.sum(axis=0)
        column_sums /= tmap_sum
        col_meta['sum'] = column_sums
        # column_filter = column_sums > 0

        # summary = summary[row_filter]
        # summary = summary[:, column_filter]
        # row_meta = row_meta.iloc[row_filter]
        # col_meta = col_meta.iloc[column_filter]

        summary /= tmap_sum

        # tmap.to_csv(prefix + '_' + str(start_time) + '_' + str(end_time) + '_transition.txt', index_label='id', sep='\t',
        #             doublequote=False, quoting=csv.QUOTE_NONE)
        # summary.to_csv(prefix + '_' + str(start_time) + '_' + str(end_time) + '_transition_summary.txt', index_label='id',
        #                sep='\t', doublequote=False, quoting=csv.QUOTE_NONE)

        wot.io.write_dataset(
            wot.Dataset(summary, row_meta, col_meta),
            args.prefix + '_' + str(args.start_time) + '_' + str(args.end_time) + '_transition_summary',
            output_format=args.format, txt_full=True)


TransitionFate.summarize_transport_map(parser.parse_args())
