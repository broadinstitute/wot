#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import wot
import pandas
import numpy as np
import sklearn.metrics.pairwise
import csv

parser = argparse.ArgumentParser(
    description='Compute transport maps between pairs of time points')

parser.add_argument('--matrix',
                    help='Gene expression tab delimited file with cells on '
                         'rows and features on columns', required=True)
parser.add_argument('--cell_growth_rates',
                    help='Two column tab delimited file without header with '
                         'cell ids and growth rates per day.',
                    required=True)
parser.add_argument('--cell_days',
                    help='Two column tab delimited file without header with '
                         'cell ids and days', required=True)
parser.add_argument('--day_pairs',
                    help='Two column tab delimited file without header with '
                         'pairs of days to compute transport maps for',
                    required=True)

parser.add_argument('--epsilon', type=float, default=0.1,
                    help='Controls the entropy of the transport map. An '
                         'extremely large entropy parameter will give a '
                         'maximally entropic transport map, and an '
                         'extremely '
                         'small entropy parameter will give a nearly '
                         'deterministic transport map (but could also '
                         'lead to '
                         'numerical instability in the algorithm')
parser.add_argument('--prefix',
                    help='Prefix for ouput file names', required=True)

parser.add_argument('--max_transport_fraction',
                    default=0.4,
                    help='The maximum fraction of cells at time t that are '
                         'transported to time t + 1',
                    type=float)
parser.add_argument('--min_transport_fraction',
                    default=0.05,
                    help='The minimum fraction of cells at time t that are '
                         'transported to time t + 1',
                    type=float)
parser.add_argument('--lambda1', default=1,
                    help='Regularization parameter that controls the '
                         'fidelity '
                         'of the constraints on p',
                    type=float)
parser.add_argument('--lambda2', default=1,
                    help='Regularization parameter that controls the '
                         'fidelity '
                         'of the constraints on q',
                    type=float)
parser.add_argument('--scaling_iter', default=250,
                    help='Number of scaling iterations', type=int)
parser.add_argument('--min_growth_fit', type=float, default=0.9)
parser.add_argument('--l0_max', type=float, default=100)
parser.add_argument('--clusters',
                    help='Two column tab delimited file without header with '
                         'cell id and cluster id. Used to summarize transport '
                         'maps.',
                    required=False)
parser.add_argument('--cluster_details', action='store_true',
                    help='Save cluster details when clusters is specified.')
parser.add_argument('--subsample_cells', help='Number of cells to sample '
                                              'without '
                                              'replacement.', type=int,
                    action='append')
# parser.add_argument('--subsample_genes', help='Number of genes to sample '
#                                               'without '
#                                               'replacement.', type=int,
#                     action='append')
parser.add_argument('--subsample_iter', help='Number of subsample iterations '
                                             'to perform',
                    type=int, default=0)
parser.add_argument('--no_save', action='store_true',
                    help='Do not save transport maps.')
parser.add_argument('--compress', action='store_true',
                    help='gzip output files')
parser.add_argument('--verbose', action='store_true',
                    help='Print progress information')

args = parser.parse_args()

# cells on rows, features on columns
gene_expression = pandas.read_table(args.matrix, index_col=0,
                                    quoting=csv.QUOTE_NONE)

day_pairs = pandas.read_table(args.day_pairs, header=None, names=['t1', 't2'],
                              index_col=False, quoting=csv.QUOTE_NONE)
days_data_frame = pandas.read_table(args.cell_days, index_col=0, header=None,
                                    names=['day'], quoting=csv.QUOTE_NONE)

cell_growth_rates = pandas.read_table(args.cell_growth_rates, index_col=0,
                                      header=None, names=['cell_growth_rate'],
                                      quoting=csv.QUOTE_NONE)

gene_expression = gene_expression.join(cell_growth_rates).join(days_data_frame)
# cell_growth_rates = cell_growth_rates.align(gene_expression, copy=False,
#                                             join='right')
# days_data_frame = days_data_frame.align(gene_expression, copy=False,
#                                         join='right')

fields_to_drop_for_distance = [days_data_frame.columns[0],
                               cell_growth_rates.columns[0]]

group_by_day = gene_expression.groupby(days_data_frame.columns[0])
if args.verbose:
    print('Computing ' + str(day_pairs.shape[0]) + ' transport map' + 's' if
          day_pairs.shape[0] > 1 else '')

cluster_transport_maps = []
resample = False
subsample_cells = [0]
subsample_genes = [0]
if args.clusters is not None:
    clusters = pandas.read_table(args.clusters, index_col=0, header=None,
                                 names=['cluster'], quoting=csv.QUOTE_NONE)
    grouped_by_cluster = clusters.groupby(clusters.columns[0], axis=0)
    cluster_ids = list(grouped_by_cluster.groups.keys())
    if args.subsample_iter > 0:
        if args.subsample_cells is None:
            print('subsample_cells required when '
                  'subsample_iter > 0')
            exit(1)
        resample = True

        subsample_writer = open(args.prefix + '_subsample_summary.txt', 'w')
        # if args.subsample_genes is not None:
        #     subsample_genes = args.subsample_genes
        #     subsample_writer.write('ngenes\t')
        if args.subsample_cells is not None:
            subsample_cells = args.subsample_cells
            subsample_writer.write('ncells\t')

        subsample_writer.write('t1' + '\t' + 't2' + '\t' + 'standard '
                                                           'deviation mean' +
                               '\n')
column_cell_ids_by_time = []
all_cell_ids = set()

for day_index in range(day_pairs.shape[0]):
    t1 = day_pairs.iloc[day_index, 0]
    t2 = day_pairs.iloc[day_index, 1]
    m1 = group_by_day.get_group(t1)
    m2 = group_by_day.get_group(t2)

    delta_t = t2 - t1
    result = None
    if not resample:
        if args.verbose:
            print(
                'Computing transport map from ' + str(
                    t1) + ' to ' + str(
                    t2))
        unnormalized_cost_matrix = sklearn.metrics.pairwise.pairwise_distances(
            m1.drop(fields_to_drop_for_distance, axis=1),
            Y=m2.drop(fields_to_drop_for_distance, axis=1),
            metric='sqeuclidean')

        cost_matrix = unnormalized_cost_matrix / np.median(
            unnormalized_cost_matrix)
        growth_rate = m1.cell_growth_rate.values
        result = wot.optimal_transport(cost_matrix=cost_matrix,
                                       growth_rate=growth_rate,
                                       delta_days=delta_t,
                                       max_transport_fraction=args.max_transport_fraction,
                                       min_transport_fraction=args.min_transport_fraction,
                                       min_growth_fit=args.min_growth_fit,
                                       l0_max=args.l0_max, lambda1=args.lambda1,
                                       lambda2=args.lambda2,
                                       epsilon=args.epsilon,
                                       scaling_iter=args.scaling_iter)

        transport_map = pandas.DataFrame(result['transport'], index=m1.index,
                                         columns=m2.index)
        if args.verbose:
            print('Done computing transport map')

        if args.clusters is not None:
            cluster_transport_map = wot.transport_map_by_cluster(
                transport_map, grouped_by_cluster, cluster_ids)
            all_cell_ids.update(transport_map.columns)
            all_cell_ids.update(transport_map.index)
            column_cell_ids_by_time.append(transport_map.columns)
            if args.verbose:
                print('Summarized transport map by cluster')
            if args.cluster_details:
                if args.verbose:
                    print('Saving cluster transport map')
                cluster_transport_map.to_csv(
                    args.prefix + '_cluster_' + str(t1) + '_' + str(
                        t2) + '.txt' + ('.gz' if
                    args.compress else ''),
                    index_label="id",
                    sep='\t',
                    compression='gzip' if args.compress
                    else None)
            cluster_transport_maps.append(cluster_transport_map)

        # save the tranport map
        if not args.no_save:
            if args.verbose:
                print('Saving transport map')
            transport_map.to_csv(args.prefix + '_' + str(t1) + '_' + str(
                t2) + '.txt' + ('.gz' if
            args.compress else ''), index_label='id',
                                 sep='\t',
                                 compression='gzip' if
                                 args.compress else None, doublequote=False,
                                 quoting=csv.QUOTE_NONE)

    if resample:
        rnd = np.random.RandomState()
        for ncells in subsample_cells:
            if args.verbose:
                print(
                    'Resampling using ' + str(ncells) + ' cells')
            subsampled_maps = []
            for subsample in range(args.subsample_iter):
                if args.verbose:
                    print('Subsample iteration ' + str(subsample + 1))
                m1_sample = m1.sample(n=ncells, replace=False, axis=0,
                                      random_state=rnd)
                m2_sample = m2.sample(n=ncells, replace=False, axis=0,
                                      random_state=rnd)
                unnormalized_cost_matrix = \
                    sklearn.metrics.pairwise.pairwise_distances(
                        m1_sample.drop(fields_to_drop_for_distance, axis=1),
                        Y=m2_sample.drop(fields_to_drop_for_distance, axis=1),
                        metric='sqeuclidean')
                cost_matrix = unnormalized_cost_matrix / np.median(
                    unnormalized_cost_matrix)
                growth_rate = m1_sample.cell_growth_rate.values
                subsampled_result = wot.optimal_transport(
                    cost_matrix=cost_matrix,
                    growth_rate=growth_rate,
                    delta_days=delta_t,
                    max_transport_fraction=args.max_transport_fraction,
                    min_transport_fraction=args.min_transport_fraction,
                    min_growth_fit=args.min_growth_fit,
                    l0_max=args.l0_max,
                    lambda1=result['lambda1'] if result is not None else
                    args.lambda1,
                    lambda2=result['lambda2'] if result is not None else
                    args.lambda2,
                    epsilon=result['epsilon'] if result is not None else
                    args.epsilon,
                    scaling_iter=args.scaling_iter)

                subsampled_transport_map = pandas.DataFrame(
                    subsampled_result['transport'],
                    index=m1_sample.index,
                    columns=m2_sample.index)
                subsampled_maps.append(wot.transport_map_by_cluster(
                    subsampled_transport_map, grouped_by_cluster, cluster_ids))

            cluster_shape = subsampled_maps[0].shape
            vals = np.zeros(
                (cluster_shape[0], cluster_shape[1], len(subsampled_maps)))
            for subsample_i in range(len(subsampled_maps)):
                subsampled_map = subsampled_maps[subsample_i]
                for i in range(cluster_shape[0]):
                    for j in range(cluster_shape[1]):
                        vals[i, j, subsample_i] = subsampled_map.iloc[i, j]

            stdevs = np.zeros(cluster_shape[0] * cluster_shape[1])
            counter = 0
            for i in range(cluster_shape[0]):
                for j in range(cluster_shape[1]):
                    stdevs[counter] = np.sqrt(np.var(vals[i, j]))
                    counter += 1
            mean_stdev = np.mean(stdevs)
            subsample_writer.write(str(ncells) +
                                   '\t' + str(t1) + '\t' + str(t2) + '\t' + str(
                mean_stdev) + '\n')
            subsample_writer.flush()

if resample:
    subsample_writer.close()
if not resample and args.clusters is not None:
    if args.verbose:
        print('Saving summarized transport map')
    weights = wot.get_weights(all_cell_ids, column_cell_ids_by_time,
                              grouped_by_cluster, cluster_ids)
    cluster_weights_by_time = weights['cluster_weights_by_time']
    combined_cluster_map = wot.transport_maps_by_time(
        cluster_transport_maps,
        cluster_weights_by_time)
    combined_cluster_map.to_csv(
        args.prefix + '_cluster_summary.txt' + ('.gz' if
        args.compress else ''),
        index_label="id",
        sep='\t',
        compression='gzip' if args.compress
        else None)
