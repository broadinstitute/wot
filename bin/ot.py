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
parser.add_argument('--compress', action='store_true',
                    help='gzip output files')
parser.add_argument('--verbose', action='store_true',
                    help='Print progress information')

args = parser.parse_args()

# cells on rows, features on columns
gene_expression = pandas.read_table(args.matrix, index_col=0)

day_pairs = pandas.read_table(args.day_pairs, header=None, names=['t1', 't2'],
                              index_col=False)
days_data_frame = pandas.read_table(args.cell_days, index_col=0, header=None,
                                    names=['day'])

cell_growth_rates = pandas.read_table(args.cell_growth_rates, index_col=0,
                                      header=None, names=['cell_growth_rate'])

gene_expression = gene_expression.join(cell_growth_rates).join(days_data_frame)

fields_to_drop_for_distance = [days_data_frame.columns[0],
                               cell_growth_rates.columns[0]]

group_by_day = gene_expression.groupby(days_data_frame.columns[0])

if args.verbose:
    print('Computing ' + str(day_pairs.shape[0]) + ' transport maps')

for i in range(day_pairs.shape[0]):
    # consecutive days only
    t1 = day_pairs.iloc[i, 0]
    t2 = day_pairs.iloc[i, 1]
    m1 = group_by_day.get_group(t1)
    m2 = group_by_day.get_group(t2)
    if args.verbose:
        print(
            'Computing transport map from ' + str(
                t1) + ' to ' + str(
                t2))
    delta_t = t2 - t1
    cost_matrix = sklearn.metrics.pairwise.pairwise_distances(
        m1.drop(fields_to_drop_for_distance, axis=1),
        Y=m2.drop(fields_to_drop_for_distance, axis=1),
        metric='sqeuclidean')
    cost_matrix = cost_matrix / np.median(cost_matrix)
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

    transport = pandas.DataFrame(result['transport'], index=m1.index,
                                 columns=m2.index)
    if args.verbose:
        print('Done computing transport map')
        print('Saving transport map')

    transport.to_csv(args.prefix + '_' + str(t1) + '_' + str(
        t2) + '.txt' + ('.gz' if
                        args.compress else ''), index_label='id',
                     sep='\t',
                     compression='gzip' if
                     args.compress else None, doublequote=False,
                     quoting=csv.QUOTE_NONE)
