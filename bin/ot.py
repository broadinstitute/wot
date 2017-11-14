#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import wot
import pandas
import numpy as np
import sklearn.metrics.pairwise

parser = argparse.ArgumentParser(
    description='Compute transport maps between time points')

parser.add_argument('--expression_file',
                    help='Gene expression tab delimited file with cells on '
                         'rows and '
                         'features on columns',
                    required=True)
parser.add_argument('--growth_file',
                    help='Two column tab delimited file without header with '
                         'cell ids and '
                         'relative growth scores. If not specified, uniform '
                         'growth is assumed.',
                    required=False)
parser.add_argument('--days_file',
                    help='Two column tab delimited file without header with '
                         'cell ids and '
                         'days',
                    required=True)
parser.add_argument('--growth_ratio_file',
                    help='Three column tab delimited file without header with '
                         'timepoint 1, '
                         'timepoint 2, '
                         'and growth ratio. Over 1 day, a cell '
                         'in the more proliferative '
                         'group '
                         'is expected to produce growth_ratio times as '
                         'many '
                         'offspring as a cell in the non-proliferative '
                         'group', required=True)
parser.add_argument('--epsilon', type=float,
                    help='Controls the entropy of the transport map. An '
                         'extremely large entropy parameter will give a '
                         'maximally entropic transport map, and an '
                         'extremely '
                         'small entropy parameter will give a nearly '
                         'deterministic transport map (but could also '
                         'lead to '
                         'numerical instability in the algorithm')
parser.add_argument('--prefix',
                    help='Prefix for ouput file names',
                    required=True)

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
gene_expression = pandas.read_table(args.expression_file, index_col=0)

days_frame = pandas.read_table(args.days_file, index_col=0, header=None,
                               names=['day'])

growth_ratios = pandas.read_table(args.growth_ratio_file,
                                  header=None,
                                  names=['t1', 't2', 'growth_ratio'])
if args.growth_file is not None:
    growth_scores = pandas.read_table(args.growth_file, index_col=0,
                                      header=None,
                                      names=['growth_score'])
else:
    growth_scores = pandas.DataFrame(index=days_frame.index,
                                     columns=['growth_score'],
                                     data=1)

gene_expression = gene_expression.join(growth_scores).join(days_frame)

growth_score_field_name = growth_scores.columns[0]
fields_to_drop_for_distance = [days_frame.columns[0], growth_score_field_name]

group_by_day = gene_expression.groupby(days_frame.columns[0])

if args.verbose:
    print('Computing ' + str(growth_ratios.shape[0]) + ' transport maps...')

for i in range(growth_ratios.shape[0]):
    # consecutive days only
    t1 = growth_ratios.iloc[i, 0]
    t2 = growth_ratios.iloc[i, 1]
    growth_ratio = growth_ratios.iloc[i, 2]
    m1 = group_by_day.get_group(t1)
    m2 = group_by_day.get_group(t2)
    if args.verbose:
        print(
            'Computing transport map from ' + str(
                t1) + ' to ' + str(
                t2) + ' with growth ratio ' + str(growth_ratio))
    delta_t = t2 - t1
    if args.verbose:
        print('delta_t' + str(delta_t))
    cost_matrix = sklearn.metrics.pairwise.pairwise_distances(
        m1.drop(fields_to_drop_for_distance, axis=1),
        Y=m2.drop(fields_to_drop_for_distance, axis=1),
        metric='sqeuclidean')
    cost_matrix = cost_matrix / np.median(cost_matrix)
    if args.verbose:
        print('Computed cost matrix.')
    growth_rate = m1.growth_score.values
    result = wot.optimal_transport(cost_matrix, growth_rate,
                                   delta_days=delta_t,
                                   max_transport_fraction=args.max_transport_fraction,
                                   min_transport_fraction=args.min_transport_fraction,
                                   min_growth_fit=args.min_growth_fit,
                                   l0_max=args.l0_max, lambda1=args.lambda1,
                                   lambda2=args.lambda2,
                                   epsilon=args.epsilon,
                                   growth_ratio=growth_ratio,
                                   scaling_iter=args.scaling_iter)
    transport = pandas.DataFrame(result['transport'], index=m1.index,
                                 columns=m2.index)
    if args.verbose:
        print('Computed transport map.')
    if args.verbose:
        print('Saving transport map')
    transport.to_csv(args.prefix + '_' + str(t1) + '_' + str(
        t2) + '.txt' + ('.gz' if
                        args.compress else ''), index_label='id',
                     sep='\t',
                     compression='gzip' if
                     args.compress else None)
