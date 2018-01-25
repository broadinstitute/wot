#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import wot.ot
import pandas as pd
import numpy as np
import sklearn.metrics.pairwise
import csv
import ot as pot

# from gslrandom import PyRNG, multinomial
# def coupling_sampler(Lineage, nf=1e-3, s=1, threads=1, nmin=10):
#     Pairs = [[] for _ in range(s)]
#     for lineage in Lineage:
#         if len(lineage) > 0:
#             l = lineage / lineage.sum(0)
#             l = l / l.sum()
#             l = l.flatten()
#             sd = np.exp(scipy.stats.entropy(l))
#             n = max(nmin, int(sd * nf)) * np.ones(s, dtype=np.uint32)
#             P = np.ones((s, len(l)), dtype=np.uint32)
#             l_tile = np.tile(l, (s, 1))
#             rngs = [PyRNG(np.random.randint(2 ** 16)) for _ in range(threads)]
#             P = multinomial(rngs, n, l_tile, P)
#             # Too slow: P = np.random.multinomial(n,l,size=s)
#             for i in range(s):
#                 pairs = np.nonzero(P[i].reshape(lineage.shape))
#                 Pairs[i].append(pairs)
#             del P, l_tile
#         else:
#             for i in range(s):
#                 Pairs[i].append([])
#     return Pairs

# more than 10 to the minus 8
# return indices and

tm_thresh = 1e-8


def sample_randomly(exp1, exp2, tm):
    p = np.ones(exp1.shape[0] * exp2.shape[0])
    p = p / p.sum()
    l = tm / tm.sum(axis=0)
    l = l.flatten()
    qq = np.where(l > tm_thresh)
    pairs = np.random.multinomial(len(qq[0]), p, size=1)
    pairs = np.nonzero(pairs.reshape(exp1.shape[0], exp2.shape[0]))
    return exp1[pairs[0]], exp2[pairs[1]]


# def sample_from_transport_map(exp1, exp2, transport_map):
#     tm = transport_map / transport_map.sum(axis=0)
#     l = tm.flatten()
#     l = l / l.sum()
#
#     pairs = np.random.multinomial(args.npairs, l, size=1)
#     pairs = np.nonzero(pairs.reshape(exp1.shape[0], exp2.shape[0]))
#     return exp1[pairs[0]], exp2[pairs[1]]

def sample_from_transport_map(exp1, exp2, tm):
    l = tm / tm.sum(axis=0)
    l = l.flatten()
    # l = l / l.sum()
    reshaped_tmap = l.reshape(exp1.shape[0], exp2.shape[0])
    q = np.where(reshaped_tmap > tm_thresh)
    qq = np.where(l > tm_thresh)
    return exp1[q[0]], exp2[q[1]], l[qq]


def split_in_two(n):
    indices = np.random.choice(n, int(n * 0.5))
    indices_c = np.zeros(n, dtype=bool)
    indices_c[indices] = True
    indices_c = np.invert(indices_c)
    return indices, indices_c


def point_cloud_distance(c1, c2, a=None, b=None):
    cloud_distances = sklearn.metrics.pairwise.pairwise_distances(c1, Y=c2, metric='sqeuclidean')
    cloud_distances = cloud_distances / np.median(cloud_distances)

    if a is None:
        a = np.ones((cloud_distances.shape[0]), dtype=np.float64) / cloud_distances.shape[0]
    else:
        a = a / a.sum()
    if b is None:
        b = np.ones((cloud_distances.shape[1]), dtype=np.float64) / cloud_distances.shape[1]
    else:
        b = b / b.sum()
    return np.sqrt(pot.emd2(a, b, cloud_distances, numItermax=max(10000000, c1.shape[0] * c2.shape[0])))


parser = argparse.ArgumentParser(
    description='Compute transport maps between pairs of time points')

parser.add_argument('--matrix',
                    help='Gene expression tab delimited file with cells on '
                         'rows and features on columns', required=True)

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

# parser.add_argument('--subsample_genes', help='Number of genes to sample '
#                                               'without '
#                                               'replacement.', type=int,
#                     action='append')

parser.add_argument('--no_save', action='store_true',
                    help='Do not save transport maps.')
parser.add_argument('--compress', action='store_true',
                    help='gzip output files')

parser.add_argument('--epsilon_adjust', help='Scaling factor to adjust epsilon',
                    type=float, default=1.1)
parser.add_argument('--lambda_adjust', help='Scaling factor to adjust lambda',
                    type=float, default=1.5)

parser.add_argument('--beta_max', type=float, default=1.7)
parser.add_argument('--beta_center', type=float, default=0.25)
parser.add_argument('--delta_max', type=float, default=1.7)

parser.add_argument('--numItermax', type=int, default=100, help='For sinkhorn_epsilon solver')
parser.add_argument('--epsilon0', type=float, default=1, help='For sinkhorn_epsilon solver')
parser.add_argument('--numInnerItermax', type=int, default=10, help='For sinkhorn_epsilon solver')
parser.add_argument('--tau', type=float, default=100000, help='For sinkhorn_epsilon solver')
parser.add_argument('--stopThr', type=float, default=1e-10, help='For sinkhorn_epsilon solver')

growth_rate_group = parser.add_mutually_exclusive_group(required=True)
growth_rate_group.add_argument('--gene_set_scores',
                               help='File containing "Cell.cycle" and "Apoptosis" scores')
growth_rate_group.add_argument('--cell_growth_rates',
                               help='Two column tab delimited file without '
                                    'header with '
                                    'cell ids and growth rates per day.')
parser.add_argument('--diagonal', help='Diagonal scaling matrix')
parser.add_argument('--power', help='Diagonal scaling power', type=float)

parser.add_argument('--subsample_iter', help='Number of subsample iterations to perform', type=int, default=0)
parser.add_argument('--solver',
                    help='Solver to use when computing transport maps. One of epsilon, sinkhorn_epsilon, unregularized',
                    choices=['epsilon', 'sinkhorn_epsilon', 'unregularized'])
parser.add_argument('--t_interpolate', help='Interpolation fraction between two time points', type=float)

parser.add_argument('--verbose', action='store_true',
                    help='Print progress information')
args = parser.parse_args()
eigenvals = None
solver = args.solver
print(solver)
if args.diagonal is not None:
    eigenvals = np.loadtxt(args.diagonal, delimiter='\n')
if eigenvals is not None and args.power is not None:
    eigenvals = np.power(eigenvals, args.power)

# cells on rows, features on columns
gene_expression = pd.read_table(args.matrix, index_col=0,
                                quoting=csv.QUOTE_NONE, engine='python',
                                sep=None)

day_pairs = pd.read_table(args.day_pairs, header=None, names=['t1', 't2'],
                          index_col=False, quoting=csv.QUOTE_NONE,
                          engine='python', sep=None, dtype=np.float32)
days_data_frame = pd.read_table(args.cell_days, index_col=0, header=None,
                                names=['day'], quoting=csv.QUOTE_NONE,
                                engine='python', sep=None,
                                dtype={'day': np.float32})

gene_set_scores = None
gene_set_sigmas = None

if eigenvals is not None:
    gene_expression = gene_expression.dot(np.diag(eigenvals))

params_writer = None
if solver is 'epsilon':
    params_writer = open(args.prefix + '_params.txt', 'w')
    params_writer.write('t1' + '\t' + 't2' + '\t' + 'epsilon' + '\t' + 'lambda1' + '\t' + 'lambda2'
                                                                                          '\n')

if args.gene_set_scores is not None:
    gene_set_scores = pd.read_table(args.gene_set_scores, index_col=0,
                                    quoting=csv.QUOTE_NONE, engine='python',
                                    sep=None)
    apoptosis = gene_set_scores['Apoptosis']
    proliferation = gene_set_scores['Cell.cycle']
    g = wot.ot.compute_growth_scores(proliferation.values, apoptosis.values, beta_max=args.beta_max,
                                     beta_center=args.beta_center,
                                     delta_max=args.delta_max)
    cell_growth_rates = pd.DataFrame(index=gene_set_scores.index,
                                     data={'cell_growth_rate': g})

else:
    cell_growth_rates = pd.read_table(args.cell_growth_rates, index_col=0,
                                      header=None, names=['cell_growth_rate'],
                                      quoting=csv.QUOTE_NONE, engine='python',
                                      sep=None)
fields_to_drop_for_distance = [days_data_frame.columns[0],
                               cell_growth_rates.columns[0]]

gene_expression = gene_expression.join(cell_growth_rates).join(days_data_frame)
group_by_day = gene_expression.groupby(days_data_frame.columns[0])
if args.verbose:
    print('Computing ' + str(day_pairs.shape[0]) + ' transport map' + 's' if
          day_pairs.shape[0] > 1 else '')

cluster_transport_maps = []
resample = False
total_cluster_size = None
if args.clusters is not None:
    clusters = pd.read_table(args.clusters, index_col=0, header=None,
                             names=['cluster'], quoting=csv.QUOTE_NONE,
                             engine='python',
                             sep=None)
    clusters = clusters.align(gene_expression, join='right', axis=0,
                              copy=False)[0]
    grouped_by_cluster = clusters.groupby(clusters.columns[0], axis=0)
    cluster_ids = list(grouped_by_cluster.groups.keys())
    total_cluster_size = wot.ot.get_column_weights(clusters.index,
                                                   grouped_by_cluster,
                                                   cluster_ids)

subsample_writer = None
if args.t_interpolate is not None:
    resample = True
    subsample_writer = open(args.prefix + '_subsample_summary.txt', 'w')
    subsample_writer.write('t1' + '\t' + 't2' + '\t' + 't_interpolate' + '\t' + 'distance' + '\t' + 'type' + '\n')
column_cell_ids_by_time = []
all_cell_ids = set()

for day_index in range(day_pairs.shape[0]):
    t1 = day_pairs.iloc[day_index, 0]
    t2 = day_pairs.iloc[day_index, 1]
    if group_by_day.groups.get(t1) is None or group_by_day.groups.get(
            t2) is None:
        print('skipping transport map from ' + str(t1) + ' to ' + str(t2))
        continue
    m1 = group_by_day.get_group(t1)
    m2 = group_by_day.get_group(t2)
    delta_t = t2 - t1
    if args.verbose:
        print(
            'Computing transport map from ' + str(
                t1) + ' to ' + str(
                t2) + '...', end='')

    c = sklearn.metrics.pairwise.pairwise_distances(m1.drop(fields_to_drop_for_distance, axis=1).values,
                                                    m2.drop(fields_to_drop_for_distance, axis=1).values,
                                                    metric='sqeuclidean')
    c = c / np.median(c)

    result = wot.ot.optimal_transport(cost_matrix=c,
                                      growth_rate=m1[
                                          cell_growth_rates.columns[0]].values,
                                      delta_days=delta_t,
                                      max_transport_fraction=args.max_transport_fraction,
                                      min_transport_fraction=args.min_transport_fraction,
                                      min_growth_fit=args.min_growth_fit,
                                      l0_max=args.l0_max, lambda1=args.lambda1,
                                      lambda2=args.lambda2,
                                      epsilon=args.epsilon,
                                      scaling_iter=args.scaling_iter,
                                      epsilon_adjust=args.epsilon_adjust,
                                      lambda_adjust=args.lambda_adjust,
                                      numItermax=args.numItermax,
                                      epsilon0=args.epsilon0,
                                      numInnerItermax=args.numInnerItermax, tau=args.tau, stopThr=args.stopThr,
                                      solver=solver)
    if solver is 'epsilon':
        params_writer.write(
            str(t1) + '\t' + str(t2) + '\t' + str(result['epsilon']) + '\t' + str(
                result['lambda1']) + '\t' + str(
                result['lambda2']) + '\n')
    transport_map = pd.DataFrame(result['transport'], index=m1.index,
                                 columns=m2.index)
    if args.verbose:
        print('done')

    if args.clusters is not None:
        cluster_transport_map = wot.ot.transport_map_by_cluster(
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

    if resample:  # resample and optionally perturb parameters
        rnd = np.random.RandomState(123125)


        def write_point_cloud_distance(point_cloud1, point_cloud2, weights1, weights2, point_cloud1_name,
                                       point_cloud2_name):
            if args.verbose:
                print('Computing distance between ' + point_cloud1_name + ' and ' + point_cloud2_name + '...', end='')
            subsample_writer.write(
                str(t1) + '\t' + str(t2) + '\t' + str(args.t_interpolate) + '\t' + str(
                    point_cloud_distance(point_cloud1,
                                         point_cloud2, weights1,
                                         weights2)) + '\t' + point_cloud1_name + ' vs ' + point_cloud2_name + '\n')
            subsample_writer.flush()
            if args.verbose:
                print('done')


        inferred_time = t1 + (t2 - t1) * args.t_interpolate
        actual_mtx = group_by_day.get_group(inferred_time)
        actual_mtx = actual_mtx.drop(fields_to_drop_for_distance, axis=1).values

        m1_mtx = m1.drop(fields_to_drop_for_distance, axis=1).values
        m2_mtx = m2.drop(fields_to_drop_for_distance, axis=1).values

        m1_subset, m2_subset, m1_m2_subset_weights = sample_from_transport_map(m1_mtx, m2_mtx,
                                                                               result['transport'])
        m1_random_subset, m2_random_subset = sample_randomly(m1_mtx, m2_mtx, result['transport'])
        inferred = m1_subset + args.t_interpolate * (m2_subset - m1_subset)

        random_inferred = m1_random_subset + args.t_interpolate * (m2_random_subset - m1_random_subset)

        point_clouds = [[m1_mtx, None], [m2_mtx, None], [actual_mtx, None], [inferred, m1_m2_subset_weights],
                        [random_inferred, None]]
        point_cloud_names = ['P0', 'P1', 'P' + str(args.t_interpolate), 'I' + str(args.t_interpolate),
                             'R' + str(args.t_interpolate)]
        for ix in range(1, len(point_cloud_names)):
            for iy in range(ix):
                write_point_cloud_distance(point_clouds[ix][0], point_clouds[iy][0], point_clouds[ix][1],
                                           point_clouds[iy][1],
                                           point_cloud_names[ix], point_cloud_names[iy])

            # self point cloud distances by splitting point clouds in 2
        for point_cloud_idx in [0, 1, 2, 4]:
            cloud = point_clouds[point_cloud_idx][0]
            split1, split2 = split_in_two(cloud.shape[0])
            write_point_cloud_distance(cloud[split1], cloud[split2], None, None,
                                       point_cloud_names[point_cloud_idx], point_cloud_names[point_cloud_idx])

        for subsample_iter in range(args.subsample_iter):
            if args.verbose:
                print('Subsample iteration ' + str(subsample_iter + 1))

            # compare pairs
            interpolated_matrices = []
            interpolated_weights = []
            m1_indices_ = split_in_two(m1.shape[0])
            m2_indices_ = split_in_two(m2.shape[0])
            for split_iter in range(2):
                m1_sample = m1.iloc[m1_indices_[split_iter]]
                m2_sample = m2.iloc[m2_indices_[split_iter]]
                # if gene_set_scores is not None and gene_set_sigmas is not None:
                #     _gene_set_scores = gene_set_scores.loc[m1.index]
                #
                #     apoptosis = _gene_set_scores['Apoptosis'] + np.random.normal(0,
                #                                                                  sigma,
                #                                            o                      _gene_set_scores.shape[
                #                                                                      0])
                #     proliferation = _gene_set_scores[
                #                         'Proliferation'] + np.random.normal(
                #         0,
                #         sigma,
                #         _gene_set_scores.shape[0])
                #     g = wot.ot.compute_growth_scores(proliferation.values,
                #                                      apoptosis.values)
                m1_mtx = m1_sample.drop(fields_to_drop_for_distance, axis=1).values
                m2_mtx = m2_sample.drop(fields_to_drop_for_distance, axis=1).values
                c = sklearn.metrics.pairwise.pairwise_distances(m1_mtx, Y=m2_mtx, metric='sqeuclidean')
                c = c / np.median(c)

                perturbed_result = wot.ot.optimal_transport(
                    cost_matrix=c,
                    growth_rate=m1_sample[
                        cell_growth_rates.columns[0]].values,
                    delta_days=delta_t,
                    max_transport_fraction=args.max_transport_fraction,
                    min_transport_fraction=args.min_transport_fraction,
                    min_growth_fit=args.min_growth_fit,
                    l0_max=args.l0_max, lambda1=args.lambda1,
                    lambda2=args.lambda2,
                    epsilon=args.epsilon,
                    scaling_iter=args.scaling_iter,
                    epsilon_adjust=args.epsilon_adjust,
                    lambda_adjust=args.lambda_adjust, numItermax=args.numItermax,
                    epsilon0=args.epsilon0,
                    numInnerItermax=args.numInnerItermax, tau=args.tau, stopThr=args.stopThr, solver=solver)
                perturbed_transport = perturbed_result['transport']

                m1_mtx, m2_mtx, m1_m2_weights = sample_from_transport_map(m1_mtx, m2_mtx, perturbed_transport)

                interpolated_matrices.append(m1_mtx + args.t_interpolate * (m2_mtx - m1_mtx))
                interpolated_weights.append(m1_m2_weights)
            write_point_cloud_distance(interpolated_matrices[0], interpolated_matrices[1], None, None, 'I\'', 'I\'')
            for i in range(len(interpolated_matrices)):
                write_point_cloud_distance(interpolated_matrices[i], actual_mtx, interpolated_weights[i], None,
                                           'P' + str(args.t_interpolate), 'I\'')

        # cluster_shape = subsampled_maps[0].shape
        # vals = np.zeros(
        #     (cluster_shape[0], cluster_shape[1], len(subsampled_maps)))
        # for subsample_i in range(len(subsampled_maps)):
        #     subsampled_map = subsampled_maps[subsample_i]
        #     for i in range(cluster_shape[0]):
        #         for j in range(cluster_shape[1]):
        #             vals[i, j, subsample_i] = subsampled_map.iloc[i, j]
        #
        # stdevs = np.zeros(cluster_shape[0] * cluster_shape[1])
        # counter = 0
        # for i in range(cluster_shape[0]):
        #     for j in range(cluster_shape[1]):
        #         stdevs[counter] = np.sqrt(np.var(vals[i, j]))
        #         counter += 1
        # mean_stdev = np.mean(stdevs)
        # subsample_writer.write(
        #     str(m1.shape[0]) + '\t' + str(m2.shape[0]) + '\t' + str(n) +
        #     '\t' + str(t1) + '\t' + str(t2) + '\t' + str(
        #         mean_stdev) + '\n')
        # subsample_writer.flush()

    if subsample_writer is not None:
        subsample_writer.close()

    if params_writer is not None:
        params_writer.close()
    if not resample and args.clusters is not None:
        if args.verbose:
            print('Saving summarized transport map')
        weights = wot.ot.get_weights(all_cell_ids, column_cell_ids_by_time,
                                     grouped_by_cluster, cluster_ids)
        cluster_weights_by_time = weights['cluster_weights_by_time']
        combined_cluster_map = wot.ot.transport_maps_by_time(
            cluster_transport_maps,
            cluster_weights_by_time)
        combined_cluster_map.to_csv(
            args.prefix + '_cluster_summary.txt' + ('.gz' if
            args.compress else ''),
            index_label="id",
            sep='\t',
            compression='gzip' if args.compress
            else None)
