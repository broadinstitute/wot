#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import wot.ot
import pandas as pd
import numpy as np
import sklearn.metrics.pairwise
import csv
import ot as pot
import os
import io


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


def sample_randomly(exp1, exp2, tm, g):
    if args.npairs is None or args.npairs <= 0:
        l = tm / tm.sum(axis=0)
        l = l.flatten()
        qq = np.where(l > args.tm_thresh)
        n = len(qq[0])
    else:
        n = args.npairs

    # column_sums = np.sum(tm, axis=0)
    # row_sums = np.sum(tm, axis=1)
    # p = np.outer(row_sums, column_sums)

    p = g
    q = np.ones(exp2.shape[0]) * np.average(g)
    p = np.outer(p, q)
    p = p.flatten()
    p = p / p.sum()
    pairs = np.random.multinomial(n, p, size=1)
    pairs = np.nonzero(pairs.reshape(exp1.shape[0], exp2.shape[0]))
    return exp1[pairs[0]], exp2[pairs[1]]


def sample_uniformly(exp1, exp2, tm):
    if args.npairs is None or args.npairs <= 0:
        l = tm / tm.sum(axis=0)
        l = l.flatten()
        qq = np.where(l > args.tm_thresh)
        n = len(qq[0])
    else:
        n = args.npairs

    p = np.ones(exp1.shape[0] * exp2.shape[0])
    p = p / p.sum()
    pairs = np.random.multinomial(n, p, size=1)
    pairs = np.nonzero(pairs.reshape(exp1.shape[0], exp2.shape[0]))
    return exp1[pairs[0]], exp2[pairs[1]]


def sample_from_transport_map(exp1, exp2, tm):
    if args.npairs is None or args.npairs <= 0:
        l = tm / tm.sum(axis=0)
        l = l.flatten()
        # l = l / l.sum()
        reshaped_tmap = l.reshape(exp1.shape[0], exp2.shape[0])
        q = np.where(reshaped_tmap > args.tm_thresh)
        qq = np.where(l > args.tm_thresh)
        return exp1[q[0]], exp2[q[1]], l[qq]
    else:
        l = tm / tm.sum(axis=0)
        l = l.flatten()
        l = l / l.sum()
        pairs = np.random.multinomial(args.npairs, l, size=1)
        pairs = np.nonzero(pairs.reshape(exp1.shape[0], exp2.shape[0]))
        return exp1[pairs[0]], exp2[pairs[1]], None


def split_in_two(n):
    indices = np.random.choice(n, int(n * 0.5))
    indices_c = np.zeros(n, dtype=bool)
    indices_c[indices] = True
    indices_c = np.invert(indices_c)
    return indices, indices_c


def point_cloud_distance(c1, c2, a=None, b=None):
    cloud_distances = sklearn.metrics.pairwise.pairwise_distances(c1, Y=c2, metric='sqeuclidean')

    if a is None:
        a = np.ones((cloud_distances.shape[0]), dtype=np.float64) / cloud_distances.shape[0]
    else:
        a = a / a.sum()
    if b is None:
        b = np.ones((cloud_distances.shape[1]), dtype=np.float64) / cloud_distances.shape[1]
    else:
        b = b / b.sum()
    return np.sqrt(pot.emd2(a, b, cloud_distances, numItermax=10000000))


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

parser.add_argument('--covariate',
                    help='Two column tab delimited file without header with '
                         'cell ids and covariate value')

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

parser.add_argument('--no_save', action='store_true',
                    help='Do not save transport maps.')
parser.add_argument('--compress', action='store_true',
                    help='gzip output files')

parser.add_argument('--epsilon_adjust', help='Scaling factor to adjust epsilon for floating_epsilon solver',
                    type=float, default=1.1)
parser.add_argument('--lambda_adjust', help='Scaling factor to adjust lambda for floating_epsilon solver',
                    type=float, default=1.5)

parser.add_argument('--beta_max', type=float, default=1.7, help='Growth function parameter')
parser.add_argument('--beta_center', type=float, default=0.25, help='Growth function parameter')
parser.add_argument('--delta_max', type=float, default=1.7, help='Growth function parameter')

parser.add_argument('--numItermax', type=int, default=100, help='For sinkhorn_epsilon solver')
parser.add_argument('--epsilon0', type=float, default=1, help='For sinkhorn_epsilon and unbalanced solvers')
parser.add_argument('--numInnerItermax', type=int, default=10, help='For sinkhorn_epsilon and unbalanced solvers')
parser.add_argument('--tau', type=float, default=100000, help='For sinkhorn_epsilon and unbalanced solvers')
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

parser.add_argument('--subsample_iter', help='Number of subsample iterations to perform', type=int, default=2)
parser.add_argument('--resample_iter', help='Number of resample iterations to perform', type=int, default=5)

parser.add_argument('--solver',
                    help='Solver to use when computing transport maps. One of unbalanced, floating_epsilon, '
                         'sinkhorn_epsilon, unregularized',
                    choices=['epsilon', 'sinkhorn_epsilon', 'unbalanced', 'unregularized'], default='sinkhorn_epsilon')
parser.add_argument('--t_interpolate', help='Interpolation fraction between two time points', type=float)
parser.add_argument('--verbose', action='store_true',
                    help='Print progress information')
parser.add_argument('--quick', action='store_true')
parser.add_argument('--npairs', type=int)
parser.add_argument('--tm_thresh', type=float, default=1e-9)
args = parser.parse_args()
eigenvals = None
solver = args.solver

if args.diagonal is not None:
    eigenvals = np.loadtxt(args.diagonal, delimiter='\n')
if eigenvals is not None and args.power is not None:
    eigenvals = np.power(eigenvals, args.power)

# cells on rows, features on columns
gene_expression = pd.read_table(args.matrix, index_col=0,
                                quoting=csv.QUOTE_NONE, engine='python',
                                sep=None)

if not os.path.isfile(args.day_pairs):
    day_pairs = pd.read_table(io.StringIO(args.day_pairs), header=None, names=['t1', 't2'],
                              index_col=False, lineterminator=';', sep=',', dtype=np.float32)
else:
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
if solver is 'floating_epsilon':
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
fields_to_drop_for_distance = [days_data_frame.columns[0], cell_growth_rates.columns[0]]
gene_expression = gene_expression.join(cell_growth_rates).join(days_data_frame)

if args.covariate is not None:
    covariate_df = pd.read_table(args.covariate, index_col=0,
                                 header=None, names=['covariate'],
                                 quoting=csv.QUOTE_NONE, engine='python',
                                 sep=None)
    import itertools

    unique_covariates = pd.unique(covariate_df[covariate_df.columns[0]].values)
    covariate_pairs = list(itertools.product(unique_covariates, unique_covariates))
    args.resample_iter = 1
    gene_expression = gene_expression.join(covariate_df)
    fields_to_drop_for_distance.append(covariate_df.columns[0])
else:
    covariate_pairs = [['', '']]

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
    subsample_writer.write(
        't1' + '\t' + 't2' + '\t' + 't_interpolate' + '\t' + 'distance' + '\t' + 'pairs' + '\t' + 'epsilon' + '\t'
        + 'beta_max' + '\t' + 'delta_max' + '\t' + 'lambda' + '\n')
column_cell_ids_by_time = []
all_cell_ids = set()

for day_index in range(day_pairs.shape[0]):
    t1 = day_pairs.iloc[day_index, 0]
    t2 = day_pairs.iloc[day_index, 1]
    if group_by_day.groups.get(t1) is None or group_by_day.groups.get(
            t2) is None:
        print('skipping transport map from ' + str(t1) + ' to ' + str(t2))
        continue
    p0 = group_by_day.get_group(t1)
    p1 = group_by_day.get_group(t2)
    delta_t = t2 - t1
    if args.verbose:
        print(
            'Computing transport map from ' + str(
                t1) + ' to ' + str(
                t2) + '...', end='')

    cost_matrix = sklearn.metrics.pairwise.pairwise_distances(p0.drop(fields_to_drop_for_distance, axis=1).values,
                                                              p1.drop(fields_to_drop_for_distance, axis=1).values,
                                                              metric='sqeuclidean')
    cost_matrix = cost_matrix / np.median(cost_matrix)
    result = wot.ot.optimal_transport(cost_matrix=cost_matrix,
                                      growth_rate=p0[
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
    if solver is 'floating_epsilon':
        params_writer.write(
            str(t1) + '\t' + str(t2) + '\t' + str(result['epsilon']) + '\t' + str(
                result['lambda1']) + '\t' + str(
                result['lambda2']) + '\n')
    transport_map = pd.DataFrame(result['transport'], index=p0.index,
                                 columns=p1.index)
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

    if resample:
        def write_point_cloud_distance(point_cloud1, point_cloud2, weights1, weights2, point_cloud1_name,
                                       point_cloud2_name):
            if args.verbose:
                print('Computing distance between ' + point_cloud1_name + ' and ' + point_cloud2_name + '...', end='')

            subsample_writer.write(
                str(t1)
                + '\t' + str(t2)
                + '\t' + str(args.t_interpolate)
                + '\t' + str(point_cloud_distance(point_cloud1, point_cloud2, weights1, weights2))
                + '\t' + point_cloud1_name + ' vs ' + point_cloud2_name
                + '\t' + str(args.epsilon)
                + '\t' + str(args.beta_max)
                + '\t' + str(args.delta_max)
                + '\t' + str(args.lambda1)
                + '\n')
            subsample_writer.flush()
            if args.verbose:
                print('done')


        inferred_time = t1 + (t2 - t1) * args.t_interpolate
        p_0_5 = group_by_day.get_group(inferred_time)

        p0_p1_map = result['transport']
        m1 = p0
        m2 = p1

        m1_mtx = m1.drop(fields_to_drop_for_distance, axis=1).values
        m2_mtx = m2.drop(fields_to_drop_for_distance, axis=1).values

        for resample_index in range(args.resample_iter):
            m1_subset, m2_subset, m1_m2_subset_weights = sample_from_transport_map(m1_mtx, m2_mtx, p0_p1_map)
            inferred = m1_subset + args.t_interpolate * (m2_subset - m1_subset)

            m1_random_subset, m2_random_subset = sample_randomly(m1_mtx, m2_mtx, p0_p1_map, m1[
                cell_growth_rates.columns[0]].values ** delta_t)
            random_inferred = m1_random_subset + args.t_interpolate * (m2_random_subset - m1_random_subset)

            m1_uniform_random_subset, m2_uniform_random_subset = sample_uniformly(m1_mtx, m2_mtx, p0_p1_map)
            uniform_random_inferred = m1_uniform_random_subset + args.t_interpolate * (
                    m2_uniform_random_subset - m1_uniform_random_subset)

            point_clouds = [{'m': p_0_5.drop(fields_to_drop_for_distance, axis=1).values, 'weights': None,
                             'name': 'P' + str(args.t_interpolate)},
                            {'m': inferred, 'weights': m1_m2_subset_weights, 'name': 'I' + str(args.t_interpolate)},
                            {'m': random_inferred, 'weights': None, 'name': 'R' + str(args.t_interpolate)},
                            {'m': uniform_random_inferred, 'weights': None, 'name': 'RU' + str(args.t_interpolate)}]

            # self distances
            # D(P0.5, P0.5) split according to condition instead of in half
            # D(R0.5, R0.5) split randomly in two or all covariate splits
            # D(I0.5, I0.5) = split randomly in two or all covariate splits

            # pairs I0.5, R0.5 and P0.5 and maybe(P0 and P1)
            if not args.quick:
                point_clouds.append({'m': p0, 'weights': None, 'name': 'P0'});
                point_clouds.append({'m': p1, 'weights': None, 'name': 'P1'});
            for i in range(1, len(point_clouds)):
                pc1 = point_clouds[i]
                for j in range(i):
                    pc2 = point_clouds[j]
                    write_point_cloud_distance(pc1['m'], pc2['m'],
                                               pc1['weights'],
                                               pc1['weights'],
                                               pc1['name'], pc2['name'])

            # P0.5 to P0.5 self point cloud distance
            if not args.covariate:
                p_0_5_mtx = p_0_5.drop(fields_to_drop_for_distance, axis=1).values
                split1, split2 = split_in_two(p_0_5_mtx.shape[0])
                write_point_cloud_distance(p_0_5_mtx[split1], p_0_5_mtx[split2], None, None,
                                           'P' + str(args.t_interpolate),
                                           'P' + str(args.t_interpolate))
            else:
                seen = {}  # keep track of [c1,c2] and [c2, c1] as they are redundant
                for covariate_pair in covariate_pairs:
                    key = list(covariate_pair)
                    key.sort()
                    key = tuple(key)
                    if covariate_pair[0] != covariate_pair[1] and seen.get(key) is None:
                        seen[key] = True
                        subset1 = p_0_5[p_0_5[covariate_df.columns[0]] == covariate_pair[0]]
                        subset2 = p_0_5[p_0_5[covariate_df.columns[0]] == covariate_pair[1]]
                        write_point_cloud_distance(subset1.drop(fields_to_drop_for_distance, axis=1).values,
                                                   subset2.drop(fields_to_drop_for_distance, axis=1).values, None, None,
                                                   'P' + str(args.t_interpolate) + '_' + str(covariate_pair[0]),
                                                   'P' + str(args.t_interpolate) + '_' + str(covariate_pair[1]))

            # I' to I' self distance and R to R self distance
            self_distance_iters = args.subsample_iter if not args.covariate else len(covariate_pairs)
            interpolated_point_clouds = []
            random_point_clouds = []
            for self_distance_iter in range(self_distance_iters):
                if args.verbose:
                    print('Computing sampled transport map...', end='')
                if args.covariate:
                    covariate_pair = covariate_pairs[self_distance_iter]  # pairs of [[C1, C2],[C2, C2], ...]
                    m1_indices_list = [np.where(m1[covariate_df.columns[0]] == covariate_pair[0])[0]]
                    m2_indices_list = [np.where(m2[covariate_df.columns[0]] == covariate_pair[1])[0]]
                else:
                    m1_indices_list = split_in_two(m1.shape[0])
                    m2_indices_list = split_in_two(m2.shape[0])

                for inner_index in range(len(m1_indices_list)):
                    m1_sample = m1.iloc[m1_indices_list[inner_index]]
                    m2_sample = m2.iloc[m2_indices_list[inner_index]]
                    m1_mtx_sample = m1_sample.drop(fields_to_drop_for_distance, axis=1).values
                    m2_mtx_sample = m2_sample.drop(fields_to_drop_for_distance, axis=1).values
                    c = sklearn.metrics.pairwise.pairwise_distances(m1_mtx_sample, Y=m2_mtx_sample,
                                                                    metric='sqeuclidean')
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
                    if args.verbose:
                        print('done')
                    perturbed_transport = perturbed_result['transport']
                    name_suffix = ('_' + str(covariate_pair[0]) + '_' + str(
                        covariate_pair[1])) if args.covariate is not None else ''

                    coupling_sample = sample_from_transport_map(m1_mtx_sample,
                                                                m2_mtx_sample,
                                                                perturbed_transport)
                    interpolated_point_clouds.append({
                        'm': coupling_sample[0] + args.t_interpolate * (coupling_sample[1] - coupling_sample[0]),
                        'weights': coupling_sample[2], 'name': 'I\'' + name_suffix})
                    random_coupling_sample = sample_randomly(m1_mtx_sample,
                                                             m2_mtx_sample,
                                                             perturbed_transport,
                                                             m1_sample[
                                                                 cell_growth_rates.columns[
                                                                     0]].values ** delta_t)

                    random_point_clouds.append({
                        'm': random_coupling_sample[0] + args.t_interpolate * (
                                random_coupling_sample[1] - random_coupling_sample[0]),
                        'weights': None, 'name': 'R' + str(args.t_interpolate) + name_suffix})

            # compare all pairs of Is and Rs

            for i in range(1, len(interpolated_point_clouds)):
                for j in range(i):
                    write_point_cloud_distance(interpolated_point_clouds[i]['m'], interpolated_point_clouds[j]['m'],
                                               interpolated_point_clouds[i]['weights'],
                                               interpolated_point_clouds[j]['weights'],
                                               interpolated_point_clouds[i]['name'],
                                               interpolated_point_clouds[j]['name'])
                    write_point_cloud_distance(random_point_clouds[i]['m'], random_point_clouds[j]['m'],
                                               random_point_clouds[i]['weights'],
                                               random_point_clouds[j]['weights'],
                                               random_point_clouds[i]['name'],
                                               random_point_clouds[j]['name'])
                # I' vs. P0.5
                # for i in range(len(interpolated_point_clouds)):
                #     write_point_cloud_distance(interpolated_point_clouds[i]['m'], p_0_5,
                #                                interpolated_point_clouds[i]['weights'], None,
                #                                'P' + str(args.t_interpolate), i_name)

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
