#!/usr/bin/env python
# -*- coding: utf-8 -*-

import wot.ot
import wot.io
import pandas as pd
import numpy as np
import sklearn.metrics.pairwise
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
        return {'pc0': exp1[pairs[0]], 'pc1': exp2[pairs[1]], 'indices0': pairs[0], 'indices1': pairs[1],
                'weights': None}


def split_in_two(n):
    indices = np.random.choice(n, int(n * 0.5))
    indices_c = np.zeros(n, dtype=bool)
    indices_c[indices] = True
    indices_c = np.invert(indices_c)
    return indices, indices_c


def merge_point_clouds(point_cloud0, point_cloud1):
    # pc0_weights = point_cloud0['weights']
    # pc1_weights = point_cloud1['weights']
    name = point_cloud0['name'] + '+' + point_cloud1['name']
    point_cloud0 = point_cloud0['m']
    point_cloud1 = point_cloud1['m']
    # if pc0_weights is None:
    pc0_weights = np.ones((point_cloud0.shape[0]), dtype=np.float64) / point_cloud0.shape[0]
    # if pc1_weights is None:
    pc1_weights = np.ones((point_cloud1.shape[0]), dtype=np.float64) / point_cloud1.shape[0]
    return {'m': np.concatenate((point_cloud0, point_cloud1)),
            'weights': np.concatenate((pc0_weights, pc1_weights)) / 2,
            'name': name}


def write_point_cloud_distance(point_cloud1, point_cloud2, weights1, weights2, point_cloud1_name,
                               point_cloud2_name, t0, t1):
    subsample_writer.write(
        str(t0)
        + '\t' + str(t1)
        + '\t' + str(point_cloud_distance(point_cloud1, point_cloud2, weights1, weights2))
        + '\t' + point_cloud1_name
        + '\t' + point_cloud2_name
        + '\t' + str(args.epsilon)
        + '\t' + str(args.lambda1)
        + '\t' + str(args.power)
        + '\n')
    subsample_writer.flush()


def compute_self_distance(pc):
    if args.verbose:
        print('Computing self distance for ' + pc['name'] + '...', end='')

    for k in range(args.resample_iter):
        split_indices1, split_indices2 = split_in_two(pc['m'].shape[0])
        write_point_cloud_distance(point_cloud1=pc['m'][split_indices1],
                                   point_cloud2=pc['m'][split_indices2],
                                   weights1=None, weights2=None, point_cloud1_name=pc['name'],
                                   point_cloud2_name=pc['name'], t0=pc['t'],
                                   t1=pc['t'])
    if args.verbose:
        print('done')


def compute_distances(pc0, pc1):
    if args.verbose:
        print('Computing distances between ' + pc0['name'] + ' and ' + pc1['name'] + '...', end='')
    merged = merge_point_clouds(pc0, pc1)
    merged_mtx = merged['m']
    # merge, then split in 2

    for k in range(args.resample_iter):
        split_indices1, split_indices2 = split_in_two(merged['m'].shape[0])
        write_point_cloud_distance(point_cloud1=merged_mtx[split_indices1],
                                   point_cloud2=merged_mtx[split_indices2],
                                   weights1=None, weights2=None, point_cloud1_name=merged['name'],
                                   point_cloud2_name=merged['name'], t0=pc0['t'],
                                   t1=pc1['t'])

    # full point cloud
    write_point_cloud_distance(point_cloud1=pc0['m'],
                               point_cloud2=pc1['m'],
                               weights1=None, weights2=None,
                               point_cloud1_name=pc0['name'],
                               point_cloud2_name=pc1['name'],
                               t0=pc0['t'],
                               t1=pc1['t'])
    # take % from each cloud
    for k in range(args.resample_iter):
        indices_i = np.random.choice(pc0['m'].shape[0],
                                     int(pc0['m'].shape[0] * subsample_percent))
        indices_j = np.random.choice(pc1['m'].shape[0],
                                     int(pc1['m'].shape[0] * subsample_percent))
        write_point_cloud_distance(point_cloud1=pc0['m'][indices_i],
                                   point_cloud2=pc1['m'][indices_j],
                                   weights1=None, weights2=None,
                                   point_cloud1_name=pc0['name'] + '_sub',
                                   point_cloud2_name=pc1['name'] + '_sub',
                                   t0=pc0['t'],
                                   t1=pc1['t'])
    if args.verbose:
        print('done')


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


parser = wot.ot.OptimalTransportHelper.create_base_parser('Compute point cloud distances')
# parser.add_argument('--covariate',
#                     help='Two column tab delimited file without header with '
#                          'cell ids and covariate value')
# parser.add_argument('--quick', action='store_true')
parser.add_argument('--resample_iter', help='Number of resample iterations to perform', type=int, default=10)
parser.add_argument('--subsample_percent', help='Percent to subsample from a point cloud', type=float, default=80)
parser.add_argument('--npairs', type=int)
parser.add_argument('--t_interpolate', help='Interpolation fraction between two time points', type=float)
parser.add_argument('--no_i', action='store_true', help='Do not include interpolated point clouds in computations')
parser.add_argument('--no_p', action='store_true', help='Do not include non-interpolated point clouds in computations')
parser.add_argument('--save', action='store_true', help='Save interpolated point clouds')
args = parser.parse_args()
ot_helper = wot.ot.OptimalTransportHelper(args)
covariate_df = None
# if args.covariate is not None:
#     covariate_df = pd.read_table(args.covariate, index_col=0,
#                                  header=None, names=['covariate'],
#                                  quoting=csv.QUOTE_NONE, engine='python',
#                                  sep=None)
#     import itertools
#
#     unique_covariates = pd.unique(covariate_df[covariate_df.columns[0]].values)
#     covariate_pairs = list(itertools.product(unique_covariates, unique_covariates))
#     args.resample_iter = 1
#
# else:
#     covariate_pairs = [['', '']]

ot_helper = wot.ot.OptimalTransportHelper(args, join=None if covariate_df is None else [covariate_df])

subsample_writer = open(args.prefix + '_subsample_summary.txt', 'w')
subsample_writer.write(
    't0'
    + '\t' + 't1'
    + '\t' + 'distance'
    + '\t' + 'pair0'
    + '\t' + 'pair1'
    + '\t' + 'epsilon'
    + '\t' + 'lambda'
    + '\t' + 'power' + '\n')

fields_to_drop_for_distance = ot_helper.fields_to_drop_for_distance
computed = {}

subsample_percent = args.subsample_percent / 100.0
t_interpolate_s = str(args.t_interpolate)
group_by_day = ot_helper.group_by_day


def callback(cb_args):
    t0 = cb_args['t0']
    t1 = cb_args['t1']
    p0 = cb_args['p0']
    p1 = cb_args['p1']
    inferred_time = t0 + (t1 - t0) * args.t_interpolate
    p_0_5 = group_by_day.get_group(inferred_time)
    p_0_5_mtx = p_0_5.drop(fields_to_drop_for_distance, axis=1).values
    p0_mtx = p0.drop(fields_to_drop_for_distance, axis=1).values
    p1_mtx = p1.drop(fields_to_drop_for_distance, axis=1).values
    static_point_clouds = list()
    static_point_clouds.append({'m': p0_mtx, 'weights': None, 'name': 'P0', 't': t0})
    static_point_clouds.append({'m': p_0_5_mtx, 'weights': None, 'name': 'P' + t_interpolate_s, 't': inferred_time})
    static_point_clouds.append({'m': p1_mtx, 'weights': None, 'name': 'P1', 't': t1})
    if not args.no_p:
        for i in range(len(static_point_clouds)):
            name = str(static_point_clouds[i]['t'])
            if computed.get(name) is None:
                computed[name] = True
                compute_self_distance(static_point_clouds[i])
        for i in range(1, len(static_point_clouds)):
            for j in range(i):
                name = str(static_point_clouds[i]['t']) + str(
                    static_point_clouds[j]['t'])  # avoid duplicate computations
                if computed.get(name) is None:
                    computed[name] = True
                    compute_distances(static_point_clouds[j], static_point_clouds[i])

    result = cb_args['result']
    if result is not None:
        p0_p1_map = result['transport']

        tm_sample = sample_from_transport_map(p0_mtx, p1_mtx, p0_p1_map)
        tm_subset0 = tm_sample['pc0']
        tm_subset1 = tm_sample['pc1']
        p0_m1_subset_weights = tm_sample['weights']
        # p0_m1_subset_weights =

        inferred = tm_subset0 + args.t_interpolate * (tm_subset1 - tm_subset0)
        if args.save:
            wot.io.write_dataset(wot.Dataset(inferred, pd.DataFrame(
                index=p0.iloc[tm_sample['indices0']].index + ';' + p1.iloc[tm_sample['indices1']].index),
                                             pd.DataFrame(
                                                 index=pd.RangeIndex(start=0, stop=inferred.shape[1], step=1))),
                                 args.prefix + '_I_' + str(t0) + '_' + str(t1) + '.txt')
        # m1_random_subset, m2_random_subset = sample_randomly(p0_mtx, p1_mtx, p0_p1_map, p0[
        #     cell_growth_rates.columns[0]].values ** delta_t)
        # random_inferred = m1_random_subset + args.t_interpolate * (m2_random_subset - m1_random_subset)

        # m1_uniform_random_subset, m2_uniform_random_subset = sample_uniformly(p0_mtx, p1_mtx, p0_p1_map)
        # uniform_random_inferred = m1_uniform_random_subset + args.t_interpolate * (
        #         m2_uniform_random_subset - m1_uniform_random_subset)

        I = {'m': inferred, 'weights': p0_m1_subset_weights, 'name': 'I' + t_interpolate_s, 't': inferred_time}
        compute_self_distance(I)
        for i in range(len(static_point_clouds)):
            compute_distances(static_point_clouds[i], I)
            # {'m': random_inferred, 'weights': None, 'name': 'R' + t_interpolate_s}
            # {'m': uniform_random_inferred, 'weights': None, 'name': 'RU' + t_interpolate_s}

            # pairs I0.5, R0.5 and P0.5 and maybe(P0, P1, P0+P1, I+P0.5)

            # point_clouds.append(
            #     {'m': np.concatenate((p0_mtx, p1_mtx)),
            #      'weights': np.concatenate((np.ones((p0_mtx.shape[0]), dtype=np.float64) / p0_mtx.shape[
            #          0], np.ones((p1_mtx.shape[0]), dtype=np.float64) / p1_mtx.shape[
            #                                     0])) / 2,
            #      'name': 'P0+P1'})
            # point_clouds.append(
            #     {'m': np.concatenate((inferred, p_0_5_mtx)),
            #      'weights': np.concatenate((np.ones((inferred.shape[0]), dtype=np.float64) / inferred.shape[
            #          0], np.ones((p_0_5_mtx.shape[0]), dtype=np.float64) / p_0_5_mtx.shape[0])) / 2,
            #      'name': 'I' + t_interpolate_s + '+P' + t_interpolate_s})
            #
            # point_clouds.append(
            #     {'m': np.concatenate((p0_mtx, p_0_5_mtx)),
            #      'weights': np.concatenate((np.ones((p0_mtx.shape[0]), dtype=np.float64) / p0_mtx.shape[
            #          0], np.ones((p_0_5_mtx.shape[0]), dtype=np.float64) / p_0_5_mtx.shape[
            #                                     0])) / 2,
            #      'name': 'P0+P' + t_interpolate_s})
            # point_clouds.append(
            #     {'m': np.concatenate((p1_mtx, p_0_5_mtx)),
            #      'weights': np.concatenate((np.ones((p1_mtx.shape[0]), dtype=np.float64) / p1_mtx.shape[
            #          0], np.ones((p_0_5_mtx.shape[0]), dtype=np.float64) / p_0_5_mtx.shape[
            #                                     0])) / 2,
            #      'name': 'P1+P' + t_interpolate_s})

            # p0_p5_p1_weights = np.concatenate((np.ones((p0_mtx.shape[0]), dtype=np.float64) / p0_mtx.shape[0],
            #                                    np.ones((p1_mtx.shape[0]), dtype=np.float64) / p1_mtx.shape[0])) / 2
            # p0_p5_p1_weights = np.concatenate(
            #     (p0_p5_p1_weights, np.ones((p_0_5_mtx.shape[0]), dtype=np.float64) / p_0_5_mtx.shape[0])) / 2
            # point_clouds.append(
            #     {'m': np.concatenate((p0_mtx, p1_mtx, p_0_5_mtx)),
            #      'weights': p0_p5_p1_weights,
            #      'name': 'P0+P' + t_interpolate_s + '+P1'})

            # pairs of point cloud distances
            #
            # for i in range(1, len(point_clouds)):
            #     pc1 = point_clouds[i]
            #     for j in range(i):
            #         pc2 = point_clouds[j]
            #         write_point_cloud_distance(pc1['m'], pc2['m'],
            #                                    pc1['weights'],
            #                                    pc2['weights'],
            #                                    pc1['name'], pc2['name'])
            #
            # # self point cloud distances
            # if False:
            #     if not args.covariate:
            #         if args.quick:
            #             split1, split2 = split_in_two(p_0_5_mtx.shape[0])
            #             write_point_cloud_distance(p_0_5_mtx[split1], p_0_5_mtx[split2], None, None,
            #                                        'P' + t_interpolate_s,
            #                                        'P' + t_interpolate_s)
            #         else:
            #             for i in range(len(point_clouds)):
            #                 pc = point_clouds[i]
            #                 split1, split2 = split_in_two(pc['m'].shape[0])
            #                 write_point_cloud_distance(pc['m'][split1], pc['m'][split2], None, None,
            #                                            pc['name'], pc['name'])
            #
            #     else:
            #         # I, ICiCj (x4), ICi-, I-Cj  (ICi- means Interpolate between just P0Ci and all of P1)
            #         seen = {}  # keep track of [c1,c2] and [c2, c1] as they are redundant
            #         # self P0.5 distance split by covariate
            #         for covariate_pair in covariate_pairs:
            #             key = list(covariate_pair)
            #             key.sort()
            #             key = tuple(key)
            #             if covariate_pair[0] != covariate_pair[1] and seen.get(key) is None:
            #                 seen[key] = True
            #                 subset1 = p_0_5[p_0_5[covariate_df.columns[0]] == covariate_pair[0]]
            #                 subset2 = p_0_5[p_0_5[covariate_df.columns[0]] == covariate_pair[1]]
            #                 write_point_cloud_distance(subset1.drop(fields_to_drop_for_distance, axis=1).values,
            #                                            subset2.drop(fields_to_drop_for_distance, axis=1).values, None,
            #                                            None,
            #                                            'P' + t_interpolate_s + '_' + str(covariate_pair[0]),
            #                                            'P' + t_interpolate_s + '_' + str(covariate_pair[1]))
            #         # P0 split by covariate to P1
            #
            # if False:
            #     # I' to I' self distance and R to R self distance
            #     self_distance_iters = args.subsample_iter if not args.covariate else len(covariate_pairs)
            #     interpolated_point_clouds = []
            #     random_point_clouds = []
            #     for self_distance_iter in range(self_distance_iters):
            #         if args.verbose:
            #             print('Computing sampled transport map...', end='')
            #         if args.covariate:
            #             covariate_pair = covariate_pairs[self_distance_iter]  # pairs of [[C1, C2],[C2, C2], ...]
            #             m1_indices_list = [np.where(p0[covariate_df.columns[0]] == covariate_pair[0])[0]]
            #             m2_indices_list = [np.where(p1[covariate_df.columns[0]] == covariate_pair[1])[0]]
            #         else:
            #             m1_indices_list = split_in_two(p0.shape[0])
            #             m2_indices_list = split_in_two(p1.shape[0])
            #
            #         for inner_index in range(len(m1_indices_list)):
            #             m1_sample = p0.iloc[m1_indices_list[inner_index]]
            #             m2_sample = p1.iloc[m2_indices_list[inner_index]]
            #             m1_mtx_sample = m1_sample.drop(fields_to_drop_for_distance, axis=1).values
            #             m2_mtx_sample = m2_sample.drop(fields_to_drop_for_distance, axis=1).values
            #             c = sklearn.metrics.pairwise.pairwise_distances(m1_mtx_sample, Y=m2_mtx_sample,
            #                                                             metric='sqeuclidean')
            #             c = c / np.median(c)
            #
            #             perturbed_result = wot.ot.optimal_transport(
            #                 cost_matrix=c,
            #                 growth_rate=m1_sample[
            #                     cell_growth_rates.columns[0]].values,
            #                 delta_days=delta_t,
            #                 max_transport_fraction=args.max_transport_fraction,
            #                 min_transport_fraction=args.min_transport_fraction,
            #                 min_growth_fit=args.min_growth_fit,
            #                 l0_max=args.l0_max, lambda1=args.lambda1,
            #                 lambda2=args.lambda2,
            #                 epsilon=args.epsilon,
            #                 scaling_iter=args.scaling_iter,
            #                 epsilon_adjust=args.epsilon_adjust,
            #                 lambda_adjust=args.lambda_adjust, numItermax=args.numItermax,
            #                 epsilon0=args.epsilon0,
            #                 numInnerItermax=args.numInnerItermax, tau=args.tau, stopThr=args.stopThr, solver=solver)
            #             if args.verbose:
            #                 print('done')
            #             perturbed_transport = perturbed_result['transport']
            #             name_suffix = ('_' + str(covariate_pair[0]) + '_' + str(
            #                 covariate_pair[1])) if args.covariate is not None else ''
            #
            #             coupling_sample = sample_from_transport_map(m1_mtx_sample,
            #                                                         m2_mtx_sample,
            #                                                         perturbed_transport)
            #             interpolated_point_clouds.append({
            #                 'm': coupling_sample[0] + args.t_interpolate * (coupling_sample[1] - coupling_sample[0]),
            #                 'weights': coupling_sample[2], 'name': 'I\'' + name_suffix})
            #             random_coupling_sample = sample_randomly(m1_mtx_sample,
            #                                                      m2_mtx_sample,
            #                                                      perturbed_transport,
            #                                                      m1_sample[
            #                                                          cell_growth_rates.columns[
            #                                                              0]].values ** delta_t)
            #
            #             random_point_clouds.append({
            #                 'm': random_coupling_sample[0] + args.t_interpolate * (
            #                         random_coupling_sample[1] - random_coupling_sample[0]),
            #                 'weights': None, 'name': 'R' + t_interpolate_s + name_suffix})
            #         if not args.covariate:
            #             # compare 50-50 splits of Is and Rs
            #             for i in range(1, len(interpolated_point_clouds)):
            #                 for j in range(i):
            #                     write_point_cloud_distance(interpolated_point_clouds[i]['m'],
            #                                                interpolated_point_clouds[j]['m'],
            #                                                interpolated_point_clouds[i]['weights'],
            #                                                interpolated_point_clouds[j]['weights'],
            #                                                interpolated_point_clouds[i]['name'],
            #                                                interpolated_point_clouds[j]['name'])
            #                     write_point_cloud_distance(random_point_clouds[i]['m'], random_point_clouds[j]['m'],
            #                                                random_point_clouds[i]['weights'],
            #                                                random_point_clouds[j]['weights'],
            #                                                random_point_clouds[i]['name'],
            #                                                random_point_clouds[j]['name'])
            #             interpolated_point_clouds = []
            #             random_point_clouds = []
            #     # compare all covariate pairs of Is and Rs
            #     if args.covariate:
            #         for i in range(1, len(interpolated_point_clouds)):
            #             for j in range(i):
            #                 write_point_cloud_distance(interpolated_point_clouds[i]['m'],
            #                                            interpolated_point_clouds[j]['m'],
            #                                            interpolated_point_clouds[i]['weights'],
            #                                            interpolated_point_clouds[j]['weights'],
            #                                            interpolated_point_clouds[i]['name'],
            #                                            interpolated_point_clouds[j]['name'])
            #                 write_point_cloud_distance(random_point_clouds[i]['m'], random_point_clouds[j]['m'],
            #                                            random_point_clouds[i]['weights'],
            #                                            random_point_clouds[j]['weights'],
            #                                            random_point_clouds[i]['name'],
            #                                            random_point_clouds[j]['name'])
            # I' vs. P0.5
            # for i in range(len(interpolated_point_clouds)):
            #     write_point_cloud_distance(interpolated_point_clouds[i]['m'], p_0_5,
            #                                interpolated_point_clouds[i]['weights'], None,
            #                                'P' + t_interpolate_s, i_name)


if not args.no_i:
    ot_helper.compute_transport_maps(callback)
else:
    for day_index in range(ot_helper.day_pairs.shape[0]):
        t0 = ot_helper.day_pairs.iloc[day_index, 0]
        t1 = ot_helper.day_pairs.iloc[day_index, 1]
        if ot_helper.group_by_day.groups.get(t0) is None or ot_helper.group_by_day.groups.get(
                t1) is None:
            print('skipping transport map from ' + str(t0) + ' to ' + str(t1))
            continue
        p0 = ot_helper.group_by_day.get_group(t0)
        p1 = ot_helper.group_by_day.get_group(t1)
        callback({'t0': t0, 't1': t1, 'p0': p0, 'p1': p1, 'result': None})

subsample_writer.close()
