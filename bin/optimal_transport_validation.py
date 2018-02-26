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
    probabilties = tm / np.power(tm.sum(axis=0), 1.0 - args.t_interpolate)
    probabilties = probabilties.flatten()
    if args.npairs is None or args.npairs <= 0:
        # l = l / l.sum()
        reshaped_tmap = probabilties.reshape(exp1.shape[0], exp2.shape[0])
        q = np.where(reshaped_tmap > args.tm_thresh)
        qq = np.where(probabilties > args.tm_thresh)
        return exp1[q[0]], exp2[q[1]], probabilties[qq]
    else:
        probabilties = probabilties / probabilties.sum()
        pairs = np.random.multinomial(args.npairs, probabilties, size=1)
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
                               point_cloud2_name, t0, t1, interval_start, interval_end):
    subsample_writer.write(
        str(interval_start)
        + '\t' + str(interval_end)
        + '\t' + str(t0)
        + '\t' + str(t1)
        + '\t' + str(point_cloud_distance(point_cloud1, point_cloud2, weights1, weights2))
        + '\t' + point_cloud1_name
        + '\t' + point_cloud2_name
        + '\t' + str(args.epsilon)
        + '\t' + str(args.lambda1)
        + '\t' + str(args.power)
        + '\t' + str(args.beta_min)
        + '\n')
    subsample_writer.flush()


def compute_self_distance(pc, interval_start, interval_end):
    if args.verbose:
        print('Computing self distance for ' + pc['name'] + ' ' + str(pc['t']) + '...', end='')

    for k in range(args.resample_iter):
        split_indices1, split_indices2 = split_in_two(pc['m'].shape[0])
        write_point_cloud_distance(point_cloud1=pc['m'][split_indices1],
                                   point_cloud2=pc['m'][split_indices2],
                                   weights1=None, weights2=None, point_cloud1_name=pc['name'],
                                   point_cloud2_name=pc['name'], t0=pc['t'],
                                   t1=pc['t'], interval_start=interval_start, interval_end=interval_end)
    if args.verbose:
        print('done')


def compute_covariate_self_distance(pc, interval_start, interval_end):
    if args.verbose:
        print('Computing self distance for ' + pc['name'] + ' ' + str(pc['t']) + '...', end='')
    # Self-dist(A, "batch-mode")  = Distance(A1,A2).
    df = pc['df']

    for i in range(1, len(ot_helper.unique_covariates)):
        cv0 = ot_helper.unique_covariates[i]
        p0 = df if cv0 is None else df[df[ot_helper.covariate_df.columns[0]] == cv0]
        for j in range(i):
            cv1 = ot_helper.unique_covariates[j]
            # order doesn't matter for self distances
            p1 = df if cv1 is None else df[df[ot_helper.covariate_df.columns[0]] == cv1]
            write_point_cloud_distance(point_cloud1=p0.drop(fields_to_drop_for_distance, axis=1).values,
                                       point_cloud2=p1.drop(fields_to_drop_for_distance, axis=1).values,
                                       weights1=None, weights2=None, point_cloud1_name=pc['name'] + '_' + str(cv0),
                                       point_cloud2_name=pc['name'] + '_' + str(cv1), t0=pc['t'],
                                       t1=pc['t'], interval_start=interval_start, interval_end=interval_end)
    if args.verbose:
        print('done')


def compute_covariate_distances(pc0, pc1, interval_start, interval_end):
    df0 = pc0['df']
    df1 = pc1['df']

    # In batch mode, when we compute Distance(A,B), we don't do random 80% subsamples.
    # D(A1,B1), D(A2,B1), D(A1,B2), D(A2,B2)

    # And in batch mode when we compute self-distance for A+B, we don't do random splits of A+B.
    # Instead, we do:
    # D(A1+B1,A2+B2),  D(A1+B2, A2+B1)

    if args.verbose:
        print('Computing distances between ' + pc0['name'] + ' ' + str(pc0['t']) + ' and ' + pc1['name'] + ' ' + str(
            pc1['t']) + '...', end='')

    # merge all pairs
    for covariate_pair1 in ot_helper.covariate_pairs:
        a_cv0 = covariate_pair1[0]
        a_cv1 = covariate_pair1[1]
        a_p0 = df0 if a_cv0 is None else df0[df0[ot_helper.covariate_df.columns[0]] == a_cv0]
        a_p1 = df1 if a_cv1 is None else df1[df1[ot_helper.covariate_df.columns[0]] == a_cv1]

        A = merge_point_clouds(
            {'m': a_p0.drop(fields_to_drop_for_distance, axis=1).values, 'weights': None,
             'name': '', 't': t0},
            {'m': a_p1.drop(fields_to_drop_for_distance, axis=1).values, 'weights': None,
             'name': '', 't': t0})
        for covariate_pair2 in ot_helper.covariate_pairs:
            b_cv0 = covariate_pair2[0]
            b_cv1 = covariate_pair2[1]
            if b_cv0 == a_cv0 and b_cv1 == a_cv1:
                continue
            p0_2 = df0 if b_cv0 is None else df0[df0[ot_helper.covariate_df.columns[0]] == b_cv0]
            p1_2 = df1 if b_cv1 is None else df1[df1[ot_helper.covariate_df.columns[0]] == b_cv1]
            B = merge_point_clouds(
                {'m': p0_2.drop(fields_to_drop_for_distance, axis=1).values, 'weights': None,
                 'name': '', 't': t0},
                {'m': p1_2.drop(fields_to_drop_for_distance, axis=1).values, 'weights': None,
                 'name': '', 't': t0})
            write_point_cloud_distance(point_cloud1=A['m'],
                                       point_cloud2=B['m'],
                                       weights1=None, weights2=None,
                                       point_cloud1_name=pc0['name'] + '_' + str(a_cv0) + '+' + pc1['name'] + '_' + str(
                                           a_cv1),
                                       point_cloud2_name=pc0['name'] + '_' + str(b_cv0) + '+' + pc1['name'] + '_' + str(
                                           b_cv1), t0=str(pc0['t']) + '_' + str(pc1['t']),
                                       t1=str(pc0['t']) + '_' + str(pc1['t']), interval_start=interval_start,
                                       interval_end=interval_end)

    for covariate_pair in ot_helper.covariate_pairs:
        # do all pairs of distances of pc0 and pc1 including covariate
        cv0 = covariate_pair[0]
        cv1 = covariate_pair[1]
        p0 = df0 if cv0 is None else df0[df0[ot_helper.covariate_df.columns[0]] == cv0]
        p1 = df1 if cv1 is None else df1[df1[ot_helper.covariate_df.columns[0]] == cv1]
        write_point_cloud_distance(point_cloud1=p0.drop(fields_to_drop_for_distance, axis=1).values,
                                   point_cloud2=p1.drop(fields_to_drop_for_distance, axis=1).values,
                                   weights1=None, weights2=None,
                                   point_cloud1_name=pc0['name'] + '_' + str(cv0),
                                   point_cloud2_name=pc1['name'] + '_' + str(cv1),
                                   t0=pc0['t'],
                                   t1=pc1['t'], interval_start=interval_start, interval_end=interval_end)
    if args.verbose:
        print('done')


def compute_distances(pc0, pc1, interval_start, interval_end):
    if args.verbose:
        print('Computing distances between ' + pc0['name'] + ' ' + str(pc0['t']) + ' and ' + pc1['name'] + ' ' + str(
            pc1['t']) + '...', end='')
    merged = merge_point_clouds(pc0, pc1)
    merged_mtx = merged['m']
    # merge, then split in 2
    for k in range(args.resample_iter):
        split_indices1, split_indices2 = split_in_two(merged['m'].shape[0])
        write_point_cloud_distance(point_cloud1=merged_mtx[split_indices1],
                                   point_cloud2=merged_mtx[split_indices2],
                                   weights1=None, weights2=None, point_cloud1_name=merged['name'],
                                   point_cloud2_name=merged['name'], t0=pc0['t'],
                                   t1=pc1['t'], interval_start=interval_start, interval_end=interval_end)

    # full point cloud
    write_point_cloud_distance(point_cloud1=pc0['m'],
                               point_cloud2=pc1['m'],
                               weights1=None, weights2=None,
                               point_cloud1_name=pc0['name'],
                               point_cloud2_name=pc1['name'],
                               t0=pc0['t'],
                               t1=pc1['t'], interval_start=interval_start, interval_end=interval_end)
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
                                   t1=pc1['t'], interval_start=interval_start, interval_end=interval_end)
    if args.verbose:
        print('done')


def point_cloud_distance(c1, c2, a=None, b=None):
    if ot_helper.eigenvals is not None:
        c1 = c1.dot(ot_helper.eigenvals)
        c2 = c2.dot(ot_helper.eigenvals)
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

parser.add_argument('--resample_iter', help='Number of resample iterations to perform', type=int, default=10)
parser.add_argument('--subsample_percent', help='Percent to subsample from a point cloud', type=float, default=80)
parser.add_argument('--npairs', type=int, default=10000)
parser.add_argument('--t_interpolate', help='Interpolation fraction between two time points', type=float)
parser.add_argument('--no_i', action='store_true', help='Do not include interpolated point clouds in computations')
parser.add_argument('--no_p', action='store_true', help='Do not include non-interpolated point clouds in computations')
parser.add_argument('--save', action='store_true', help='Save interpolated point clouds')
parser.add_argument('--quick', action='store_true', help='Compute I0.5 vs P0.5, P0 vs P0.5, and P1 vs P0.5 only')

args = parser.parse_args()
ot_helper = wot.ot.OptimalTransportHelper(args)

ot_helper = wot.ot.OptimalTransportHelper(args)
subsample_writer = None
if not args.no_i or not args.no_p:
    subsample_writer = open(args.prefix + '_subsample_summary.txt', 'w')
    subsample_writer.write(
        'interval_start'
        + '\t' + 'interval_end'
        + '\t' + 't0'
        + '\t' + 't1'
        + '\t' + 'distance'
        + '\t' + 'pair0'
        + '\t' + 'pair1'
        + '\t' + 'epsilon'
        + '\t' + 'lambda'
        + '\t' + 'power'
        + '\t' + 'beta_min' +
        '\n')

fields_to_drop_for_distance = ot_helper.fields_to_drop_for_distance
computed = {}  # avoid duplicate computations

subsample_percent = args.subsample_percent / 100.0
t_interpolate_s = str(args.t_interpolate)
group_by_day = ot_helper.group_by_day


def transport_map_callback(cb_args):
    t0 = cb_args['t0']
    t1 = cb_args['t1']

    inferred_time = t0 + (t1 - t0) * args.t_interpolate
    p0_5 = group_by_day.get_group(inferred_time)
    p0 = group_by_day.get_group(t0)
    p1 = group_by_day.get_group(t1)

    point_clouds_to_compare_with_i = list()
    if not args.quick:
        point_clouds_to_compare_with_i.append(
            {'m': p0.drop(fields_to_drop_for_distance, axis=1).values, 'weights': None, 'name': 'P0', 't': t0})
        point_clouds_to_compare_with_i.append(
            {'m': p1.drop(fields_to_drop_for_distance, axis=1).values, 'weights': None, 'name': 'P1', 't': t1})

    point_clouds_to_compare_with_i.append(
        {'m': p0_5.drop(fields_to_drop_for_distance, axis=1).values, 'weights': None, 'name': 'P' + t_interpolate_s,
         't': inferred_time})

    self_distance_function = compute_self_distance
    compute_distances_function = compute_distances
    point_clouds = point_clouds_to_compare_with_i

    if ot_helper.covariate_df is not None:
        self_distance_function = compute_covariate_self_distance
        compute_distances_function = compute_covariate_distances
        dfs = list()
        if not args.quick:
            dfs.append({'df': p0, 'weights': None, 'name': 'P0', 't': t0})
            dfs.append({'df': p1, 'weights': None, 'name': 'P1', 't': t1})
        dfs.append({'df': p0_5, 'weights': None, 'name': 'P' + t_interpolate_s, 't': inferred_time})

        point_clouds = dfs
        if not args.no_i:
            for cloud in dfs:
                df = cloud['df']
                for i in range(len(ot_helper.unique_covariates)):
                    cv = ot_helper.unique_covariates[i]
                    p = df if cv is None else df[df[ot_helper.covariate_df.columns[0]] == cv]
                    if p.shape[0] > 0:
                        point_clouds_to_compare_with_i.append(
                            {'m': p.drop(fields_to_drop_for_distance, axis=1), 'weights': None,
                             'name': cloud['name'] + '_' + str(cv), 't': cloud['t']})

    if not args.no_p:
        for i in range(len(point_clouds)):
            name = str(point_clouds[i]['t'])
            if computed.get(name) is None:
                # computed[name] = True Do duplicate computes for ease of file parsing
                self_distance_function(point_clouds[i], t0, t1)

        for i in range(1, len(point_clouds)):
            for j in range(i):
                name = str(point_clouds[i]['t']) + '_' + str(
                    point_clouds[j]['t'])
                if computed.get(name) is None:
                    # computed[name] = True
                    compute_distances_function(point_clouds[j], point_clouds[i], t0, t1)

    transport_result = cb_args['result']
    if transport_result is not None:
        p0_p1_map = transport_result['transport']

        tm_sample = sample_from_transport_map(cb_args['p0'], cb_args['p1'], p0_p1_map)
        pc0 = tm_sample['pc0']
        pc1 = tm_sample['pc1']
        p0_m1_subset_weights = tm_sample['weights']

        inferred = pc0 + args.t_interpolate * (pc1 - pc0)
        I = {'m': inferred, 'weights': p0_m1_subset_weights, 'name': 'I' + t_interpolate_s + '_' + cb_args['name'],
             't': inferred_time}
        if not args.no_i:
            compute_self_distance(I, t0, t1)
            if args.covariate is not None:
                for i in range(len(point_clouds_to_compare_with_i)):
                    write_point_cloud_distance(point_cloud1=point_clouds_to_compare_with_i[i]['m'],
                                               point_cloud2=I['m'],
                                               weights1=None, weights2=None,
                                               point_cloud1_name=point_clouds_to_compare_with_i[i]['name'],
                                               point_cloud2_name=I['name'],
                                               t0=point_clouds_to_compare_with_i[i]['t'],
                                               t1=I['t'], interval_start=t0, interval_end=t1)
            else:
                for i in range(len(point_clouds_to_compare_with_i)):
                    compute_distances(point_clouds_to_compare_with_i[i], I, t0, t1)
        if args.save:
            inferred_row_meta = pd.DataFrame(
                index=cb_args['df0'].iloc[tm_sample['indices0']].index + ';' + cb_args['df1'].iloc[
                    tm_sample['indices1']].index)
            # save inferred matrix
            wot.io.write_dataset(wot.Dataset(inferred, inferred_row_meta,
                                             pd.DataFrame(
                                                 index=p0_5.drop(fields_to_drop_for_distance, axis=1).columns)),
                                 args.prefix + '_I_' + str(inferred_time) + '.txt')

            # save actual matrix
            # wot.io.write_dataset(wot.Dataset(p_0_5_mtx, pd.DataFrame(index=p0_5.index),
            #                                  pd.DataFrame(index=p0_5.columns)),
            #                      args.prefix + '_P_' + str(inferred_time) + '.txt')
            #
            # I_P0_5_distance = sklearn.metrics.pairwise.pairwise_distances(inferred, Y=p_0_5_mtx, metric='sqeuclidean')
            # save cost matrix
            # wot.io.write_dataset(wot.Dataset(I_P0_5_distance, inferred_row_meta, pd.DataFrame(index=p0_5.index)),
            #                      args.prefix + '_cost_I' + str(inferred_time) + '_P_' + str(
            #                          inferred_time) + '.txt')
            # coupling = pot.emd(np.ones((I_P0_5_distance.shape[0]), dtype=np.float64) / I_P0_5_distance.shape[0],
            #                    np.ones((I_P0_5_distance.shape[1]), dtype=np.float64) / I_P0_5_distance.shape[1],
            #                    I_P0_5_distance, numItermax=10000000)
            # save coupling
            # wot.io.write_dataset(wot.Dataset(coupling, inferred_row_meta, pd.DataFrame(index=p0_5.index)),
            #                      args.prefix + '_transport_I' + str(inferred_time) + '_P_' + str(
            #                          inferred_time) + '.txt')

            # m1_random_subset, m2_random_subset = sample_randomly(p0_mtx, p1_mtx, p0_p1_map, p0[
            #     cell_growth_rates.columns[0]].values ** delta_t)
            # random_inferred = m1_random_subset + args.t_interpolate * (m2_random_subset - m1_random_subset)

            # m1_uniform_random_subset, m2_uniform_random_subset = sample_uniformly(p0_mtx, p1_mtx, p0_p1_map)
            # uniform_random_inferred = m1_uniform_random_subset + args.t_interpolate * (
            #         m2_uniform_random_subset - m1_uniform_random_subset)
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


if not args.no_i or args.save:
    ot_helper.compute_transport_maps(transport_map_callback)
else:
    for day_index in range(ot_helper.day_pairs.shape[0]):
        t0 = ot_helper.day_pairs.iloc[day_index, 0]
        t1 = ot_helper.day_pairs.iloc[day_index, 1]
        if ot_helper.group_by_day.groups.get(t0) is None or ot_helper.group_by_day.groups.get(
                t1) is None:
            print('skipping transport map from ' + str(t0) + ' to ' + str(t1))
            continue
        p0_full = ot_helper.group_by_day.get_group(t0)
        p1_full = ot_helper.group_by_day.get_group(t1)
        transport_map_callback({'t0': t0, 't1': t1, 'p0': p0_full.drop(fields_to_drop_for_distance, axis=1).values,
                                'p1': p1_full.drop(fields_to_drop_for_distance, axis=1).values,
                                'result': None, 'name': None})

if subsample_writer is not None:
    subsample_writer.close()
