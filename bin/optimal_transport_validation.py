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
    probabilties = p / p.sum()
    # pairs = np.random.multinomial(n, p, size=1)
    # pairs = np.nonzero(pairs.reshape(exp1.shape[0], exp2.shape[0]))

    s = np.random.multinomial(n, probabilties, size=1)
    reshaped_s = s.reshape(exp1.shape[0], exp2.shape[0])
    pairs = np.nonzero(reshaped_s)
    weights = reshaped_s[pairs]
    return {'pc0': exp1[pairs[0]], 'pc1': exp2[pairs[1]],
            'indices0': pairs[0],
            'indices1': pairs[1],
            'weights': weights}


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
        # pairs = np.random.multinomial(args.npairs, probabilties, size=1)
        # pairs = np.nonzero(pairs.reshape(exp1.shape[0], exp2.shape[0]))
        #  return {'pc0': exp1[pairs[0]], 'pc1': exp2[pairs[1]], 'indices0': pairs[0], 'indices1': pairs[1], 'weights': None}
        s = np.random.multinomial(args.npairs, probabilties, size=1)
        reshaped_s = s.reshape(exp1.shape[0], exp2.shape[0])
        pairs = np.nonzero(reshaped_s)
        weights = reshaped_s[pairs]
        return {'pc0': exp1[pairs[0]], 'pc1': exp2[pairs[1]],
                'indices0': pairs[0],
                'indices1': pairs[1],
                'weights': weights}


def split_in_two(n):
    indices = np.random.choice(n, int(n * 0.5))
    indices_c = np.zeros(n, dtype=bool)
    indices_c[indices] = True
    indices_c = np.invert(indices_c)
    return indices, indices_c


def write_point_cloud_distance(point_cloud1, point_cloud2, weights1, weights2, point_cloud1_name,
                               point_cloud2_name, t0, t1, interval_start, interval_end, interval_start_n,
                               interval_middle_n, interval_end_n):
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
        + '\t' + str(args.delta_min)
        + '\t' + str(args.beta_max)
        + '\t' + str(interval_start_n)
        + '\t' + str(interval_middle_n)
        + '\t' + str(interval_end_n)
        + '\n')
    subsample_writer.flush()


def point_cloud_s(p, q, interval_start, interval_end, interval_start_n,
                  interval_middle_n, interval_end_n):
    # Given a pair of point clouds: P,Q.
    # Compute random splits P_A, P_B  Q_A, Q_B.
    # Compute distances
    # D(P_A,P_B)
    # D(P_A, Q_A)
    # D(Q_A,Q_B)
    # D(P_B, Q_B)

    write_point_cloud_distance(point_cloud1=p['m'],
                               point_cloud2=q['m'],
                               weights1=p['weights'], weights2=q['weights'],
                               point_cloud1_name=p['name'],
                               point_cloud2_name=q['name'],
                               t0=p['t'],
                               t1=q['t'], interval_start=interval_start, interval_end=interval_end,
                               interval_start_n=interval_start_n, interval_middle_n=interval_middle_n,
                               interval_end_n=interval_end_n)

    for k in range(args.resample_iter):
        p_indices_a, p_indices_b = split_in_two(p['m'].shape[0])
        q_indices_a, q_indices_b = split_in_two(q['m'].shape[0])

        pairs = [{'p': p, 'indices': p_indices_a}, {'p': p, 'indices': p_indices_b}, {'p': q, 'indices': q_indices_a},
                 {'p': q, 'indices': q_indices_b}]

        for i in range(1, len(pairs)):
            p1 = pairs[i]
            weights1 = p1['p']['weights']
            indices1 = p1['indices']
            if weights1 is not None:
                weights1 = weights1[indices1]
            for j in range(i):
                p2 = pairs[j]
                weights2 = p2['p']['weights']
                indices2 = p2['indices']
                if weights2 is not None:
                    weights2 = weights2[indices2]
                write_point_cloud_distance(point_cloud1=p1['p']['m'][indices1],
                                           point_cloud2=p2['p']['m'][indices2],
                                           weights1=weights1, weights2=weights2,
                                           point_cloud1_name=p1['p']['name'] + '_sub',
                                           point_cloud2_name=p2['p']['name'] + '_sub',
                                           t0=p1['p']['t'], t1=p2['p']['t'], interval_start=interval_start,
                                           interval_end=interval_end,
                                           interval_start_n=interval_start_n, interval_middle_n=interval_middle_n,
                                           interval_end_n=interval_end_n)


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
parser.add_argument('--npairs', type=int, default=10000)
parser.add_argument('--t_interpolate', help='Interpolation fraction between two time points', type=float)
parser.add_argument('--save', action='store_true', help='Save interpolated point clouds')
parser.add_argument('--save_transport', action='store_true', help='Save transport maps')

args = parser.parse_args()
covariate_df = None
unique_covariates = None
if args.covariate is not None:
    covariate_df = pd.read_table(args.covariate, index_col=0,
                                 header=None, names=['covariate'],
                                 engine='python',
                                 sep=None)

    import itertools

    unique_covariates = list(pd.unique(covariate_df[covariate_df.columns[0]].values))
    # unique_covariates.append(None)
    covariate_pairs = list(itertools.product(unique_covariates, unique_covariates))
ot_helper = wot.ot.OptimalTransportHelper(args)
pair_names = None
t_interpolate_s = str(args.t_interpolate)
pair_names = [['P' + t_interpolate_s, 'R' + t_interpolate_s], ['P' + t_interpolate_s, 'I' + t_interpolate_s],
              ['I' + t_interpolate_s, 'R' + t_interpolate_s]]

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
    + '\t' + 'beta_min'
    + '\t' + 'delta_min'
    + '\t' + 'beta_max'
    + '\t' + 'interval_start_n'
    + '\t' + 'interval_middle_n'
    + '\t' + 'interval_end_n'
    + '\n')

fields_to_drop_for_distance = ot_helper.fields_to_drop_for_distance

group_by_day = ot_helper.group_by_day


def transport_map_callback(cb_args):
    t0 = cb_args['t0']
    t1 = cb_args['t1']

    inferred_time = t0 + (t1 - t0) * args.t_interpolate
    p0 = group_by_day.get_group(t0)
    p0_5 = group_by_day.get_group(inferred_time)
    p1 = group_by_day.get_group(t1)
    interval_start_n = p0.shape[0]
    interval_middle_n = p0_5.shape[0]
    interval_end_n = p1.shape[0]
    point_clouds = list()

    point_clouds.append(
        {'m': p0.drop(fields_to_drop_for_distance, axis=1).values, 'weights': None, 'name': 'P0', 't': t0})
    point_clouds.append(
        {'m': p1.drop(fields_to_drop_for_distance, axis=1).values, 'weights': None, 'name': 'P1', 't': t1})
    point_clouds.append(
        {'m': p0_5.drop(fields_to_drop_for_distance, axis=1).values, 'weights': None, 'name': 'P' + t_interpolate_s,
         't': inferred_time})

    transport_result = cb_args['result']
    if transport_result is not None:
        p0_p1_map = transport_result['transport']
        if args.save_transport:
            transport_map = pd.DataFrame(p0_p1_map, index=cb_args['df0'].index, columns=cb_args['df1'].index)
            transport_map.to_csv(args.prefix + '_' + str(cb_args['t0']) + '_' + str(cb_args['t1']) + '.txt',
                                 index_label='id',
                                 sep='\t',
                                 doublequote=False)

        tm_sample = sample_from_transport_map(cb_args['p0'], cb_args['p1'], p0_p1_map)
        pc0 = tm_sample['pc0']
        pc1 = tm_sample['pc1']
        p0_m1_subset_weights = tm_sample['weights']
        inferred = pc0 + args.t_interpolate * (pc1 - pc0)
        point_clouds.append(
            {'m': inferred, 'weights': p0_m1_subset_weights, 'name': 'I' + t_interpolate_s,
             't': inferred_time})

        if args.save:
            inferred_row_meta = pd.DataFrame(
                index=cb_args['df0'].iloc[tm_sample['indices0']].index + ';' + cb_args['df1'].iloc[
                    tm_sample['indices1']].index)
            # save inferred matrix
            wot.io.write_dataset(wot.Dataset(inferred, inferred_row_meta,
                                             pd.DataFrame(
                                                 index=p0_5.drop(fields_to_drop_for_distance, axis=1).columns)),
                                 args.prefix + '_I_' + str(inferred_time) + '.txt')

        random_sample = sample_randomly(cb_args['p0'], cb_args['p1'], p0_p1_map,
                                        p0[ot_helper.cell_growth_rates.columns[0]].values ** (args.t_interpolate - t0))
        pc0 = random_sample['pc0']
        pc1 = random_sample['pc1']
        p0_m1_subset_weights = random_sample['weights']
        inferred = pc0 + args.t_interpolate * (pc1 - pc0)
        point_clouds.append(
            {'m': inferred, 'weights': p0_m1_subset_weights, 'name': 'R' + t_interpolate_s,
             't': inferred_time})
        if args.save:
            inferred_row_meta = pd.DataFrame(
                index=cb_args['df0'].iloc[random_sample['indices0']].index + ';' + cb_args['df1'].iloc[
                    random_sample['indices1']].index)
            wot.io.write_dataset(wot.Dataset(inferred, inferred_row_meta,
                                             pd.DataFrame(
                                                 index=p0_5.drop(fields_to_drop_for_distance, axis=1).columns)),
                                 args.prefix + '_random_' + str(inferred_time) + '.txt')

    if covariate_df is not None:
        batch_names = []
        p0_5_cv = p0_5.copy(False).join(covariate_df)
        fields_to_drop_for_distance_cv = list(fields_to_drop_for_distance)
        fields_to_drop_for_distance_cv.append(covariate_df.columns[0])
        for i in range(len(unique_covariates)):
            cv = unique_covariates[i]
            p = p0_5_cv[p0_5_cv[covariate_df.columns[0]] == cv]

            if p.shape[0] > 0:
                name = 'P' + t_interpolate_s + '_' + str(cv)
                point_clouds.append(
                    {'m': p.drop(fields_to_drop_for_distance_cv, axis=1).values, 'weights': None,
                     'name': name, 't': inferred_time})
                batch_names.append(name)
        for i in range(1, len(batch_names)):
            for j in range(i):
                pair_names.append([batch_names[i], batch_names[j]])

    def get_cloud(cloud_name):
        for i in range(len(point_clouds)):
            if point_clouds[i]['name'] == cloud_name:
                return point_clouds[i]
        raise ValueError(cloud_name + ' not found')

    for pair in pair_names:
        point_cloud_s(get_cloud(pair[0]), get_cloud(pair[1]), t0, t1,
                      interval_start_n=interval_start_n,
                      interval_middle_n=interval_middle_n, interval_end_n=interval_end_n)


ot_helper.compute_transport_maps(transport_map_callback)
subsample_writer.close()
