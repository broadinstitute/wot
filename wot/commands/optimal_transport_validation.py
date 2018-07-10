#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import wot.io
import wot.ot


def main(argv):
    parser = wot.ot.OptimalTransportHelper.create_base_parser('Compute point cloud distances')
    parser.add_argument('--t_interpolate', help='Interpolation fraction between two time points', type=float,
                        required=True)
    parser.add_argument('--covariate',
                        help='Two column file with headers "id" and "covariate" indicating cell ids and covariate value')
    parser.add_argument('--resample_iter', help='Number of resample iterations to perform', type=int, default=0)
    parser.add_argument('--npairs', type=int, default=10000)
    parser.add_argument('--save_interpolated', action='store_true', help='Save interpolated point clouds')
    parser.add_argument('--save_transport', action='store_true', help='Save transport maps')
    parser.add_argument('--output_param', action='append', help='Parameters to save in output file')
    parser.add_argument('--progress', action='store_true', help='Print a progress bar while computing')

    args = parser.parse_args(argv)
    args_dict = vars(args)
    parameters_to_write = args.output_param
    if parameters_to_write is None:
        parameters_to_write = []
    covariate_df = None
    unique_covariates = None
    covariate_pairs = None
    if args.covariate is not None:
        covariate_df = pd.read_table(args.covariate, index_col='id', engine='python', sep=None)

        import itertools

        unique_covariates = list(pd.unique(covariate_df[covariate_df.columns[0]].values))
        # unique_covariates.append(None)
        covariate_pairs = list(itertools.product(unique_covariates, unique_covariates))
    ot_helper = wot.ot.OptimalTransportHelper(args, covariate_df=covariate_df, covariate_pairs=covariate_pairs)
    t_interpolate_s = str(args.t_interpolate)
    all_pair_names = [['P' + t_interpolate_s, 'R' + t_interpolate_s], ['P' + t_interpolate_s, 'I' + t_interpolate_s],
                      ['I' + t_interpolate_s, 'R' + t_interpolate_s], ['P0', 'P' + t_interpolate_s],
                      ['P1', 'P' + t_interpolate_s]]

    subsample_writer = open(args.out + '_validation_summary.txt', 'w')
    subsample_writer.write(
        'interval_start'
        + '\t' + 'interval_end'
        + '\t' + 't0'
        + '\t' + 't1'
        + '\t' + 'cv0'
        + '\t' + 'cv1'
        + '\t' + 'distance'
        + '\t' + 'pair0'
        + '\t' + 'pair1'
        + '\t' + 'interval_start_n'
        + '\t' + 'interval_middle_n'
        + '\t' + 'interval_end_n')
    for param in parameters_to_write:
        subsample_writer.write('\t' + param)
    subsample_writer.write('\n')

    def transport_map_callback(cb_args):
        t0 = cb_args['t0']
        t1 = cb_args['t1']
        cv0 = cb_args['cv0']
        cv1 = cb_args['cv1']

        inferred_time = t0 + (t1 - t0) * args.t_interpolate
        p0 = cb_args['P0']  # with cv
        p1 = cb_args['P1']  # with cv
        p0_5_full = cb_args['P0.5']  # full

        interval_start_n = p0.x.shape[0]
        interval_middle_n = p0_5_full.x.shape[0]
        interval_end_n = p1.x.shape[0]
        point_clouds = list()

        point_clouds.append(
            {'m': p0.x, 'weights': None, 'name': 'P0', 't': t0, 'suffix': cb_args['P0_suffix']})
        point_clouds.append(
            {'m': p1.x, 'weights': None, 'name': 'P1', 't': t1, 'suffix': cb_args['P1_suffix']})
        point_clouds.append(
            {'m': p0_5_full.x, 'weights': None, 'name': 'P' + t_interpolate_s, 't': inferred_time,
             'suffix': ''})  # p0_5 is full

        transport_result = cb_args['result']
        if transport_result is not None:
            p0_p1_map = transport_result['transport']
            if args.save_transport:
                transport_map = pd.DataFrame(p0_p1_map, index=cb_args['df0'].index, columns=cb_args['df1'].index)
                transport_map.to_csv(args.out + '_' + str(cb_args['t0']) + '_' + str(cb_args['t1']) + '.txt',
                                     index_label='id',
                                     sep='\t',
                                     doublequote=False)

            tm_sample = wot.ot.sample_from_transport_map(cb_args['P0'].x, cb_args['P1'].x, p0_p1_map, args.npairs,
                                                         args.t_interpolate)
            pc0 = tm_sample['pc0']
            pc1 = tm_sample['pc1']
            p0_m1_subset_weights = tm_sample['weights']
            inferred = pc0 + args.t_interpolate * (pc1 - pc0)
            point_clouds.append(
                {'m': inferred, 'weights': p0_m1_subset_weights, 'name': 'I' + t_interpolate_s,
                 't': inferred_time, 'suffix': cb_args['P0_suffix'] + cb_args['P1_suffix']})

            if args.save_interpolated:
                inferred_row_meta = pd.DataFrame(
                    index=cb_args['df0'].iloc[tm_sample['indices0']].index + ';' + cb_args['df1'].iloc[
                        tm_sample['indices1']].index)
                # save inferred matrix
                wot.io.write_dataset(wot.Dataset(inferred, inferred_row_meta, p0.col_meta),
                                     args.out + '_I_' + str(inferred_time) + '.txt')

            # random sample
            random_sample = wot.ot.sample_randomly(cb_args['P0'].x, cb_args['P1'].x, p0_p1_map,
                                                   p0.row_meta[ot_helper.cell_growth_rates.columns[0]].values ** (
                                                           args.t_interpolate - t0), args.npairs)
            pc0 = random_sample['pc0']
            pc1 = random_sample['pc1']
            p0_m1_subset_weights = random_sample['weights']
            inferred = pc0 + args.t_interpolate * (pc1 - pc0)
            point_clouds.append(
                {'m': inferred, 'weights': p0_m1_subset_weights, 'name': 'R' + t_interpolate_s,
                 't': inferred_time, 'suffix': cb_args['P0_suffix'] + cb_args['P1_suffix']})
            if args.save_interpolated:
                inferred_row_meta = pd.DataFrame(
                    index=cb_args['df0'].iloc[random_sample['indices0']].index + ';' + cb_args['df1'].iloc[
                        random_sample['indices1']].index)
                wot.io.write_dataset(wot.Dataset(inferred, inferred_row_meta, p0.col_meta),
                                     args.out + '_random_' + str(inferred_time) + '.txt')

        pair_names = all_pair_names.copy()

        def get_cloud(cloud_name):
            for i in range(len(point_clouds)):
                if point_clouds[i]['name'] == cloud_name:
                    return point_clouds[i]
            raise ValueError(cloud_name + ' not found')

        if covariate_df is not None:
            batch_names = []
            p0_5_full_copy = wot.Dataset(p0_5_full.x, p0_5_full.row_meta.copy(False), p0_5_full.col_meta)
            # add covariates for P0.5
            for i in range(len(unique_covariates)):
                cv = unique_covariates[i]
                p0_5_cv_filter = np.where(p0_5_full_copy.row_meta[covariate_df.columns[0]] == cv)[0]
                p0_5_cv = wot.Dataset(p0_5_full_copy.x[p0_5_cv_filter], p0_5_full.row_meta.iloc[p0_5_cv_filter],
                                      p0_5_full_copy.col_meta)
                if p0_5_cv.x.shape[0] > 0:
                    name = 'P' + t_interpolate_s + '_cv-' + str(cv)
                    point_clouds.append(
                        {'m': p0_5_cv.x, 'weights': None,
                         'name': name, 't': inferred_time, 'suffix': ''})
                    batch_names.append(name)
            pair_names_i = []
            pair_names_r = []
            for i in range(len(batch_names)):
                pair_names_i.append(['I' + t_interpolate_s, batch_names[i]])
                pair_names_r.append(['R' + t_interpolate_s, batch_names[i]])
            for i in range(1, len(batch_names)):  # batch splits
                for j in range(i):
                    pair_names_i.append([batch_names[i], batch_names[j]])
                    pair_names_r.append([batch_names[i], batch_names[j]])
            pair_names = pair_names + pair_names_i + pair_names_r

        for pair in pair_names:
            if args.progress:
                wot.io.output_progress(cb_args['call_count'] / cb_args['total_call_count'])
            point_cloud_s(get_cloud(pair[0]), get_cloud(pair[1]), t0, t1,
                          interval_start_n=interval_start_n,
                          interval_middle_n=interval_middle_n, interval_end_n=interval_end_n,
                          resample_iter=args.resample_iter, cv0=cv0, cv1=cv1)

    def point_cloud_s(p, q, interval_start, interval_end, interval_start_n,
                      interval_middle_n, interval_end_n, resample_iter, cv0, cv1):
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
                                   point_cloud1_name=p['name'] + p['suffix'],
                                   point_cloud2_name=q['name'] + q['suffix'],
                                   t0=p['t'],
                                   t1=q['t'], interval_start=interval_start, interval_end=interval_end,
                                   interval_start_n=interval_start_n, interval_middle_n=interval_middle_n,
                                   interval_end_n=interval_end_n, cv0=cv0, cv1=cv1)

        def dist(p1, p2):

            weights1 = p1['p']['weights']
            indices1 = p1['indices']
            if weights1 is not None:
                weights1 = weights1[indices1]

            weights2 = p2['p']['weights']
            indices2 = p2['indices']
            if weights2 is not None:
                weights2 = weights2[indices2]
            write_point_cloud_distance(point_cloud1=p1['p']['m'][indices1],
                                       point_cloud2=p2['p']['m'][indices2],
                                       weights1=weights1, weights2=weights2,
                                       point_cloud1_name=p1['p']['name'] + p1['suffix'],
                                       point_cloud2_name=p2['p']['name'] + p2['suffix'],
                                       t0=p1['p']['t'], t1=p2['p']['t'], interval_start=interval_start,
                                       interval_end=interval_end,
                                       interval_start_n=interval_start_n, interval_middle_n=interval_middle_n,
                                       interval_end_n=interval_end_n, cv0=cv0, cv1=cv1)

        for k in range(resample_iter):
            p_indices_a, p_indices_b = wot.ot.split_in_two(p['m'].shape[0])
            q_indices_a, q_indices_b = wot.ot.split_in_two(q['m'].shape[0])

            pairs = [{'p': p, 'indices': p_indices_a, 'suffix': '_split-A' + p['suffix']},
                     {'p': p, 'indices': p_indices_b, 'suffix': '_split-B' + p['suffix']},
                     {'p': q, 'indices': q_indices_a, 'suffix': '_split-A' + q['suffix']},
                     {'p': q, 'indices': q_indices_b, 'suffix': '_split-B' + q['suffix']}]

            dist(pairs[0], pairs[1])
            dist(pairs[0], pairs[2])
            dist(pairs[2], pairs[3])
            dist(pairs[1], pairs[3])

    def write_point_cloud_distance(point_cloud1, point_cloud2, weights1, weights2, point_cloud1_name,
                                   point_cloud2_name, t0, t1, interval_start, interval_end, interval_start_n,
                                   interval_middle_n, interval_end_n, cv0, cv1):
        subsample_writer.write(
            str(interval_start)
            + '\t' + str(interval_end)
            + '\t' + str(t0)
            + '\t' + str(t1)
            + '\t' + str(cv0)
            + '\t' + str(cv1)
            + '\t' + str(
                wot.ot.point_cloud_distance(point_cloud1, point_cloud2, weights1, weights2, ot_helper.eigenvals))
            + '\t' + point_cloud1_name
            + '\t' + point_cloud2_name
            + '\t' + str(interval_start_n)
            + '\t' + str(interval_middle_n)
            + '\t' + str(interval_end_n))
        for param in parameters_to_write:
            subsample_writer.write('\t' + str(args_dict[param]))
        subsample_writer.write('\n')
        subsample_writer.flush()

    if args.progress:
        wot.io.output_progress(0)
    ot_helper.compute_transport_maps(transport_map_callback)
    subsample_writer.close()
