#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
from itertools import product

import numpy as np
import pandas as pd
import scipy.sparse

import wot
import wot.io
import wot.ot


def compute_validation_summary(ot_model, interp_pattern=(0.5, 1), save_interpolated=False, interp_size=10000):
    """
    Compute the validation summary for the given OTModel

    Parameters
    ----------
    ot_model : wot.OTModel
        The OTModel to validate
    interp_pattern : (float, float), optional, default: (0.5,1)
        The interpolation pattern : (x, y) will compute transport maps from t to t+y and interpolate at t+x
    save_interpolated : bool, optional, default: False
        Wether to save or discard the interpolated and random point clouds

    Returns
    -------
    validation_summary : pandas.DataFrame
        The validation summary
    """

    times = np.array(ot_model.timepoints)
    day_pairs = {}
    day_pairs_triplets = []
    tol = 0.01

    for i in range(len(times)):
        t0 = times[i]
        t05_value = t0 + interp_pattern[0]
        t1_value = t0 + interp_pattern[1]
        t1_indices = np.where(np.isclose(times, t1_value, tol))[0]
        t05_indices = np.where(np.isclose(times, t05_value, tol))[0]

        if len(t1_indices) is 1 and len(t05_indices) is 1:
            t1 = times[t1_indices[0]]
            t05 = times[t05_indices[0]]
            day_pairs[(t0, t1)] = {}
            day_pairs_triplets.append((t0, t05, t1))

    ot_model.day_pairs = day_pairs
    if 'covariate' not in ot_model.matrix.row_meta.columns:
        print('Warning-no covariate specified.')
        wot.add_cell_metadata(ot_model.matrix, 'covariate', 0)

    ot_model.compute_all_transport_maps(with_covariates=True)
    # Now validate
    summary = []
    summary_columns = ['interval_start', 'interval_mid', 'interval_end', 't0', 't1', 'cv0', 'cv1', 'pair0', 'pair1',
                       'distance']
    local_pca = ot_model.ot_config['local_pca']
    tmap_model = wot.model.TransportMapModel.from_directory(os.path.join(ot_model.tmap_dir, ot_model.tmap_prefix), True)

    for triplet in day_pairs_triplets:
        t0, t05, t1 = triplet
        interp_frac = (t05 - t0) / (t1 - t0)
        p0_ds = ot_model.matrix.where(day=t0)
        p05_ds = ot_model.matrix.where(day=t05)
        p1_ds = ot_model.matrix.where(day=t1)

        if local_pca > 0:
            matrices = list()
            matrices.append(p0_ds.x if not scipy.sparse.isspmatrix(p0_ds.x) else p0_ds.x.toarray())
            matrices.append(p1_ds.x if not scipy.sparse.isspmatrix(p1_ds.x) else p1_ds.x.toarray())
            p0_ds.x, p1_ds.x, pca, mean_shift = wot.ot.compute_pca(p0_ds.x, p1_ds.x, local_pca)
            eigenvals = np.diag(pca.singular_values_)
            U = np.vstack(matrices).T.dot(pca.components_.T).dot(np.diag(1 / pca.singular_values_))
            y = p05_ds.x - mean_shift
            p05_ds.x = np.diag(1 / pca.singular_values_).dot(U.T.dot(y.T)).T
            p0_ds.col_meta = pd.DataFrame(index=pd.RangeIndex(start=0, stop=local_pca, step=1))
            p1_ds.col_meta = pd.DataFrame(index=pd.RangeIndex(start=0, stop=local_pca, step=1))
            p05_ds.col_meta = pd.DataFrame(index=pd.RangeIndex(start=0, stop=local_pca, step=1))

        tmap = tmap_model.get_transport_map(t0, t1)
        i05 = wot.ot.interpolate_with_ot(p0_ds.x, p1_ds.x, tmap.x, interp_frac, interp_size)
        r05_with_growth = wot.ot.interpolate_randomly_with_growth(p0_ds.x, p1_ds.x, interp_frac, interp_size,
                                                                  p0_ds.row_meta['cell_growth_rate'].values ** (
                                                                      interp_frac))
        r05_no_growth = wot.ot.interpolate_randomly(p0_ds.x, p1_ds.x, interp_frac, interp_size)

        def update_full_summary(pop, t, name):
            dist = wot.ot.earth_mover_distance(pop, p05_ds.x, eigenvals if local_pca > 0 else None)
            summary.append([t0, t05, t1, t, t05, 'full', 'full', name, 'P', dist])

        update_full_summary(i05, t05, 'I')
        update_full_summary(r05_with_growth, t05, 'Rg')
        update_full_summary(r05_no_growth, t05, 'R')
        update_full_summary(p0_ds.x, t0, 'F')
        update_full_summary(p1_ds.x, t1, 'L')

        p0 = p0_ds.split_by('covariate')
        p05 = p05_ds.split_by('covariate')
        p1 = p1_ds.split_by('covariate')
        for cv0, cv1 in product(p0.keys(), p1.keys()):
            tmap = tmap_model.get_transport_map(t0, t1, covariate=(cv0, cv1))
            # interp_size = (len(p0[cv0]) + len(p1[cv1])) / 2
            # pca, mean = wot.ot.get_pca(local_pca, p0[cv0].x, p1[cv1].x)
            # p0_x = wot.ot.pca_transform(pca, mean, p0[cv0].x)
            # p1_x = wot.ot.pca_transform(pca, mean, p1[cv1].x)
            p0_x = p0[cv0].x
            p1_x = p1[cv1].x
            i05 = wot.ot.interpolate_with_ot(p0_x, p1_x, tmap.x, interp_frac, interp_size)
            r05_with_growth = wot.ot.interpolate_randomly_with_growth(p0_x, p1_x, interp_frac, interp_size,
                                                                      p0[cv0].row_meta['cell_growth_rate'].values ** (
                                                                          interp_frac))
            r05_no_growth = wot.ot.interpolate_randomly(p0_x, p1_x, interp_frac, interp_size)

            if save_interpolated:
                prefix = os.path.join(ot_model.tmap_dir, ot_model.tmap_prefix)
                prefix += '_{}_{}_cv{}_cv{}'.format(t0, t1, cv0, cv1)
                wot.io.write_dataset(wot.dataset_from_x(i05),
                                     prefix + '_interp.txt')
                # wot.io.write_dataset(wot.dataset_from_x(r05),
                #                      prefix + '_random.txt')

            for cv05 in p05.keys():
                # p05_x = wot.ot.pca_transform(pca, mean, p05[cv05].x)
                p05_x = p05[cv05].x

                def update_summary(pop, t, name):
                    name_05 = 'P_cv{}'.format(cv05)
                    dist = wot.ot.earth_mover_distance(pop, p05_x, eigenvals if local_pca > 0 else None)
                    summary.append([t0, t05, t1, t, t05, cv0, cv1, name, name_05, dist])

                if cv0 == cv1:
                    update_summary(p0_x, t0, 'F_cv{}'.format(cv0))
                    update_summary(p1_x, t1, 'L_cv{}'.format(cv1))
                if cv0 == cv1 and cv0 < cv05:
                    # p05_cv0_x = wot.ot.pca_transform(pca, mean, p05[cv0].x)
                    update_summary(p05[cv0].x, t05, 'P_cv{}'.format(cv0))

                update_summary(i05, t05, 'I_cv{}_cv{}'.format(cv0, cv1))
                update_summary(r05_with_growth, t05, 'Rg_cv{}_cv{}'.format(cv0, cv1))
                update_summary(r05_no_growth, t05, 'R_cv{}_cv{}'.format(cv0, cv1))

    return pd.DataFrame(summary, columns=summary_columns)


def group_ot_validation_summary(df):
    df = df.copy()
    df['time'] = (df['interval_start'] + df['interval_end']) / 2
    df['type'] = df['pair0'].astype(str).str.split('_').str.get(0)
    # full_df = df[df['cv0'] == 'full']
    # full_df.set_index(['time', 'type'], inplace=True)
    # full_df = full_df.rename(columns={'distance': 'mean'})['mean']
    cv_df = df[df['cv0'] != 'full']

    cv_agg = cv_df.groupby(['time', 'type'])['distance'].agg([np.mean, np.std])
    # cv_agg.update(full_df)
    # mean and std from covariates only CVs
    return cv_agg


def main(argv):
    parser = argparse.ArgumentParser(description='Compute a validation summary')
    wot.commands.add_model_arguments(parser)
    wot.commands.add_ot_parameters_arguments(parser)
    parser.add_argument('--covariate', help='Covariate values for each cell')
    parser.add_argument('--save_interpolated', type=bool, default=False,
                        help='Save interpolated and random point clouds')
    parser.add_argument('--interp_pattern', default='0.5,1',
                        help='The interpolation pattern "x,y" will compute transport maps from time t to t+y and interpolate at t+x')
    parser.add_argument('--out', default='./tmaps_val',
                        help='Prefix for output file names')
    parser.add_argument('--interp_size', default=10000, type=int)
    args = parser.parse_args(argv)
    ot_model = wot.initialize_ot_model(args.matrix, args.cell_days,
                                       tmap_out=args.out,
                                       local_pca=args.local_pca,
                                       growth_iters=args.growth_iters,
                                       epsilon=args.epsilon,
                                       lambda1=args.lambda1,
                                       lambda2=args.lambda2,
                                       max_threads=args.max_threads,
                                       epsilon0=args.epsilon0,
                                       tau=args.tau,
                                       day_pairs=args.config,
                                       cell_day_filter=args.cell_day_filter,
                                       cell_growth_rates=args.cell_growth_rates,
                                       gene_filter=args.gene_filter,
                                       cell_filter=args.cell_filter,
                                       sampling_bias=args.sampling_bias,
                                       scaling_iter=args.scaling_iter,
                                       inner_iter_max=args.inner_iter_max,
                                       covariate=args.covariate
                                       )
    summary = compute_validation_summary(ot_model,
                                         interp_pattern=[float(x) for x in args.interp_pattern.split(',')],
                                         save_interpolated=args.save_interpolated, interp_size=args.interp_size)

    summary.to_csv(os.path.join(ot_model.tmap_dir, ot_model.tmap_prefix + '_validation_summary.txt'), sep='\t',
                   index=False)
    summary_stats = group_ot_validation_summary(summary)
    summary_stats.to_csv(os.path.join(ot_model.tmap_dir, ot_model.tmap_prefix + '_validation_summary_stats.txt'),
                         sep="\t", )
    wot.graphics.plot_ot_validation_summary(summary_stats, os.path.join(ot_model.tmap_dir,
                                                                        ot_model.tmap_prefix + '_validation_summary.png'))
