#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os

import anndata
import itertools
import numpy as np
import pandas as pd
import scipy.sparse

import wot
import wot.graphics
import wot.io
import wot.ot


def compute_validation_summary(ot_model, day_pairs_triplets, save_interpolated=False, compute_full_distances=False,
                               interp_size=10000):
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
    day_pairs = {}
    for triplet in day_pairs_triplets:
        day_pairs[(triplet[0], triplet[2])] = {}
    ot_model.day_pairs = day_pairs
    if 'covariate' not in ot_model.matrix.obs.columns:
        print('Warning-no covariate specified.')
        wot.add_cell_metadata(ot_model.matrix, 'covariate', 0)

    ot_model.compute_all_transport_maps(with_covariates=True)
    if compute_full_distances:
        ot_model.compute_all_transport_maps()
    # Now validate
    summary = []
    # 'P': ["#e41a1c", "between real batches"],
    # 'I': ["#377eb8", "between interpolated and real"],
    # 'F': ["#4daf4a", "between first and real"],
    # 'L': ["#984ea3", "between last and real"],
    # 'R': ["#ff7f00", "between random (no growth) and real"],
    # 'Rg': ["#ffff33", "between random (with growth) and real"]
    summary_columns = ['interval_start', 'interval_mid', 'interval_end', 't0', 't1', 'cv0', 'cv1', 'pair0',
                       'pair1', 'distance']

    local_pca = ot_model.ot_config['local_pca']
    tmap_model = wot.tmap.TransportMapModel.from_directory(os.path.join(ot_model.tmap_dir, ot_model.tmap_prefix), True)

    for triplet in day_pairs_triplets:
        t0, t05, t1 = triplet
        interp_frac = (t05 - t0) / (t1 - t0)

        p0_ds = ot_model.matrix[ot_model.matrix.obs['day'] == float(t0), :]
        p05_ds = ot_model.matrix[ot_model.matrix.obs['day'] == float(t05), :]
        p1_ds = ot_model.matrix[ot_model.matrix.obs['day'] == float(t1), :]

        if local_pca > 0:
            matrices = list()
            matrices.append(p0_ds.X if not scipy.sparse.isspmatrix(p0_ds.X) else p0_ds.X.toarray())
            matrices.append(p1_ds.X if not scipy.sparse.isspmatrix(p1_ds.X) else p1_ds.X.toarray())
            p0_pca, p1_pca, pca, mean_shift = wot.ot.compute_pca(p0_ds.X, p1_ds.X, local_pca)
            p0_ds = anndata.AnnData(p0_pca, obs=p0_ds.obs,
                                    var=pd.DataFrame(index=pd.RangeIndex(start=0, stop=local_pca, step=1)))
            p1_ds = anndata.AnnData(p1_pca, obs=p1_ds.obs,
                                    var=pd.DataFrame(index=pd.RangeIndex(start=0, stop=local_pca, step=1)))

            eigenvals = np.diag(pca.singular_values_)
            U = np.vstack(matrices).T.dot(pca.components_.T).dot(np.diag(1 / pca.singular_values_))
            y = p05_ds.X - mean_shift

            p05_ds = anndata.AnnData(np.diag(1 / pca.singular_values_).dot(U.T.dot(y.T)).T, obs=p05_ds.obs,
                                     var=pd.DataFrame(index=pd.RangeIndex(start=0, stop=local_pca, step=1)))

        if compute_full_distances:
            tmap = tmap_model.get_transport_map(t0, t1)
            i05 = wot.ot.interpolate_with_ot(p0_ds.X, p1_ds.X, tmap.X, interp_frac, interp_size)
            r05_with_growth = wot.ot.interpolate_randomly_with_growth(p0_ds.X, p1_ds.X, interp_frac, interp_size,
                                                                      p0_ds.obs['cell_growth_rate'].values ** (
                                                                          interp_frac))
            r05_no_growth = wot.ot.interpolate_randomly(p0_ds.X, p1_ds.X, interp_frac, interp_size)

            def update_full_summary(pop, t, name):
                dist = wot.ot.earth_mover_distance(pop, p05_ds.X, eigenvals if local_pca > 0 else None)
                summary.append([t0, t05, t1, t, t05, 'full', 'full', name, 'P', dist])

            update_full_summary(i05, t05, 'I')
            update_full_summary(r05_with_growth, t05, 'Rg')
            update_full_summary(r05_no_growth, t05, 'R')
            update_full_summary(p0_ds.X, t0, 'F')
            update_full_summary(p1_ds.X, t1, 'L')

        p0 = wot.split_anndata(p0_ds, 'covariate')
        p05 = wot.split_anndata(p05_ds, 'covariate')
        p1 = wot.split_anndata(p1_ds, 'covariate')
        for cv05 in p05.keys():
            p05_x = p05[cv05].X
            seen_first = set()
            seen_last = set()

            def distance_to_p05(pop, t, name, cv):
                dist = wot.ot.earth_mover_distance(pop, p05_x, eigenvals if local_pca > 0 else None)

                summary.append([t0, t05, t1, t, t05, cv, cv05, name, 'P_cv{}'.format(cv05), dist])

            # p05_x = wot.ot.pca_transform(pca, mean, p05[cv05].X)

            for cv05_2 in p05.keys():
                if cv05_2 != cv05:
                    distance_to_p05(p05[cv05_2].X, t05, 'P_cv{}'.format(cv05_2), cv05_2)

            for cv0, cv1 in itertools.product(p0.keys(), p1.keys()):
                tmap = tmap_model.get_transport_map(t0, t1, covariate=(cv0, cv1))
                if tmap is None:
                    # no data for combination of day and covariate
                    continue
                # interp_size = (len(p0[cv0]) + len(p1[cv1])) / 2
                # pca, mean = wot.ot.get_pca(local_pca, p0[cv0].X, p1[cv1].X)
                # p0_x = wot.ot.pca_transform(pca, mean, p0[cv0].X)
                # p1_x = wot.ot.pca_transform(pca, mean, p1[cv1].X)
                p0_x = p0[cv0].X
                p1_x = p1[cv1].X
                i05 = wot.ot.interpolate_with_ot(p0_x, p1_x, tmap.X, interp_frac, interp_size)
                r05_with_growth = wot.ot.interpolate_randomly_with_growth(p0_x, p1_x, interp_frac, interp_size,
                                                                          p0[cv0].obs['cell_growth_rate'].values ** (
                                                                              interp_frac))
                r05_no_growth = wot.ot.interpolate_randomly(p0_x, p1_x, interp_frac, interp_size)
                distance_to_p05(i05, t05, 'I_cv{}_cv{}'.format(cv0, cv1), (cv0, cv1))
                distance_to_p05(r05_with_growth, t05, 'Rg_cv{}_cv{}'.format(cv0, cv1), (cv0, cv1))
                distance_to_p05(r05_no_growth, t05, 'R_cv{}_cv{}'.format(cv0, cv1), (cv0, cv1))

                if cv0 == cv05 and cv0 not in seen_first:
                    seen_first.update(cv0)
                    distance_to_p05(p0_x, t0, 'F_cv{}'.format(cv0), cv0)
                if cv1 == cv05 and cv1 not in seen_last:
                    seen_last.update(cv1)
                    distance_to_p05(p1_x, t1, 'L_cv{}'.format(cv1), cv1)

                if save_interpolated:
                    prefix = os.path.join(ot_model.tmap_dir, ot_model.tmap_prefix)
                    prefix += '_{}_{}_cv{}_cv{}'.format(t0, t1, cv0, cv1)
                    wot.io.write_dataset(wot.dataset_from_x(i05),
                                         prefix + '_interp.txt')
                    # wot.io.write_dataset(wot.dataset_from_x(r05),
                    #                      prefix + '_random.txt')

    return pd.DataFrame(summary, columns=summary_columns)


def group_ot_validation_summary(df):
    df = df.copy()

    df['type'] = df['pair0'].astype(str).str.split('_').str.get(0)
    # full_df = df[df['cv0'] == 'full']
    # full_df.set_index(['time', 'type'], inplace=True)
    # full_df = full_df.rename(columns={'distance': 'mean'})['mean']
    cv_df = df[df['cv0'].astype(str) != 'full']
    cv_agg = cv_df.groupby(['interval_mid', 'type'])['distance'].agg([np.mean, np.std])
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

    parser.add_argument('--day_triplets',
                        help='Three column file without a header containing start time, interpolation time, and end time')
    parser.add_argument('--out', default='tmaps_val',
                        help='Prefix for output file names')
    parser.add_argument('--interp_size', default=10000, type=int)
    args = parser.parse_args(argv)

    ot_model = wot.ot.initialize_ot_model(args.matrix, args.cell_days,
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
                                          force=args.force,
                                          ncells=args.ncells,
                                          ncounts=args.ncounts,
                                          covariate=args.covariate,
                                          transpose=args.transpose
                                          )
    day_pairs_triplets = []
    if args.day_triplets is None:
        unique_times = np.array(ot_model.timepoints)
        for i in range(len(unique_times) - 2):
            t0 = unique_times[i]
            t05 = unique_times[i + 1]
            t1 = unique_times[i + 2]
            day_pairs_triplets.append((t0, t05, t1))
    else:
        day_triplets_df = pd.read_table(args.day_triplets)
        unique_times = np.array(ot_model.timepoints)

        for i in range(len(day_triplets_df.shape[0])):
            t0 = unique_times[np.abs(unique_times - day_triplets_df.iloc[i, 0]).argmin()]

            t05 = unique_times[np.abs(unique_times - day_triplets_df.iloc[i, 1]).argmin()]
            t1 = unique_times[np.abs(unique_times - day_triplets_df.iloc[i, 2]).argmin()]

            day_pairs_triplets.append((t0, t05, t1))

    summary = compute_validation_summary(ot_model,
                                         day_pairs_triplets=day_pairs_triplets,
                                         save_interpolated=args.save_interpolated,
                                         interp_size=args.interp_size)

    summary.to_csv(os.path.join(ot_model.tmap_dir, ot_model.tmap_prefix + '_validation_summary.txt'), sep='\t',
                   index=False)
    summary_stats = group_ot_validation_summary(summary)
    summary_stats.to_csv(os.path.join(ot_model.tmap_dir, ot_model.tmap_prefix + '_validation_summary_stats.txt'),
                         sep="\t", )
    wot.graphics.plot_ot_validation_summary(summary_stats, os.path.join(ot_model.tmap_dir,
                                                                        ot_model.tmap_prefix + '_validation_summary.png'))
