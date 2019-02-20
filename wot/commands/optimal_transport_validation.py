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


def compute_validation_summary(ot_model, day_pairs_triplets=None, save_interpolated=False,
                               interp_size=10000, compute_full_distances=False):
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
    if day_pairs_triplets is None:
        day_pairs_triplets = []
        unique_times = np.array(ot_model.timepoints)
        for i in range(len(unique_times) - 2):
            t0 = unique_times[i]
            t05 = unique_times[i + 1]
            t1 = unique_times[i + 2]
            day_pairs_triplets.append((t0, t05, t1))

    day_pairs = {}
    for triplet in day_pairs_triplets:
        day_pairs[(triplet[0], triplet[2])] = {}
    ot_model.day_pairs = day_pairs
    has_covariate = 'covariate' in ot_model.matrix.obs.columns
    if not has_covariate and not compute_full_distances:
        print('No covariate specified. Please provide a covariate or compute full distances')
        exit(1)

    if has_covariate:
        ot_model.compute_all_transport_maps(with_covariates=True)
    if compute_full_distances:
        ot_model.compute_all_transport_maps()
    # Now validate
    summary_list = []
    # 'P': ["#e41a1c", "between real batches"],
    # 'I': ["#377eb8", "between interpolated and real batch"],
    # 'F': ["#4daf4a", "between first and real "],
    # 'L': ["#984ea3", "between last and real"],
    # 'R': ["#ff7f00", "between random (no growth) and real"],
    # 'Rg': ["#ffff33", "between random (with growth) and real"]

    local_pca = ot_model.ot_config['local_pca']
    tmap_model = wot.tmap.TransportMapModel.from_directory(os.path.join(ot_model.tmap_dir, ot_model.tmap_prefix), True)
    if compute_full_distances:
        tmap_model_full = wot.tmap.TransportMapModel.from_directory(
            os.path.join(ot_model.tmap_dir, ot_model.tmap_prefix))

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
            tmap_full = tmap_model_full.get_transport_map(t0, t1)

            def update_full_summary(pop, t, name, pop2=p05_ds.X):
                dist = wot.ot.earth_mover_distance(pop, pop2, eigenvals if local_pca > 0 else None)
                summary_list.append(
                    {'interval_start': t0,
                     'interval_mid': t05,
                     'interval_end': t1,
                     't0': t,
                     't1': t05,
                     'cv0': '',
                     'cv1': '',
                     'name': name,
                     'distance': dist,
                     'full': True})

            r05_with_growth = wot.ot.interpolate_randomly_with_growth(p0_ds.X, p1_ds.X, interp_frac, interp_size,
                                                                      p0_ds.obs['cell_growth_rate'].values ** (
                                                                          interp_frac))

            r05_no_growth = wot.ot.interpolate_randomly(p0_ds.X, p1_ds.X, interp_frac, interp_size)
            try:
                i05 = wot.ot.interpolate_with_ot(p0_ds.X, p1_ds.X, tmap_full.X, interp_frac,
                                                 interp_size)  # TODO handle downsampling cells case
                update_full_summary(i05, t05, 'I')
            except ValueError:
                pass

            update_full_summary(r05_with_growth, t05, 'Rg')
            update_full_summary(r05_no_growth, t05, 'R')
            update_full_summary(p0_ds.X, t0, 'F')
            update_full_summary(p1_ds.X, t1, 'L')
            update_full_summary(p0_ds.X, t1, 'A', p1_ds.X)

        if not has_covariate:
            continue
        p0 = wot.split_anndata(p0_ds, 'covariate')
        p05 = wot.split_anndata(p05_ds, 'covariate')
        p1 = wot.split_anndata(p1_ds, 'covariate')
        for cv05 in p05.keys():
            p05_x = p05[cv05].X
            seen_first = set()
            seen_last = set()

            def distance_to_p05(pop, t, name, cv):
                dist = wot.ot.earth_mover_distance(pop, p05_x, eigenvals if local_pca > 0 else None)
                summary_list.append(
                    {'interval_start': t0,
                     'interval_mid': t05,
                     'interval_end': t1,
                     't0': t,
                     't1': t05,
                     'cv0': cv,
                     'cv1': cv05,
                     'name': name,
                     'distance': dist,
                     'full': False})

            # p05_x = wot.ot.pca_transform(pca, mean, p05[cv05].X)

            for cv05_2 in p05.keys():  # distance between batches
                if cv05_2 != cv05:
                    distance_to_p05(p05[cv05_2].X, t05, 'P', cv05_2)

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
                distance_to_p05(i05, t05, 'I', (cv0, cv1))
                distance_to_p05(r05_with_growth, t05, 'Rg', (cv0, cv1))
                distance_to_p05(r05_no_growth, t05, 'R', (cv0, cv1))

                if cv0 == cv05 and cv0 not in seen_first:
                    seen_first.add(cv0)
                    distance_to_p05(p0_x, t0, 'F', cv0)
                if cv1 == cv05 and cv1 not in seen_last:
                    seen_last.add(cv1)
                    distance_to_p05(p1_x, t1, 'L', cv1)

                if save_interpolated:
                    prefix = os.path.join(ot_model.tmap_dir, ot_model.tmap_prefix)
                    prefix += '_{}_{}_cv{}_cv{}'.format(t0, t1, cv0, cv1)
                    wot.io.write_dataset(wot.dataset_from_x(i05),
                                         prefix + '_interp.txt')
                    # wot.io.write_dataset(wot.dataset_from_x(r05),
                    #                      prefix + '_random.txt')

    return pd.DataFrame(summary_list)


def main(argv):
    parser = argparse.ArgumentParser(description='Compute a validation summary')
    wot.commands.add_model_arguments(parser)
    wot.commands.add_ot_parameters_arguments(parser)
    parser.add_argument('--covariate', help='Covariate values for each cell')
    parser.add_argument('--save_interpolated', type=bool, default=False,
                        help='Save interpolated and random point clouds')
    parser.add_argument('--full_distances', action='store_true',
                        help='Compute full distances')
    parser.add_argument('--day_triplets',
                        help='Three column file without a header containing start time, interpolation time, and end time')
    parser.add_argument('--out', default='tmaps_val',
                        help='Prefix for output file names')
    parser.add_argument('--interp_size', default=10000, type=int)
    args = parser.parse_args(argv)
    ot_model = wot.commands.initialize_ot_model_from_args(args)
    day_pairs_triplets = None
    if args.day_triplets is not None:
        day_pairs_triplets = []
        if os.path.isfile(args.day_triplets):
            day_triplets_df = pd.read_csv(args.day_triplets, engine='python', sep=None)
        else:
            triplets = args.day_triplets.split(';')
            array_of_arrays = []
            for triplet in triplets:
                tokens = triplet.split(',')
                array_of_arrays.append([float(tokens[0].strip()), float(tokens[1].strip()), float(tokens[2].strip())])
            day_triplets_df = pd.DataFrame(array_of_arrays)
        unique_times = np.array(ot_model.timepoints)  # find closest to actual time

        for i in range(day_triplets_df.shape[0]):
            t0 = unique_times[np.abs(unique_times - day_triplets_df.iloc[i, 0]).argmin()]

            t05 = unique_times[np.abs(unique_times - day_triplets_df.iloc[i, 1]).argmin()]
            t1 = unique_times[np.abs(unique_times - day_triplets_df.iloc[i, 2]).argmin()]

            day_pairs_triplets.append((t0, t05, t1))

    summary = compute_validation_summary(ot_model,
                                         day_pairs_triplets=day_pairs_triplets,
                                         save_interpolated=args.save_interpolated,
                                         interp_size=args.interp_size,
                                         compute_full_distances=args.full_distances)

    summary.to_csv(os.path.join(ot_model.tmap_dir, ot_model.tmap_prefix + '_validation_summary.txt'), sep='\t',
                   index=False)

    if args.covariate is not None:
        summary_stats = summary[summary['full'] == False]
        summary_stats = summary_stats.groupby(['interval_mid', 'name'])['distance'].agg([np.mean, np.std])
        summary_stats.to_csv(os.path.join(ot_model.tmap_dir, ot_model.tmap_prefix + '_cv_validation_summary_stats.txt'),
                             sep="\t", )
        wot.graphics.plot_ot_validation_summary(summary_stats, os.path.join(ot_model.tmap_dir,
                                                                            ot_model.tmap_prefix + '_cv_validation_summary.png'))

    if args.full_distances:
        summary_stats = summary[summary['full']]
        summary_stats = summary_stats.groupby(['interval_mid', 'name'])['distance'].agg([np.mean, np.std])
        summary_stats.to_csv(
            os.path.join(ot_model.tmap_dir, ot_model.tmap_prefix + '_full_validation_summary_stats.txt'),
            sep="\t", )
        wot.graphics.plot_ot_validation_summary(summary_stats, os.path.join(ot_model.tmap_dir,
                                                                            ot_model.tmap_prefix + '_full_validation_summary.png'))
