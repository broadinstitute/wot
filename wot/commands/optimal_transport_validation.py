#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import pandas as pd
import wot.io
import wot.ot
import wot
import os
from itertools import product
import time


def compute_validation_summary(ot_model, interp_pattern=(1, 2), save_interpolated=False):
    """
    Compute the validation summary for the given OTModel

    Parameters
    ----------
    ot_model : wot.OTModel
        The OTModel to validate
    interp_pattern : (int, int), optional, default: (1,2)
        The interpolation pattern : (x, y) will compute transport maps t[i] to t[i+y] and interpolate at t[i+x]
    save_interpolated : bool, optional, default: False
        Wether to save or discard the interpolated and random point clouds

    Returns
    -------
    validation_summary : pandas.DataFrame
        The validation summary
    """
    i_mid, i_last = interp_pattern
    times = ot_model.timepoints
    # Skip a timepoint and validate using the skipped timepoint
    ot_model.day_pairs = {(times[i], times[i + i_last]): {} for i in range(len(times) - i_last)}
    if 'covariate' not in ot_model.matrix.row_meta.columns:
        print('Warning-no covariate specified.')
        wot.add_cell_metadata(ot_model.matrix, 'covariate', 0)
    ot_model.compute_all_transport_maps(with_covariates=True)
    # Now validate
    summary = []
    local_pca = ot_model.ot_config['local_pca']
    tmap_model = wot.model.TransportMapModel.from_directory(os.path.join(ot_model.tmap_dir, ot_model.tmap_prefix), True)

    for t_cur in range(len(times) - i_last):
        start_time = time.time()
        emd_time = 0

        t0, t05, t1 = times[t_cur], times[t_cur + i_mid], times[t_cur + i_last]
        p0 = ot_model.matrix.where(day=t0).split_by('covariate')
        p05 = ot_model.matrix.where(day=t05).split_by('covariate')
        p1 = ot_model.matrix.where(day=t1).split_by('covariate')
        interp_frac = (t05 - t0) / (t1 - t0)

        for cv0, cv1 in product(p0.keys(), p1.keys()):
            tmap = tmap_model.get_transport_map(t0, t1, covariate=(cv0, cv1))
            interp_size = (len(p0[cv0]) + len(p1[cv1])) // 2
            pca = wot.ot.get_pca(local_pca, p0[cv0].x, p1[cv1].x)
            p0_x = wot.ot.pca_transform(pca, p0[cv0].x)
            p1_x = wot.ot.pca_transform(pca, p1[cv1].x)
            i05 = wot.ot.interpolate_with_ot(p0_x, p1_x, tmap.x, interp_frac, interp_size)
            r05 = wot.ot.interpolate_randomly(p0_x, p1_x, interp_frac, interp_size)

            if save_interpolated:
                prefix = os.path.join(ot_model.tmap_dir, ot_model.tmap_prefix)
                prefix += '_{}_{}_cv{}_cv{}'.format(t0, t1, cv0, cv1)
                wot.io.write_dataset(wot.dataset_from_x(i05),
                                     prefix + '_interp.txt')
                wot.io.write_dataset(wot.dataset_from_x(r05),
                                     prefix + '_random.txt')

            for cv05 in p05.keys():
                p05_x = wot.ot.pca_transform(pca, p05[cv05].x)

                def update_summary(pop, t, name):
                    name_05 = 'P_cv{}'.format(cv05)
                    emd_tmp = time.time()
                    dist = wot.ot.earth_mover_distance(pop, p05_x)
                    summary.append([t0, t05, t1, t, t05, cv0, cv1, name, name_05, dist])
                    return time.time() - emd_tmp

                if cv0 == cv1:
                    emd_time += update_summary(p0_x, t0, 'F_cv{}'.format(cv0))
                    emd_time += update_summary(p1_x, t1, 'L_cv{}'.format(cv1))
                if cv0 == cv1 and cv0 < cv05:
                    p05_cv0_x = wot.ot.pca_transform(pca, p05[cv0].x)
                    emd_time += update_summary(p05_cv0_x, t05, 'P_cv{}'.format(cv0))

                emd_time += update_summary(i05, t05, 'I_cv{}_cv{}'.format(cv0, cv1))
                emd_time += update_summary(r05, t05, 'R_cv{}_cv{}'.format(cv0, cv1))

        total_time = time.time() - start_time
        wot.io.verbose("Processed ({}, {}, {}) for {} covariate pairs, in {:.2f}s ({:.2f}% on EMD)" \
                       .format(t0, t05, t1, len(p0) * len(p1), total_time, 100 * emd_time / total_time))

    # Post-process summary to make it a pandas DataFrame with proper column names
    cols = ['interval_start', 'interval_mid', 'interval_end', 't0', 't1', 'cv0', 'cv1', 'pair0', 'pair1', 'distance']
    return pd.DataFrame(summary, columns=cols)


def main(argv):
    parser = argparse.ArgumentParser(description='Compute a validation summary')
    wot.commands.add_model_arguments(parser)
    wot.commands.add_ot_parameters_arguments(parser)
    parser.add_argument('--covariate', help='Covariate values for each cell')
    parser.add_argument('--save_interpolated', type=bool, default=False,
                        help='Save interpolated and random point clouds')
    parser.add_argument('--interp_pattern', default='1,2',
                        help='The interpolation pattern. "x,y" will compute transport from time t[i] to t[i+y] and interpolate at t[i+x]')
    parser.add_argument('--out', default='./tmaps_val',
                        help='Prefix for output file names')
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
                                       cell_growth_rates=args.cell_growth_rates,
                                       gene_filter=args.gene_filter,
                                       cell_filter=args.cell_filter,
                                       sampling_bias=args.sampling_bias,
                                       scaling_iter=args.scaling_iter,
                                       covariate=args.covariate
                                       )
    summary = compute_validation_summary(ot_model,
                                         interp_pattern=(int(x) for x in args.interp_pattern.split(',')),
                                         save_interpolated=args.save_interpolated)

    summary.to_csv(os.path.join(ot_model.tmap_dir, ot_model.tmap_prefix + '_validation_summary.txt'), sep='\t',
                   index=False)
    res = wot.graphics.group_ot_validation_summary(summary,
                                                   os.path.join(ot_model.tmap_dir,
                                                                ot_model.tmap_prefix + '_validation_summary_stats.txt'))
    wot.graphics.plot_ot_validation_summary(res, os.path.join(ot_model.tmap_dir,
                                                              ot_model.tmap_prefix + '_validation_summary.png'))
