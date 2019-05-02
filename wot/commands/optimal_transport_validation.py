#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import logging
import os

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

import wot


def main(argv):
    parser = argparse.ArgumentParser(description='Compute a validation summary')
    wot.commands.add_ot_parameters_arguments(parser)
    parser.add_argument('--covariate', help='Covariate values for each cell')
    # parser.add_argument('--save_interpolated', type=bool, default=False,
    #                     help='Save interpolated and random point clouds')
    parser.add_argument('--full_distances', action='store_true',
                        help='Compute full distances')
    parser.add_argument('--day_triplets',
                        help='Three column file without a header containing start time, interpolation time, and end time')
    parser.add_argument('--out', default='tmaps_val',
                        help='Prefix for output file names')
    parser.add_argument('--interp_size', default=10000, type=int)
    parser.add_argument('--covariate_field',
                        help='Field name in covariate file that contains covariate',
                        default='covariate')
    args = parser.parse_args(argv)
    if args.verbose:
        logger = logging.getLogger('wot')
        logger.setLevel(logging.DEBUG)
        logger.addHandler(logging.StreamHandler())
    ot_model = wot.commands.initialize_ot_model_from_args(args)
    tmap_dir, tmap_prefix = os.path.split(args.out) if args.out is not None else (None, None)
    tmap_prefix = tmap_prefix or "tmaps"
    tmap_dir = tmap_dir or '.'

    day_triplets = None
    if args.day_triplets is not None:
        day_triplets = []
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

            day_triplets.append((t0, t05, t1))

    # tmap_dir=tmap_dir, tmap_prefix=tmap_prefix,
    # no_overwrite=args.no_overwrite, output_file_format=args.format,
    summary = wot.ot.compute_validation_summary(ot_model,
                                                day_triplets=day_triplets,
                                                interp_size=args.interp_size,
                                                compute_full_distances=args.full_distances)

    summary.to_csv(os.path.join(tmap_dir, tmap_prefix + '_validation_summary.txt'), sep='\t',
                   index=False)

    summary_stats = summary[summary['full'] == False]
    summary_stats = summary_stats.groupby(['interval_mid', 'name'])['distance'].agg([np.mean, np.std])
    summary_stats.to_csv(os.path.join(tmap_dir, tmap_prefix + '_cv_validation_summary_stats.txt'),
                         sep="\t", )
    wot.graphics.plot_ot_validation_summary_stats(summary_stats)
    plt.savefig(os.path.join(tmap_dir, tmap_prefix + '_cv_validation_summary.png'))

    wot.graphics.plot.plot_ot_validation_ratio(summary_stats, os.path.join(tmap_dir,
                                                                           tmap_prefix + '_cv_validation_summary_ratio.png'))
    if args.full_distances:
        summary_stats = summary[summary['full']]
        summary_stats = summary_stats.groupby(['interval_mid', 'name'])['distance'].agg([np.mean, np.std])
        summary_stats.to_csv(
            os.path.join(tmap_dir, tmap_prefix + '_full_validation_summary_stats.txt'),
            sep="\t", )
        wot.graphics.plot_ot_validation_summary_stats(summary_stats)
        plt.savefig(os.path.join(tmap_dir, tmap_prefix + '_full_validation_summary.png'))
