#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import logging

import pandas as pd

import wot
import wot.ot

logger = logging.getLogger('wot')


def main(argv):
    parser = argparse.ArgumentParser(
        'Computes distance between trajectories across time')
    parser.add_argument('--out', help='Prefix for output file names', default='wot-trajectory')
    parser.add_argument('--matrix', help=wot.commands.MATRIX_HELP, required=True)
    parser.add_argument('--distance_metric', help='Distance metric (earth mover\'s distance or total variation)',
                        choices=['emd', 'total_variation'], default='emd')
    parser.add_argument('--trajectory', help='One or more trajectory datasets as produced by the trajectory tool',
                        action='append')
    parser.add_argument('--compare',
                        help='If "match", compare trajectories with the same name. ' + 'If "all", compare all pairs. '
                             + 'If "within" compare within a trajectory. If a trajectory name, compare to the specified trajectory',
                        default='within')
    parser.add_argument('--local_pca', type=int, default=30,
                        help='Convert day matrix to local PCA coordinates.'
                             'Set to 0 to disable')
    parser.add_argument('--cell_days', help=wot.commands.CELL_DAYS_HELP, required=True)
    parser.add_argument('--plot', help='Plot results', action='store_true')
    parser.add_argument('--cell_filter', help='File with one cell id per line to include')
    parser.add_argument('--gene_filter',
                        help='File with one gene id per line to use for computing'
                             'cost matrices (e.g. variable genes)')
    parser.add_argument('--cell_day_filter',
                        help='Comma separated list of days to include (e.g. 12,14,16)', type=str)
    # parser.add_argument('--covariate',
    #                     help='Covariate (batch) values for each cell. Used to compute batch to batch distance within a timepoint.')
    parser.add_argument('--cell_days_field', help='Field name in cell_days file that contains cell days',
                        default='day')

    # parser.add_argument('--covariate_field',
    #                     help='Field name in covariate file that contains covariate',
    #                     default='covariate')
    parser.add_argument('--verbose', help='Print progress', action='store_true')
    args = parser.parse_args(argv)
    if args.verbose:
        logger.setLevel(logging.DEBUG)
        logger.addHandler(logging.StreamHandler())
    day_filter = args.cell_day_filter

    compare = args.compare
    cell_days_field = args.cell_days_field
    local_pca = args.local_pca
    days_df = pd.read_csv(args.cell_days, sep=None, index_col='id', engine='python')
    adata = wot.io.read_dataset(args.matrix, obs_filter=args.cell_filter,
                                var_filter=args.gene_filter)
    days = None
    if day_filter is not None:
        days = [float(day) for day in day_filter.split(',')] if type(day_filter) == str else day_filter

    trajectory_files = args.trajectory
    # batch_field_name = args.covariate_field

    distance_metric = args.distance_metric

    trajectory_datasets = []
    for f in trajectory_files:
        trajectory_ds = wot.io.read_dataset(f)
        if len(trajectory_files) > 1:
            trajectory_ds.var.index = trajectory_ds.var.index + '/' + wot.io.get_filename_and_extension(f)[0]
        trajectory_ds.obs = trajectory_ds.obs.join(days_df)
        if days is not None:
            trajectory_ds = trajectory_ds[trajectory_ds.obs[cell_days_field].isin(days)]
        trajectory_datasets.append(trajectory_ds)

    df = wot.tmap.trajectory_divergence(adata, trajectory_datasets, cell_days_field=cell_days_field,
                                        local_pca=local_pca,
                                        compare=compare)
    df.to_csv(args.out + '.csv', index=False)
    if args.plot:
        import matplotlib.pyplot as plt
        wot.tmap.plot_trajectory_divergence(df)
        plt.savefig(args.out + '.png')
