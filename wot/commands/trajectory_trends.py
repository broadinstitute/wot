#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

from matplotlib import pyplot as plt

import wot.io
import wot.tmap


def main(argv):
    parser = argparse.ArgumentParser(description='Generate mean expression profiles for ' \
                                                 'ancestors and descendants of each trajectory')
    parser.add_argument('--matrix', help=wot.commands.MATRIX_HELP, required=True)
    parser.add_argument('--trajectory', help='Trajectory dataset as produced by the trajectory tool', required=True)
    parser.add_argument('--out', help='Prefix for output file names', default='trends')
    parser.add_argument('--plot', help='Generate plots for each trajectory', action='store_true')
    parser.add_argument('--cell_days', help=wot.commands.CELL_DAYS_HELP, required=True)
    parser.add_argument('--format', help=wot.commands.FORMAT_HELP, choices=wot.commands.FORMAT_CHOICES, default='txt')
    parser.add_argument('--gene_filter',
                        help='File with one gene id per line or comma separated string of list of genes to include from the matrix')
    parser.add_argument('--cell_days_field', help='Field name in cell_days file that contains cell days',
                        default='day')

    args = parser.parse_args(argv)

    cell_days_field = args.cell_days_field
    trajectory_ds = wot.io.read_dataset(args.trajectory, obs=args.cell_days)
    matrix = wot.io.read_dataset(args.matrix, var_filter=args.gene_filter)
    results = wot.tmap.trajectory_trends_from_trajectory(trajectory_ds, matrix, day_field=cell_days_field)
    # output genes on columns, time on rows, one file per trajectory

    genes = matrix.var.index
    for j in range(len(results)):  # each trajectory
        mean = results[j]

        mean.obs.index = mean.obs.index.astype('category')
        # variance.obs.index = variance.obs.index.astype('category')
        trajectory_name = trajectory_ds.var.index.values[j]
        basename = args.out + '_' + trajectory_name
        wot.io.write_dataset(mean, basename + '.mean', args.format)
        # wot.io.write_dataset(variance, basename + '.variance', args.format)

        if args.plot:
            plt.figure(figsize=(5, 5))

            # stds = np.sqrt(variance.X)
            timepoints = mean.obs.index.values.astype(float)

            for i in range(mean.shape[1]):  # each gene

                mean_i = mean[:, i].X
                # std_i = stds[:, i] if len(stds.shape) > 1 else stds[i]
                plt.plot(timepoints, mean_i, label=genes[i])
                # pyplot.fill_between(timepoints, mean_i - std_i,
                #                     mean_i + std_i, alpha=.4)

            plt.xlabel("Day")
            plt.ylabel("Gene expression")
            plt.legend(loc='best')
            plt.savefig(basename + '.png')
