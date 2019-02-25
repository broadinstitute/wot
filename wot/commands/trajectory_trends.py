#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import numpy as np
from matplotlib import pyplot

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
    args = parser.parse_args(argv)

    trajectory_ds = wot.io.read_dataset(args.trajectory)
    wot.io.add_row_metadata_to_dataset(dataset=trajectory_ds, days_path=args.cell_days)
    matrix = wot.io.read_dataset(args.matrix)

    results = wot.tmap.compute_trajectory_trends_from_trajectory(trajectory_ds, matrix)
    # output genes on columns, time on rows, one file per trajectory

    genes = matrix.var.index
    for j in range(len(results)):  # each trajectory
        mean, variance = results[j]
        trajectory_name = trajectory_ds.var.index.values[j]
        basename = args.out + '_' + trajectory_name
        wot.io.write_dataset(mean, basename + '.mean')
        wot.io.write_dataset(variance, basename + '.variance')

        if args.plot:
            pyplot.figure(figsize=(5, 5))

            stds = np.sqrt(variance.X)
            timepoints = mean.obs.index.values

            for i in range(mean.shape[1]):  # each gene
                pyplot.plot(timepoints, mean.X[:, i], label=genes[i])
                pyplot.fill_between(timepoints, mean.X[:, i] - stds[:, i],
                                    mean.X[:, i] + stds[:, i], alpha=.5)

            pyplot.xlabel("Day")
            pyplot.ylabel("Gene expression")
            pyplot.legend()
            pyplot.savefig(basename + '.png')
