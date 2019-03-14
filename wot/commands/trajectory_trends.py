#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os

import numpy as np
import pandas as pd
import wot.io
import wot.tmap
from matplotlib import pyplot


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
    args = parser.parse_args(argv)

    trajectory_ds = wot.io.read_dataset(args.trajectory)
    wot.io.add_row_metadata_to_dataset(dataset=trajectory_ds, days_path=args.cell_days)
    matrix = wot.io.read_dataset(args.matrix)
    if args.gene_filter is not None:
        if os.path.isfile(args.gene_filter):
            gene_ids = pd.read_csv(args.gene_filter, index_col=0, header=None) \
                .index.values
        else:
            specified_list = args.gene_filter.split(',')
            s = set()
            for f in specified_list:
                s.add(f.strip())
            gene_ids = matrix.var.index.intersection(s)
        col_indices = matrix.var.index.isin(gene_ids)
        if np.sum(col_indices) is 0:
            raise ValueError('No genes passed the gene filter')
        matrix = matrix[:, col_indices].copy()

    results = wot.tmap.compute_trajectory_trends_from_trajectory(trajectory_ds, matrix)
    # output genes on columns, time on rows, one file per trajectory

    genes = matrix.var.index
    for j in range(len(results)):  # each trajectory
        mean, variance = results[j]

        mean.obs.index = mean.obs.index.astype('category')
        variance.obs.index = variance.obs.index.astype('category')
        trajectory_name = trajectory_ds.var.index.values[j]
        basename = args.out + '_' + trajectory_name
        wot.io.write_dataset(mean, basename + '.mean', args.format)
        wot.io.write_dataset(variance, basename + '.variance', args.format)

        if args.plot:
            pyplot.figure(figsize=(5, 5))

            stds = np.sqrt(variance.X)
            timepoints = mean.obs.index.values.astype(float)

            for i in range(mean.shape[1]):  # each gene

                mean_i = mean[:, i].X
                # std_i = stds[:, i] if len(stds.shape) > 1 else stds[i]
                pyplot.plot(timepoints, mean_i, label=genes[i])
                # pyplot.fill_between(timepoints, mean_i - std_i,
                #                     mean_i + std_i, alpha=.4)

            pyplot.xlabel("Day")
            pyplot.ylabel("Gene expression")
            pyplot.legend()
            pyplot.savefig(basename + '.png')
