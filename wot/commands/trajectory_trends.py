#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import wot.io
import wot.tmap


def main(argv):
    parser = argparse.ArgumentParser(description='Generate mean expression profiles for ' \
                                                 'ancestors and descendants of each trajectory')
    parser.add_argument('--matrix', help=wot.commands.MATRIX_HELP, required=True)
    parser.add_argument('--trajectory', help='Trajectory dataset as produced by the trajectory tool', required=True)
    parser.add_argument('--out', help='Prefix for output file names', default='trends')
    parser.add_argument('--cell_days', help=wot.commands.CELL_DAYS_HELP, required=True)
    args = parser.parse_args(argv)
    # tmap_model = wot.tmap.TransportMapModel.from_directory(args.tmap)
    # cell_sets = wot.io.read_sets(args.cell_set, as_dict=True)
    # populations = tmap_model.population_from_cell_sets(cell_sets, at_time=args.time)

    trajectory_ds = wot.io.read_dataset(args.trajectory)
    wot.io.add_row_metadata_to_dataset(dataset=trajectory_ds, days_path=args.cell_days)
    matrix = wot.io.read_dataset(args.matrix)

    results = wot.tmap.compute_trajectory_trends_from_trajectory(trajectory_ds, matrix)
    # output genes on columns, time on rows, one file per trajectory

    for j in range(len(results)):
        mean, variance = results[j]
        trajectory_name = trajectory_ds.var.index.values[j]
        basename = args.out + '_' + trajectory_name
        wot.io.write_dataset(mean, basename + '.mean')
        wot.io.write_dataset(variance, basename + '.variance')
