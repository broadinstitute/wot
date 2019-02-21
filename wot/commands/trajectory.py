#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import numpy as np

import wot.io


def main(argv):
    parser = argparse.ArgumentParser(
        'Generate ancestors and descendants of cell sets generated at the given time and computes divergence between trajectories across time')
    parser.add_argument('--tmap', help=wot.commands.TMAP_HELP, required=True)
    parser.add_argument('--cell_set', help=wot.commands.CELL_SET_HELP, required=True)
    parser.add_argument('--day', help='Day to consider', required=True, type=float)
    parser.add_argument('--out', help='Output file name', default='wot_trajectory')
    parser.add_argument('--format', help='Output trajectory matrix file format', default='txt')
    args = parser.parse_args(argv)
    tmap_model = wot.tmap.TransportMapModel.from_directory(args.tmap)
    cell_sets = wot.io.read_sets(args.cell_set, as_dict=True)
    populations = tmap_model.population_from_cell_sets(cell_sets, at_time=args.day)
    trajectory_ds = tmap_model.compute_trajectories(populations)
    # for each timepoint, compute all pairwise distances

    # dataset has cells on rows and cell sets on columns
    wot.io.write_dataset(trajectory_ds, args.out, args.format)

    unique_days = list(set(trajectory_ds.obs['day']))
    pair_to_divergenes = {}
    for i in range(trajectory_ds.shape[1]):
        for j in range(i):
            pair_to_divergenes[trajectory_ds.var.index[i], trajectory_ds.var.index[j]] = []

    # create table with time on columns and trajectory pairs on rows
    for day in unique_days:
        trajectory_at_day = trajectory_ds[trajectory_ds.obs['day'] == day]
        for i in range(trajectory_ds.shape[1]):
            for j in range(i):
                divergence = 1.0 - 0.5 * np.sum(np.abs(trajectory_at_day.X[:, i] - trajectory_at_day.X[:, j]))
                pair_to_divergenes[trajectory_ds.var.index[i], trajectory_ds.var.index[j]].append(divergence)
    # pair_names = []
    # pair_divergences = []
    # for pair in pair_to_divergenes:
    #     pair_names.append(pair[0] + ' vs. ' + pair[1])
    #     pair_divergences.append(pair_to_divergenes[pair])
    # pd.DataFrame(index=pair_names, data=pair_divergences, columns=unique_days).to_csv(args.out + '_divergence.txt',
    #                                                                                   sep='\t')

    with open(args.out + '_divergence.txt', 'w') as f:
        f.write('pair' + '\t' + 'time' + '\t' + 'divergence' + '\n')
        for pair in pair_to_divergenes:
            divergenes = pair_to_divergenes[pair]
            for i in range(len(divergenes)):
                f.write(pair[0] + ' vs. ' + pair[1])
                f.write('\t')
                f.write(str(unique_days[i]))
                f.write('\t')
                f.write(str(divergenes[i]))
                f.write('\n')
