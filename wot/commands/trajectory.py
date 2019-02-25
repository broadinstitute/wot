#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import numpy as np
import pandas as pd
from matplotlib import pyplot

import wot.io


def main(argv):
    parser = argparse.ArgumentParser(
        'Generate ancestors and descendants of cell sets generated at the given time and computes divergence between trajectories across time')
    parser.add_argument('--tmap', help=wot.commands.TMAP_HELP, required=True)
    parser.add_argument('--cell_set', help=wot.commands.CELL_SET_HELP, required=True)
    parser.add_argument('--day', help='Day to consider', required=True, type=float)
    parser.add_argument('--out', help='Prefix for output file names', default='trajectory')
    parser.add_argument('--format', help='Output trajectory matrix file format', default='txt')
    parser.add_argument('--embedding', help='Optional file with id, x, y used to plot trajectory probabilities')
    parser.add_argument('--plot_divergence', help='Whether to plot divergence over time', action='store_true')
    args = parser.parse_args(argv)
    tmap_model = wot.tmap.TransportMapModel.from_directory(args.tmap)
    cell_sets = wot.io.read_sets(args.cell_set, as_dict=True)
    populations = tmap_model.population_from_cell_sets(cell_sets, at_time=args.day)
    trajectory_ds = tmap_model.compute_trajectories(populations)
    # for each timepoint, compute all pairwise distances

    # dataset has cells on rows and cell sets on columns
    wot.io.write_dataset(trajectory_ds, args.out, args.format)
    if args.embedding:
        full_embedding_df = pd.read_csv(args.embedding, sep=None, engine='python', index_col='id')
        for j in range(trajectory_ds.shape[1]):
            color_df = pd.DataFrame(index=trajectory_ds.obs.index, data={'color': trajectory_ds.X[:, j]})
            embedding_df = full_embedding_df.copy().join(color_df)
            figure = pyplot.figure(figsize=(10, 10))
            pyplot.axis('off')
            pyplot.tight_layout(pad=0)
            pyplot.scatter(embedding_df['x'], embedding_df['y'], c=embedding_df['color'],
                           s=0.8, marker=',', edgecolors='none', cmap='viridis')
            pyplot.colorbar()
            pyplot.title(str(trajectory_ds.var.index[j]))
            figure.savefig(args.out + '_' + str(trajectory_ds.var.index[j]) + '.png')
        # plot probabilties on embedding

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
                divergence = 0.5 * np.sum(np.abs(trajectory_at_day.X[:, i] - trajectory_at_day.X[:, j]))
                pair_to_divergenes[trajectory_ds.var.index[i], trajectory_ds.var.index[j]].append(divergence)

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

    if args.plot_divergence:
        divergence_df = pd.read_csv(args.out + '_divergence.txt', sep='\t')
        pyplot.figure(figsize=(10, 10))

        pyplot.xlabel("Day")
        pyplot.ylabel("Divergence")

        for p, d in divergence_df.groupby('pair'):
            day = np.asarray(d['time'])
            divergence = np.asarray(d['divergence'])
            pyplot.plot(day, divergence, label=p)
        pyplot.legend()
        pyplot.savefig(args.out + '_divergence.png')
        # plot all pairs over time
