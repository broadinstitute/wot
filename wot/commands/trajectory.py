#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import numpy as np
import pandas as pd
import wot.io
from matplotlib import pyplot


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

    # dataset has cells on rows and cell sets on columns
    wot.io.write_dataset(trajectory_ds, args.out + '_trajectory', args.format)
    if args.embedding:
        nbins = 500
        full_embedding_df = pd.read_csv(args.embedding, sep=None, engine='python', index_col='id')
        xrange = full_embedding_df['x'].min(), full_embedding_df['x'].max()
        yrange = full_embedding_df['y'].min(), full_embedding_df['y'].max()
        full_embedding_df['x'] = np.floor(
            np.interp(full_embedding_df['x'], [xrange[0], xrange[1]], [0, nbins - 1])).astype(int)
        full_embedding_df['y'] = np.floor(
            np.interp(full_embedding_df['y'], [yrange[0], yrange[1]], [0, nbins - 1])).astype(int)
        population_list = list(populations.values())
        for j in range(trajectory_ds.shape[1]):
            color_df = pd.DataFrame(index=trajectory_ds.obs.index, data={'color': trajectory_ds[:, j].X})
            embedding_df = color_df.join(full_embedding_df)
            figure = pyplot.figure(figsize=(10, 10))
            pyplot.axis('off')
            pyplot.tight_layout()

            pyplot.scatter(full_embedding_df['x'], full_embedding_df['y'], c='#f0f0f0',
                           s=4, marker=',', edgecolors='none', alpha=0.8)
            summed_df = embedding_df.groupby(['x', 'y'], as_index=False).agg('sum')

            pyplot.scatter(summed_df['x'], summed_df['y'], c=summed_df['color'],
                           s=6, marker=',', edgecolors='none', cmap='viridis_r', alpha=1)
            pyplot.colorbar()
            ncells = (population_list[j].p > 0).sum()
            pyplot.title('{}, day {}, {}/{} cells'.format(trajectory_ds.var.index[j], args.day, ncells,
                                                          len(population_list[j].p)))
            figure.savefig(args.out + '_' + str(trajectory_ds.var.index[j]) + '_trajectory.png')
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
                divergence = 0.5 * np.sum(np.abs(trajectory_at_day[:, i].X - trajectory_at_day[:, j].X))
                pair_to_divergenes[trajectory_ds.var.index[i], trajectory_ds.var.index[j]].append(divergence)

    if len(populations) > 1:
        with open(args.out + '_trajectory_divergence.txt', 'w') as f:
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
            divergence_df = pd.read_csv(args.out + '__trajectory_divergence.txt', sep='\t')
            pyplot.figure(figsize=(10, 10))

            pyplot.xlabel("Day")
            pyplot.ylabel("Divergence")

            for p, d in divergence_df.groupby('pair'):
                day = np.asarray(d['time'])
                divergence = np.asarray(d['divergence'])
                pyplot.plot(day, divergence, label=p)
            pyplot.legend()
            pyplot.savefig(args.out + '_trajectory_divergence.png')
            # plot all pairs over time
