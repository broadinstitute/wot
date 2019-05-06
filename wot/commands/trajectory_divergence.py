#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import numpy as np
import pandas as pd
import scipy.sparse
import sklearn
import sklearn.decomposition

import wot
import wot.ot


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
                        help='Compare across trajectories when more than one trajectory is supplied. If value is "match" only compare trajectories with the same name. If "all", compare all pairs',
                        choices=['match', 'all'], required=True)
    parser.add_argument('--local_pca', type=int, default=30,
                        help='Convert day matrix to local PCA coordinates.'
                             'Set to 0 to disable')
    parser.add_argument('--cell_days', help=wot.commands.CELL_DAYS_HELP, required=True)
    parser.add_argument('--plot', help='Plot results', action='store_true')
    parser.add_argument('--cell_filter', help='File with one cell id per line to include')
    parser.add_argument('--gene_filter',
                        help='File with one gene id per line to use for computing'
                             'cost matrices (e.g. variable genes)')
    parser.add_argument('--covariate',
                        help='Covariate (batch) values for each cell. Used to compute batch to batch distance within a timepoint.')
    parser.add_argument('--cell_days_field', help='Field name in cell_days file that contains cell days',
                        default='day')
    parser.add_argument('--covariate_field',
                        help='Field name in covariate file that contains covariate',
                        default='covariate')
    args = parser.parse_args(argv)
    compare = args.compare
    cell_days_field = args.cell_days_field
    local_pca = args.local_pca
    expression_matrix = wot.io.read_dataset(args.matrix)
    trajectory_files = args.trajectory
    batch_field_name = args.covariate_field
    wot.io.add_row_metadata_to_dataset(expression_matrix, days=args.cell_days, covariate=args.covariate)

    distance_metric = args.distance_metric

    trajectory_names = []

    for f in trajectory_files:
        trajectory_ds = wot.io.read_dataset(f)
        if len(trajectory_files) > 1:
            trajectory_ds.var.index = trajectory_ds.var.index + '/' + wot.io.get_filename_and_extension(f)[0]
        expression_matrix.obs = expression_matrix.obs.join(
            pd.DataFrame(index=trajectory_ds.obs.index, data=trajectory_ds.X, columns=trajectory_ds.var.index))
        trajectory_names += list(trajectory_ds.var.index)

    if args.cell_filter is not None:
        expression_matrix = expression_matrix[
            expression_matrix.obs.index.isin(wot.io.read_sets(args.cell_filter).obs.index.values)].copy()

    if args.gene_filter is not None:
        expression_matrix = expression_matrix[:, expression_matrix.var.index.isin(
            wot.io.read_sets(args.gene_filter).obs.index.values)].copy()

    output = open(args.out + '.txt', 'wt')
    output.write('trajectory_1' + '\t' + 'trajectory_2' + '\t' + 'day' + '\t' + 'distance\n')

    if batch_field_name is not None:
        batch_output = open(args.out + '_batch.txt', 'wt')
        batch_output.write('covariate_1' + '\t' + 'covariate_2' + '\t' + 'day' + '\t' + 'distance' + '\n')

    if compare == 'all':  # all pairs
        base_trajectory_names_to_trajectory_names = {'': trajectory_names}
    elif compare == 'match':  # all matching names
        base_trajectory_names_to_trajectory_names = {}
        for i in range(len(trajectory_names)):
            full_trajectory_name = trajectory_names[i]
            base_trajectory_name = full_trajectory_name[0:full_trajectory_name.rindex('/')]
            names = base_trajectory_names_to_trajectory_names.get(base_trajectory_name)
            if names is None:
                names = []
                base_trajectory_names_to_trajectory_names[base_trajectory_name] = names
            names.append(full_trajectory_name)
    for day, group in expression_matrix.obs.groupby(cell_days_field):
        expression_matrix_day = expression_matrix[group.index]
        x = expression_matrix_day.X
        if scipy.sparse.isspmatrix(x):
            x = x.toarray()

        eigenvals = None
        if local_pca > 0 and local_pca < x.shape[1]:
            mean_shift = x.mean(axis=0)
            x = x - mean_shift
            n_components = min(local_pca, x.shape[0])  # n_components must be <= ncells
            pca = sklearn.decomposition.PCA(n_components=n_components, random_state=58951)
            pca.fit(x.T)
            x = pca.components_.T
            eigenvals = np.diag(pca.singular_values_)
        if batch_field_name is not None:
            # compute distances between all pairs of batches
            adata_day = expression_matrix[group.index]
            batches = adata_day.obs[batch_field_name].unique()
            for batch_index1 in range(len(batches)):
                batch1_indices = np.where(adata_day.obs[batch_field_name] == batches[batch_index1])
                for batch_index2 in range(batch_index1):
                    batch2_indices = np.where(adata_day.obs[batch_field_name] == batches[batch_index2])
                    print('{} vs. {}, day {}'.format(batches[batch_index1], batches[batch_index2], day))
                    d = wot.ot.earth_mover_distance(x[batch1_indices], x[batch2_indices], eigenvals=eigenvals)
                    batch_output.write('{}\t{}\t{}\t{}\n'.format(batches[batch_index1], batches[batch_index2], day, d))

        for base_trajectory_name in base_trajectory_names_to_trajectory_names:
            trajectory_names = base_trajectory_names_to_trajectory_names[base_trajectory_name]
            for i in range(len(trajectory_names)):
                for j in range(i):
                    print('{} vs {}, day {}'.format(trajectory_names[j], trajectory_names[i], day))
                    trajectory_ds_i = expression_matrix_day.obs[trajectory_names[i]].values.astype(np.float64)
                    trajectory_ds_j = expression_matrix_day.obs[trajectory_names[j]].values.astype(np.float64)

                    if distance_metric == 'total_variation':
                        q = np.where((np.isnan(trajectory_ds_i) == False) & (np.isnan(trajectory_ds_j) == False))
                        d = 0.5 * np.sum(
                            np.abs(trajectory_ds_i[q] - trajectory_ds_j[q]))
                    else:
                        trajectory_ds_i_filter = np.where(np.isnan(trajectory_ds_i) == False)
                        trajectory_ds_j_filter = np.where(np.isnan(trajectory_ds_j) == False)
                        d = wot.ot.earth_mover_distance(x[trajectory_ds_i_filter], x[trajectory_ds_j_filter],
                                                        eigenvals=eigenvals,
                                                        weights1=trajectory_ds_i[trajectory_ds_i_filter],
                                                        weights2=trajectory_ds_j[trajectory_ds_j_filter])
                    output.write(
                        '{}\t{}\t{}\t{}\n'.format(trajectory_names[j], trajectory_names[i], day, d))
                    output.flush()

    output.close()
    if batch_field_name is not None:
        batch_output.close()

    if args.plot:
        import matplotlib.pyplot as plt
        df = pd.read_csv(args.out + '.txt', sep='\t')
        df['name'] = df['trajectory_1'].str.split('/').str.get(0) + ' vs. ' + df['trajectory_2'].str.split('/').str.get(
            0)
        plt.clf()
        plt.figure(figsize=(10, 10))
        plt.xlabel("Day")
        plt.ylabel("Distance")
        for p, d in df.groupby('name'):
            plt.plot(d['day'], d['distance'], '-o', label=p)

        if batch_field_name is not None:
            df = pd.read_csv(args.out + '_batch.txt', sep='\t')
            df = df.groupby('day', as_index=False).mean()
            plt.plot(df['day'], df['distance'], '-o', label='covariate')
        plt.legend(loc='best')
        plt.savefig(args.out + '.png')
        plt.close()
