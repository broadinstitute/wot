#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging

import anndata
import numpy as np
import pandas as pd
import scipy.sparse
import sklearn
import sklearn.decomposition

import wot

logger = logging.getLogger('wot')


def trajectory_divergence(adata: anndata.AnnData, trajectory_datasets: [anndata.AnnData], cell_days_field: str = 'day',
                          local_pca: int = 30, distance_metric: str = 'emd', compare: str = 'all') -> pd.DataFrame:
    """
    Multi-dimensional linear interpolation.
    Returns the N-dimensional piecewise linear interpolant to a function
    with given values at discrete data-points.

    Parameters
    ----------
    adata :
        The expression matrix.
    trajectory_datasets :
        List of trajectory datasets.
    cell_days_field
        Cell days field in trajectory datasets.
    local_pca :
        Number of PCA components or 0 to disable
    distance_metric :
        emd (earth movers distance) or total_variation.
    compare :
        within, match, all, or trajectory name

    Returns
    -------
    df :
        A dataframe with the columns name1, name2, day1, day2, distance
    """

    adata = anndata.AnnData(adata.X, adata.obs.copy(), adata.var)
    has_covariate = False
    batch_output = []  # name1, name2, covariate_1, covariate_2, day1, day2, distance
    # if covariate_field is not None and covariate_field in adata.obs:
    #     has_covariate = True
    trajectory_names = []
    if not isinstance(trajectory_datasets, list) and not isinstance(trajectory_datasets, tuple):
        trajectory_datasets = [trajectory_datasets]

    # add trajectory dataset to obs
    for trajectory_ds in trajectory_datasets:
        trajectory_names += list(trajectory_ds.var.index)
        adata.obs = adata.obs.join(
            pd.DataFrame(index=trajectory_ds.obs.index, data=trajectory_ds.X, columns=trajectory_ds.var.index))
    # add cell days from trajectory dataset to obs
    for trajectory_ds in trajectory_datasets:
        adata.obs = adata.obs.combine_first(trajectory_ds.obs[[cell_days_field]])

    unique_days = np.array(sorted(trajectory_ds.obs[cell_days_field].unique().astype(float)))
    unique_days = unique_days[np.isnan(unique_days) == False]
    logger.info('{} days'.format(len(unique_days)))
    comparisons = wot.tmap.generate_comparisons(comparison_names=trajectory_names, compare=compare,
                                                days=unique_days,
                                                delta_days=0, reference_day='end')
    output = []  # name1, name2, day1, day2, distance

    for comparison in comparisons:
        names = comparison[0]
        days = comparison[1]
        name1 = names[0]
        name2 = names[1]
        day1 = days[0]
        day2 = days[1]
        adata_both_days = adata[((adata.obs[cell_days_field] == day1) | (adata.obs[cell_days_field] == day2))
                                & ((np.isnan(adata.obs[name1]) == False) | (
                np.isnan(adata.obs[name2]) == False))]
        x = adata_both_days.X
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

        # if has_covariate:
        #     # compute distances between all pairs of batches
        #     batches = adata_day.obs[covariate_field].unique()
        #     for batch_index1 in range(len(batches)):
        #         batch1_indices = np.where(adata_day.obs[covariate_field] == batches[batch_index1])
        #         for batch_index2 in range(batch_index1):
        #             batch2_indices = np.where(adata_day.obs[covariate_field] == batches[batch_index2])
        #             logger.info('{} vs. {}, day {}'.format(batches[batch_index1], batches[batch_index2], day))
        #             d = wot.ot.earth_mover_distance(x[batch1_indices], x[batch2_indices], eigenvals=eigenvals)
        #             batch_output.write('{}\t{}\t{}\t{}\n'.format(batches[batch_index1], batches[batch_index2], day, d))

        trajectory1_filter = adata_both_days.obs[cell_days_field] == day1
        trajectory2_filter = adata_both_days.obs[cell_days_field] == day2
        trajectory1 = adata_both_days[trajectory1_filter].obs[name1].values
        trajectory2 = adata_both_days[trajectory2_filter].obs[name2].values

        if distance_metric == 'total_variation':
            q = np.where((np.isnan(trajectory1) == False) & (np.isnan(trajectory2) == False))
            distance = 0.5 * np.sum(
                np.abs(trajectory1[q] - trajectory2[q]))
        else:
            distance = wot.ot.earth_mover_distance(x[trajectory1_filter], x[trajectory2_filter],
                                                   eigenvals=eigenvals,
                                                   weights1=trajectory1,
                                                   weights2=trajectory2)
        logger.info(
            '{} vs {}, day {}, day {}, {} cells, distance {}'.format(name1, name2, day1, day2, adata_both_days.shape[0],
                                                                     distance))
        output.append([name1, name2, day1, day2, distance])
    return pd.DataFrame(data=output, columns=['name1', 'name2', 'day1', 'day2', 'distance'])


def plot_trajectory_divergence(df):
    import matplotlib.pyplot as plt
    df['name'] = df['name1'].str.split('/').str.get(0) + ' vs. ' + df['name2'].str.split('/').str.get(
        0)
    plt.clf()
    plt.figure(figsize=(10, 10))
    plt.xlabel("Day")
    plt.ylabel("Distance")
    for p, d in df.groupby('name'):
        plt.plot(d['day2'], d['distance'], '-o', label=p)

    # if covariate_field is not None:
    #     df = pd.read_csv(args.out + '_batch.txt', sep='\t')
    #     df = df.groupby('day', as_index=False).mean()
    #     plt.plot(df['day'], df['distance'], '-o', label='covariate')
    plt.legend(loc='best')
