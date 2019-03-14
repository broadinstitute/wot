#!/usr/bin/env python
# -*- coding: utf-8 -*-

import anndata
import numpy as np
import pandas as pd
import scipy.sparse
import wot.tmap


def trajectory_similarity_score(p1, p2):
    return 1.0 - 0.5 * np.sum(np.abs(p1 - p2))


def trajectory_similarities(trajectory_ds):
    """
    Computes the similarity for all pairs of trajectories across time.

    Parameters
    ----------
    trajectory_ds : anndata.AnnData
       anndata.AnnData returned by wot.tmap.TransportModel.compute_trajectories

    Returns
    -------
    distances : dict
       A dict that maps names of trajectory pairs to a dict containing 'similarity' and 'time'.
       Each element in the list is a dict containing similarity and time
    """
    # group by time

    distances = {}
    split_by_day_dict = wot.split_anndata(trajectory_ds, 'day')
    split_by_day_keys = list(split_by_day_dict.keys())
    # for each pair of trajectories
    for i in range(1, trajectory_ds.X.shape[1]):
        for j in range(i):
            similarities = np.zeros(len(split_by_day_dict))
            times = np.zeros(len(split_by_day_dict))
            for k in range(len(split_by_day_keys)):
                split_by_day_ds = split_by_day_dict[split_by_day_keys[k]]
                similarities[k] = trajectory_similarity_score(split_by_day_ds.X[:, i], split_by_day_ds.X[:, j])
                times[k] = split_by_day_ds.obs.day.values[0]

            distances[(trajectory_ds.var.index.values[i], trajectory_ds.var.index.values[j])] = {
                'similarity': similarities, 'time': times}
    return distances


def compute_trajectory_trends_from_trajectory(trajectory_ds, expression_ds):
    """
    Computes the mean and variance of each gene over time for the given trajectories

    Parameters
    ----------
    trajectory_ds : anndata.AnnData
       anndata.AnnData returned by wot.tmap.TransportModel.compute_trajectories
    ds : anndata.AnnData
        Dataset used to compute mean and variance

    Returns
    -------
    results : list
        The list of mean and variance datasets, one dataset per trajectory
        The dataset has time on the rows and genes on the columns
    """

    # align gene expression matrix with trajectory matrix
    ds_indices = trajectory_ds.obs.index.get_indexer_for(expression_ds.obs.index)
    ds_indices = ds_indices[ds_indices != -1]
    if len(ds_indices) != trajectory_ds.X.shape[0]:
        raise ValueError('Dataset does not match transport map')
    expression_ds = expression_ds[ds_indices].copy()
    timepoints = []
    mean_list = []
    variance_list = []
    for j in range(trajectory_ds.shape[1]):
        mean_list.append(None)
        variance_list.append(None)

    for day, group in trajectory_ds.obs.groupby('day'):
        timepoints.append(day)
        cell_indices_at_day = trajectory_ds.obs.index.get_indexer_for(group.index)
        trajectory_weights = trajectory_ds[cell_indices_at_day].X
        expression_values = expression_ds[cell_indices_at_day].X
        if scipy.sparse.isspmatrix(expression_values):
            expression_values = expression_values.toarray()
        # if inverse_log:
        #     expression_values = np.expm1(expression_values)

        for j in range(trajectory_ds.shape[1]):  # each trajectory
            weights_per_cell = trajectory_weights[:, j] if len(trajectory_weights.shape) > 1 else trajectory_weights

            mean = np.average(expression_values, weights=weights_per_cell, axis=0)
            var = np.average((expression_values - mean) ** 2, weights=weights_per_cell, axis=0)

            if mean_list[j] is None:
                mean_list[j] = mean.T
                variance_list[j] = var.T
            else:
                mean_list[j] = np.vstack((mean_list[j], mean.T))
                variance_list[j] = np.vstack((variance_list[j], var.T))

    results = []
    for j in range(len(mean_list)):
        mean_ds = anndata.AnnData(mean_list[j], pd.DataFrame(index=timepoints), expression_ds.var.copy())
        variance_ds = anndata.AnnData(variance_list[j], pd.DataFrame(index=timepoints), expression_ds.var.copy())
        results.append((mean_ds, variance_ds))

    return results


def compute_trajectory_trends(tmap_model, *populations):
    """
    Computes the mean and variance of each gene over time for the given populations

    Parameters
    ----------
    tmap_model : wot.TransportMapModel
        The TransportMapModel used to find ancestors and descendants of the population
    *populations : wot.Population
        The target populations

    Returns
    -------
    timepoints : 1-D array
        The list of timepoints indexing the other two return values
    means : ndarray
        The list of the means of each gene at each timepoint
    variances : ndarray
        The list of the variances of each gene at each timepoint

    Notes
    -----
    If only one population is given, means and variances will have two dimensions, otherwise three
    """
    initial_populations = populations
    timepoints = []
    traj, variances = [], []

    def update(head, populations):
        x = 0 if head else len(traj)
        m, v = tmap_model.population_mean_and_variance(*populations)
        timepoints.insert(x, wot.tmap.unique_timepoint(*populations))
        traj.insert(x, m)
        variances.insert(x, v)

    update(True, populations)
    while tmap_model.can_pull_back(*populations):
        populations = tmap_model.pull_back(*populations, as_list=True)
        update(True, populations)
    populations = initial_populations
    while tmap_model.can_push_forward(*populations):
        populations = tmap_model.push_forward(*populations, as_list=True)
        update(False, populations)

    def unpack(arr):
        arr = np.asarray(arr)
        if arr.ndim == 3:
            # rearrange dimensions when more than one population is passed
            arr = [arr[:, i, :] for i in range(arr.shape[1])]
        return np.asarray(arr) if len(arr) > 1 else arr[0]

    return timepoints, unpack(traj), unpack(variances)
