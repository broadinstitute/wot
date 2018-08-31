#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import wot.model
import scipy


def trajectory_similarity_score(p1, p2):
    return 1.0 - 0.5 * np.sum(np.abs(p1 - p2))


def trajectory_similarities(trajectory_ds):
    """
    Computes the similarity for all pairs of trajectories across time.

    Parameters
    ----------
    trajectory_ds : wot.Dataset
       Dataset returned by wot.ot.compute_trajectories

    Returns
    -------
    distances : dict
       A dict that maps names of trajectory pairs to a dict containing 'similarity' and 'time'.
       Each element in the list is a dict containing similarity and time
    """
    # group by time

    distances = {}
    split_by_day_dict = trajectory_ds.split_by('day')
    split_by_day_keys = list(split_by_day_dict.keys())
    # for each pair of trajectories
    for i in range(1, trajectory_ds.x.shape[1]):
        for j in range(i):
            similarities = np.zeros(len(split_by_day_dict))
            times = np.zeros(len(split_by_day_dict))
            for k in range(len(split_by_day_keys)):
                split_by_day_ds = split_by_day_dict[split_by_day_keys[k]]
                similarities[k] = trajectory_similarity_score(split_by_day_ds.x[:, i], split_by_day_ds.x[:, j])
                times[k] = split_by_day_ds.row_meta.day.values[0]

            distances[(trajectory_ds.col_meta.index.values[i], trajectory_ds.col_meta.index.values[j])] = {
                'similarity': similarities, 'time': times}
    return distances


def compute_trajectory_trends_from_trajectory(trajectory_ds, ds):
    """
    Computes the mean and variance of each gene over time for the given trajectories

    Parameters
    ----------
    trajectory_ds : wot.Dataset
       Dataset returned by wot.ot.compute_trajectories
    ds : wot.Dataset
        Dataset used to compute mean and variance

    Returns
    -------
    results : list
        The list of mean and variance datasets, one dataset per trajectory
        The dataset has time on the rows and genes on the columns
    """

    ds_indices = trajectory_ds.row_meta.index.get_indexer_for(ds.row_meta.index)

    if (ds_indices[ds_indices == -1]).sum() > 0:
        raise ValueError('Cell ids not found')
    ds = wot.Dataset(ds.x[ds_indices], ds.row_meta.iloc[ds_indices], ds.col_meta)
    timepoints = []
    mean_list = []
    variance_list = []
    for j in range(trajectory_ds.x.shape[1]):
        mean_list.append(None)
        variance_list.append(None)

    for day, group in trajectory_ds.row_meta.groupby('day'):
        timepoints.append(day)
        indices = trajectory_ds.row_meta.index.get_indexer_for(group.index)
        p = trajectory_ds.x[indices]
        values = ds.x[indices]
        if scipy.sparse.isspmatrix(values):
            values = values.toarray()
        for j in range(trajectory_ds.x.shape[1]):  # each trajectory
            mean = np.average(values, weights=p[:, j], axis=0)
            var = np.average((values - mean) ** 2, weights=p[:, j], axis=0)

            if mean_list[j] is None:
                mean_list[j] = mean.T
                variance_list[j] = var.T
            else:
                mean_list[j] = np.vstack((mean_list[j], mean.T))
                variance_list[j] = np.vstack((variance_list[j], var.T))
    row_meta = pd.DataFrame(index=timepoints)
    results = []
    for j in range(len(variance_list)):
        mean_ds = wot.Dataset(mean_list[j], row_meta, ds.col_meta)
        variance_ds = wot.Dataset(variance_list[j], row_meta, ds.col_meta)
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
        timepoints.insert(x, wot.model.unique_timepoint(*populations))
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


def compute_trajectories(tmap_model, population_dict):
    """
    Computes a trajectory for each population

    Parameters
    ----------
    tmap_model : wot.TransportMapModel
        The TransportMapModel used to find ancestors and descendants of the population
    *population_dict : dict of str: wot.Population
        The target populations such as ones from tmap_model.population_from_cell_sets

    Returns
    -------
    trajectories : wot.Dataset
        Rows : all cells, Columns : populations index. At point (i, j) : the probability that cell i is an
        ancestor/descendant of population j
    """
    trajectories = []
    populations = population_dict.values()
    population_names = list(population_dict.keys())
    initial_populations = populations

    # timepoints = []

    def update(head, populations):
        idx = 0 if head else len(trajectories)
        # timepoints.insert(idx, wot.model.unique_timepoint(*populations))
        trajectories.insert(idx, np.array([pop.p for pop in populations]).T)

    update(True, populations)
    while tmap_model.can_pull_back(*populations):
        populations = tmap_model.pull_back(*populations, as_list=True)
        update(True, populations)
    populations = initial_populations
    while tmap_model.can_push_forward(*populations):
        populations = tmap_model.push_forward(*populations, as_list=True)
        update(False, populations)

    # # list of trajectories. Each trajectory is a list of wot.Datasets split by time
    # name_to_list = {}
    # start = 0
    # for j in range(len(population_names)):
    #     name_to_list[population_names[j]] = []
    # for i in range(len(timepoints)):
    #     t = timepoints[i]
    #     probs = trajectories[i]
    #     ncells = probs.shape[0]
    #     row_meta = tmap_model.matrix.row_meta.iloc[np.arange(start, start + ncells)].copy()
    #     row_meta['day'] = t
    #
    #     for j in range(probs.shape[1]):
    #         name_to_list[population_names[j]].append(
    #             wot.Dataset(x=probs[:, j], row_meta=row_meta, col_meta=pd.DataFrame(index=[population_names[j]])))
    #     start += ncells
    # return list(name_to_list.values())
    return wot.Dataset(x=np.concatenate(trajectories), row_meta=tmap_model.meta,
                       col_meta=pd.DataFrame(index=population_names))
