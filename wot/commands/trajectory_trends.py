#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os

import h5py
import numpy as np
import pandas as pd
import scipy
import wot.io
import wot.ot

def compute_trajectory_trends(ot_model, *populations):
    """
    Computes the mean and variance of each gene over time for the given populations

    Parameters
    ----------
    ot_model : wot.OTModel
        The OTModel used to find ancestors and descendants of the population
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

    Examples
    --------
    >>> timepoints, means, variances = compute_trajectory_trends(ot_model, pop1, pop2, pop3)
    >>> means[i][j][k] # -> the mean value of the ancestors of population i at time j for gene k
    """
    initial_populations = populations
    timepoints = []
    traj, variances = [], []
    def update(head, populations):
        x = 0 if head else len(traj)
        m, v = ot_model.population_mean_and_variance(*populations)
        timepoints.insert(x, wot.model.unique_timepoint(*populations))
        traj.insert(x, m); variances.insert(x, v)

    update(True, populations)
    while ot_model.can_pull_back(*populations):
        populations = ot_model.pull_back(*populations, as_list=True)
        update(True, populations)
    populations = initial_populations
    while ot_model.can_push_forward(*populations):
        populations = ot_model.push_forward(*populations, as_list=True)
        update(False, populations)

    def unpack(arr):
        arr = np.asarray(arr)
        arr = [ arr[:,i,:] for i in range(arr.shape[1]) ]
        return arr if len(arr) > 1 else arr[0]
    return timepoints, unpack(traj), unpack(variances)

def main(argv):
    parser = argparse.ArgumentParser(description='Generate mean expression profiles for '\
            'ancestors and descendants of each cell set at the given timepoint')
    parser.add_argument('--matrix', help=wot.commands.MATRIX_HELP, required=True)
    parser.add_argument('--cell_days', help=wot.commands.CELL_DAYS_HELP, required=True)
    parser.add_argument('--tmap', help=wot.commands.TMAP_HELP, required=True)
    parser.add_argument('--cell_set', help=wot.commands.CELL_SET_HELP, required=True)
    parser.add_argument('--time', help='Timepoint to consider', required=True)
    parser.add_argument('--out', help='Prefix for output file names', default='wot_trajectory')

    args = parser.parse_args(argv)

    ot_model = wot.load_ot_model(args.matrix, args.cell_days, args.tmap)
    cell_sets = wot.io.read_cell_sets(args.cell_set)
    populations = ot_model.population_from_cell_sets(cell_sets, at_time=args.time)

    if len(populations) == 0:
        raise ValueError("No cells from the given cell sets are present at that time")

    timepoints, trends, variances = compute_trajectory_trends(ot_model, *populations.values())

    row_meta = pd.DataFrame([], index=timepoints, columns=[])
    col_meta = ot_model.matrix.col_meta.copy()
    keys = populations.keys()
    for i in range(len(trends)):
        cs_name = keys[i]
        res = wot.Dataset(trends[i], row_meta, col_meta)
        # TODO: write the variances to a different file if a flag is passed
        wot.io.write_dataset(res, args.out + '_' + cs_name, output_format='txt', txt_full=False)
