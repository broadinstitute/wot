#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import numpy as np
import pandas as pd
import wot
import wot.io

def compute_ancestor_census(ot_model, cset_matrix, *populations):
    """
    Computes the census for the populations (for both ancestors and descendants).

    Parameters
    ----------
    ot_model : wot.OTMolde
        The OTModel used to find ancestors and descendants of the population
    cset_matrix : wot.Dataset
        The cell set matrix, cells as rows, cell sets as columns. 1s denote membership.
    *populations : wot.Population
        The target populations
    """
    initial_populations = populations
    timepoints = []
    census = []

    def update(head, populations):
        x = 0 if head else len(census)
        timepoints.insert(x, wot.model.unique_timepoint(*populations))
        census.insert(x, ot_model.population_census(cset_matrix, *populations))

    update(True, populations)
    while ot_model.can_pull_back(*populations):
        populations = ot_model.pull_back(*populations, as_list=True)
        update(True, populations)
    populations = initial_populations
    while ot_model.can_push_forward(*populations):
        populations = ot_model.push_forward(*populations, as_list=True)
        update(False, populations)
    census = np.asarray(census)
    if census.ndim == 3:
        # rearrange dimensions when more than one population is passed
        census = np.asarray([census[:,i,:] for i in range(census.shape[1])])
    return timepoints, census

def main(argv):
    parser = argparse.ArgumentParser(
            description='Generate ancestor census for each time point given an initial cell set')
    parser.add_argument('--matrix', help=wot.commands.MATRIX_HELP, required=True)
    parser.add_argument('--tmap', help=wot.commands.TMAP_HELP, required=True)
    parser.add_argument('--cell_days', help=wot.commands.CELL_DAYS_HELP, required=True)
    parser.add_argument('--cell_set', help=wot.commands.CELL_SET_HELP, required=True)
    parser.add_argument('--time', help='The starting timepoint at which to consider the cell sets', required=True)
    parser.add_argument('--out', help='Output files prefix', default='census')

    args = parser.parse_args(argv)

    ot_model = wot.load_ot_model(args.matrix, args.cell_days, args.tmap)
    cell_sets_matrix = wot.io.read_gene_sets(args.cell_set)
    cell_sets = wot.io.read_cell_sets(args.cell_set)
    populations = ot_model.population_from_cell_sets(cell_sets, at_time=args.time)

    timepoints, census = compute_ancestor_census(ot_model, cell_sets_matrix, *populations.values())

    row_meta = pd.DataFrame([], index=timepoints, columns=[])
    keys = populations.keys()
    for i in range(len(census)):
        cs_name = keys[i]
        res = wot.Dataset(census[i], row_meta, cell_sets_matrix.col_meta)
        wot.io.write_dataset(res, args.out + '_' + cs_name, output_format='txt', txt_full=False)
