#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import numpy as np
import pandas as pd
import wot
import wot.io

def ancestor_census(ot_model, cset_matrix, *populations):
    timepoints = [wot.model.unique_timepoint(*populations)]
    census = []
    census.append(ot_model.population_census(cset_matrix, *populations))
    while ot_model.can_pull_back(*populations):
        populations = ot_model.pull_back(*populations)
        timepoints.append(wot.model.unique_timepoint(*populations))
        census.append(ot_model.population_census(cset_matrix, *populations))
    census = np.asarray(census)[::-1,:,:]
    return timepoints[::-1], [ census[:,i,:] for i in range(census.shape[1]) ]

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
    cell_sets = wot.io.read_cell_sets(args.cell_set)
    cell_sets_matrix = wot.io.read_gene_sets(args.cell_set)
    keys = list(cell_sets.keys())
    populations = ot_model.population_from_ids(*[ cell_sets[name] for name in keys ], at_time = float(args.time))
    # Get rid of empty populations : just ignore them
    keys = [ keys[i] for i in range(len(keys)) if populations[i] is not None ]
    populations = [ p for p in populations if p is not None ]

    timepoints, census = ancestor_census(ot_model, cell_sets_matrix, *populations)

    row_meta = pd.DataFrame([], index=timepoints, columns=[])
    for i in range(len(census)):
        cs_name = keys[i]
        res = wot.Dataset(census[i], row_meta, cell_sets_matrix.col_meta)
        wot.io.write_dataset(res, args.out + '_' + cs_name, output_format='txt', txt_full=False)
