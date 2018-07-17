#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import numpy as np
import pandas as pd
import wot
import wot.io

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

    core = wot.initialize_core(args.matrix, args.cell_days, transport_maps_directory = args.tmap)
    cell_sets = wot.io.read_cell_sets(args.cell_set)
    keys = list(cell_sets.keys())
    populations = core.population_from_ids(*[ cell_sets[name] for name in keys ], at_time = float(args.time))
    cell_sets_matrix = wot.io.read_gene_sets(args.cell_set)

    timepoints = [float(args.time)]
    ancestor_census = []
    ancestor_census.append(core.population_census(cell_sets_matrix, *populations))
    while core.can_pull_back(*populations):
        populations = core.pull_back(*populations)
        timepoints.append(wot.core.unique_timepoint(*populations))
        ancestor_census.append(core.population_census(cell_sets_matrix, *populations))
    ancestor_census = np.asarray(ancestor_census)

    row_meta = pd.DataFrame([], index=timepoints[::-1], columns=[])
    for i in range(ancestor_census.shape[1]):
        cs_name = keys[i]
        res = wot.Dataset(ancestor_census[::-1, i, :], row_meta, cell_sets_matrix.col_meta)
        wot.io.write_dataset(res, args.out + '_' + cs_name, output_format='txt', txt_full=False)
