#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import wot.io


def main(argv):
    parser = argparse.ArgumentParser(
        description='Generate a transition table from one cell set to another cell set')
    parser.add_argument('--tmap', help=wot.commands.TMAP_HELP, required=True)
    parser.add_argument('--cell_set', help=wot.commands.CELL_SET_HELP, required=True, action='append')
    parser.add_argument('--start_time',
                        help='The start time for the cell sets to compute the transitions to cell sets at end_time',
                        required=True, type=float)
    parser.add_argument('--end_time', help='The end time', required=True, type=float)
    parser.add_argument('--out', help='Prefix for ouput file.')
    parser.add_argument('--format', help=wot.commands.FORMAT_HELP, default='h5ad', choices=wot.commands.FORMAT_CHOICES)

    args = parser.parse_args(argv)

    tmap_model = wot.tmap.TransportMapModel.from_directory(args.tmap)
    cell_sets = wot.io.read_sets(args.cell_set, as_dict=True)
    start_populations = tmap_model.population_from_cell_sets(cell_sets, at_time=args.start_time)
    end_populations = tmap_model.population_from_cell_sets(cell_sets, at_time=args.end_time)
    ds = tmap_model.transition_table(start_populations, end_populations)
    wot.io.write_dataset(ds, args.out + '_' + str(args.start_time) + '_' + str(args.end_time) + '_transition_table',
                         output_format=args.format)
