#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import wot.io


def main(argv):
    parser = argparse.ArgumentParser('Generate ancestors and descendants of cell sets at the given time')
    parser.add_argument('--tmap', help=wot.commands.TMAP_HELP, required=True)
    parser.add_argument('--cell_set', help=wot.commands.CELL_SET_HELP, required=True)
    parser.add_argument('--time', help='Timepoint to consider', required=True)
    parser.add_argument('--out', help='Output file name', default='wot_trajectory.txt')
    args = parser.parse_args(argv)
    tmap_model = wot.tmap.TransportMapModel.from_directory(args.tmap)
    cell_sets = wot.io.read_sets(args.cell_set, as_dict=True)
    populations = tmap_model.population_from_cell_sets(cell_sets, at_time=args.time)
    trajectories = tmap_model.compute_trajectories(populations)
    wot.io.write_dataset(trajectories, args.out)
