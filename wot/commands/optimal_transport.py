#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import argparse
import wot.io
import wot.ot
import wot.commands

def main(argv):
    parser = argparse.ArgumentParser('Compute transport maps between pairs of time points')
    wot.commands.add_model_arguments(parser)
    wot.commands.add_ot_parameters_arguments(parser)
    parser.add_argument('--day_pairs',
            help='Two column file without header with '
            'pairs of days to compute transport maps for')
    parser.add_argument('--out', default='./tmaps',
            help='Prefix for ouput file names')

    args = parser.parse_args(argv)

    # TODO: add support for the following arguments :
    # '--cell_growth_rates'
    # '--gene_filter'
    # '--cell_filter'
    # '--ncells'
    # '--ncounts'
    # '--numInnerItermax'

    tmap_dir, tmap_prefix = os.path.split(args.out)
    ot_model = wot.initialize_ot_model(args.matrix, args.cell_days,
            transport_maps_directory = tmap_dir,
            transport_maps_prefix = tmap_prefix,
            local_pca = args.local_pca,
            growth_iters = args.growth_iters,
            epsilon = args.epsilon,
            lambda1 = args.lambda1,
            lambda2 = args.lambda2,
            scaling_iter = args.scaling_iter,
            epsilon0 = args.epsilon0,
            tau = args.tau
            )

    if args.day_pairs is None:
        day_pairs = None
    else:
        day_pairs = wot.io.read_day_pairs(args.day_pairs).values

    ot_model.compute_all_transport_maps(day_pairs, force = True)
