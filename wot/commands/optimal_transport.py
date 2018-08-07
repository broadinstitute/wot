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
    parser.add_argument('--out', default='./tmaps',
            help='Prefix for ouput file names')

    args = parser.parse_args(argv)

    # TODO: add support for the following arguments :
    # '--gene_filter'
    # '--cell_filter'
    # '--ncells'
    # '--ncounts'

    tmap_dir, tmap_prefix = os.path.split(args.out)
    ot_model = wot.initialize_ot_model(args.matrix, args.cell_days,
            tmap_dir = tmap_dir,
            tmap_prefix = tmap_prefix,
            local_pca = args.local_pca,
            growth_iters = args.growth_iters,
            epsilon = args.epsilon,
            lambda1 = args.lambda1,
            lambda2 = args.lambda2,
            max_iter = args.max_iter,
            max_threads = args.max_threads,
            epsilon0 = args.epsilon0,
            tau = args.tau,
            day_pairs = args.config,
            tolerance = args.tolerance,
            batch_size = args.batch_size,
            cell_growth_rates = args.cell_growth_rates,
            )

    ot_model.compute_all_transport_maps(force = True)
