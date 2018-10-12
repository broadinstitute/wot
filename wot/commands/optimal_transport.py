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
                        help='Prefix for output file names')
    # parser.add_argument('--format', default='loom', help='Transport map file format.',
    #                     choices=wot.commands.FORMAT_CHOICES)
    args = parser.parse_args(argv)
    ot_model = wot.initialize_ot_model(args.matrix, args.cell_days,
                                       tmap_out=args.out,
                                       local_pca=args.local_pca,
                                       growth_iters=args.growth_iters,
                                       epsilon=args.epsilon,
                                       lambda1=args.lambda1,
                                       lambda2=args.lambda2,
                                       max_threads=args.max_threads,
                                       epsilon0=args.epsilon0,
                                       tau=args.tau,
                                       day_pairs=args.config,
                                       cell_day_filter=args.cell_day_filter,
                                       cell_growth_rates=args.cell_growth_rates,
                                       gene_filter=args.gene_filter,
                                       cell_filter=args.cell_filter,
                                       # output_file_format=args.format,
                                       sampling_bias=args.sampling_bias,
                                       scaling_iter=args.scaling_iter,
                                       inner_iter_max=args.inner_iter_max,
                                       force=args.force,
                                       ncells=args.ncells,
                                       ncounts=args.ncounts
                                       )
    ot_model.compute_all_transport_maps()
