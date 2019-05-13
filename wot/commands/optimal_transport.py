#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import logging

import wot.commands
import wot.io
import wot.ot


def main(argv):
    parser = argparse.ArgumentParser('Compute transport maps between pairs of time points')
    wot.commands.add_ot_parameters_arguments(parser)
    parser.add_argument('--format', help='Output file format', default='h5ad', choices=['h5ad', 'loom'])
    parser.add_argument('--no_overwrite', help='Do not overwrite existing transport maps if they exist',
                        action='store_true')
    parser.add_argument('--out', default='./tmaps',
                        help='Prefix for output file names')
    args = parser.parse_args(argv)
    if args.verbose:
        logger = logging.getLogger('wot')
        logger.setLevel(logging.DEBUG)
        logger.addHandler(logging.StreamHandler())
    ot_model = wot.commands.initialize_ot_model_from_args(args)
    ot_model.compute_all_transport_maps(overwrite=not args.no_overwrite, output_file_format=args.format,
                                        tmap_out=args.out)
