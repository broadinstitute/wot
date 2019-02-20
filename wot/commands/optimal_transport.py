#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import wot.commands
import wot.io
import wot.ot


def main(argv):
    parser = argparse.ArgumentParser('Compute transport maps between pairs of time points')
    wot.commands.add_model_arguments(parser)
    wot.commands.add_ot_parameters_arguments(parser)
    parser.add_argument('--out', default='./tmaps',
                        help='Prefix for output file names')
    args = parser.parse_args(argv)
    ot_model = wot.commands.initialize_ot_model_from_args(args)
    ot_model.compute_all_transport_maps()
