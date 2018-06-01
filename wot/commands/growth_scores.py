#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import wot.ot
import pandas as pd
import numpy as np
import sklearn.metrics.pairwise
import csv
import ot as pot
import os
import io


def main(argv):
    parser = argparse.ArgumentParser(
        description='Convert gene set scores to growth scores')

    parser.add_argument('--gene_set_scores', required=True, help='File containing "Cell.cycle" and "Apoptosis" scores')
    parser.add_argument('--beta_max', type=float, default=1.7, help='Growth function parameter')
    parser.add_argument('--beta_center', type=float, default=0.25, help='Growth function parameter')
    parser.add_argument('--delta_max', type=float, default=1.7, help='Growth function parameter')
    parser.add_argument('--output', help='Output file name', required=True)
    args = parser.parse_args(argv)

    gene_set_scores = pd.read_table(args.gene_set_scores, index_col=0,
                                    quoting=csv.QUOTE_NONE, engine='python',
                                    sep=None)
    apoptosis = gene_set_scores['Apoptosis']
    proliferation = gene_set_scores['Cell.cycle']
    g = wot.ot.compute_growth_scores(proliferation.values, apoptosis.values, beta_max=args.beta_max,
                                     beta_center=args.beta_center,
                                     delta_max=args.delta_max)
    cell_growth_rates = pd.DataFrame(index=gene_set_scores.index,
                                     data={'cell_growth_rate': g})

    cell_growth_rates.to_csv(args.output, sep='\t', header=False)
