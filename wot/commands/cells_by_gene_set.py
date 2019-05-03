#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import sys
import wot.io


def main(argv):
    parser = argparse.ArgumentParser(description='Generate cell sets from gene set scores')
    parser.add_argument('--score',
                        help='Gene sets scores generated from the gene_set_scores command', required=True,
                        action='append')
    parser.add_argument('--out', help='Output file name prefix', required=True)
    parser.add_argument('--quantile', default=99,
                        help='Quantile for cells to be considered a member of a cell set', type=float)

    args = parser.parse_args(argv)
    quantile = args.quantile
    cell_set_name_to_ids = {}
    for path in args.score:
        ds = wot.io.read_dataset(path)
        percentiles = np.percentile(ds.X, axis=0, q=quantile)
        if len(percentiles.shape) == 0:
            percentiles = [percentiles]

        for j in range(ds.shape[1]):  # each gene set
            x = ds[:, j].X
            selected = x >= percentiles[j]
            cell_ids = ds[selected].obs.index
            if len(cell_ids) > 0:
                cell_set_name_to_ids[ds.var.index[j]] = cell_ids

    if len(cell_set_name_to_ids) == 0:
        sys.exit('No cells found')

    wot.io.write_sets(cell_set_name_to_ids, args.out)
