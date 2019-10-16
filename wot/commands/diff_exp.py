#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import logging
import sys

import anndata
import wot.io

logger = logging.getLogger('wot')


def create_parser():
    parser = argparse.ArgumentParser(
        description='Compute differentially expressed genes from the output of the fate tool')
    parser.add_argument('--matrix', help=wot.commands.MATRIX_HELP, required=True)
    parser.add_argument('--fate', help='Fate dataset produced by the fate tool', required=True)
    parser.add_argument('--cell_days', help=wot.commands.CELL_DAYS_HELP)
    parser.add_argument('--out', help='Output file name', default='wot_diff_exp.csv')
    # parser.add_argument('--compare',
    #                     help='If "match", compare fates with the same name. ' + 'If "all", compare all pairs. '
    #                          + 'If "within" compare within a fate. If a fate name, compare to the specified fate',
    #                     default='all')
    # parser.add_argument('--delta',
    #                     help='Delta days to compare sampled expression matrix against within a fate. If not specified all comparison are done against the first day.',
    #                     type=float)
    parser.add_argument('--cell_days_field', help='Field name in cell_days file that contains cell days',
        default='day')
    parser.add_argument('--cell_day_filter',
        help='Comma separated list of days to include (e.g. 12,14,16)', type=str)
    parser.add_argument('--gene_filter',
        help='File with one gene id per line')
    parser.add_argument('--verbose', help='Print progress', action='store_true')
    return parser


def main(args):
    if args.verbose:
        logger.setLevel(logging.DEBUG)
        logger.addHandler(logging.StreamHandler())
    compare = 'all'  # all pairs args.compare
    fate_files = [args.fate]
    delta_days = 0  # args.delta
    expression_file = args.matrix
    cell_days_file = args.cell_days
    cell_days_field = args.cell_days_field
    adata = wot.io.read_dataset(expression_file, var_filter=args.gene_filter)
    day_filter = args.cell_day_filter

    if delta_days is None:
        delta_days = 0
    delta_days = abs(delta_days)
    if cell_days_file is not None:
        wot.io.add_row_metadata_to_dataset(dataset=adata, days=cell_days_file)
    if day_filter is not None:
        days = [float(day) for day in day_filter.split(',')] if type(day_filter) == str else day_filter
        adata = adata[adata.obs[cell_days_field].isin(days)]
        adata = anndata.AnnData(adata.X, adata.obs.copy(), adata.var)
    if adata.shape[1] is 0:
        sys.exit('Expression matrix has 0 genes')
    fate_datasets = []
    for f in fate_files:
        fate_ds = wot.io.read_dataset(f)
        if len(fate_files) > 1:
            fate_ds.var.index = fate_ds.var.index + '/' + wot.io.get_filename_and_extension(f)[0]
        # adata.obs = adata.obs.join(
        #     pd.DataFrame(index=fate_ds.obs.index, data=fate_ds.X, columns=fate_ds.var.index))
        fate_datasets.append(fate_ds)
    diff_exp_results = wot.tmap.diff_exp(adata=adata, fate_datasets=fate_datasets, cell_days_field=cell_days_field,
        compare=compare, delta_days=delta_days)
    diff_exp_results.to_csv(args.out, header=True, index_label='id')
