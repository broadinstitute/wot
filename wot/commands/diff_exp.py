#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import logging
import sys

import anndata

import wot.io

logger = logging.getLogger('wot')


def main(argv):
    parser = argparse.ArgumentParser(
        description='Compute differentially expressed genes from the output of the fate tool')
    parser.add_argument('--matrix', help=wot.commands.MATRIX_HELP, required=True)
    parser.add_argument('--fates', help='One or more fate datasets as produced by the fate tool',
                        action='append')
    parser.add_argument('--cell_days', help=wot.commands.CELL_DAYS_HELP, required=True)
    parser.add_argument('--compare',
                        help='If "match", compare fates with the same name. ' + 'If "all", compare all pairs. '
                             + 'If "within" compare within a fate. If a fate name, compare to the specified fate',
                        default='all')
    parser.add_argument('--delta',
                        help='Delta days to compare sampled expression matrix against within a fate. If not specified all comparison are done against the first day.',
                        type=float)
    parser.add_argument('--nperm',
                        help='Number of permutations', type=int)
    parser.add_argument('--smooth_p_values',
                        help='Smooth p-values', action='store_true')
    parser.add_argument('--fold_change', type=float, default=0.25,
                        help='Limit permutations to genes which show at least X-fold difference (log-scale) between the two groups of cells.')
    parser.add_argument('--cell_days_field', help='Field name in cell_days file that contains cell days',
                        default='day')
    parser.add_argument('--cell_day_filter',
                        help='Comma separated list of days to include (e.g. 12,14,16)', type=str)
    parser.add_argument('--gene_filter',
                        help='File with one gene id per line')
    parser.add_argument('--verbose', help='Print progress', action='store_true')
    args = parser.parse_args(argv)
    if args.verbose:
        logger.setLevel(logging.DEBUG)
        logger.addHandler(logging.StreamHandler())
    compare = args.compare
    fate_files = args.fates
    smooth_p_values = args.smooth_p_values
    expression_file = args.matrix
    delta_days = args.delta
    cell_days_file = args.cell_days
    nperm = args.nperm
    min_fold_change = args.fold_change
    cell_days_field = args.cell_days_field
    adata = wot.io.read_dataset(expression_file, var_filter=args.gene_filter)
    day_filter = args.cell_day_filter

    if delta_days is None:
        delta_days = 0
    delta_days = abs(delta_days)
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
                                         nperm=nperm, min_fold_change=min_fold_change, smooth_p_values=smooth_p_values,
                                         compare=compare, delta_days=delta_days)
    for name in diff_exp_results:
        df = diff_exp_results[name]
        df.to_csv(name + '.csv', header=True, index_label='id')
