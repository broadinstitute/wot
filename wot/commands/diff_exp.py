#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import logging
import sys

import anndata
import numpy as np
import pandas as pd
import scipy.sparse
import statsmodels.stats.multitest

import wot.io

logger = logging.getLogger('wot')


class DiffExp:

    def __init__(self, expression_matrix, delta_days, fate_names,
                 compare, nperm, min_fold_change, day_field, smooth_p_values):
        self.expression_matrix = expression_matrix
        self.delta_days = delta_days
        self.fate_names = fate_names
        if len(fate_names) == 1:
            compare = 'within'
        self.compare = compare
        self.nperm = nperm
        self.min_fold_change = min_fold_change
        self.features = expression_matrix.var.index
        self.day_field = day_field
        days = np.array(sorted(expression_matrix.obs[day_field].unique().astype(float)))
        self.days = days[np.isnan(days) == False]
        self.smooth_p_values = smooth_p_values

    def add_stats(self, expression_values, weights, df, suffix):
        # expression_values = np.expm1(expression_values)
        mean = np.average(expression_values, weights=weights, axis=0)
        fraction_expressed = weights.dot(expression_values > 0)
        # variance = np.average((expression_values - mean) ** 2, weights=weights, axis=0)
        # variance = np.log1p(variance)
        # mean = np.log1p(mean)
        if 'mean{}'.format(suffix) not in df:
            return df.join(pd.DataFrame(index=self.features,
                                        data={
                                            'mean{}'.format(suffix): mean,
                                            'fraction_expressed{}'.format(suffix): fraction_expressed
                                        }))
        return df

    def get_expression_and_weights(self, day, fate_name):
        ds = self.expression_matrix[
            (self.expression_matrix.obs[self.day_field] == day) & (
                    False == self.expression_matrix.obs[fate_name].isna())]
        weights = ds.obs[fate_name].values
        expression_values = ds.X
        if scipy.sparse.isspmatrix(expression_values):
            expression_values = expression_values.toarray()
        weights = weights / weights.sum()
        return expression_values, weights

    def do_comparison(self, expression_values1, weights1, day1, expression_values2, weights2, day2):
        # expression_values1 = np.expm1(expression_values1)
        # expression_values2 = np.expm1(expression_values2)
        mean1 = np.average(expression_values1, weights=weights1, axis=0)
        mean2 = np.average(expression_values2, weights=weights2, axis=0)
        # variance1 = np.average((expression_values1 - mean1) ** 2, weights=weights1, axis=0)
        # variance2 = np.average((expression_values2 - mean2) ** 2, weights=weights2, axis=0)
        # fold_change = np.log1p(mean1) - np.log1p(mean2)
        observed = (mean1 - mean2)
        suffix = "_{}_{}".format(day1, day2)
        results = pd.DataFrame(index=self.features, data={'fold_change' + suffix: observed})

        if self.nperm is not None and self.nperm > 0:
            genes_use = np.abs(observed) >= self.min_fold_change
            if genes_use.sum() > 0:
                expression_values1 = expression_values1[:, genes_use]
                expression_values2 = expression_values2[:, genes_use]
                observed_use = observed[genes_use]
                weights1 = weights1.copy()
                weights2 = weights2.copy()
                p = np.zeros(expression_values2.shape[1])

                for i in range(self.nperm):
                    np.random.shuffle(weights1)
                    np.random.shuffle(weights2)
                    mean1 = np.average(expression_values1, weights=weights1, axis=0)
                    mean2 = np.average(expression_values2, weights=weights2, axis=0)

                    # permuted_fold_change = np.log1p(mean1) - np.log1p(mean2)
                    permuted = (mean1 - mean2)
                    p[(permuted >= observed_use)] += 1
                # 2-sided p-value
                k = p
                if self.smooth_p_values:
                    p = (p + 1) / (self.nperm + 2)
                else:
                    p = p / self.nperm
                one_minus_p = 1.0 - p;
                expr = one_minus_p < p
                p[expr] = 1 - p[expr]
                p *= 2
                fdr = statsmodels.stats.multitest.multipletests(p)[1]
                results = results.join(
                    pd.DataFrame(index=self.features[genes_use], data={'p_value' + suffix: p,
                                                                       'fdr' + suffix: fdr,
                                                                       'k' + suffix: k}))
        return results

    def execute(self):
        comparisons = wot.tmap.generate_comparisons(comparison_names=self.fate_names, compare=self.compare,
                                                    days=self.days,
                                                    delta_days=self.delta_days, reference_day='start')

        current_name1 = None
        current_name2 = None
        df = None
        for comparison in comparisons:
            names = comparison[0]
            days = comparison[1]
            name1 = names[0]
            name2 = names[1]
            day1 = days[0]
            day2 = days[1]
            if current_name1 != name1 or current_name2 != name2:
                if df is not None:
                    df.to_csv('{}_{}.tsv'.format(current_name1, current_name2).replace('/', '-'),
                              sep='\t', header=True, index_label='id')

                df = pd.DataFrame(index=self.features)
                current_name1 = name1
                current_name2 = name2

            logger.info('{} vs {}, day {}, day {}'.format(name1, name2, day1, day2))
            values1, weights1 = self.get_expression_and_weights(day1, name1)
            values2, weights2 = self.get_expression_and_weights(day2, name2)

            df = self.add_stats(values1, weights1, df, '_{}_{}'.format(name1, day1))

            df = self.add_stats(values2, weights2, df, '_{}_{}'.format(name2, day2))

            df = df.join(self.do_comparison(values1, weights1, day1, values2, weights2, day2))

        if df is not None:
            df.to_csv('{}_{}.tsv'.format(current_name1, current_name2).replace('/', '-'),
                      sep='\t', header=True, index_label='id')


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
    expression_file = args.matrix
    delta_days = args.delta
    cell_days_file = args.cell_days
    nperm = args.nperm
    min_fold_change = args.fold_change
    cell_days_field = args.cell_days_field
    expression_matrix = wot.io.read_dataset(expression_file, var_filter=args.gene_filter)
    day_filter = args.cell_day_filter

    if delta_days is None:
        delta_days = 0
    delta_days = abs(delta_days)
    wot.io.add_row_metadata_to_dataset(dataset=expression_matrix, days=cell_days_file)
    if day_filter is not None:
        days = [float(day) for day in day_filter.split(',')] if type(day_filter) == str else day_filter
        expression_matrix = expression_matrix[expression_matrix.obs[cell_days_field].isin(days)]
        expression_matrix = anndata.AnnData(expression_matrix.X, expression_matrix.obs.copy(), expression_matrix.var)
    if expression_matrix.shape[1] is 0:
        sys.exit('Expression matrix has 0 genes')
    fate_names = []
    for f in fate_files:
        fate_ds = wot.io.read_dataset(f)
        if len(fate_files) > 1:
            fate_ds.var.index = fate_ds.var.index + '/' + wot.io.get_filename_and_extension(f)[0]
        expression_matrix.obs = expression_matrix.obs.join(
            pd.DataFrame(index=fate_ds.obs.index, data=fate_ds.X, columns=fate_ds.var.index))
        fate_names += list(fate_ds.var.index)
    d = DiffExp(expression_matrix=expression_matrix, delta_days=delta_days, fate_names=fate_names,
                compare=compare, nperm=nperm, min_fold_change=min_fold_change, day_field=cell_days_field,
                smooth_p_values=args.smooth_p_values)
    d.execute()
