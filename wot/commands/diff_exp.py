#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import numpy as np
import pandas as pd
import scipy.sparse
import statsmodels.stats.multitest

import wot.io


class DiffExp:

    def __init__(self, expression_matrix, delta_days, trajectory_names,
                 compare, nperm, min_fold_change):
        self.expression_matrix = expression_matrix
        self.delta_days = delta_days
        self.trajectory_names = trajectory_names
        if len(trajectory_names) == 1:
            compare = 'within'
        self.compare = compare
        self.nperm = nperm
        self.min_fold_change = min_fold_change
        self.features = expression_matrix.var.index
        days = np.array(sorted(expression_matrix.obs['day'].unique().astype(float)))
        self.days = days[np.isnan(days) == False]

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

    def get_expression_and_weights(self, day, trajectory_name):
        ds = self.expression_matrix[
            (self.expression_matrix.obs['day'] == day) & (False == self.expression_matrix.obs[trajectory_name].isna())]
        weights = ds.obs[trajectory_name].values
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
        suffix = "_{}_{}".format(day1, day2) if day1 != day2 else '_{}'.format(day1)
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
                p = (p + 1) / (self.nperm + 2)
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
        if self.compare != 'within':
            if self.compare == 'all':
                base_trajectory_names_to_trajectory_names = {'': self.trajectory_names}
            elif self.compare == 'match':
                base_trajectory_names_to_trajectory_names = {}
                for i in range(len(self.trajectory_names)):
                    full_trajectory_name = self.trajectory_names[i]
                    base_trajectory_name = full_trajectory_name[0:full_trajectory_name.rindex('/')]
                    names = base_trajectory_names_to_trajectory_names.get(base_trajectory_name)
                    if names is None:
                        names = []
                        base_trajectory_names_to_trajectory_names[base_trajectory_name] = names
                    names.append(full_trajectory_name)

            for base_trajectory_name in base_trajectory_names_to_trajectory_names:
                trajectory_names = base_trajectory_names_to_trajectory_names[base_trajectory_name]
                for i in range(len(trajectory_names)):
                    for j in range(i):
                        df = pd.DataFrame(index=self.features)

                        for day_index in range(len(self.days)):
                            day = self.days[day_index]
                            print('{} vs {}, day {}'.format(trajectory_names[j], trajectory_names[i], day))
                            values1, weights1 = self.get_expression_and_weights(day, trajectory_names[j])
                            values2, weights2 = self.get_expression_and_weights(day, trajectory_names[i])

                            df = self.add_stats(values1, weights1, df, '_{}_{}'.format(trajectory_names[j], day))
                            df = self.add_stats(values2, weights2, df, '_{}_{}'.format(trajectory_names[i], day))
                            df = df.join(self.do_comparison(values1, weights1, day, values2, weights2, day))

                        df.to_csv('{}_{}.tsv'.format(trajectory_names[j], trajectory_names[i]).replace('/', '-'),
                                  sep='\t', header=True)



        else:
            # within

            for name in self.trajectory_names:
                df = pd.DataFrame(index=self.features)
                for day_index in range(1, len(self.days)):
                    day2 = self.days[day_index]
                    if self.delta_days > 0:
                        requested_day = day2 - self.delta_days
                        day1 = self.days[np.abs(self.days - requested_day).argmin()]
                        if day1 == day2 or np.abs(
                                day1 - day2 - self.delta_days) > 0.1:  # too big or small a gap
                            continue
                    else:
                        day1 = self.days[0]
                    print('{}, day {} vs day {}'.format(name, day1, day2))
                    values1, weights1 = self.get_expression_and_weights(day1, name)
                    values2, weights2 = self.get_expression_and_weights(day2, name)
                    df = self.add_stats(values1, weights1, df, '_{}'.format(day1))
                    df = self.add_stats(values2, weights2, df, '_{}'.format(day2))
                    df = df.join(self.do_comparison(values1, weights1, day1, values2, weights2, day2))

                df.to_csv(name + '.tsv', sep='\t', header=True)


def main(argv):
    parser = argparse.ArgumentParser(
        description='Compute differentially expressed genes from the output of the trajectory tool')
    parser.add_argument('--matrix', help=wot.commands.MATRIX_HELP, required=True)
    parser.add_argument('--trajectory', help='One or more trajectory datasets as produced by the trajectory tool',
                        action='append')
    parser.add_argument('--cell_days', help=wot.commands.CELL_DAYS_HELP, required=True)
    parser.add_argument('--compare',
                        help='Compare across trajectories when more than one trajectory is supplied. If "match" compare trajectories with the same name. If "all", compare all pairs. If "within" compare within a trajectory.',
                        choices=['within', 'match', 'all'], default='within')
    parser.add_argument('--delta',
                        help='Delta days to compare sampled expression matrix against within a trajectory. If not specified all comparison are done against the first day.',
                        type=float)
    parser.add_argument('--nperm',
                        help='Number of permutations', type=int)
    parser.add_argument('--fold_change', type=float, default=0.25,
                        help='Limit permutations to genes which show at least X-fold difference (log-scale) between the two groups of cells.')

    args = parser.parse_args(argv)
    compare = args.compare
    trajectory_files = args.trajectory
    expression_file = args.matrix
    delta_days = args.delta
    cell_days_file = args.cell_days
    nperm = args.nperm
    min_fold_change = args.fold_change
    expression_matrix = wot.io.read_dataset(expression_file)

    if delta_days is None:
        delta_days = 0
    delta_days = abs(delta_days)
    wot.io.add_row_metadata_to_dataset(dataset=expression_matrix, days_path=cell_days_file)
    trajectory_names = []
    for f in trajectory_files:
        trajectory_ds = wot.io.read_dataset(f)
        if len(trajectory_files) > 1:
            trajectory_ds.var.index = trajectory_ds.var.index + '/' + wot.io.get_filename_and_extension(f)[0]
        expression_matrix.obs = expression_matrix.obs.join(
            pd.DataFrame(index=trajectory_ds.obs.index, data=trajectory_ds.X, columns=trajectory_ds.var.index))
        trajectory_names += list(trajectory_ds.var.index)
    d = DiffExp(expression_matrix=expression_matrix, delta_days=delta_days, trajectory_names=trajectory_names,
                compare=compare, nperm=nperm, min_fold_change=min_fold_change)
    d.execute()
