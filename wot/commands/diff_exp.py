#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import numpy as np
import pandas as pd
import scipy.sparse
import wot.io
from statsmodels.stats.multitest import multipletests


def main(argv):
    parser = argparse.ArgumentParser(
        description='Compute differentially expressed genes from the output of the trajectory tool')
    parser.add_argument('--matrix', help=wot.commands.MATRIX_HELP, required=True)
    parser.add_argument('--trajectory', help='One or more trajectory datasets as produced by the trajectory tool',
                        action='append')
    parser.add_argument('--out', help='Prefix for output file names', default='enrichment')
    parser.add_argument('--cell_days', help=wot.commands.CELL_DAYS_HELP, required=True)
    parser.add_argument('--ncells', help='Number of cells to sample at each day', type=int, default=4000)
    parser.add_argument('--between', help='Compare across trajectories when more than one trajectory is supplied',
                        action='store_true')
    parser.add_argument('--delta',
                        help='Delta days to compare sampled expression matrix against within a trajectory. If not specified all comparison are done against the first day.',
                        type=float)

    args = parser.parse_args(argv)
    ncells = args.ncells
    trajectory_files = args.trajectory
    expression_file = args.matrix
    delta_days = args.delta
    cell_days_file = args.cell_days

    expression_matrix = wot.io.read_dataset(expression_file)

    if delta_days is None:
        delta_days = 0
    delta_days = abs(delta_days)
    wot.io.add_row_metadata_to_dataset(dataset=expression_matrix, days_path=cell_days_file)
    trajectory_names = []
    for f in trajectory_files:
        trajectory_ds = wot.io.read_dataset(f)
        if len(trajectory_files) > 1:
            trajectory_ds.var.index = trajectory_ds.var.index + wot.io.get_filename_and_extension(f)[0]
        expression_matrix.obs = expression_matrix.obs.join(
            pd.DataFrame(index=trajectory_ds.obs.index, data=trajectory_ds.X, columns=trajectory_ds.var.index))
        trajectory_names += list(trajectory_ds.var.index)

    def get_sampled_expression(day, trajectory_name, invert=False):
        ds = expression_matrix[
            (expression_matrix.obs['day'] == day) & (False == expression_matrix.obs[trajectory_name].isna())]
        weights = ds.obs[trajectory_name].values
        expression_values = ds.X
        if scipy.sparse.isspmatrix(expression_values):
            expression_values = expression_values.toarray()
        if invert:
            weights = 1 - weights
        weights = weights / weights.sum()
        sampled_indices = np.random.choice(np.arange(ds.shape[0]), size=ncells, replace=True, p=weights)
        sampled_expression = expression_values[sampled_indices]
        return sampled_expression

    def add_stats(values, df, suffix):
        if 'fraction_expressed{}'.format(suffix) not in df:
            return df.join(pd.DataFrame(index=features,
                                        data={'fraction_expressed{}'.format(suffix): (
                                                (values > 0).sum(axis=0) / values.shape[0]),
                                            'mean{}'.format(suffix): values.mean(axis=0),
                                            # 'min{}'.format(suffix): values.min(axis=0),
                                            # 'max{}'.format(suffix): values.max(axis=0)
                                        }))
        return df

    def do_comparison(values, compare_to_values, day, compare_to_day):

        # adata = anndata.AnnData(X=values1, var=pd.DataFrame(index=features)).concatenate(
        #     anndata.AnnData(X=values2, var=pd.DataFrame(index=features)))
        # # method 'logreg', 't-test', 'wilcoxon', 't-test_overestim_var'
        # # sc.tl.rank_genes_groups(adata, 'batch', method='t-test_overestim_var', n_genes=len(features), groups=['0'],
        # #                         reference='1')
        # sc.tl.rank_genes_groups(adata, 'batch', method='logreg')
        # sorted_names = adata.uns['rank_genes_groups']['names']
        # scores = adata.uns['rank_genes_groups']['scores']
        # logfoldchanges = adata.uns['rank_genes_groups']['logfoldchanges']
        # pvals = adata.uns['rank_genes_groups']['pvals']
        # pvals_adj = adata.uns['rank_genes_groups']['pvals_adj']
        # sc.pl.rank_genes_groups_violin(adata, n_genes=10, save='{}_{}'.format(day, compare_to_day), show=False)
        # # sc.pl.rank_genes_groups_violin(adata, groups='0', n_genes=10, save=True)

        # scores_wilcoxon = np.zeros(expression_matrix.shape[1])
        # pvals_ks = np.zeros(values.shape[1])
        # scores_ks = np.zeros(values.shape[1])
        scores_t, pvals_t = scipy.stats.ttest_ind(values, compare_to_values, equal_var=False, axis=0)
        # pvals_t[np.isnan(pvals_t)] = 1
        # for i in range(values.shape[1]):  # each gene
        #     # Wilcoxon rank-sum statistic
        #     # scores_wilcoxon[i], pvals_wilcoxon[i] = scipy.stats.ranksums(values[:, i], y=compare_to_values[:, i])
        #     scores_ks[i], pvals_ks[i] = scipy.stats.ks_2samp(values[:, i], compare_to_values[:, i])

        mean1 = values.mean(axis=0)
        mean2 = compare_to_values.mean(axis=0)
        mean1[mean1 == 0] = 1e-9  # set 0s to small value
        mean2[mean2 == 0] = 1e-9  # set 0s to small value
        foldchanges = mean1 / mean2
        mean_difference = mean1 - mean2
        fraction_expressed_difference = ((values > 0).sum(axis=0) / values.shape[0]) - (
                (compare_to_values > 0).sum(axis=0) / compare_to_values.shape[0])
        return pd.DataFrame(index=features,
                            data={
                                'p_value_ttest_day_{}_vs_day_{}'.format(day, compare_to_day): pvals_t,
                                'fdr_ttest_day_{}_vs_day_{}'.format(day, compare_to_day):
                                    multipletests(pvals_t, alpha=0.1, method='fdr_bh')[1],
                                'ttest_{}_vs_day_{}'.format(day, compare_to_day): scores_t,
                                'fold_change_{}_vs_day_{}'.format(day, compare_to_day): foldchanges,
                                'mean_difference_{}_vs_day_{}'.format(day, compare_to_day): mean_difference,
                                'fraction_expressed_difference_{}_vs_day_{}'.format(day,
                                                                                    compare_to_day): fraction_expressed_difference
                            })

    features = expression_matrix.var.index
    days = np.array(sorted(expression_matrix.obs['day'].unique().astype(float)))
    days = days[np.isnan(days) == False]

    if args.between and len(trajectory_names) > 1:
        for trajectory_index1 in range(len(trajectory_names)):
            for trajectory_index2 in range(trajectory_index1):
                df = pd.DataFrame(index=features)
                for day_index in range(len(days)):
                    day = days[day_index]
                    values1 = get_sampled_expression(day, trajectory_names[trajectory_index1])
                    values2 = get_sampled_expression(day, trajectory_names[trajectory_index2])
                    df = add_stats(values1, df, '_{}_day_{}'.format(trajectory_names[trajectory_index1], day))
                    df = add_stats(values2, df, '_{}_day_{}'.format(trajectory_names[trajectory_index2], day))
                    df = df.join(do_comparison(values2, values1, day, day))

                    print('{} vs {}, day {}'.format(trajectory_names[trajectory_index2],
                                                    trajectory_names[trajectory_index1], day))
                df.to_csv('{}_{}.tsv'.format(trajectory_names[trajectory_index2], trajectory_names[trajectory_index1]),
                          sep='\t', header=True)

    else:
        # within
        for name in trajectory_names:
            day_zero_values = get_sampled_expression(days[0], name)  # (cells, genes)
            df = pd.DataFrame(index=features)

            for day_index in range(1, len(days)):
                day = days[day_index]
                compare_to_day = days[0]
                compare_to_values = day_zero_values
                if delta_days > 0:
                    requested_day = day - delta_days
                    compare_to_day = days[np.abs(days - requested_day).argmin()]
                    if compare_to_day == day or np.abs(
                            compare_to_day - day - delta_days) > 0.1:  # too big or small a gap
                        continue
                    compare_to_values = get_sampled_expression(compare_to_day, name)
                print('{}, day {} vs day {}'.format(name, day, compare_to_day))
                sampled_values = get_sampled_expression(day, name)
                df = df.join(do_comparison(sampled_values, compare_to_values, day, compare_to_day))
                df = add_stats(sampled_values, df, '_day_{}'.format(day))
                df = add_stats(compare_to_values, df, '_day_{}'.format(compare_to_day))
            df.to_csv(name + '.tsv', sep='\t', header=True)
