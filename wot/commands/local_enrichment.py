#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import numpy as np
import pandas as pd

import wot.io


def get_scores(ds1, ds2, ds1_time_index, ds2_time_index, score_function):
    scores = np.zeros(shape=ds1.shape[1])
    for feature_idx in range(ds1.shape[1]):  # each feature
        m1 = ds1[ds1_time_index, feature_idx].X
        m2 = ds2[ds2_time_index, feature_idx].X
        v1 = ds1.variance[ds1_time_index, feature_idx]
        v2 = ds2.variance[ds2_time_index, feature_idx]
        scores[feature_idx] = score_function(m1, m2, v1, v2)
    return scores


def main(argv):
    parser = argparse.ArgumentParser(
        description='Compute differentially expressed genes from output of trajectory_trends. Outputs a ranked list for each comparison.')
    parser.add_argument('--matrix1', help='Gene expression matrix with timepoints on rows, features on columns',
                        required=True)
    parser.add_argument('--matrix2', help='Gene expression matrix with timepoints on rows, features on columns')
    parser.add_argument('--variance1', help='Variance matrix with timepoints on rows, features on columns',
                        required=True)
    parser.add_argument('--variance2', help='Variance matrix with timepoints on rows, features on columns')
    parser.add_argument('--score',
                        help='Method to compute differential gene expression score. Choices are signal to noise, mean difference, t-test, and fold change',
                        choices=['s2n', 'mean_difference', 'fold_change'], default='s2n')
    parser.add_argument('--comparisons',
                        help='Comparisons to generate ranked lists for. By default, for one matrix signatures are created with respect to the first timepoint. For two matrices for all matching timepoints.')

    args = parser.parse_args(argv)
    # dataset has time on rows, genes on columns
    ds1 = wot.io.read_dataset(args.matrix1)
    ds1.variance = wot.io.read_dataset(args.variance1).X
    ds2 = None
    if args.matrix2 is not None:
        ds2 = wot.io.read_dataset(args.matrix2)
        ds2_variance = wot.io.read_dataset(args.variance2).X
        ds_indices = ds1.var.index.get_indexer_for(ds2.var.index)

        if (ds_indices[ds_indices == -1]).sum() > 0:
            raise ValueError('Unable to align datasets')
        ds2 = ds2[:, ds_indices]
        ds2.variance = ds2_variance[:, ds_indices]

    def s2n(m1, m2, v1, v2, *varargs):
        denom = (np.sqrt(v1) + np.sqrt(v2))
        s2n = 0.0 if denom == 0.0 else ((m1 - m2) / denom)
        return s2n

    def fold_change(m1, m2, *varargs):
        return m1 / m2

    def mean_difference(m1, m2, *varargs):
        return m1 - m2

    score_function = locals()[args.score]
    if args.comparisons is not None:
        comparisons = pd.read_csv(args.comparisons, header=None, index_col=False, engine='python', sep=None)
        for comparison_idx in range(comparisons.shape[0]):
            i = np.where(ds1.obs.index.values == comparisons.iloc[comparison_idx, 0])[0][0]
            j = np.where(ds2.obs.index.values == comparisons.iloc[comparison_idx, 1])[0][0]
            name = str(ds1.obs.index.values[i]) + '_' + str(ds2.obs.index.values[j])
            scores = get_scores(ds1, ds2, i, j, score_function)
            pd.DataFrame(index=ds1.var.index.str.upper().values, data={'scores': scores}).to_csv(
                name + '.rnk', sep='\t', header=False)
    else:
        if ds2 is not None:
            for i in range(ds1.shape[0]):  # each time
                scores = get_scores(ds1, ds2, i, i, score_function)
                name = str(ds1.obs.index.values[i])
                pd.DataFrame(index=ds1.var.index.str.upper().values, data={'scores': scores}).to_csv(
                    name + '.rnk', sep='\t', header=False)
        else:
            reference_index = 0
            for i in range(1, ds1.shape[0]):
                scores = get_scores(ds1, ds1, i, reference_index, score_function)
                name = str(ds1.obs.index.values[i]) + '_' + str(ds1.obs.index.values[reference_index])
                pd.DataFrame(index=ds1.var.index.str.upper().values, data={'scores': scores}).to_csv(
                    name + '.rnk', sep='\t',
                    header=False)

    # if args.gsea is not None and len(args.gsea) > 0:
    #     urls = []
    #     for c in args.gsea:
    #         if c.find('.') == -1:
    #             c += '.all'
    #         urls.append('gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/' + c.lower() + '.v6.1.symbols.gmt')
    #     gmx = ','.join(urls)
    #     for name in names:
    #         classpath = pkg_resources.resource_filename('wot', 'commands/resources/gsea-3.0.jar')
    #         gsea_args = ['java', '-Djava.awt.headless=true', '-cp', classpath, '-Xmx4g', 'xtools.gsea.GseaPreranked',
    #                      '-gmx', gmx,
    #                      '-norm', 'meandiv', '-nperm', '1000', '-rnk', name + '.rnk', '-scoring_scheme', 'weighted',
    #                      '-rpt_label', name, '-create_svgs', 'false',
    #                      '-make_sets', 'true', '-plot_top_x', '50', '-rnd_seed', 'timestamp', '-set_max', '500',
    #                      '-set_min', '15', '-out',
    #                      name, '-gui', 'false']
    #         subprocess.check_call(gsea_args)
