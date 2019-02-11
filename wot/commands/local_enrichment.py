#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import subprocess

import anndata
import numpy as np
import pandas as pd
import pkg_resources

import wot.io


def get_scores(ds1, ds2, ds1_time_index, ds2_time_index, score_function):
    scores = np.zeros(shape=ds1.X.shape[1])
    for feature_idx in range(ds1.X.shape[1]):  # each feature
        m1 = ds1.X[ds1_time_index, feature_idx]
        m2 = ds2.X[ds2_time_index, feature_idx]
        v1 = ds1.variance[ds1_time_index, feature_idx]
        v2 = ds2.variance[ds2_time_index, feature_idx]

        n1 = ds1.obs.iloc[ds1_time_index]
        n2 = ds2.obs.iloc[ds2_time_index]
        scores[feature_idx] = score_function(m1, m2, v1, v2, n1, n1)
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
                        choices=['s2n', 'mean_difference', 'fold_change', 't_test'])
    parser.add_argument('--gsea',
                        help='Run GSEA on the specified MSigDB collections (http://software.broadinstitute.org/gsea/msigdb/collections.jsp). H (hallmark gene sets), C1 (positional gene sets), C2 (curated gene sets), C3 (motif gene sets), C4 (computational gene sets), C5 (GO gene sets), C6 (oncogenic signatures), C7 (immunologic signatures)',
                        action='append')
    parser.add_argument('--comparisons',
                        help='Comparisons to generate ranked lists for. By default, for one matrix signatures are created for all consecutive timepoints. For two matrices for all matching timepoints.')

    args = parser.parse_args(argv)
    # dataset has time on rows, genes on columns
    ds1 = wot.io.read_dataset(args.matrix1)
    ds1.variance = wot.io.read_dataset(args.variance1).X
    ds2 = ds1
    if args.matrix2 is not None:
        ds2 = wot.io.read_dataset(args.matrix2)
        ds2.variance = wot.io.read_dataset(args.variance2).X
        ds_indices = ds1.var.index.get_indexer_for(ds2.var.index)

        if (ds_indices[ds_indices == -1]).sum() > 0:
            raise ValueError('Unable to align datasets')
        ds2 = anndata.AnnData(ds2.X[:, ds_indices], obs=ds2.obs, var=ds2.var.iloc[ds_indices])
        ds2.variance = ds2.variance[:, ds_indices]

    def s2n(m1, m2, v1, v2, *varargs):
        denom = (np.sqrt(v1) + np.sqrt(v2))
        s2n = 0.0 if denom == 0.0 else ((m1 - m2) / denom)
        return s2n

    def fold_change(m1, m2, *varargs):
        return m1 / m2

    def mean_difference(m1, m2, *varargs):
        return m1 - m2

    def t_test(m1, m2, v1, v2, n1, n2):
        return (m1 - m2) / np.sqrt((v1 / n1) + (v2 / n2))

    score_function = locals()[args.score]
    names = []
    if args.comparisons is not None:
        comparisons = pd.read_csv(args.comparisons, header=None, index_col=False, engine='python', sep=None)
        for comparison_idx in range(comparisons.shape[0]):
            i = np.where(ds1.obs.index.values == comparisons.iloc[comparison_idx, 0])[0][0]
            j = np.where(ds2.obs.index.values == comparisons.iloc[comparison_idx, 1])[0][0]
            name = str(ds1.obs.index.values[i]) + '_' + str(ds2.obs.index.values[j])
            scores = get_scores(ds1, ds2, i, j, score_function)
            names.append(name)
            pd.DataFrame(index=ds1.var.index.str.upper().values, data={'scores': scores}).to_csv(
                name + '.rnk', sep='\t', header=False)
    else:
        if ds2 is not None:
            for i in range(ds1.X.shape[0]):  # each time
                scores = get_scores(ds1, ds2, i, i, score_function)
                name = str(ds1.obs.index.values[i])
                names.append(name)
                pd.DataFrame(index=ds1.var.index.str.upper().values, data={'scores': scores}).to_csv(
                    name + '.rnk', sep='\t', header=False)
        else:
            for i in range(1, ds1.X.shape[0]):
                scores = get_scores(ds1, ds1, i - 1, i, score_function)
                name = str(ds1.obs.index.values[i - 1]) + '_' + str(ds1.obs.index.values[i])
                names.append(name)
                pd.DataFrame(index=ds1.var.index.str.upper().values, data={'scores': scores}).to_csv(
                    name + '.rnk', sep='\t',
                    header=False)

    if args.gsea is not None and len(args.gsea) > 0:
        urls = []
        for c in args.gsea:
            if c.find('.') == -1:
                c += '.all'
            urls.append('gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/' + c.lower() + '.v6.1.symbols.gmt')
        gmx = ','.join(urls)
        for name in names:
            classpath = pkg_resources.resource_filename('wot', 'commands/resources/gsea-3.0.jar')
            gsea_args = ['java', '-Djava.awt.headless=true', '-cp', classpath, '-Xmx4g', 'xtools.gsea.GseaPreranked',
                         '-gmx', gmx,
                         '-norm', 'meandiv', '-nperm', '1000', '-rnk', name + '.rnk', '-scoring_scheme', 'weighted',
                         '-rpt_label', name, '-create_svgs', 'false',
                         '-make_sets', 'true', '-plot_top_x', '50', '-rnd_seed', 'timestamp', '-set_max', '500',
                         '-set_min', '15', '-out',
                         name, '-gui', 'false']
            subprocess.check_call(gsea_args)
