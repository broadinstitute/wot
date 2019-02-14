#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import subprocess
import unittest

import anndata
import numpy as np
import pandas as pd
import scipy.sparse

import wot


class TestGeneSetScores(unittest.TestCase):

    def test_score_gene_set_command(self):
        subprocess.call(args=['wot', 'gene_set_scores',
                              '--matrix',
                              os.path.abspath(
                                  'inputs/score_gene_sets/matrix.txt'),
                              '--gene_sets', os.path.abspath(
                'inputs/score_gene_sets/gene_sets.gmx'),
                              '--out', 'test_gene_set_test_output',
                              '--method', 'mean', '--format', 'txt'],
                        cwd=os.getcwd(),
                        stderr=subprocess.STDOUT)
        set_names = ['s1', 's2', 's3']
        scores = np.array([[1, 0, 1.5], [4, 0, 4.5]])
        for i in range(len(set_names)):
            output_file = 'test_gene_set_test_output_' + set_names[i] + '.txt'
            output = pd.read_csv(output_file, index_col=0, sep='\t')
            np.testing.assert_array_equal(output[set_names[i]].values, scores[:, i])
            os.remove(output_file)

    def test_p_value1(self):
        ds = anndata.AnnData(X=np.array([[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]]), obs=None, var=None)
        gs = anndata.AnnData(X=np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0]], dtype=np.uint8).T, obs=None, var=None)
        result = wot.score_gene_sets(ds=ds, gs=gs, method=None, permutations=10,
                                     random_state=1234)
        np.testing.assert_array_equal(result['p_value'][0], 11.0 / 12.0)

    # def test_p_value2(self):
    #     ds = anndata.AnnData(X=np.array([[1, 2, 3], [4, 5, 6]]), obs=None, var=None)
    #     gs = anndata.AnnData(X=np.array([[1, 1, 1]], dtype=np.uint8).T, obs=None, var=None)
    #     result = wot.score_gene_sets(ds=ds, gs=gs, method=None, permutations=100, smooth_p_values=False,
    #                                  random_state=1234)

    def test_score_gene_sets_sparse_ds(self):
        ds = anndata.AnnData(
            X=scipy.sparse.csr_matrix(np.array([[0, 0, 0, 4, 5, 6, 7, 8, 9, 10]])),
            obs=None,
            var=None)

        gs = anndata.AnnData(X=np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 1]], dtype=np.uint8).T, obs=None, var=None)
        result = wot.score_gene_sets(ds=ds, gs=gs, method=None, permutations=100,
                                     random_state=1234,
                                     smooth_p_values=False)
        self.assertEqual(result['k'][0], 100)

    def test_score_gene_sets_sparse_ds_zscore(self):
        ds = anndata.AnnData(
            X=scipy.sparse.csr_matrix(np.array([[0, 0, 0, 4, 5, 6, 7, 8, 9, 10], [1, 2, 3, 5, 6, 7, 9, 9, 19, 11]])),
            obs=None,
            var=None)

        gs = anndata.AnnData(X=np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 1]], dtype=np.uint8).T, obs=None, var=None)
        result = wot.score_gene_sets(ds=ds, gs=gs, method='mean_z_score', permutations=100,
                                     random_state=1234,
                                     smooth_p_values=False)

        self.assertEqual(result['k'][0], 100)

    def test_score_gene_sets_sparse_gs(self):
        ds = anndata.AnnData(
            X=np.array([[0, 0, 0, 4, 5, 6, 7, 8, 9, 10]]),
            obs=None,
            var=None)

        gs = anndata.AnnData(X=scipy.sparse.csr_matrix(np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 1]], dtype=np.uint8).T),
                             obs=None, var=None)
        result = wot.score_gene_sets(ds=ds, gs=gs, method=None,
                                     permutations=100,
                                     random_state=1234,
                                     smooth_p_values=False)

        self.assertEqual(result['k'][0], 100)

    def test_score_gene_sets_basic(self):
        ds = anndata.AnnData(X=np.array([[1.0, 2.0, 3, 0], [4, 5, 6.0, 0]]),
                             obs=pd.DataFrame(
                                 index=['c1', 'c2']),
                             var=pd.DataFrame(
                                 index=['g1', 'g2',
                                        'g3', 'g4']))

        gs = anndata.AnnData(X=np.array([[1, 0, 1], [0, 0, 1], [0, 0, 0], [0, 1, 0]], dtype=np.uint8),
                             obs=pd.DataFrame(
                                 index=['g1', 'g2', 'g3', 'g4']),
                             var=pd.DataFrame(
                                 index=['s1', 's2', 's3']))

        expected = np.array([[1, 0, 1.5], [4, 0, 4.5]])

        for j in range(gs.X.shape[1]):
            result = wot.score_gene_sets(ds=ds,
                                         gs=anndata.AnnData(gs.X[:, [j]], gs.obs, gs.var.iloc[[j]]), method=None,
                                         permutations=10)
            np.testing.assert_array_equal(result['score'], expected[:, j])


if __name__ == '__main__':
    unittest.main()
