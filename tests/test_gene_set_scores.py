#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import subprocess
import unittest

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
            output = pd.read_table(output_file, index_col=0)
            np.testing.assert_array_equal(output[set_names[i]].values, scores[:, i])

            os.remove(output_file)

    def test_score_gene_sets_drop(self):
        ds = wot.Dataset(x=np.array([[1.0, 2, 3, 4, 5, 6, 7, 8, 9, 10]]), row_meta=None,
                         col_meta=None)

        gs = wot.Dataset(x=np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 1]], dtype=np.uint8).T, row_meta=None, col_meta=None)
        result = wot.score_gene_sets(dataset_to_score=ds, gs=gs, method=None, permutations=100, nbins=1,
                                     random_state=1234, drop_frequency=100, drop_p_value_threshold=1,
                                     smooth_p_values=False)
        self.assertTrue(result['p_value'][0] == 8.0 / 100.0)

    def test_p_value1(self):
        ds = wot.Dataset(x=np.array([[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]]), row_meta=None, col_meta=None)
        gs = wot.Dataset(x=np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0]], dtype=np.uint8).T, row_meta=None, col_meta=None)
        result = wot.score_gene_sets(dataset_to_score=ds, gs=gs, method=None, permutations=10, nbins=1,
                                     drop_frequency=0)
        np.testing.assert_array_equal(result['p_value'][0], 11.0 / 12.0)

    def test_score_gene_sets_no_drop(self):
        ds = wot.Dataset(x=np.array([[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]]), row_meta=None, col_meta=None)

        gs = wot.Dataset(x=np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 1]], dtype=np.uint8).T, row_meta=None, col_meta=None)
        result = wot.score_gene_sets(dataset_to_score=ds, gs=gs, method=None, permutations=1000, nbins=1,
                                     random_state=1234, drop_frequency=0)
        self.assertTrue(result['p_value'][0] == 108.0 / 1002.0)

    def test_score_gene_sets_sparse_matrix(self):
        ds = wot.Dataset(
            x=scipy.sparse.csr_matrix(np.array([[0, 0, 0, 4, 5, 6, 7, 8, 9, 10]])),
            row_meta=None,
            col_meta=None)

        gs = wot.Dataset(x=np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 1]], dtype=np.uint8).T, row_meta=None, col_meta=None)
        result = wot.score_gene_sets(dataset_to_score=ds, gs=gs, method=None, permutations=100, nbins=1,
                                     random_state=1234, drop_frequency=100, drop_p_value_threshold=1,
                                     smooth_p_values=False)

        self.assertTrue(result['p_value'][0] == 8.0 / 100.0)

    def test_score_gene_sets_basic(self):
        ds = wot.Dataset(x=np.array([[1.0, 2.0, 3, 0],
                                     [4, 5, 6.0, 0]]),
                         row_meta=pd.DataFrame(
                             index=['c1', 'c2']),
                         col_meta=pd.DataFrame(
                             index=['g1', 'g2',
                                    'g3', 'g4']))

        gs = wot.Dataset(x=np.array([[1, 0, 1],
                                     [0, 0, 1],
                                     [0, 0, 0],
                                     [0, 1, 0]
                                     ],
                                    dtype=np.uint8),
                         row_meta=pd.DataFrame(
                             index=['g1', 'g2', 'g3', 'g4']),
                         col_meta=pd.DataFrame(
                             index=['s1', 's2', 's3']))

        expected = np.array([[1, 0, 1.5], [4, 0, 4.5]])
        scores = []
        for j in range(gs.x.shape[1]):
            result = wot.score_gene_sets(dataset_to_score=ds,
                                         gs=wot.Dataset(gs.x[:, [j]], gs.row_meta, gs.col_meta.iloc[[j]]), method=None,
                                         permutations=10, nbins=1)
            scores.append(result['score'])
        np.testing.assert_array_equal(np.hstack(scores), expected)


if __name__ == '__main__':
    unittest.main()
