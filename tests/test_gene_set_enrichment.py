#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
import numpy as np
import pandas as pd
import subprocess
import os
import wot


class TestGeneSetEnrichment(unittest.TestCase):

    def test_score_gene_set_command(self):
        subprocess.call(args=['python', os.path.abspath(
            '../bin/gene_set.py'),
                              '--matrix',
                              os.path.abspath(
                                  'inputs/score_gene_sets/matrix.txt'),
                              '--gene_sets', os.path.abspath(
                'inputs/score_gene_sets/gene_sets.gmx'),
                              '--prefix', 'gene_set_test_output',
                              '--no_zscore'],
                        cwd=os.getcwd(),
                        stderr=subprocess.STDOUT)
        output_file = 'gene_set_test_output.txt'
        output = pd.read_table(output_file, index_col=0)
        np.testing.assert_array_equal(output.values,
                                      np.array([[1, 0.5, 1.5], [4, 2,
                                                                4.5]]))

        os.remove(output_file)

    def test_score_gene_sets(self):
        ds = wot.Dataset(x=np.array([[1, 2, 3, 0],
                                     [4, 5, 6, 0]]),
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
        result = wot.score_gene_sets(ds=ds, gs=gs, z_score_ds=False)
        np.testing.assert_array_equal(result.x,
                                      np.array([[1, 0, 1.5], [4, 0,
                                                              4.5]]))


if __name__ == '__main__':
    unittest.main()
