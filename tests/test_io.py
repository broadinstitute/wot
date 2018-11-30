import os
import unittest

import anndata
import numpy as np
import pandas as pd

import wot.io


class TestIO(unittest.TestCase):

    def test_read_grp(self):
        gs = wot.io.read_sets(os.path.abspath('inputs/io/test.grp'))
        self.assertTrue(np.sum(gs.X) == 4)

    def test_read_grp_subset(self):
        gs = wot.io.read_sets(os.path.abspath('inputs/io/test.grp'), ['a', 'e', 'c', 'f'])
        self.assertTrue(gs.X[0, 0] == 1 and gs.X[2, 0] == 1)

    def test_read_gmx(self):
        gs = wot.io.read_sets(os.path.abspath(
            'inputs/io/test_gene_sets.gmx'))
        np.testing.assert_array_equal(gs.X,
                                      np.array([[1, 0], [0, 1], [1, 0], [0,
                                                                         1],
                                                [1, 0]]))

    def test_read_gmt(self):
        gs = wot.io.read_sets(
            os.path.abspath('inputs/io/msigdb.v6.1.symbols.gmt'))

        set_id_filter = gs.var.index.values == \
                        'GO_NUCLEOTIDE_TRANSMEMBRANE_TRANSPORT'
        d = gs.X[:, set_id_filter]
        gene_id_filter = np.where(d == 1)
        ids = gs.obs.index.values[gene_id_filter[0]]
        expected_ids = (
            'SLC35B2', 'SLC35E3', 'SLC35D1', 'SLC35B3', 'SLC35A2', 'SLC35D2',
            'SLC25A42', 'SLC35B1', 'SLC35A3', 'SLC35B4', 'SLC25A33', 'SLC25A17')
        self.assertTrue(len(ids) == len(expected_ids))
        for id in expected_ids:
            index = np.where(ids == id)
            self.assertTrue(len(index) == 1)
            self.assertTrue(index[0] >= 0)

    def test_read_gmt_write_gmx(self):
        gs = wot.io.read_sets(
            os.path.abspath('inputs/io/msigdb.v6.1.symbols.gmt'))

        set_id_filter = gs.var.index.values == \
                        'GO_NUCLEOTIDE_TRANSMEMBRANE_TRANSPORT'
        d = gs.X[:, set_id_filter]
        gene_id_filter = np.where(d == 1)
        ids = gs.obs.index.values[gene_id_filter[0]]
        expected_ids = (
            'SLC35B2', 'SLC35E3', 'SLC35D1', 'SLC35B3', 'SLC35A2', 'SLC35D2',
            'SLC25A42', 'SLC35B1', 'SLC35A3', 'SLC35B4', 'SLC25A33', 'SLC25A17')
        self.assertTrue(len(ids) == len(expected_ids))
        for id in expected_ids:
            index = np.where(ids == id)
            self.assertTrue(len(index) == 1)
            self.assertTrue(index[0] >= 0)

    def test_read_gmt_gmx_order(self):
        expected_ids_array = [['a', 'b'], ['c', 'd']]
        inputs = ['inputs/io/test_gene_sets2.gmx', 'inputs/io/test_gene_sets2.gmt'];
        for input in inputs:

            for j in range(len(expected_ids_array)):
                gs = wot.io.read_sets(
                    os.path.abspath(input), feature_ids=expected_ids_array[j])
                gs_ids = gs.obs.index.values
                expected_ids = expected_ids_array[j]

                for i in range(len(expected_ids)):
                    self.assertTrue(expected_ids[i] == gs_ids[i])

    def test_mtx_to_gct(self):
        ds = wot.io.read_dataset(
            'inputs/io/filtered_gene_bc_matrices/hg19/matrix.mtx')

        wot.io.write_dataset(ds, 'test.gct', 'gct')
        ds2 = wot.io.read_dataset('test.gct')
        np.testing.assert_array_equal(ds.X.toarray(),
                                      ds2.X)
        pd.testing.assert_frame_equal(
            ds.obs,
            ds2.obs)
        pd.testing.assert_frame_equal(
            ds.var,
            ds2.var)

    def test_mtx_to_loom(self):
        ds = wot.io.read_dataset(
            'inputs/io/filtered_gene_bc_matrices/hg19/matrix.mtx')
        ds.obs.index.rename('id', inplace=True)
        ds.var.index.rename('id', inplace=True)
        wot.io.write_dataset(ds, 'test.loom', 'loom')
        ds2 = wot.io.read_dataset('test.loom')

        np.testing.assert_array_equal(ds.X.toarray(),
                                      ds2.X.toarray())
        pd.testing.assert_frame_equal(
            ds.obs,
            ds2.obs)
        pd.testing.assert_frame_equal(
            ds.var,
            ds2.var)

    def test_read_mtx(self):
        ds2 = wot.io.read_dataset(
            'inputs/io/filtered_gene_bc_matrices/hg19/matrix.mtx')
        self.assertTrue(ds2.X[2248, 30957] == 10)

        self.assertTrue(ds2.var.index.values[39] == 'ENSG00000272512')
        self.assertTrue(ds2.var['symbol'].values[39] == 'RP11-54O7.17')

    def test_read_loom(self):
        ds = anndata.AnnData(X=np.array([[1, 2, 3, 0],
                                         [4, 5, 6, 0]]),
                             obs=pd.DataFrame(
                                 index=['c1', 'c2']),
                             var=pd.DataFrame(
                                 index=['g1', 'g2',
                                        'g3', 'g4']))
        ds.obs.index.rename('id', inplace=True)
        ds.var.index.rename('id', inplace=True)
        wot.io.write_dataset(ds, 'test.loom', 'loom')
        ds2 = wot.io.read_dataset('test.loom')

        np.testing.assert_array_equal(ds.X,
                                      ds2.X)
        pd.testing.assert_frame_equal(
            ds.obs,
            ds2.obs)
        pd.testing.assert_frame_equal(
            ds.var,
            ds2.var)
