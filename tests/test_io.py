import unittest
import wot.io
import numpy as np
import os
import pandas as pd
import dask.array as da
import dask.dataframe as dd
import dask


class TestIO(unittest.TestCase):

    def test_read_gmx(self):
        gs = wot.io.read_gene_sets(os.path.abspath(
            'inputs/io/test_gene_sets.gmx'))
        np.testing.assert_array_equal(gs.x,
                                      np.array([[1, 0], [0, 1], [1, 0], [0,
                                                                         1],
                                                [1, 0]]))

    def test_read_gmt(self):
        gs = wot.io.read_gene_sets(
            os.path.abspath('inputs/io/msigdb.v6.1.symbols.gmt'))

        set_id_filter = gs.col_meta.index.values == \
                        'GO_NUCLEOTIDE_TRANSMEMBRANE_TRANSPORT'
        d = gs.x[:, set_id_filter]
        gene_id_filter = np.where(d == 1)
        ids = gs.row_meta.index.values[gene_id_filter[0]]
        expected_ids = (
            'SLC35B2', 'SLC35E3', 'SLC35D1', 'SLC35B3', 'SLC35A2', 'SLC35D2',
            'SLC25A42', 'SLC35B1', 'SLC35A3', 'SLC35B4', 'SLC25A33', 'SLC25A17')
        self.assertTrue(len(ids) == len(expected_ids))
        for id in expected_ids:
            index = np.where(ids == id)
            self.assertTrue(len(index) == 1)
            self.assertTrue(index[0] >= 0)

    def test_mtx_to_loom(self):
        ds = wot.io.read_dataset(
            'inputs/io/filtered_gene_bc_matrices/hg19/matrix.mtx')
        ds.row_meta.index.rename('id', inplace=True)
        ds.col_meta.index.rename('id', inplace=True)
        wot.io.write_dataset(ds, 'test.loom', 'loom')
        ds2 = wot.io.read_dataset('test.loom')
        ds2.row_meta.set_index('id', inplace=True)
        ds2.col_meta.set_index('id', inplace=True)
        np.testing.assert_array_equal(ds.x.todense(),
                                      ds2.x)
        pd.testing.assert_frame_equal(
            ds.row_meta,
            ds2.row_meta)
        pd.testing.assert_frame_equal(
            ds.col_meta,
            ds2.col_meta)

    def test_read_gmt_dataset_order(self):
        sets = (('ENSG00000248473', 'ENSG00000272431'), ('ENSG00000181323', 'ENSG00000266824'))
        ds = wot.io.read_dataset(
            'inputs/io/filtered_gene_bc_matrices/hg19/matrix.mtx')
        gs = wot.io.read_gene_sets(
            os.path.abspath('inputs/io/test_gene_sets2.gmx'), feature_ids=ds.col_meta.index.values)
        for j in range(len(sets)):
            expected_ids = sets[j]
            indices = np.isin(gs.row_meta.index.values, expected_ids)
            gs_x = gs.x[indices, [j]]
            gs_row_meta = gs.row_meta.iloc[indices]
            self.assertTrue(np.sum(gs_x) == 2)
            ids = gs_row_meta.index.values

            for id in expected_ids:
                index = np.where(ids == id)
                self.assertTrue(len(index) == 1)
                self.assertTrue(index[0] >= 0)

    def test_read_mtx(self):
        ds2 = wot.io.read_dataset(
            'inputs/io/filtered_gene_bc_matrices/hg19/matrix.mtx')
        ds2.x, ds2.row_meta, ds2.col_meta = dask.compute(ds2.x, ds2.row_meta,
                                                         ds2.col_meta)
        self.assertTrue(ds2.x[2248, 30957] == 10)

        self.assertTrue(ds2.col_meta.index.values[39] == 'ENSG00000272512')
        self.assertTrue(ds2.col_meta['Symbol'].values[39] == 'RP11-54O7.17')

    def test_read_loom(self):
        ds = wot.Dataset(x=np.array([[1, 2, 3, 0],
                                     [4, 5, 6, 0]]),
                         row_meta=pd.DataFrame(
                             index=['c1', 'c2']),
                         col_meta=pd.DataFrame(
                             index=['g1', 'g2',
                                    'g3', 'g4']))
        ds.row_meta.index.rename('id', inplace=True)
        ds.col_meta.index.rename('id', inplace=True)
        wot.io.write_dataset(ds, 'test.loom', 'loom')
        ds2 = wot.io.read_dataset('test.loom', use_dask=True)
        ds2.x, ds2.row_meta, ds2.col_meta = dask.compute(ds2.x, ds2.row_meta,
                                                         ds2.col_meta)
        ds2.row_meta.set_index('id', inplace=True)
        ds2.col_meta.set_index('id', inplace=True)
        np.testing.assert_array_equal(ds.x,
                                      ds2.x)
        pd.testing.assert_frame_equal(
            ds.row_meta,
            ds2.row_meta)
        pd.testing.assert_frame_equal(
            ds.col_meta,
            ds2.col_meta)

    def test_read_gmt_dask(self):
        gs = wot.io.read_gene_sets(
            os.path.abspath('inputs/io/msigdb.v6.1.symbols.gmt'))
        gs = wot.Dataset(x=da.from_array(gs.x, chunks=(1000, 1000)),
                         row_meta=dd.from_pandas(gs.row_meta, sort=False,
                                                 npartitions=4),
                         col_meta=dd.from_pandas(gs.col_meta, npartitions=4,
                                                 sort=False))

        set_id_filter = gs.col_meta.index.values == \
                        'GO_NUCLEOTIDE_TRANSMEMBRANE_TRANSPORT'
        set_id_filter = set_id_filter.compute()
        d = gs.x[:, set_id_filter]
        gene_id_filter = da.where(d == 1)
        ids = gs.row_meta.index.values.compute()[gene_id_filter[0]]
        expected_ids = (
            'SLC35B2', 'SLC35E3', 'SLC35D1', 'SLC35B3', 'SLC35A2', 'SLC35D2',
            'SLC25A42', 'SLC35B1', 'SLC35A3', 'SLC35B4', 'SLC25A33', 'SLC25A17')
        self.assertTrue(len(ids) == len(expected_ids))
        for id in expected_ids:
            index = np.where(ids == id)
            self.assertTrue(len(index) == 1)
            self.assertTrue(index[0] >= 0)
