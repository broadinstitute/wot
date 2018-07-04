# -*- coding: utf-8 -*-
import pandas as pd


class Dataset:
    """
       A Dataset consists of a 2-d matrix x, row metadata, and column metadata.
       Args:
           x (ndarray|scipy.spmatrix): 2-d matrix
           row_meta (pd.DataFrame) Row metadata
           col_meta (pd.DataFrame) Column metadata
       """

    def __init__(self, x, row_meta, col_meta):
        self.x = x
        self.row_meta = row_meta
        self.col_meta = col_meta
        self.layers = {}
        if type(row_meta) == pd.DataFrame and x.shape[0] != row_meta.shape[0]:
            raise Exception('Row dimensions do not match: ' + str(x.shape[0]) +
                            '!=' + str(row_meta.shape[0]))
        if type(col_meta) == pd.DataFrame and x.shape[1] != col_meta.shape[0]:
            raise Exception(
                'Column dimensions do not match: ' + str(x.shape[1]) +
                '!=' + str(col_meta.shape[0]))

    def transpose(self):
        self.x = self.x.transpose()
        tmp = self.row_meta
        self.row_meta = self.col_meta
        self.col_meta = tmp

    def __len__(self):
        return len(self.x)
