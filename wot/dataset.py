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
        row_duplicates = row_meta.index.duplicated()
        if any(row_duplicates):
            raise Exception('Duplicates in row indices : ',
                    list(row_meta.index[row_duplicates])[:25],
                    '(list may be truncated)')
        col_duplicates = col_meta.index.duplicated()
        if any(col_duplicates):
            raise Exception('Duplicates in column indices : ',
                    list(col_meta.index[col_duplicates])[:25],
                    '(list may be truncated)')
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

    def where(self, **kwargs):
        """
        Get subset of cells that match required metadata

        Parameters
        ----------
        **kwargs : dict
            The dictionnary of  metadata to match

        Returns
        -------
        result : wot.Dataset
            The subset of cells

        Examples
        --------
        >>> ds.where(day=1.2) # -> only cells that have a metadata for 'day' of 1.2
        >>> ds.where(day=1, covariate=3) # -> only cells at day 1 that have a covariate value of 3
        """
        if not kwargs:
            return wot.Dataset(self.x, self.row_meta.copy(), self.col_meta.copy())

        for key in kwargs:
            if key not in self.row_meta.columns:
                raise ValueError("No such metadata in dataset : \"{}\"".format(key))
        query = lambda s : all(s[k] == kwargs[k] for k in kwargs)

        all_ids = self.row_meta.index[self.row_meta.apply(query, axis=1)]
        all_indices = self.row_meta.index.get_indexer_for(all_ids)
        return Dataset(self.x[all_indices], self.row_meta.iloc[all_indices].copy(), self.col_meta.copy())
