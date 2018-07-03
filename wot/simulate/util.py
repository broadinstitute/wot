# -*- coding: utf-8 -*-

import wot.io

def load_dataset(matrix_file, days=None, covariate=None, gene_sets=None):
    ds = wot.io.read_dataset('matrix.txt')
    if days is not None:
        wot.io.incorporate_days_information_in_datatset(ds, days)
    if covariate is not None:
        pass
    if gene_sets is not None:
        pass
    return ds
