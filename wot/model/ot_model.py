# -*- coding: utf-8 -*-

import itertools
import os
from multiprocessing import Process

import pandas as pd
import numpy as np
import sklearn
import scipy

import wot.io
import wot.model


class OTModel:
    """
    The OTModel computes transport maps.

    Parameters
    ----------
    matrix : wot.Dataset
        The gene expression matrix for this OTModel. Matrix must have the row meta data field 'day'.
    transport_maps_directory : str
        Path to the transport map directory, where transport maps are written.
    transport_maps_prefix : str, optional
        Prefix to use for the transport maps. This can highly speed up
        initialization if the directory is filled with other non-tmap files,
        and allows to have several transport maps not overriding each other.
        If None, all files named `{prefix}_{t0}_{t1}.{extension}` will be
        considered as transport maps.
        The default prefix for transport maps is 'tmaps'
    max_threads : int, optional
        Maximum number of threads to use when computing transport maps
    **kwargs : dict
        Dictionnary of parameters. Will be inserted as is into OT configuration.
    """

    def __init__(self, matrix, tmap_out, max_threads=None, **kwargs):
        tmap_dir, tmap_prefix = os.path.split(tmap_out) if tmap_out is not None else (None, None)
        self.matrix = matrix
        self.tmap_dir = tmap_dir or '.'
        self.tmap_prefix = tmap_prefix or "tmaps"
        self.day_pairs = wot.model.parse_configuration(kwargs.pop('day_pairs', None))
        cell_filter = kwargs.pop('cell_filter', None)
        gene_filter = kwargs.pop('gene_filter', None)
        self.output_file_format = kwargs.pop('output_file_format', 'loom')
        if gene_filter is not None:
            if os.path.isfile(gene_filter):
                gene_ids = pd.read_table(gene_filter, index_col=0, header=None) \
                    .index.values
            else:
                import re
                expr = re.compile(gene_filter)
                gene_ids = [e for e in self.matrix.col_meta.index.values if expr.match(e)]
            col_indices = self.matrix.col_meta.index.isin(gene_ids)
            self.matrix = wot.Dataset(self.matrix.x[:, col_indices],
                                      self.matrix.row_meta, self.matrix.col_meta[col_indices].copy(False))
            wot.io.verbose('Successfuly applied gene_filter: "{}"'.format(gene_filter))
        if cell_filter is not None:
            if os.path.isfile(cell_filter):
                cell_ids = pd.read_table(cell_filter, index_col=0, header=None) \
                    .index.values
            else:
                import re
                expr = re.compile(cell_filter)
                cell_ids = [e for e in self.matrix.row_meta.index.values if expr.match(e)]
            row_indices = self.matrix.row_meta.index.isin(cell_ids)
            self.matrix = wot.Dataset(self.matrix.x[row_indices, :],
                                      self.matrix.row_meta[row_indices].copy(False), self.matrix.col_meta)
            wot.io.verbose('Successfuly applied cell_filter: "{}"'.format(cell_filter))
        self.timepoints = sorted(set(self.matrix.row_meta['day']))
        wot.io.verbose(len(self.timepoints), "timepoints loaded :", self.timepoints)

        if max_threads is None or max_threads == 0:
            try:
                max_usable_cores = len(os.sched_getaffinity(0))
            except Exception:
                import multiprocessing
                max_usable_cores = multiprocessing.cpu_count()
            if kwargs.pop('fast', False):
                wot.io.verbose("Fast mode. Using all but one core")
                self.max_threads = max_usable_cores - 1
            else:
                self.max_threads = 1
            wot.io.verbose("Argument max_threads not set. Using " + str(self.max_threads))
        else:
            self.max_threads = max_threads
        wot.io.verbose("Using", self.max_threads, "thread(s) at most")
        if self.max_threads > 1:
            wot.io.verbose("Warning : Multiple threads are being used. Time estimates will be inaccurate")

        self.ot_config = {}
        for k in kwargs.keys():
            self.ot_config[k] = kwargs[k]
        local_pca = self.get_ot_config()['local_pca']
        if local_pca > self.matrix.x.shape[1]:
            print("Warning : local_pca set to {}, above gene count of {}. Disabling PCA" \
                  .format(local_pca, self.matrix.x.shape[1]))
            self.ot_config['local_pca'] = 0
        if 'day' not in self.matrix.row_meta.columns:
            raise ValueError("Days information not available for matrix")
        if any(self.matrix.row_meta['day'].isnull()):
            query = self.matrix.row_meta['day'].isnull()
            faulty = list(self.matrix.row_meta.index[query])
            raise ValueError("Days information missing for cells : {}".format(faulty))

    def get_ot_config(self):
        """
        Get valid parameters for the Optimal Transport computation.

        Returns
        -------
        ot_config : dict
            Dictionnary of valid parameters, or defaults if unspecified.
        """
        # WARNING: Any value in ot_config that does not appear in the following dict will be ignored
        ot_defaults = {
            'epsilon': .05, 'lambda1': 1, 'lambda2': 50,
            'epsilon0': 1, 'tau': 1e4,
            'growth_iters': 3, 'batch_size': 50,
            'local_pca': 30, 'max_iter': 1e7,
            'tolerance': 1e-2,
        }
        # TODO: support ncells and ncounts
        config = self.ot_config

        return {x: config[x] if x in config else ot_defaults[x] for x in ot_defaults}

    def get_covariate_pairs(self):
        """Get all covariate pairs in the dataset"""
        if 'covariate' not in self.matrix.row_meta.columns:
            raise ValueError("Covariate value not available in dataset")
        from itertools import product
        covariate = sorted(set(self.matrix.row_meta['covariate']))
        return product(covariate, covariate)

    def compute_all_transport_maps(self, with_covariates=False):
        """
        Computes all required transport maps.

        Parameters
        ----------
        force : bool, optional, default : False
            Force recomputation of each transport map, after config update for instance.
        with_covariates : bool, optional, default : False
            Compute all covariate-restricted transport maps as well

        Returns
        -------
        None
            Only computes and saves all transport maps, does not return them.
        """
        t = self.timepoints
        day_pairs = self.day_pairs
        if day_pairs is None:
            day_pairs = [(t[i], t[i + 1]) for i in range(len(t) - 1)]

        if with_covariates:
            day_pairs = [(*d, c) for d, c in itertools.product(day_pairs, self.get_covariate_pairs())]

        # if not force:
        #     if with_covariates:
        #         day_pairs = [(t0, t1, cv) for t0, t1, cv in day_pairs
        #                      if self.cov_tmaps.get((t0, t1, *cv), None) is None]
        #     else:
        #         day_pairs = [x for x in day_pairs if self.tmaps.get(x, None) is None]

        m = self.max_threads

        if not day_pairs:
            print('No day pairs')
            return

        if m > 1:
            procs = []
            for x in day_pairs:
                p = Process(target=self.compute_transport_map, args=(*x,))
                procs.append(p)

            for i in range(len(procs) + m):
                if i >= m:
                    procs[i - m].join()
                if i < len(procs):
                    procs[i].start()
        else:
            for x in day_pairs:
                self.compute_transport_map(*x)

    def compute_transport_map(self, t0, t1, covariate=None):
        """
        Computes the transport map from time t0 to time t1

        Parameters
        ----------
        t0 : float
            Source timepoint for the transport map
        t1 : float
            Destination timepoint for the transport map
        covariate : None or (int, int)
            The covariate restriction on cells from t0 and t1. None to skip

        Returns
        -------
        None
            Only computes and saves the transport maps, does not return it.

        Raises
        ------
        ValueError
            If the OTModel was initialized with day_pairs and the given pair is not present.
        """
        path = self.tmap_prefix
        wot.io.verbose("Computing tmap ({},{})".format(t0, t1))
        # If day_pairs is not None, its configuration takes precedence
        if self.day_pairs is not None:
            if (t0, t1) not in self.day_pairs:
                raise ValueError("Transport map ({},{}) is not present in day_pairs".format(t0, t1))
            local_config = self.day_pairs[(t0, t1)]
        else:
            local_config = {}

        config = {**self.get_ot_config(), **local_config, 't0': t0, 't1': t1, 'covariate': covariate}
        tmap = OTModel.compute_single_transport_map(self.matrix, config)
        if covariate is None:
            path += "_{}_{}".format(t0, t1)

        else:
            path += "_{}_{}_cv{}_cv{}".format(t0, t1, *covariate)

        wot.io.write_dataset(tmap, os.path.join(self.tmap_dir, path),
                             output_format=self.output_file_format)
        wot.io.verbose("Created tmap ({}, {}) : {}".format(t0, t1, path))

    @staticmethod
    def compute_default_cost_matrix(a, b, eigenvals=None):
        if eigenvals is not None:
            a = a.dot(eigenvals)
            b = b.dot(eigenvals)

        cost_matrix = sklearn.metrics.pairwise.pairwise_distances(a.toarray() if scipy.sparse.isspmatrix(a) else a,
                                                                  b.toarray() if scipy.sparse.isspmatrix(b) else b,
                                                                  metric='sqeuclidean')
        cost_matrix = cost_matrix / np.median(cost_matrix)
        return cost_matrix

    @staticmethod
    def compute_single_transport_map(ds, config):
        """
        Computes a single transport map

        Parameters
        ----------
        ds : wot.Dataset
            The gene expression matrix to consider.
            It is assumed to have a valid day column for each cell.
        config : dict
            Configuration to use for all parameters for the couplings :
            - t0, t1
            - lambda1, lambda2, epsilon, g
        """
        t0 = config.pop('t0', None)
        t1 = config.pop('t1', None)
        if t0 is None or t1 is None:
            raise ValueError("config must have both t0 and t1, indicating target timepoints")

        covariate = config.pop('covariate', None)
        if covariate is None:
            p0 = ds.where(day=float(t0))
            p1 = ds.where(day=float(t1))
        else:
            p0 = ds.where(day=float(t0), covariate=covariate[0])
            p1 = ds.where(day=float(t1), covariate=covariate[1])

        if 'cell_growth_rate' in p0.row_meta.columns:
            config['g'] = np.asarray(p0.row_meta['cell_growth_rate'].values)
        if 'pp' in p0.row_meta.columns:
            config['pp'] = np.asarray(p0.row_meta['pp'].values)
        if 'pp' in p1.row_meta.columns:
            config['qq'] = np.asarray(p1.row_meta['pp'].values)

        local_pca = config.pop('local_pca', None)
        if local_pca is not None and local_pca > 0:
            pca = wot.ot.get_pca(local_pca, p0.x, p1.x)
            p0_x = wot.ot.pca_transform(pca, p0.x)
            p1_x = wot.ot.pca_transform(pca, p1.x)
        else:
            p0_x = p0.x
            p1_x = p1.x

        C = OTModel.compute_default_cost_matrix(p0_x, p1_x)
        if config.get('g') is None:
            config['g'] = np.ones(C.shape[0])
        tmap = wot.ot.transport_stablev1_learnGrowth(C, **config)
        return wot.Dataset(tmap, p0.row_meta.copy(), p1.row_meta.copy())
