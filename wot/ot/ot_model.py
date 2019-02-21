# -*- coding: utf-8 -*-

import os

import anndata
import itertools
import numpy as np
import pandas as pd
import scipy
import sklearn

import wot.io
import wot.ot


class OTModel:
    """
    The OTModel computes transport maps.

    Parameters
    ----------
    matrix : anndata.AnnData
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
        self.day_pairs = wot.ot.parse_configuration(kwargs.pop('config', None))

        cell_filter = kwargs.pop('cell_filter', None)
        gene_filter = kwargs.pop('gene_filter', None)
        day_filter = kwargs.pop('cell_day_filter', None)
        ncounts = kwargs.pop('ncounts', None)
        ncells = kwargs.pop('ncells', None)
        self.force = kwargs.pop('force', False)
        self.output_file_format = kwargs.pop('format', 'loom')
        if gene_filter is not None:
            if os.path.isfile(gene_filter):
                gene_ids = pd.read_csv(gene_filter, index_col=0, header=None) \
                    .index.values
            else:
                import re
                expr = re.compile(gene_filter)
                gene_ids = [e for e in self.matrix.var.index.values if expr.match(e)]
            col_indices = self.matrix.var.index.isin(gene_ids)
            if np.sum(col_indices) is 0:
                raise ValueError('No genes passed the gene filter')
            self.matrix = anndata.AnnData(self.matrix.X[:, col_indices],
                                          self.matrix.obs, self.matrix.var[col_indices].copy(False))
            wot.io.verbose('Successfuly applied gene_filter: "{}"'.format(gene_filter))
        if cell_filter is not None:
            if os.path.isfile(cell_filter):
                cell_ids = pd.read_csv(cell_filter, index_col=0, header=None) \
                    .index.values
            else:
                import re
                expr = re.compile(cell_filter)
                cell_ids = [e for e in self.matrix.obs.index.values if expr.match(e)]
            row_indices = self.matrix.obs.index.isin(cell_ids)
            if np.sum(row_indices) is 0:
                raise ValueError('No cells passed the cell filter')
            self.matrix = anndata.AnnData(self.matrix.X[row_indices, :],
                                          self.matrix.obs[row_indices].copy(False), self.matrix.var)

            wot.io.verbose('Successfuly applied cell_filter: "{}"'.format(cell_filter))
        if day_filter is not None:
            days = day_filter.split(',')
            row_indices = self.matrix.obs['day'].isin(days)
            self.matrix = anndata.AnnData(self.matrix.X[row_indices, :],
                                          self.matrix.obs[row_indices].copy(False), self.matrix.var)

            wot.io.verbose('Successfuly applied day_filter: "{}"'.format(day_filter))
        self.timepoints = sorted(set(self.matrix.obs['day']))

        cvs = set(self.matrix.obs['covariate']) if 'covariate' in self.matrix.obs else [None]
        if ncells is not None:
            index_list = []
            for day in self.timepoints:
                day_query = self.matrix.obs['day'] == day
                for cv in cvs:
                    if cv is None:
                        indices = np.where(day_query)[0]
                    else:
                        indices = np.where(day_query & (self.matrix.obs['covariate'] == cv))[0]
                    if len(indices) > ncells:
                        np.random.shuffle(indices)
                        indices = indices[0:ncells]
                    index_list.append(indices)
            row_indices = np.concatenate(index_list)
            self.matrix = anndata.AnnData(self.matrix.X[row_indices, :],
                                          self.matrix.obs.iloc[row_indices].copy(False), self.matrix.var)
        if ncounts is not None:
            for i in range(self.matrix.X.shape[0]):
                p = self.matrix.X[i]
                if scipy.sparse.isspmatrix(p):
                    p = p.toarray()
                p = p.astype('float64')
                total = p.sum()
                if total > ncounts:
                    p /= total
                    self.matrix.X[i] = np.random.multinomial(ncounts, p, size=1)[0]

        if self.matrix.X.shape[0] is 0:
            print('No cells in matrix')
            exit(1)
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

        self.ot_config = {'local_pca': 30, 'growth_iters': 3, 'scaling_iter': 3000, 'inner_iter_max': 50,
                          'epsilon': 0.05, 'lambda1': 1, 'lambda2': 50, 'epsilon0': 1, 'tau': 10000}

        parameters_from_file = kwargs.pop('parameters', None)
        for k in kwargs.keys():
            self.ot_config[k] = kwargs[k]

        if parameters_from_file is not None:
            config_dict = wot.ot.parse_parameter_file(parameters_from_file)
            for k in config_dict.keys():
                self.ot_config[k] = config_dict[k]

        local_pca = self.ot_config['local_pca']
        if local_pca > self.matrix.X.shape[1]:
            print("Warning : local_pca set to {}, above gene count of {}. Disabling PCA" \
                  .format(local_pca, self.matrix.X.shape[1]))
            self.ot_config['local_pca'] = 0
        if 'day' not in self.matrix.obs.columns:
            raise ValueError("Days information not available for matrix")
        if any(self.matrix.obs['day'].isnull()):
            print("Days information missing for {} cells".format(self.matrix.obs['day'].isnull().sum()))
            self.matrix = self.matrix[self.matrix.obs['day'].isnull() == False]

    def get_covariate_pairs(self):
        """Get all covariate pairs in the dataset"""
        if 'covariate' not in self.matrix.obs.columns:
            raise ValueError("Covariate value not available in dataset")
        from itertools import product
        covariate = set(self.matrix.obs['covariate'])
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

        if day_pairs is None or len(day_pairs) == 0:
            day_pairs = [(t[i], t[i + 1]) for i in range(len(t) - 1)]

        if with_covariates:
            covariate_day_pairs = [(*d, c) for d, c in itertools.product(day_pairs, self.get_covariate_pairs())]
            # if type(day_pairs) is dict:
            #     day_pairs = list(day_pairs.keys())
            day_pairs = covariate_day_pairs

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
            from joblib import Parallel, delayed
            Parallel(n_jobs=m)(delayed(self.compute_transport_map)(*x) for x in day_pairs)
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
        anndata.AnnData
            The transport map from t0 to t1

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

        if covariate is None:
            path += "_{}_{}".format(t0, t1)
        else:
            path += "_{}_{}_cv{}_cv{}".format(t0, t1, *covariate)
        output_file = os.path.join(self.tmap_dir, path)
        output_file = wot.io.check_file_extension(output_file, self.output_file_format)
        if os.path.exists(output_file) and not self.force:
            wot.io.verbose('Found existing tmap at ' + output_file + '. Use --force to overwrite.')
            return wot.io.read_dataset(output_file)

        config = {**self.ot_config, **local_config, 't0': t0, 't1': t1, 'covariate': covariate}
        tmap = OTModel.compute_single_transport_map(self.matrix, config)
        if tmap is not None:
            wot.io.write_dataset(tmap, output_file, output_format=self.output_file_format)
            wot.io.verbose("Created tmap ({}, {}) : {}".format(t0, t1, path))
        return tmap

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
        Computes a single transport map.
        Note that None is returned if no data is available at the specified timepoints or covariates.

        Parameters
        ----------
        ds : anndata.AnnData
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
            p0_indices = ds.obs['day'] == float(t0)
            p1_indices = ds.obs['day'] == float(t1)
        else:
            p0_indices = (ds.obs['day'] == float(t0)) & (ds.obs['covariate'] == covariate[0])
            p1_indices = (ds.obs['day'] == float(t1)) & (ds.obs['covariate'] == covariate[1])

        if p0_indices.sum() == 0 or p1_indices.sum() == 0:
            return None
        p0 = ds[p0_indices, :]
        p1 = ds[p1_indices, :]

        if 'cell_growth_rate' in p0.obs.columns:
            config['g'] = np.asarray(p0.obs['cell_growth_rate'].values)
        if 'pp' in p0.obs.columns:
            config['pp'] = np.asarray(p0.obs['pp'].values)
        if 'pp' in p1.obs.columns:
            config['qq'] = np.asarray(p1.obs['pp'].values)

        local_pca = config.pop('local_pca', None)
        eigenvals = None
        if local_pca is not None and local_pca > 0:
            # pca, mean = wot.ot.get_pca(local_pca, p0.X, p1.X)
            # p0_x = wot.ot.pca_transform(pca, mean, p0.X)
            # p1_x = wot.ot.pca_transform(pca, mean, p1.X)
            p0_x, p1_x, pca, mean = wot.ot.compute_pca(p0.X, p1.X, local_pca)
            eigenvals = np.diag(pca.singular_values_)
        else:
            p0_x = p0.X
            p1_x = p1.X

        C = OTModel.compute_default_cost_matrix(p0_x, p1_x, eigenvals)
        if config.get('g') is None:
            config['g'] = np.ones(C.shape[0])
        delta_days = t1 - t0
        config['g'] = config['g'] ** delta_days
        tmap = wot.ot.transport_stable_learn_growth(C, **config)
        return anndata.AnnData(tmap, p0.obs.copy(), p1.obs.copy())
