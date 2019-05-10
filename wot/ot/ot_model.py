# -*- coding: utf-8 -*-

import itertools
import logging
import os

import anndata
import numpy as np
import pandas as pd
import scipy
import sklearn

import wot.io
import wot.ot

logger = logging.getLogger('wot')


class OTModel:
    """
    The OTModel computes transport maps.

    Parameters
    ----------
    matrix : anndata.AnnData
        The gene expression matrix for this OTModel.
    day_field : str, optional
        Cell day obs name
    covariate_field : str, optional
        Cell covariate obs name
    cell_growth_rate_field : str, optional
        Cell growth rate obs name
    **kwargs : dict
        Dictionary of parameters. Will be inserted as is into OT configuration.
    """

    def __init__(self, matrix, day_field='day', covariate_field='covariate',
                 growth_rate_field='cell_growth_rate', **kwargs):
        self.matrix = matrix
        self.day_field = day_field
        self.covariate_field = covariate_field
        self.cell_growth_rate_field = growth_rate_field
        self.day_pairs = wot.ot.parse_configuration(kwargs.pop('config', None))
        cell_filter = kwargs.pop('cell_filter', None)
        gene_filter = kwargs.pop('gene_filter', None)
        day_filter = kwargs.pop('cell_day_filter', None)
        ncounts = kwargs.pop('ncounts', None)
        ncells = kwargs.pop('ncells', None)
        self.matrix = wot.io.filter_adata(self.matrix, obs_filter=cell_filter, var_filter=gene_filter)
        if day_filter is not None:
            days = [float(day) for day in day_filter.split(',')] if type(day_filter) == str else day_filter
            row_indices = self.matrix.obs[self.day_field].isin(days)
            self.matrix = self.matrix[row_indices].copy()

        cvs = set(self.matrix.obs[self.covariate_field]) if self.covariate_field in self.matrix.obs else [None]
        if ncells is not None:
            index_list = []
            for day in self.timepoints:
                day_query = self.matrix.obs[self.day_field] == day
                for cv in cvs:
                    if cv is None:
                        indices = np.where(day_query)[0]
                    else:
                        indices = np.where(day_query & (self.matrix.obs[self.covariate_field] == cv))[0]
                    if len(indices) > ncells:
                        np.random.shuffle(indices)
                        indices = indices[0:ncells]
                    index_list.append(indices)
            row_indices = np.concatenate(index_list)
            self.matrix = self.matrix[row_indices]
        if ncounts is not None:
            for i in range(self.matrix.X.shape[0]):
                p = self.matrix[i].X
                if scipy.sparse.isspmatrix(p):
                    p = p.toarray()
                p = p.astype('float64')
                total = p.sum()
                if total > ncounts:
                    p /= total
                    self.matrix.X[i] = np.random.multinomial(ncounts, p, size=1)[0]

        if self.matrix.X.shape[0] is 0:
            raise ValueError('No cells in matrix')

        self.ot_config = {'local_pca': 30, 'growth_iters': 1, 'epsilon': 0.05, 'lambda1': 1, 'lambda2': 50,
                          'epsilon0': 1, 'tau': 10000, 'scaling_iter': 3000, 'inner_iter_max': 50, 'tolerance': 1e-8,
                          'max_iter': 1e7, 'batch_size': 5, 'extra_iter': 1000}
        solver = kwargs.pop('solver', 'duality_gap')
        if solver == 'fixed_iters':
            self.solver = wot.ot.transport_stablev2
        elif solver == 'duality_gap':
            self.solver = wot.ot.optimal_transport_duality_gap
        else:
            raise ValueError('Unknown solver')

        parameters_from_file = kwargs.pop('parameters', None)
        for k in kwargs.keys():
            self.ot_config[k] = kwargs[k]

        if parameters_from_file is not None:
            config_dict = wot.ot.parse_parameter_file(parameters_from_file)
            for k in config_dict.keys():
                self.ot_config[k] = config_dict[k]

        local_pca = self.ot_config['local_pca']
        if local_pca > self.matrix.X.shape[1]:
            logger.warning("local_pca set to {}, above gene count of {}. Disabling PCA" \
                           .format(local_pca, self.matrix.X.shape[1]))
            self.ot_config['local_pca'] = 0
        if self.day_field not in self.matrix.obs.columns:
            raise ValueError("Days information not available for matrix")
        if any(self.matrix.obs[self.day_field].isnull()):
            self.matrix = self.matrix[self.matrix.obs[self.day_field].isnull() == False]
        self.timepoints = sorted(set(self.matrix.obs[self.day_field]))

    def get_covariate_pairs(self):
        """Get all covariate pairs in the dataset"""
        if self.covariate_field not in self.matrix.obs.columns:
            raise ValueError("Covariate value not available in dataset")
        from itertools import product
        covariate = set(self.matrix.obs[self.covariate_field])
        return product(covariate, covariate)

    def compute_all_transport_maps(self, tmap_out='tmaps', overwrite=True, output_file_format='h5ad',
                                   with_covariates=False):
        """
        Computes all required transport maps.

        Parameters
        ----------
        tmap_out : str, optional
            Path and prefix for output transport maps
        overwrite : bool, optional
            Overwrite existing transport maps
        output_file_format: str, optional
            Transport map file format
        with_covariates : bool, optional, default : False
            Compute all covariate-restricted transport maps as well

        Returns
        -------
        None
            Only computes and saves all transport maps, does not return them.
        """

        tmap_dir, tmap_prefix = os.path.split(tmap_out) if tmap_out is not None else (None, None)
        tmap_prefix = tmap_prefix or "tmaps"
        tmap_dir = tmap_dir or '.'
        if not os.path.exists(tmap_dir):
            os.makedirs(tmap_dir)
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

        if not day_pairs:
            logger.info('No day pairs')
            return

        full_learned_growth_df = None
        save_learned_growth = self.ot_config.get('growth_iters', 1) > 1
        for day_pair in day_pairs:
            path = tmap_prefix
            if not with_covariates:
                path += "_{}_{}".format(*day_pair)
            else:
                path += "_{}_{}_cv{}_cv{}".format(*day_pair)
            output_file = os.path.join(tmap_dir, path)
            output_file = wot.io.check_file_extension(output_file, output_file_format)
            if os.path.exists(output_file) and not overwrite:
                logger.info('Found existing tmap at ' + output_file + '. ')
                continue

            tmap = self.compute_transport_map(*day_pair)
            wot.io.write_dataset(tmap, output_file, output_format=output_file_format)
            if save_learned_growth:
                learned_growth_df = tmap.obs
                full_learned_growth_df = learned_growth_df if full_learned_growth_df is None else pd.concat(
                    (full_learned_growth_df, learned_growth_df), copy=False)
        if full_learned_growth_df is not None:
            full_learned_growth_df.to_csv(os.path.join(tmap_dir, tmap_prefix + '_g.txt'), sep='\t', index_label='id')

    def compute_transport_map(self, t0, t1, covariate=None):
        """
        Computes the transport map from time t0 to time t1

        Parameters
        ----------
        t0 : float
            Source timepoint for the transport map
        t1 : float
            Destination timepoint for the transport map
        covariate : None or (str, str)
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
        if self.day_pairs is not None:
            if (t0, t1) not in self.day_pairs:
                raise ValueError("Transport map ({},{}) is not present in day_pairs".format(t0, t1))
            local_config = self.day_pairs[(t0, t1)]
        else:
            local_config = {}
        if covariate is None:
            logger.info('Computing transport map from {} to {}'.format(t0, t1))
        else:
            logger.info('Computing transport map from {} {} to {} {}'.format(t0, covariate[0], t1, covariate[1]))
        config = {**self.ot_config, **local_config, 't0': t0, 't1': t1, 'covariate': covariate}
        return self.compute_single_transport_map(config)

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

    def compute_single_transport_map(self, config):
        """
        Computes a single transport map.

        Parameters
        ----------
        config : dict
            Configuration to use for all parameters for the couplings :
            - t0, t1
            - lambda1, lambda2, epsilon, g
        """

        import gc
        gc.collect()

        t0 = config.pop('t0', None)
        t1 = config.pop('t1', None)
        if t0 is None or t1 is None:
            raise ValueError("config must have both t0 and t1, indicating target timepoints")
        ds = self.matrix
        covariate = config.pop('covariate', None)
        if covariate is None:
            p0_indices = ds.obs[self.day_field] == float(t0)
            p1_indices = ds.obs[self.day_field] == float(t1)
        else:
            p0_indices = (ds.obs[self.day_field] == float(t0)) & (ds.obs[self.covariate_field] == covariate[0])
            p1_indices = (ds.obs[self.day_field] == float(t1)) & (ds.obs[self.covariate_field] == covariate[1])

        p0 = ds[p0_indices, :]
        p1 = ds[p1_indices, :]

        if p0.shape[0] == 0:
            logger.info('No cells at {}'.format(t0))
            return None
        if p1.shape[0] == 0:
            logger.info('No cells at {}'.format(t1))
            return None

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
        config['C'] = C
        delta_days = t1 - t0

        if self.cell_growth_rate_field in p0.obs.columns:
            config['G'] = np.power(p0.obs[self.cell_growth_rate_field].values, delta_days)
        else:
            config['G'] = np.ones(C.shape[0])
        tmap, learned_growth = wot.ot.compute_transport_matrix(solver=self.solver, **config)
        learned_growth.append(tmap.sum(axis=1))
        obs_growth = {}
        for i in range(len(learned_growth)):
            g = learned_growth[i]
            g = np.power(g, 1.0 / delta_days)
            obs_growth['g' + str(i)] = g
        obs = pd.DataFrame(index=p0.obs.index, data=obs_growth)
        return anndata.AnnData(tmap, obs, pd.DataFrame(index=p1.obs.index))
