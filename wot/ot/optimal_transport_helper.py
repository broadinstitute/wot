import numpy as np
import pandas as pd
import csv
import scipy
import sklearn.metrics
import sklearn.decomposition
import argparse
import io
import wot.ot
import wot.io
import os


class OptimalTransportHelper:

    @staticmethod
    def create_base_parser(description):
        parser = argparse.ArgumentParser(
            description=description)

        parser.add_argument('--matrix',
                            help='Gene expression tab delimited file with cells on '
                                 'rows and features on columns', required=True)

        parser.add_argument('--cell_days',
                            help='Two column tab delimited file without header with '
                                 'cell ids and days', required=True)
        parser.add_argument('--day_pairs',
                            help='Two column tab delimited file without header with '
                                 'pairs of days to compute transport maps for',
                            required=True)

        parser.add_argument('--epsilon', type=float, default=0.05,
                            help='Controls the entropy of the transport map. An '
                                 'extremely large entropy parameter will give a '
                                 'maximally entropic transport map, and an '
                                 'extremely '
                                 'small entropy parameter will give a nearly '
                                 'deterministic transport map (but could also '
                                 'lead to '
                                 'numerical instability in the algorithm')
        parser.add_argument('--prefix',
                            help='Prefix for ouput file names', required=True)
        parser.add_argument('--ncells', help='Number of cells to sample from each timepoint', type=int)

        parser.add_argument('--max_transport_fraction',
                            default=0.4,
                            help='The maximum fraction of cells at time t that are '
                                 'transported to time t + 1',
                            type=float)
        parser.add_argument('--min_transport_fraction',
                            default=0.05,
                            help='The minimum fraction of cells at time t that are '
                                 'transported to time t + 1',
                            type=float)
        parser.add_argument('--lambda1',
                            help='Regularization parameter that controls the '
                                 'fidelity of the constraints on p', type=float, default=1)
        parser.add_argument('--lambda2', default=50,
                            help='Regularization parameter that controls the '
                                 'fidelity of the constraints on q', type=float)
        parser.add_argument('--scaling_iter', default=3000,
                            help='Number of scaling iterations', type=int)
        parser.add_argument('--min_growth_fit', type=float, default=0.9)
        parser.add_argument('--l0_max', type=float, default=100)

        parser.add_argument('--epsilon_adjust', help='Scaling factor to adjust epsilon for floating_epsilon solver',
                            type=float, default=1.1)
        parser.add_argument('--lambda_adjust', help='Scaling factor to adjust lambda for floating_epsilon solver',
                            type=float, default=1.5)

        parser.add_argument('--numItermax', type=int, default=100, help='For sinkhorn_epsilon solver')
        parser.add_argument('--epsilon0', type=float, default=1, help='For sinkhorn_epsilon and unbalanced solvers')
        parser.add_argument('--numInnerItermax', type=int, default=50,
                            help='For sinkhorn_epsilon and unbalanced solvers')
        parser.add_argument('--tau', type=float, default=10000, help='For sinkhorn_epsilon and unbalanced solvers')
        parser.add_argument('--stopThr', type=float, default=1e-10, help='For sinkhorn_epsilon solver')

        parser.add_argument('--beta_min', type=float, default=0.3, help='Growth function parameter')
        parser.add_argument('--delta_min', type=float, default=0.3, help='Growth function parameter')
        parser.add_argument('--beta_max', type=float, default=1.7, help='Growth function parameter')
        parser.add_argument('--delta_max', type=float, default=1.7, help='Growth function parameter')
        parser.add_argument('--beta_center', type=float, default=0.25, help='Growth function parameter')

        parser.add_argument('--growth_iters', type=int, default=3, help='Number of growth iterations')

        growth_rate_group = parser.add_mutually_exclusive_group(required=False)
        growth_rate_group.add_argument('--gene_set_scores', help='File containing "Cell.cycle" and "Apoptosis" scores')

        growth_rate_group.add_argument('--cell_growth_rates',
                                       help='Two column tab delimited file without '
                                            'header with '
                                            'cell ids and growth rates per day.')
        parser.add_argument('--diagonal', help='Diagonal scaling matrix')
        parser.add_argument('--power', help='Diagonal scaling power', type=float)
        parser.add_argument('--local_pca', help='Convert day pairs matrix to PCA coordinates', type=int)

        parser.add_argument('--solver',
                            help='Solver to use when computing transport maps. One of unbalanced, floating_epsilon, '
                                 'sinkhorn_epsilon, unregularized',
                            choices=['epsilon', 'sinkhorn_epsilon', 'unbalanced', 'unregularized'],
                            default='unbalanced')
        parser.add_argument('--cell_filter',
                            help='File with one cell id per line to include or or a python regular expression of cell ids to include')
        parser.add_argument('--verbose', action='store_true',
                            help='Print progress information')
        return parser

    def __init__(self, args, covariate_df=None, covariate_pairs=None):
        eigenvals = None
        if args.diagonal is not None:
            eigenvals = \
                pd.read_table(args.diagonal, header=None, names=['eigenvals'], index_col=False, dtype=np.float64)[
                    'eigenvals'].values
        if eigenvals is not None and args.power is not None:
            eigenvals = np.power(eigenvals, args.power)

        # cells on rows, features on columns
        ds = wot.io.read_dataset(args.matrix)

        if args.cell_filter is not None:
            prior = ds.x.shape[0]
            if not os.path.isfile(args.cell_filter):
                import re
                expr = re.compile(args.cell_filter)
                cell_ids = [elem for elem in ds.row_meta.index.values if expr.match(elem)]
            else:
                cell_ids = pd.read_table(args.cell_filter, index_col=0, header=None).index.values

            # row_indices = np.isin(ds.row_meta.index.values, cell_ids, assume_unique=True)
            row_indices = ds.row_meta.index.isin(cell_ids)
            nkeep = np.sum(row_indices)
            if args.verbose and len(cell_ids) > nkeep:
                print(str(len(cell_ids) - nkeep) + ' are in cell filter, but not in matrix')

            ds = wot.Dataset(ds.x[row_indices], ds.row_meta.iloc[row_indices], ds.col_meta)
            if args.verbose:
                print('Keeping ' + str(ds.x.shape[0]) + '/' + str(prior) + ' cells')

        if not os.path.isfile(args.day_pairs):
            day_pairs = pd.read_table(io.StringIO(args.day_pairs), header=None, names=['t0', 't1'],
                                      index_col=False, lineterminator=';', sep=',', dtype=np.float64)
        else:
            day_pairs = pd.read_table(args.day_pairs, header=None, names=['t0', 't1'],
                                      index_col=False, quoting=csv.QUOTE_NONE,
                                      engine='python', sep=None, dtype=np.float64)

        days_data_frame = pd.read_table(args.cell_days, index_col=0, header=None,
                                        names=['day'], quoting=csv.QUOTE_NONE,
                                        engine='python', sep=None,
                                        dtype={'day': np.float64})

        self.eigenvals = np.diag(eigenvals) if eigenvals is not None else None

        if args.gene_set_scores is not None:
            ext = wot.io.get_file_basename_and_extension(args.gene_set_scores)[1]
            if ext == 'loom' or ext == 'gct':
                gene_set_scores_ds = wot.io.read_dataset(args.gene_set_scores)
                apoptosis = gene_set_scores_ds[:, np.where(gene_set_scores_ds.col_meta.columns == 'Apoptosis')[0]]
                proliferation = gene_set_scores_ds[:, np.where(gene_set_scores_ds.col_meta.columns == 'Cell.cycle')[0]]
                gene_set_scores_ids = gene_set_scores_ds.row_meta.index
            else:
                gene_set_scores = pd.read_table(args.gene_set_scores, index_col=0, quoting=csv.QUOTE_NONE,
                                                engine='python', sep=None)
                apoptosis = gene_set_scores['Apoptosis'].values
                proliferation = gene_set_scores['Cell.cycle'].values
                gene_set_scores_ids = gene_set_scores.index
            g = wot.ot.compute_growth_scores(proliferation, apoptosis, beta_max=args.beta_max,
                                             beta_center=args.beta_center,
                                             delta_max=args.delta_max, delta_min=args.delta_min,
                                             beta_min=args.beta_min)
            cell_growth_rates = pd.DataFrame(index=gene_set_scores_ids, data={'cell_growth_rate': g})

        elif args.cell_growth_rates is not None:
            cell_growth_rates = pd.read_table(args.cell_growth_rates, index_col=0,
                                              header=None, names=['cell_growth_rate'],
                                              quoting=csv.QUOTE_NONE, engine='python',
                                              sep=None)
        else:
            cell_growth_rates = pd.DataFrame(index=ds.row_meta.index.values, data={'cell_growth_rate': 1})
            if args.verbose:
                print('Using growth rate of 1')
        ds.row_meta = ds.row_meta.join(cell_growth_rates).join(days_data_frame)

        if covariate_df is not None:
            self.covariate_df = covariate_df
            ds.row_meta = ds.row_meta.join(covariate_df)
            self.covariate_pairs = covariate_pairs + [[None, None]]
        else:
            self.covariate_df = None
            self.covariate_pairs = [[None, None]]

        day_to_indices = {}
        days = ds.row_meta[days_data_frame.columns[0]].values

        for i in range(len(days)):
            val = days[i]
            if val is not None:
                indices = day_to_indices.get(val)
                if indices is None:
                    indices = []
                    day_to_indices[val] = indices
                indices.append(i)

        self.ds = ds
        if args.ncells is not None:
            for day in day_to_indices:
                indices = day_to_indices[day]
                if len(indices) > args.ncells:
                    np.random.shuffle(indices)
                    indices = indices[0:args.ncells]
                    day_to_indices[day] = indices

        if day_pairs.shape[0] is 0:
            print('No day pairs found')
            exit(1)
        if args.verbose:
            print('Computing ' + str(day_pairs.shape[0]) + ' transport map' + ('s' if
                                                                               day_pairs.shape[0] > 1 else ''))
        self.day_pairs = day_pairs
        self.day_to_indices = day_to_indices
        self.cell_growth_rates = cell_growth_rates
        self.args = args
        self.t_interpolate = vars(args).get('t_interpolate')

    def compute_cost_matrix(self, a, b):
        if self.eigenvals is not None:
            a = a.dot(self.eigenvals)
            b = b.dot(self.eigenvals)

        cost_matrix = sklearn.metrics.pairwise.pairwise_distances(a.toarray() if scipy.sparse.isspmatrix(a) else a,
                                                                  b.toarray() if scipy.sparse.isspmatrix(b) else b,
                                                                  metric='sqeuclidean')
        cost_matrix = cost_matrix / np.median(cost_matrix)
        return cost_matrix

    def compute_transport_maps(self, callback):
        day_pairs = self.day_pairs
        cell_growth_rates = self.cell_growth_rates
        args = self.args
        covariate_df = self.covariate_df
        covariate_pairs = self.covariate_pairs
        day_to_indices = self.day_to_indices
        ds = self.ds
        for day_index in range(day_pairs.shape[0]):
            t0 = day_pairs.iloc[day_index, 0]
            t1 = day_pairs.iloc[day_index, 1]
            t0_indices = day_to_indices.get(t0)
            if t0_indices is None:
                print('No data for time ' + str(t0))
                continue
            t1_indices = day_to_indices.get(t1)
            if t1_indices is None:
                print('No data for time ' + str(t1))
                continue
            p0_full = wot.Dataset(ds.x[t0_indices], ds.row_meta.iloc[t0_indices], ds.col_meta)
            p1_full = wot.Dataset(ds.x[t1_indices], ds.row_meta.iloc[t1_indices], ds.col_meta)
            p0_5_full = None

            if self.t_interpolate is not None:
                t0_5 = t0 + (t1 - t0) * self.t_interpolate
                t0_5_indices = day_to_indices.get(t0_5)
                if t0_5_indices is not None:
                    p0_5_full = wot.Dataset(ds.x[t0_5_indices], ds.row_meta.iloc[t0_5_indices], ds.col_meta)
                else:
                    print('Unable to find time ' + str(t0_5) + ' - skipping.')
                    continue

            if args.local_pca is not None:
                import scipy.sparse
                matrices = list()
                matrices.append(p0_full.x if not scipy.sparse.isspmatrix(p0_full.x) else p0_full.x.toarray())
                matrices.append(p1_full.x if not scipy.sparse.isspmatrix(p1_full.x) else p1_full.x.toarray())
        #        if p0_5_full is not None:
        #            matrices.append(p0_5_full.x if not scipy.sparse.isspmatrix(p0_5_full.x) else p0_5_full.x.toarray())

                x = np.vstack(matrices)

                x = x - x.mean(axis=0)
                pca = sklearn.decomposition.PCA(n_components=args.local_pca)
                pca.fit(x.T)
                x = pca.components_.T
                p0_full = wot.Dataset(x[0:len(t0_indices)],
                                      p0_full.row_meta,
                                      pd.DataFrame(index=pd.RangeIndex(start=0, stop=args.local_pca, step=1)))

                p1_full = wot.Dataset(x[len(t0_indices):len(t0_indices) + len(t1_indices)],
                                      p1_full.row_meta,
                                      pd.DataFrame(index=pd.RangeIndex(start=0, stop=args.local_pca, step=1)))
                if p0_5_full is not None: # change here
                    U = np.vstack(matrices).T.dot(pca.components_.T).dot(np.diag(1/pca.singular_values_))
                    print(U.shape)
                    print(p0_5_full.x.shape)
                    p0_5_full = wot.Dataset(np.diag(1/pca.singular_values_).dot(U.T.dot(p0_5_full.x.T)).T, p0_5_full.row_meta,
                                    pd.DataFrame(index=pd.RangeIndex(start=0, stop=args.local_pca, step=1)))
#                    p0_5_full = wot.Dataset(x[len(t0_indices) + len(t1_indices):],
#                                            p0_5_full.row_meta,
#                                            pd.DataFrame(index=pd.RangeIndex(start=0, stop=args.local_pca, step=1)))
                self.eigenvals = np.diag(pca.singular_values_)

            delta_t = t1 - t0
            for covariate_pair in covariate_pairs:
                cv0 = covariate_pair[0]
                p0_expr = None if cv0 is None else np.where(p0_full.row_meta[covariate_df.columns[0]] == cv0)[0]
                cv1 = covariate_pair[1]
                p1_expr = None if cv1 is None else np.where(p1_full.row_meta[covariate_df.columns[0]] == cv1)[0]

                if p0_expr is None:
                    p0 = wot.Dataset(p0_full.x, p0_full.row_meta, p0_full.col_meta)
                else:
                    p0 = wot.Dataset(p0_full.x[p0_expr], p0_full.row_meta.iloc[p0_expr], p0_full.col_meta)
                if p1_expr is None:
                    p1 = wot.Dataset(p1_full.x, p1_full.row_meta, p1_full.col_meta)
                else:
                    p1 = wot.Dataset(p1_full.x[p1_expr], p1_full.row_meta.iloc[p1_expr], p1_full.col_meta)

                if args.verbose:
                    print('Computing cost matrix...', end='')

                cost_matrix = self.compute_cost_matrix(p0.x, p1.x)
                if args.verbose:
                    print('done')

                if args.verbose:
                    if covariate_df is not None:
                        print(
                            'Computing transport map from ' + str(
                                t0) + ' ' + (str(cv0) if cv0 is not None else 'full') + ' to ' + str(
                                t1) + ' ' + (str(cv1) if cv1 is not None else 'full') + '...', end='')
                    else:
                        print(
                            'Computing transport map from ' + str(
                                t0) + ' to ' + str(
                                t1) + '...', end='')
                growth_rate = p0.row_meta[cell_growth_rates.columns[0]].values
                result = wot.ot.optimal_transport(cost_matrix=cost_matrix,
                                                  growth_rate=growth_rate,
                                                  delta_days=delta_t,
                                                  max_transport_fraction=args.max_transport_fraction,
                                                  min_transport_fraction=args.min_transport_fraction,
                                                  min_growth_fit=args.min_growth_fit,
                                                  l0_max=args.l0_max, lambda1=args.lambda1,
                                                  lambda2=args.lambda2,
                                                  epsilon=args.epsilon,
                                                  scaling_iter=args.scaling_iter,
                                                  epsilon_adjust=args.epsilon_adjust,
                                                  lambda_adjust=args.lambda_adjust,
                                                  numItermax=args.numItermax,
                                                  epsilon0=args.epsilon0,
                                                  numInnerItermax=args.numInnerItermax, tau=args.tau,
                                                  stopThr=args.stopThr,
                                                  solver=args.solver, growth_iters=args.growth_iters)

                if args.verbose:
                    print('done')

                callback({'t0': t0, 't1': t1, 'result': result, 'df0': p0.row_meta, 'df1': p1.row_meta,
                          'P0': p0, 'P1': p1, 'P0.5': p0_5_full, 'g': growth_rate ** delta_t,
                          'P0_suffix': '_cv-' + str(cv0) if cv0 is not None else '',
                          'P1_suffix': '_cv-' + str(cv1) if cv1 is not None else ''})
