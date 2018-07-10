import argparse
import io
import os

import numpy as np
import pandas as pd
import scipy
import sklearn.decomposition
import sklearn.metrics
import wot.io
import wot.ot


class OptimalTransportHelper:

    @staticmethod
    def create_base_parser(description):
        parser = argparse.ArgumentParser(
            description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

        parser.add_argument('--matrix', help=wot.commands.MATRIX_HELP, required=True)
        parser.add_argument('--cell_days', help=wot.commands.CELL_DAYS_HELP, required=True)
        parser.add_argument('--day_pairs',
                            help='Two column file without header with '
                                 'pairs of days to compute transport maps for')
        parser.add_argument('--local_pca', help='Convert day pairs matrix to PCA coordinates. Set to 0 to disable',
                            type=int, default=30)
        parser.add_argument('--growth_iters', type=int, default=3,
                            help='Number of growth iterations for learning the growth rate.')
        parser.add_argument('--cell_growth_rates',
                            help='Two column file with "id" and "cell_growth_rate" headers corresponding to cell id and growth rate per day.')
        parser.add_argument('--gene_filter',
                            help='File with one gene id per line to use for computing cost matrices (e.g. variable genes)')
        parser.add_argument('--cell_filter',
                            help='File with one cell id per line to include or or a python regular expression of cell ids to include')
        parser.add_argument('--epsilon', type=float, default=0.05,
                            help='Controls the entropy of the transport map. An '
                                 'extremely large entropy parameter will give a '
                                 'maximally entropic transport map, and an '
                                 'extremely '
                                 'small entropy parameter will give a nearly '
                                 'deterministic transport map (but could also '
                                 'lead to '
                                 'numerical instability in the algorithm')

        parser.add_argument('--lambda1',
                            help='Regularization parameter that controls the '
                                 'fidelity of the constraints on p', type=float, default=1)
        parser.add_argument('--lambda2', default=50,
                            help='Regularization parameter that controls the '
                                 'fidelity of the constraints on q', type=float)
        parser.add_argument('--scaling_iter', default=3000,
                            help='Number of scaling iterations', type=int)
        parser.add_argument('--out', help='Prefix for ouput file names')
        parser.add_argument('--ncells', help='Number of cells to downsample from each timepoint', type=int)
        parser.add_argument('--ncounts', help='Number of counts to downsample from each cell', type=int)

        # parser.add_argument('--max_transport_fraction',
        #                     default=0.4,
        #                     help='The maximum fraction of cells at time t that are '
        #                          'transported to time t + 1 for floating_epsilon solver',
        #                     type=float)
        # parser.add_argument('--min_transport_fraction',
        #                     default=0.05,
        #                     help='The minimum fraction of cells at time t that are '
        #                          'transported to time t + 1 for floating_epsilon solver',
        #                     type=float)

        # parser.add_argument('--min_growth_fit', type=float, default=0.9)
        # parser.add_argument('--l0_max', type=float, default=100)

        # parser.add_argument('--epsilon_adjust', help='Scaling factor to adjust epsilon for floating_epsilon solver',
        #                     type=float, default=1.1)
        # parser.add_argument('--lambda_adjust', help='Scaling factor to adjust lambda for floating_epsilon solver',
        #                     type=float, default=1.5)

        parser.add_argument('--numItermax', type=int, default=100)
        parser.add_argument('--epsilon0', type=float, default=1, help='Warm starting value for epsilon')
        parser.add_argument('--numInnerItermax', type=int, default=50)
        parser.add_argument('--tau', type=float, default=10000)
        # parser.add_argument('--stopThr', type=float, default=1e-10, help='For sinkhorn_epsilon solver')

        # parser.add_argument('--beta_min', type=float, default=0.3, help='Growth function parameter')
        # parser.add_argument('--delta_min', type=float, default=0.3, help='Growth function parameter')
        # parser.add_argument('--beta_max', type=float, default=1.7, help='Growth function parameter')
        # parser.add_argument('--delta_max', type=float, default=1.7, help='Growth function parameter')
        # parser.add_argument('--beta_center', type=float, default=0.25, help='Growth function parameter')

        # parser.add_argument('--diagonal', help='Diagonal scaling matrix')
        # parser.add_argument('--power', help='Diagonal scaling power', type=float)

        # parser.add_argument('--solver',
        #                     help='Solver to use when computing transport maps. One of unbalanced, floating_epsilon, '
        #                          'sinkhorn_epsilon, unregularized',
        #                     choices=['epsilon', 'sinkhorn_epsilon', 'unbalanced', 'unregularized'],
        #                     default='unbalanced')

        parser.add_argument('--verbose', action='store_true',
                            help='Print progress information')
        return parser

    def __init__(self, args, covariate_df=None, covariate_pairs=None):

        eigenvals = None
        # if args.diagonal is not None:
        #     eigenvals = \
        #         pd.read_table(args.diagonal, header=None, names=['eigenvals'], index_col=False, dtype=np.float64)[
        #             'eigenvals'].values
        # if eigenvals is not None and args.power is not None:
        #     eigenvals = np.power(eigenvals, args.power)

        if args.out is None:
            args.out = wot.io.get_filename_and_extension(os.path.basename(args.matrix))[0] + '_ot'

        # cells on rows, features on columns
        ds = wot.io.read_dataset(args.matrix)
        ds = wot.io.filter_ds_from_command_line(ds, args)

        if args.ncounts is not None:
            for i in range(ds.x.shape[0]):
                p = ds.x[i]
                if scipy.sparse.isspmatrix(p):
                    p = p.toarray()
                p = p.astype('float64')
                counts_p = p.sum()
                if counts_p > args.ncounts:
                    p /= counts_p
                    ds.x[i] = np.random.multinomial(args.ncounts, p, size=1)[0]

        days_data_frame = wot.io.read_days_data_frame(args.cell_days)
        day_pairs = None
        if args.day_pairs is not None:
            if not os.path.isfile(args.day_pairs):
                day_pairs = pd.read_table(io.StringIO(args.day_pairs), header=None, names=['t0', 't1'],
                                          index_col=False, lineterminator=';', sep=',',
                                          dtype={'t0': np.float64, 't1': np.float64})
            else:
                day_pairs = pd.read_table(args.day_pairs, header=None, names=['t0', 't1'],
                                          index_col=False, engine='python', sep=None,
                                          dtype={'t0': np.float64, 't1': np.float64})

        # self.eigenvals = np.diag(eigenvals) if eigenvals is not None else None
        self.eigenvals = None

        # if args.gene_set_scores is not None:
        #     ext = wot.io.get_filename_and_extension(args.gene_set_scores)[1]
        #     if ext == 'loom' or ext == 'gct':
        #         gene_set_scores_ds = wot.io.read_dataset(args.gene_set_scores)
        #         apoptosis = gene_set_scores_ds[:, np.where(gene_set_scores_ds.col_meta.columns == 'Apoptosis')[0]]
        #         proliferation = gene_set_scores_ds[:, np.where(gene_set_scores_ds.col_meta.columns == 'Cell.cycle')[0]]
        #         gene_set_scores_ids = gene_set_scores_ds.row_meta.index
        #     else:
        #         gene_set_scores = pd.read_table(args.gene_set_scores, index_col=0, engine='python', sep=None)
        #         apoptosis = gene_set_scores['Apoptosis'].values
        #         proliferation = gene_set_scores['Cell.cycle'].values
        #         gene_set_scores_ids = gene_set_scores.index
        #     g = wot.ot.compute_growth_scores(proliferation, apoptosis, beta_max=args.beta_max,
        #                                      beta_center=args.beta_center,
        #                                      delta_max=args.delta_max, delta_min=args.delta_min,
        #                                      beta_min=args.beta_min)
        #     cell_growth_rates = pd.DataFrame(index=gene_set_scores_ids, data={'cell_growth_rate': g})

        if args.cell_growth_rates is not None:
            cell_growth_rates = pd.read_table(args.cell_growth_rates, index_col='id', engine='python', sep=None)
        else:
            cell_growth_rates = pd.DataFrame(index=ds.row_meta.index.values, data={'cell_growth_rate': 1})
            if args.verbose:
                print('Using growth rate of 1')
        ds.row_meta = ds.row_meta.join(cell_growth_rates).join(days_data_frame)
        if day_pairs is None:
            unique_days = list(set(days_data_frame['day'].values))
            unique_days.sort()
            _unique_days = []
            for i in range(len(unique_days)):
                if unique_days[i] >= 0:
                    indices = np.where(ds.row_meta['day'] == unique_days[i])[0]
                    if len(indices) > 0:
                        _unique_days.append(unique_days[i])
            unique_days = _unique_days
            pairs = []
            for i in range(len(unique_days) - 1):
                d1 = unique_days[i]
                d2 = unique_days[i + 1]
                pairs.append([d1, d2])
            day_pairs = pd.DataFrame(data=pairs, columns=['t0', 't1'])
        if day_pairs.shape[0] is 0:
            print('No day pairs found')
            exit(1)
        if covariate_df is not None:
            self.covariate_df = covariate_df
            ds.row_meta = ds.row_meta.join(covariate_df)
            self.covariate_pairs = covariate_pairs + [[None, None]]
        else:
            self.covariate_df = None
            self.covariate_pairs = [[None, None]]
        self.day_pairs = day_pairs
        self.cell_growth_rates = cell_growth_rates
        self.args = args
        self.t_interpolate = vars(args).get('t_interpolate')
        day_to_indices = {}
        self.ds = ds
        if args.ncells is not None:
            unique_cvs = set()
            if covariate_pairs is None:
                unique_cvs.add(None)
            else:
                for cv_pair in covariate_pairs:
                    unique_cvs.add(cv_pair[0])
                    unique_cvs.add(cv_pair[1])
            unique_days = set()
            for day_index in range(day_pairs.shape[0]):
                t0 = day_pairs.iloc[day_index, 0]
                t1 = day_pairs.iloc[day_index, 1]
                unique_days.add(t0)
                unique_days.add(t1)
                if self.t_interpolate is not None:
                    t0_5 = t0 + (t1 - t0) * self.t_interpolate
                    unique_days.add(t0_5)
            for day in unique_days:
                index_list = []
                day_query = ds.row_meta['day'] == day
                for cv in unique_cvs:
                    if cv is None:
                        indices = np.where(day_query)[0]
                    else:
                        indices = np.where(day_query & (ds.row_meta[covariate_df.columns[0]] == cv))[0]
                    if len(indices) > args.ncells:
                        np.random.shuffle(indices)
                        indices = indices[0:args.ncells]
                    index_list.append(indices)
                day_to_indices[day] = np.concatenate(index_list)

        else:
            days = ds.row_meta['day'].values
            for i in range(len(days)):
                val = days[i]
                if val is not None:
                    indices = day_to_indices.get(val)
                    if indices is None:
                        indices = []
                        day_to_indices[val] = indices
                    indices.append(i)

        if args.verbose:
            print('Computing ' + str(day_pairs.shape[0]) + ' transport map' + ('s' if
                                                                               day_pairs.shape[0] > 1 else ''))
        self.day_to_indices = day_to_indices

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

    def compute_cost_matrix(self, a, b):
        return OptimalTransportHelper.compute_default_cost_matrix(a, b, self.eigenvals)

    def compute_transport_maps(self, callback):
        day_pairs = self.day_pairs
        args = self.args
        covariate_df = self.covariate_df
        covariate_pairs = self.covariate_pairs
        day_to_indices = self.day_to_indices
        callback_call_count = 0
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

            if args.local_pca is not None and args.local_pca > 0:
                import scipy.sparse
                matrices = list()
                matrices.append(p0_full.x if not scipy.sparse.isspmatrix(p0_full.x) else p0_full.x.toarray())
                matrices.append(p1_full.x if not scipy.sparse.isspmatrix(p1_full.x) else p1_full.x.toarray())

                x = np.vstack(matrices)
                mean_shift = x.mean(axis=0)
                x = x - mean_shift
                pca = sklearn.decomposition.PCA(n_components=args.local_pca)
                pca.fit(x.T)
                x = pca.components_.T
                p0_full = wot.Dataset(x[0:len(t0_indices)],
                                      p0_full.row_meta,
                                      pd.DataFrame(index=pd.RangeIndex(start=0, stop=args.local_pca, step=1)))

                p1_full = wot.Dataset(x[len(t0_indices):len(t0_indices) + len(t1_indices)],
                                      p1_full.row_meta,
                                      pd.DataFrame(index=pd.RangeIndex(start=0, stop=args.local_pca, step=1)))
                if p0_5_full is not None:  # compute PCA only on local coordinates
                    U = np.vstack(matrices).T.dot(pca.components_.T).dot(np.diag(1 / pca.singular_values_))
                    y = p0_5_full.x - mean_shift
                    p0_5_full = wot.Dataset(np.diag(1 / pca.singular_values_).dot(U.T.dot(y.T)).T, p0_5_full.row_meta,
                                            pd.DataFrame(index=pd.RangeIndex(start=0, stop=args.local_pca, step=1)))
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
                    print('Computing cost matrix in ' + str(p0.x.shape[1]) + ' dimensions...', end='')

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
                growth_rate = p0.row_meta['cell_growth_rate'].values
                solver = 'unbalanced'
                result = wot.ot.optimal_transport(cost_matrix=cost_matrix,
                                                  growth_rate=growth_rate,
                                                  delta_days=delta_t,
                                                  lambda1=args.lambda1,
                                                  lambda2=args.lambda2,
                                                  epsilon=args.epsilon,
                                                  scaling_iter=args.scaling_iter,
                                                  numItermax=args.numItermax,
                                                  epsilon0=args.epsilon0,
                                                  numInnerItermax=args.numInnerItermax,
                                                  tau=args.tau,
                                                  solver=solver, growth_iters=max(1, args.growth_iters))

                if args.verbose:
                    print('done')

                callback({'t0': t0, 't1': t1, 'result': result, 'df0': p0.row_meta, 'df1': p1.row_meta,
                          'P0': p0, 'P1': p1, 'P0.5': p0_5_full, 'g': growth_rate ** delta_t,
                          'P0_suffix': '_cv-' + str(cv0) if cv0 is not None else '',
                          'P1_suffix': '_cv-' + str(cv1) if cv1 is not None else '',
                          'cv0': cv0, 'cv1': cv1,
                          'call_count': callback_call_count,
                          'total_call_count': len(day_pairs) * len(covariate_pairs) })
                callback_call_count += 1
