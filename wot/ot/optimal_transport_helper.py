import numpy as np
import pandas as pd
import csv
import sklearn.metrics
import argparse
import io
import wot.ot
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
        parser.add_argument('--lambda1', default=1,
                            help='Regularization parameter that controls the '
                                 'fidelity '
                                 'of the constraints on p',
                            type=float)
        parser.add_argument('--lambda2', default=1,
                            help='Regularization parameter that controls the '
                                 'fidelity '
                                 'of the constraints on q',
                            type=float)
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
        parser.add_argument('--delta_min', type=float, default=0.38, help='Growth function parameter')
        parser.add_argument('--beta_max', type=float, default=1.7, help='Growth function parameter')
        parser.add_argument('--beta_center', type=float, default=0.25, help='Growth function parameter')
        parser.add_argument('--delta_max', type=float, default=1.7, help='Growth function parameter')

        growth_rate_group = parser.add_mutually_exclusive_group(required=True)
        growth_rate_group.add_argument('--gene_set_scores', help='File containing "Cell.cycle" and "Apoptosis" scores')

        growth_rate_group.add_argument('--cell_growth_rates',
                                       help='Two column tab delimited file without '
                                            'header with '
                                            'cell ids and growth rates per day.')
        parser.add_argument('--diagonal', help='Diagonal scaling matrix')
        parser.add_argument('--power', help='Diagonal scaling power', type=float)
        parser.add_argument('--covariate',
                            help='Two column tab delimited file without header with '
                                 'cell ids and covariate value')
        parser.add_argument('--solver',
                            help='Solver to use when computing transport maps. One of unbalanced, floating_epsilon, '
                                 'sinkhorn_epsilon, unregularized',
                            choices=['epsilon', 'sinkhorn_epsilon', 'unbalanced', 'unregularized'],
                            default='sinkhorn_epsilon')

        parser.add_argument('--verbose', action='store_true',
                            help='Print progress information')
        return parser

    def __init__(self, args):
        eigenvals = None
        if args.diagonal is not None:
            eigenvals = np.loadtxt(args.diagonal, delimiter='\n')
        if eigenvals is not None and args.power is not None:
            eigenvals = np.power(eigenvals, args.power)

        # cells on rows, features on columns
        gene_expression = pd.read_table(args.matrix, index_col=0, quoting=csv.QUOTE_NONE, engine='python', sep=None)

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
            gene_set_scores = pd.read_table(args.gene_set_scores, index_col=0,
                                            quoting=csv.QUOTE_NONE, engine='python',
                                            sep=None)
            apoptosis = gene_set_scores['Apoptosis']
            proliferation = gene_set_scores['Cell.cycle']
            g = wot.ot.compute_growth_scores(proliferation.values, apoptosis.values, beta_max=args.beta_max,
                                             beta_center=args.beta_center,
                                             delta_max=args.delta_max, delta_min=args.delta_min,
                                             beta_min=args.beta_min)
            cell_growth_rates = pd.DataFrame(index=gene_set_scores.index,
                                             data={'cell_growth_rate': g})

        else:
            cell_growth_rates = pd.read_table(args.cell_growth_rates, index_col=0,
                                              header=None, names=['cell_growth_rate'],
                                              quoting=csv.QUOTE_NONE, engine='python',
                                              sep=None)
        fields_to_drop_for_distance = [days_data_frame.columns[0], cell_growth_rates.columns[0]]
        gene_expression = gene_expression.join(cell_growth_rates).join(days_data_frame)

        self.covariate_pairs = [[None, None]]
        self.covariate_df = None
        group_by_day = gene_expression.groupby(days_data_frame.columns[0])

        if args.ncells is not None:
            tmp = group_by_day.apply(
                lambda x: x.sample(n=args.ncells, axis=0) if x.shape[0] > args.ncells else x)
            tmp.index = tmp.index.droplevel()
            group_by_day = tmp.groupby(
                days_data_frame.columns[0])

        if args.verbose:
            print('Computing ' + str(day_pairs.shape[0]) + ' transport map' + 's' if
                  day_pairs.shape[0] > 1 else '')
        self.day_pairs = day_pairs
        self.group_by_day = group_by_day
        self.fields_to_drop_for_distance = fields_to_drop_for_distance
        self.cell_growth_rates = cell_growth_rates
        self.args = args

    def compute_cost_matrix(self, a, b):
        if self.eigenvals is not None:
            a = a.dot(self.eigenvals)
            b = b.dot(self.eigenvals)
        cost_matrix = sklearn.metrics.pairwise.pairwise_distances(a, b, metric='sqeuclidean')
        cost_matrix = cost_matrix / np.median(cost_matrix)
        return cost_matrix

    def compute_transport_maps(self, callback):
        day_pairs = self.day_pairs
        group_by_day = self.group_by_day
        fields_to_drop_for_distance = self.fields_to_drop_for_distance
        cell_growth_rates = self.cell_growth_rates
        args = self.args
        covariate_df = self.covariate_df
        covariate_pairs = self.covariate_pairs
        for day_index in range(day_pairs.shape[0]):
            t0 = day_pairs.iloc[day_index, 0]
            t1 = day_pairs.iloc[day_index, 1]
            if group_by_day.groups.get(t0) is None or group_by_day.groups.get(
                    t1) is None:
                print('skipping transport map from ' + str(t0) + ' to ' + str(t1))
                continue
            p0_full = group_by_day.get_group(t0)
            p1_full = group_by_day.get_group(t1)

            delta_t = t1 - t0
            for covariate_pair in covariate_pairs:
                cv0 = covariate_pair[0]
                p0 = p0_full if cv0 is None else p0_full[p0_full[covariate_df.columns[0]] == cv0]
                cv1 = covariate_pair[1]
                p1 = p1_full if cv1 is None else p1_full[p1_full[covariate_df.columns[0]] == cv1]
                if args.verbose:
                    print(
                        'Computing transport map from ' + str(
                            t0) + ' to ' + str(
                            t1) + '...', end='')
                cost_matrix = self.compute_cost_matrix(
                    p0.drop(fields_to_drop_for_distance, axis=1).values,
                    p1.drop(fields_to_drop_for_distance, axis=1).values)

                result = wot.ot.optimal_transport(cost_matrix=cost_matrix,
                                                  growth_rate=p0[
                                                      cell_growth_rates.columns[0]].values,
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
                                                  solver=args.solver)

                if args.verbose:
                    print('done')
                name = (str(cv0) if cv0 is not None else 'full') + '_' + (str(cv1) if cv1 is not None else 'full')
                callback({'t0': t0, 't1': t1, 'p0': p0.drop(fields_to_drop_for_distance, axis=1).values,
                          'p1': p1.drop(fields_to_drop_for_distance, axis=1).values,
                          'result': result, 'name': name, 'df0': p0, 'df1': p1})
