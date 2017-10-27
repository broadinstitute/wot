#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import wot
import pandas
import numpy as np
import sklearn.metrics.pairwise

parser = argparse.ArgumentParser(
    description="Compute transport maps between successive time points")

parser.add_argument("--expression_file",
                    help="Gene expression file with cells on rows and "
                         "features on columns",
                    required=True)
parser.add_argument("--growth_file",
                    help="Two column file without header with cell ids and "
                         "growth scores",
                    required=True)
parser.add_argument("--days_file",
                    help="Two column file without header with cell ids and "
                         "days",
                    required=True)
parser.add_argument("--epsilon", type=float,
                    help="Controls the entropy of the transport map. An "
                         "extremely large entropy parameter will give a "
                         "maximally entropic transport map, and an "
                         "extremely "
                         "small entropy parameter will give a nearly "
                         "deterministic transport map (but could also "
                         "lead to "
                         "numerical instability in the algorithm")
parser.add_argument("--prefix",
                    help="Prefix for ouput file names",
                    required=True)
parser.add_argument("--max_transport_fraction",
                    default=0.4,
                    help="The maximum fraction of cells at time t that are "
                         "transported to time t + 1",
                    type=float)
parser.add_argument("--min_transport_fraction",
                    default=0.05,
                    help="The minimum fraction of cells at time t that are "
                         "transported to time t + 1",
                    type=float)
parser.add_argument("--lambda1", default=1,
                    help="Regularization parameter that controls the "
                         "fidelity "
                         "of the constraints on p",
                    type=float)
parser.add_argument("--lambda2", default=1,
                    help="Regularization parameter that controls the "
                         "fidelity "
                         "of the constraints on q",
                    type=float)
parser.add_argument("--growth_ratio", type=float, default=2.5,
                    help="Over 1 day, a cell in the more proliferative "
                         "group "
                         "is expected to produce growth_ratio times as "
                         "many "
                         "offspring as a cell in the non-proliferative "
                         "group")
parser.add_argument("--scaling_iter", default=250,
                    help="Number of scaling iterations", type=int)
parser.add_argument("--min_growth_fit", type=float, default=0.9)
parser.add_argument("--l0_max", type=float, default=100)
parser.add_argument("--verbose", action="store_true",
                    help="Print progress information")

args = parser.parse_args()

# cells on rows, features on columns
gene_expression = pandas.read_table(args.expression_file, index_col=0)
growth_scores = pandas.read_table(args.growth_file, index_col=0, header=None,
                                names=["growth_score"])
days = pandas.read_table(args.days_file, index_col=0, header=None,
                       names=["day"])

gene_expression = gene_expression.join(growth_scores).join(days)

growth_score_field_name = growth_scores.columns[0]
day_field_name = days.columns[0]
fields_to_drop_for_distance = [day_field_name, growth_score_field_name]

group_by_day = gene_expression.groupby(day_field_name)
days = list(group_by_day.groups.keys())
days.sort()
if len(days) == 1:
    print("Only one unique day found.")
    exit(1)

if args.verbose:
    print(str(len(days)) + " unique days")
for i in range(len(days) - 1):
    # consecutive days only
    m1 = group_by_day.get_group(days[i])
    m2 = group_by_day.get_group(days[i + 1])
    if args.verbose:
        print(
            "Computing transport map from day " + str(
                days[i]) + " to day " + str(
                days[i + 1]))
    delta_t = days[i + 1] - days[i]
    cost_matrix = sklearn.metrics.pairwise.pairwise_distances(
        m1.drop(fields_to_drop_for_distance, axis=1),
        Y=m2.drop(fields_to_drop_for_distance, axis=1),
        metric="sqeuclidean")
    cost_matrix = cost_matrix / np.median(cost_matrix)
    growth_rate = m1.growth_score.values
    result = wot.optimal_transport(cost_matrix, growth_rate,
                                   delta_days=delta_t,
                                   max_transport_fraction=args.max_transport_fraction,
                                   min_transport_fraction=args.min_transport_fraction,
                                   min_growth_fit=args.min_growth_fit,
                                   l0_max=args.l0_max, lambda1=args.lambda1,
                                   lambda2=args.lambda2,
                                   epsilon=args.epsilon,
                                   growth_ratio=args.growth_ratio,
                                   scaling_iter=args.scaling_iter)
    transport = pandas.DataFrame(result["transport"], index=m1.index,
                                 columns=m2.index)
    if args.verbose:
        print("Saving transport map")
    transport.to_csv(args.prefix + "_" + str(days[i]) + "_" + str(
        days[i + 1]) + ".txt", index_label="id", sep="\t")
