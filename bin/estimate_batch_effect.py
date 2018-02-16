#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import wot.io
import ot as pot
import numpy as np
import sklearn.metrics
import pandas as pd


def point_cloud_distance2(cloud_distances, a=None, b=None):
    if a is None:
        a = np.ones((cloud_distances.shape[0]), dtype=np.float64) / cloud_distances.shape[0]
    else:
        a = a / a.sum()
    if b is None:
        b = np.ones((cloud_distances.shape[1]), dtype=np.float64) / cloud_distances.shape[1]
    else:
        b = b / b.sum()
    return np.sqrt(pot.emd2(a, b, cloud_distances, numItermax=10000000))


def point_cloud_distance(c1, c2, a=None, b=None):
    cloud_distances = sklearn.metrics.pairwise.pairwise_distances(c1, Y=c2, metric='sqeuclidean')

    if a is None:
        a = np.ones((cloud_distances.shape[0]), dtype=np.float64) / cloud_distances.shape[0]
    else:
        a = a / a.sum()
    if b is None:
        b = np.ones((cloud_distances.shape[1]), dtype=np.float64) / cloud_distances.shape[1]
    else:
        b = b / b.sum()
    return np.sqrt(pot.emd2(a, b, cloud_distances, numItermax=10000000))


def split_in_two(n):
    indices = np.random.choice(n, int(n * 0.5))
    indices.sort()
    indices_c = np.zeros(n, dtype=bool)
    indices_c[indices] = True
    indices_c = np.invert(indices_c)
    return indices, np.where(indices_c == 1)[0]


parser = argparse.ArgumentParser(
    description='Estimate batch effect')

parser.add_argument('--matrix',
                    help='Gene expression file with cells on rows and features on columns', required=True)

parser.add_argument('--verbose', action='store_true',
                    help='Print progress information')

parser.add_argument('--prefix', help='Output file name prefix', required=True)
parser.add_argument('--diagonal', help='Diagonal scaling matrix')
parser.add_argument('--power', help='Diagonal scaling power', type=float)
parser.add_argument('--cell_days',
                    help='Two column tab delimited file without header with '
                         'cell ids and days', required=True)
parser.add_argument('--cell_days_filter',
                    help='List of cell days to estimate batch effect for', action='append', type=float)
parser.add_argument('--covariate',
                    help='Two column tab delimited file without header with '
                         'cell ids and covariate value')
args = parser.parse_args()
days_df = pd.read_table(args.cell_days, index_col=0, header=None, names=['day'], engine='python', sep=None,
                        dtype={'day': np.float32})

days_df = days_df[~days_df.index.duplicated(keep='first')]

ds = wot.io.read_dataset(args.matrix)

covariate_df = pd.read_table(args.covariate, index_col=0, header=None, names=['covariate'], engine='python', sep=None)

days_df = days_df.align(ds.row_meta, join='right', axis=0, fill_value=np.nan, copy=False)[0]
covariate_df = covariate_df.align(ds.row_meta, join='right', axis=0, copy=False)[0]

days_df = days_df.dropna()
covariate_df = covariate_df.dropna()

unique_covariates = pd.unique(covariate_df[covariate_df.columns[0]].values)
unique_days = pd.unique(days_df[days_df.columns[0]].values)

if args.cell_days_filter is not None:
    unique_days = np.intersect1d(unique_days, np.array(args.cell_days_filter))
if args.diagonal is not None:
    eigenvals = np.loadtxt(args.diagonal, delimiter='\n')
if eigenvals is not None and args.power is not None:
    eigenvals = np.power(eigenvals, args.power)
if eigenvals is not None:
    ds.x = ds.x.dot(np.diag(eigenvals))

niter = 50
writer = open(args.prefix + '_batch_estimate.txt', 'w')
writer.write('time' + '\t' + 'comparison' + '\t' + 'distance' + '\n')

ncovariates = len(unique_covariates) if len(unique_covariates) > 2 else 1
for day in unique_days:
    day_indices = np.where(days_df[days_df.columns[0]] == day)[0]
    x = ds.x[day_indices]
    covariate_df_day = covariate_df.iloc[day_indices]
    all_distances = sklearn.metrics.pairwise.pairwise_distances(x, metric='sqeuclidean')
    print('Computing ' + str(day) + '...')
    for covariate_index in range(ncovariates):
        c = unique_covariates[covariate_index]
        batch_indices = np.where(covariate_df_day[covariate_df_day.columns[0]] == c)[
            0]
        other_indices = np.where(covariate_df_day[covariate_df_day.columns[0]] != c)[0]
        print('Computing ' + str(c) + '...')
        # dist = point_cloud_distance(x[batch_indices], x[batch_indices2])
        m = all_distances[batch_indices]
        m = m[:, other_indices]
        dist = point_cloud_distance2(m if m.flags['C_CONTIGUOUS'] else m.copy(order='C'))
        writer.write(str(day) + '\t' + str(c) + '\t' + str(dist) + '\n')
    writer.flush()
    for i in range(niter):
        print('Iteration ' + str(i + 1))
        indices1, indices2 = split_in_two(x.shape[0])
        m = all_distances[indices1]
        m = m[:, indices2]

        dist = point_cloud_distance2(m if m.flags['C_CONTIGUOUS'] else m.copy(order='C'))
        # dist = point_cloud_distance(x[indices1], x[indices2])

        writer.write(str(day) + '\t' + 'random' + '\t' + str(dist) + '\n')
        writer.flush()

writer.close()
