#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `wot` package."""

import unittest
import numpy as np
import scipy.stats
import sklearn.metrics
import pandas
import os.path
import wot


class TestWOT(unittest.TestCase):
    """Tests for `wot` package."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_same_distrubution(self):
        # test same distributions, diagonals along transport map should be equal
        m1 = np.random.rand(2, 3)
        m2 = m1
        result = wot.transport_stable(np.ones(m1.shape[0]),
                                      np.ones(m2.shape[0]),
                                      sklearn.metrics.pairwise.pairwise_distances(
                                          m1, Y=m2, metric='sqeuclidean'), 1, 1,
                                      0.1, 250, np.ones(m1.shape[0]))
        self.assertEqual(result[0, 0], result[1, 1])
        self.assertEqual(result[0, 1], result[1, 0])

    def test_growth_rate(self):
        # as growth rate goes up, sum of mass goes up
        m1 = np.random.rand(2, 3)
        m2 = np.random.rand(4, 3)
        g = [1, 2, 3]
        cost_matrix = sklearn.metrics.pairwise.pairwise_distances(m1, Y=m2,
                                                                  metric='sqeuclidean')
        last = -1
        for i in xrange(len(g)):
            result = wot.transport_stable(np.ones(m1.shape[0]),
                                          np.ones(m2.shape[0]), cost_matrix, 1,
                                          1, 0.1, 250,
                                          np.ones(m1.shape[0]) * g[i])
            sum = np.sum(result)
            self.assertTrue(sum > last)
            last = sum

    def test_epsilon(self):
        # as epsilon goes up, entropy goes up
        m1 = np.random.rand(2, 3)
        m2 = np.random.rand(4, 3)
        e = [0.01, 0.1, 1]
        cost_matrix = sklearn.metrics.pairwise.pairwise_distances(m1, Y=m2,
                                                                  metric='sqeuclidean')
        last = -1
        for i in xrange(len(e)):
            result = wot.transport_stable(np.ones(m1.shape[0]),
                                          np.ones(m2.shape[0]), cost_matrix, 1,
                                          1, e[i], 250,
                                          np.ones(m1.shape[0]))
            sum = np.sum([scipy.stats.entropy(r) for r in result])
            self.assertTrue(sum > last)
            last = sum

    def test_known_output(self):
        diffusion_matrix = pandas.read_csv("data/diffusion_map.csv.gz",
                                           index_col=0)  # cells on rows, diffusion components on columns
        growth_scores = pandas.read_csv("data/growth_scores.csv.gz",
                                        index_col=0)
        days = pandas.read_csv("data/days.csv.gz", index_col=0)
        precomputed_transport_map = None
        if os.path.isfile("data/day0_to_day2.txt.gz"):
            precomputed_transport_map = pandas.read_csv(
                "data/day0_to_day2.txt.gz",
                delimiter='\t')
        diffusion_matrix = diffusion_matrix.join(growth_scores).join(days)
        group_by_day = diffusion_matrix.groupby('day')
        timepoints = group_by_day.groups.keys()
        timepoints.sort()
        self.assertTrue(timepoints == [0, 2])
        max_transport_fraction = 0.4
        min_transport_fraction = 0.05
        min_growth_fit = 0.9
        l0_max = 100
        lambda1 = 1
        lambda2 = 1
        epsilon = 0.1
        growth_ratio = 2.5
        scaling_iter = 250
        expected_lambda = 1.5
        expected_epsilon = 0.01015255979947704383092865754179001669399440288543701171875

        m1 = group_by_day.get_group(0)
        m2 = group_by_day.get_group(2)
        delta_t = timepoints[1] - timepoints[0]
        cost_matrix = sklearn.metrics.pairwise.pairwise_distances(
            m1.drop(["day", "growth_score"], axis=1),
            Y=m2.drop(["day", "growth_score"], axis=1),
            metric='sqeuclidean')
        cost_matrix = cost_matrix / np.median(cost_matrix)
        growth_rate = m1.growth_score.values
        result = wot.optimal_transport(cost_matrix, growth_rate,
                                       delta_days=delta_t,
                                       max_transport_fraction=max_transport_fraction,
                                       min_transport_fraction=min_transport_fraction,
                                       min_growth_fit=min_growth_fit,
                                       l0_max=l0_max, lambda1=lambda1,
                                       lambda2=lambda2, epsilon=epsilon,
                                       growth_ratio=growth_ratio,
                                       scaling_iter=scaling_iter)
        self.assertTrue(result["epsilon"] == expected_epsilon)
        self.assertTrue(result["lambda1"] == expected_lambda)
        transport = pandas.DataFrame(result["transport"], index=m1.index,
                                     columns=m2.index)

        if precomputed_transport_map is not None:
            pandas.testing.assert_index_equal(left=transport.index,
                                              right=precomputed_transport_map.index,
                                              check_names=False)
            pandas.testing.assert_index_equal(left=transport.columns,
                                              right=precomputed_transport_map.columns,
                                              check_names=False)
            np.testing.assert_allclose(transport.values,
                                       precomputed_transport_map.values,
                                       atol=0.0002)


if __name__ == '__main__':
    unittest.main()
