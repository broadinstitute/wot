import unittest

import numpy as np
import pandas as pd
import scipy.stats
import sklearn.metrics

import wot.ot


class TestOT(unittest.TestCase):
    """Tests for `wot` package."""

    def test_same_distribution(self):
        # test same distributions, diagonals along transport map should be equal
        m1 = np.random.rand(2, 3)
        m2 = m1
        result = wot.ot.transport_stable(np.ones(m1.shape[0]),
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
        for i in range(len(g)):
            result = wot.ot.transport_stable(np.ones(m1.shape[0]),
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
        for i in range(len(e)):
            result = wot.ot.transport_stable(np.ones(m1.shape[0]),
                                             np.ones(m2.shape[0]), cost_matrix, 1,
                                             1, e[i], 250,
                                             np.ones(m1.shape[0]))
            sum = np.sum([scipy.stats.entropy(r) for r in result])
            self.assertTrue(sum > last)
            last = sum

    def test_growth_scores(self):
        scores = wot.ot.compute_growth_scores(np.array([-0.399883307]),
                                              np.array([0.006853961]))
        np.testing.assert_allclose(np.array([0.705444456674597]),
                                   scores,
                                   atol=0.000001)

    def test_trajectory(self):
        transport_maps = list()
        ncells = [4, 2, 5, 6, 3]
        for i in range(0, len(ncells) - 1):
            transport_map = pd.read_csv('inputs/transport_maps/t' + str(i) +
                                        '_t' + str(i + 1) + '.csv',
                                        index_col=0)

            self.assertTrue(transport_map.shape[0] == ncells[i])
            self.assertTrue(transport_map.shape[1] == ncells[i + 1])
            transport_maps.append({'transport_map': transport_map,
                                   't_minus_1': i, 't': i + 1})
        trajectory_id = ['c4-t3']
        result = wot.tmap.trajectory(trajectory_id, transport_maps, 3, False)
        ancestors = result['ancestors']

        # not messing up already computed ancestors
        ids = ['c1-t2', 'c2-t2',
               'c3-t2', 'c4-t2',
               'c5-t2']
        pd.testing.assert_frame_equal(
            ancestors[ancestors.index.isin(ids)],
            pd.DataFrame(
                {trajectory_id[0]: [5, 4, 3, 2, 1]},
                index=ids), check_names=False)

        # t1
        ids = ['c1-t1']
        pd.testing.assert_frame_equal(
            ancestors[ancestors.index.isin(ids)],
            pd.DataFrame(
                {trajectory_id[0]: [50]},
                index=ids), check_names=False)

        # t0
        ids = ['c3-t0']
        pd.testing.assert_frame_equal(
            ancestors[ancestors.index.isin(ids)],
            pd.DataFrame(
                {trajectory_id[0]: [1175]},
                index=ids), check_names=False)

        trajectory_id = ['c1-t1']
        result = wot.tmap.trajectory(trajectory_id, transport_maps, 1, False)
        descendants = result['descendants']
        # t3
        ids = ['c1-t3', 'c2-t3', 'c3-t3', 'c4-t3', 'c5-t3', 'c6-t3']
        pd.testing.assert_frame_equal(
            descendants[descendants.index.isin(ids)],
            pd.DataFrame(
                {trajectory_id[0]: [90, 190, 290, 50, 390, 490]},
                index=ids), check_names=False)

        # t4
        ids = ['c3-t4']
        pd.testing.assert_frame_equal(
            descendants[descendants.index.isin(ids)],
            pd.DataFrame(
                {trajectory_id[0]: [25450]},
                index=ids), check_names=False)


if __name__ == '__main__':
    unittest.main()
