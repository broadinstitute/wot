import unittest
import numpy as np
import pandas as pd
import scipy.stats
import sklearn.metrics

import subprocess
import os
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

    def test_transport_maps_by_time(self):
        clusters = pd.DataFrame([1, 1, 2, 3, 1, 1, 2, 1],
                                index=['a', 'b', 'c', 'd', 'e', 'f', 'g',
                                       'h'])
        map1 = pd.DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]],
                            index=[1, 2, 3], columns=[1, 2, 3])
        map2 = pd.DataFrame([[10, 11, 12], [13, 14, 15], [16, 17, 18]],
                            index=[1, 2, 3], columns=[1, 2, 3])
        # weighted average across time
        cluster_weights_by_time = [[0.4, 0.5, 0], [0.6, 0.5, 1]]
        result = wot.ot.transport_maps_by_time([map1, map2],
                                               cluster_weights_by_time)
        pd.testing.assert_frame_equal(
            result,
            pd.DataFrame(
                [[6.4, 6.5, 12.0],
                 [9.4, 9.5, 15.0],
                 [12.4, 12.5, 18.0],
                 ],
                index=[1, 2, 3], columns=[1, 2, 3]))

    def test_get_weights_intersection(self):
        clusters = pd.DataFrame([1, 1, 2, 3, 1, 1, 2, 1],
                                index=['a', 'b', 'c', 'd', 'e', 'f', 'g',
                                       'h'], columns=['cluster'])
        map1 = pd.DataFrame(index=['x', 'y', 'z'], columns=['a', 'b', 'c'])
        # note f, g, h are not present in map2
        map2 = pd.DataFrame(index=['a', 'b', 'c'], columns=['d', 'e'])
        # weighted average across time

        grouped_by_cluster = clusters.groupby(clusters.columns[0], axis=0)
        cluster_ids = list(grouped_by_cluster.groups.keys())
        all_cell_ids = set()
        column_cell_ids_by_time = []
        for transport_map in [map1, map2]:
            all_cell_ids.update(transport_map.columns)
            column_cell_ids_by_time.append(transport_map.columns)
            all_cell_ids.update(transport_map.index)
        result = wot.ot.get_weights(all_cell_ids, column_cell_ids_by_time,
                                    grouped_by_cluster, cluster_ids)

        self.assertTrue(
            np.array_equal(result['cluster_weights_by_time'],
                           [[2 / 3, 1 / 1, 0 / 1], [1 / 3, 0 / 1, 1 / 1]]))
        self.assertTrue(
            np.array_equal(result['cluster_size'],
                           [3, 1, 1]))

    def test_get_weights(self):
        clusters = pd.DataFrame([1, 1, 2, 3, 1, 1, 2, 1],
                                index=['a', 'b', 'c', 'd', 'e', 'f', 'g',
                                       'h'], columns=['cluster'])
        map1 = pd.DataFrame(index=['x', 'y', 'z'], columns=['a', 'b', 'c'])
        map2 = pd.DataFrame(index=['a', 'b', 'c'], columns=['d', 'e',
                                                            'f', 'g', 'h'])
        # weighted average across time
        grouped_by_cluster = clusters.groupby(clusters.columns[0], axis=0)
        cluster_ids = list(grouped_by_cluster.groups.keys())
        all_cell_ids = set()
        column_cell_ids_by_time = []
        for transport_map in [map1, map2]:
            all_cell_ids.update(transport_map.columns)
            column_cell_ids_by_time.append(transport_map.columns)
            all_cell_ids.update(transport_map.index)
        result = wot.ot.get_weights(all_cell_ids, column_cell_ids_by_time,
                                    grouped_by_cluster, cluster_ids)

        self.assertTrue(
            np.array_equal(result['cluster_weights_by_time'],
                           [[2 / 5, 1 / 2, 0 / 1], [3 / 5, 1 / 2, 1 / 1]]))
        self.assertTrue(
            np.array_equal(result['cluster_size'],
                           [5, 2, 1]))

    def test_transport_map_by_cluster(self):
        row_ids = ['a', 'b', 'c']
        column_ids = ['d', 'e', 'f', 'g', 'h'];
        clusters = pd.DataFrame([3, 1, 2, 3, 1, 1, 2, 3],
                                index=['a', 'b', 'c', 'd', 'e', 'f', 'g',
                                       'h'])
        grouped_by_cluster = clusters.groupby(clusters.columns[0], axis=0)
        cluster_ids = [1, 2, 3]
        transport_map = pd.DataFrame(
            [[1, 2, 3, 4, 5], [6, 7, 8, 9, 10], [11, 12, 13, 14, 15]],
            index=row_ids,
            columns=column_ids)

        # sum mass by cluster
        result = wot.ot.transport_map_by_cluster(transport_map, grouped_by_cluster,
                                                 cluster_ids)

        pd.testing.assert_frame_equal(
            result,
            pd.DataFrame(
                [[15, 9, 16],
                 [25, 14, 26],
                 [5, 4, 6]
                 ],
                index=cluster_ids, columns=cluster_ids))

    def test_growth_scores(self):
        scores = wot.ot.compute_growth_scores(np.array([-0.399883307]),
                                              np.array([0.006853961]))
        np.testing.assert_allclose(np.array([0.705444456674597]),
                                   scores,
                                   atol=0.000001)

    def test_ot_commmand_line(self):
        subprocess.call(args=['python', os.path.abspath('../bin/optimal_transport'),
                              '--matrix',
                              os.path.abspath(
                                  '../finalInput/dmap_2i_normalized.txt'),
                              '--cell_growth_rates', os.path.abspath(
                '../paper/growth_rates.txt'),
                              '--cell_days', os.path.abspath(
                '../paper/days.txt'),
                              '--day_pairs', os.path.abspath(
                '../paper/pairs_2i.txt'),
                              '--prefix', 'mytest',
                              '--verbose',
                              '--compress'],
                        cwd=os.getcwd(),
                        stderr=subprocess.STDOUT)
        timepoints = [0, 2, 4, 6, 8, 9, 10, 11, 12, 16]
        for timepoint in range(0, len(timepoints) - 1):
            my_transport = pd.read_table(
                'mytest_' + str(
                    timepoints[timepoint]) + '_' + str(
                    timepoints[timepoint + 1]) +
                '.txt.gz', index_col=0)

            # precomputed_transport_map = np.load(
            #     '../paper/transport_maps/2i/npy/lineage.day-' + str(
            #         timepoints[timepoint + 1]) + '.npy')
            precomputed_transport_map = precomputed_transport_map = \
                pd.read_table(
                    '../paper/transport_maps/2i/lineage.day-' + str(
                        timepoints[timepoint + 1]) +
                    '.txt', index_col=0)
            pd.testing.assert_index_equal(left=my_transport.index,
                                          right=precomputed_transport_map.index,
                                          check_names=False)
            pd.testing.assert_index_equal(left=my_transport.columns,
                                          right=precomputed_transport_map.columns,
                                          check_names=False)
            total = 0
            count = 0
            for i in range(my_transport.shape[0]):
                for j in range(my_transport.shape[1]):
                    diff = abs(precomputed_transport_map.values[i,
                                                                j] -
                               my_transport.values[
                                   i, j])
                    total += diff
                    if diff > 0.000001:
                        count += 1
            print('lineage_' + str(
                timepoints[timepoint]) + '_' + str(
                timepoints[timepoint + 1]))
            print('total diff: ' + str(total))
            print('pre total: ' + str(precomputed_transport_map.sum().sum()))
            print('my total: ' + str(my_transport.sum().sum()))
            print('count: ' + str(count) + '/' + str(
                my_transport.shape[0] * my_transport.shape[1]))

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
        result = wot.ot.trajectory(trajectory_id, transport_maps, 3, False)
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
        result = wot.ot.trajectory(trajectory_id, transport_maps, 1, False)
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

    def test_ot_known_output(self):
        gene_expression = pd.read_table('../paper/2i_dmap_20.txt',
                                        index_col=0)  # cells on rows,
        # diffusion components on columns
        growth_scores = pd.read_table('../paper/growth_scores.txt',
                                      index_col=0, header=None,
                                      names=['id', 'growth_score'])
        days = pd.read_table(
            '../paper/days.txt', header=None,
            index_col=0, names=['id', 'days'])

        gene_expression = gene_expression.join(growth_scores).join(days)
        growth_score_field_name = growth_scores.columns[0]
        day_field_name = days.columns[0]
        group_by_day = gene_expression.groupby(day_field_name)
        timepoints = list(group_by_day.groups.keys())
        timepoints.sort()
        self.assertTrue(timepoints[0] == 0)
        max_transport_fraction = 0.4
        min_transport_fraction = 0.05
        min_growth_fit = 0.9
        l0_max = 100
        lambda1 = 1
        lambda2 = 1
        epsilon = 0.1
        growth_ratio = 2.5
        scaling_iter = 250
        expected_lambda_t0_t2 = 1.5
        expected_epsilon_t0_t2 = \
            0.01015255979947704383092865754179001669399440288543701171875
        for i in range(0, len(timepoints) - 1):
            m1 = group_by_day.get_group(timepoints[i])
            m2 = group_by_day.get_group(timepoints[i + 1])
            delta_t = timepoints[i + 1] - timepoints[i]
            cost_matrix = sklearn.metrics.pairwise.pairwise_distances(
                m1.drop([day_field_name, growth_score_field_name], axis=1),
                Y=m2.drop([day_field_name, growth_score_field_name], axis=1),
                metric='sqeuclidean')
            cost_matrix = cost_matrix / np.median(cost_matrix)
            growth_rate = m1.growth_score.values
            result = wot.ot.optimal_transport(cost_matrix, growth_rate,
                                              delta_days=delta_t,
                                              max_transport_fraction=max_transport_fraction,
                                              min_transport_fraction=min_transport_fraction,
                                              min_growth_fit=min_growth_fit,
                                              l0_max=l0_max, lambda1=lambda1,
                                              lambda2=lambda2, epsilon=epsilon,
                                              growth_ratio=growth_ratio,
                                              scaling_iter=scaling_iter)
            if i == 0:
                self.assertTrue(result['epsilon'] == expected_epsilon_t0_t2)
                self.assertTrue(result['lambda1'] == expected_lambda_t0_t2)
                self.assertTrue(result['lambda2'] == expected_lambda_t0_t2)
            transport = pd.DataFrame(result['transport'], index=m1.index,
                                     columns=m2.index)

            precomputed_transport_map = pd.read_table(
                '../paper/transport_maps/lineage_' + str(timepoints[i])
                + '_' + str(timepoints[i + 1]) + '.txt', index_col=0)
            pd.testing.assert_index_equal(left=transport.index,
                                          right=precomputed_transport_map.index,
                                          check_names=False)
            pd.testing.assert_index_equal(left=transport.columns,
                                          right=precomputed_transport_map.columns,
                                          check_names=False)

            np.testing.assert_allclose(transport.values,
                                       precomputed_transport_map.values,
                                       atol=0.0004)


if __name__ == '__main__':
    unittest.main()
