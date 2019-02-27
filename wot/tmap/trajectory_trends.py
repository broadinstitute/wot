import anndata
import numpy as np
import scipy

import wot


class TrajectoryTrends:

    @staticmethod
    def __weighted_average(weights, ds, value_transform=None):
        if weights is not None:
            weights = weights / np.sum(weights)
        values = ds.X
        if value_transform is not None:
            values = value_transform(values)
        if scipy.sparse.isspmatrix(values):
            mean = values.T @ weights
            variance = np.average((values - mean).A1 ** 2, weights=weights, axis=0)
        else:
            mean = np.average(values, weights=weights, axis=0)
            variance = np.average((values - mean) ** 2, weights=weights, axis=0)
        return {'mean': mean, 'variance': variance}

    @staticmethod
    def compute_dataset_name_to_trends(trajectory_results, unaligned_datasets, dataset_names, value_transform=None):
        """

        Args:
            trajectory_results (list): Results from trajectory_for_cell_sets_at_time_t
            unaligned_datasets (list): List of datasets that can be unaligned with transport map matrices to compute the mean and variance
            dataset_names (list): List of dataset names
            value_transform (function) A function to transform the values in the dataset that takes a numpy array and returns a numpy array of the same shape.
        Return:
            Dict that maps dataset name to trends. A trend is a dict of mean, variance, times, ncells, cell_set, and features. Means and variance have time on rows and features on columns.
        """
        cell_set_name_to_trajectories = wot.ot.Trajectory.group_trajectories_by_cell_set(trajectory_results)
        dataset_name_to_trends = {}
        for ds_index in range(len(unaligned_datasets)):
            unaligned_ds = unaligned_datasets[ds_index]
            trends = []
            for cell_set_name in cell_set_name_to_trajectories:
                means = []
                variances = []
                times = []
                ncells = []
                trajectories = cell_set_name_to_trajectories[cell_set_name]

                # trajectories are sorted by time
                for trajectory in trajectories:
                    p = trajectory['p']
                    cell_ids = trajectory['cell_ids']
                    times.append(trajectory['t'])
                    ncells.append(len(cell_ids))
                    # align dataset with cell_ids
                    ds_order = unaligned_ds.obs.index.get_indexer_for(cell_ids)
                    ds_order = ds_order[ds_order != -1]
                    aligned_dataset = unaligned_ds[ds_order]
                    mean_and_variance = TrajectoryTrends.__weighted_average(weights=p, ds=aligned_dataset,
                                                                            value_transform=value_transform)
                    means.append(mean_and_variance['mean'])
                    variances.append(mean_and_variance['variance'])

                mean = np.array(means)
                variance = np.array(variances)
                trends.append(
                    {'mean': mean, 'variance': variance, 'times': times, 'ncells': ncells, 'cell_set': cell_set_name,
                     'features': unaligned_ds.var.index.values})
            dataset_name_to_trends[dataset_names[ds_index]] = trends
        return dataset_name_to_trends
