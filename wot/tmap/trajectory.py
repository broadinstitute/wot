#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import scipy.stats

import wot.io


class Trajectory:

    @staticmethod
    def group_trajectories_by_cell_set(trajectories):
        cell_set_name_to_trajectories = {}
        for trajectory_result in trajectories:
            cell_set = trajectory_result['cell_set']
            results = cell_set_name_to_trajectories.get(cell_set)
            if results is None:
                results = []
                cell_set_name_to_trajectories[cell_set] = results
            results.append(trajectory_result)

        for cell_set_name in cell_set_name_to_trajectories:
            trajectories = cell_set_name_to_trajectories[cell_set_name]
            trajectories.sort(key=lambda x: x['t'])

        return cell_set_name_to_trajectories

    @staticmethod
    def ancestry_similarity_score(ancestor_dist1, ancestor_dist2):
        return 1.0 - 0.5 * np.sum(np.abs(ancestor_dist1 - ancestor_dist2))

    @staticmethod
    def ancestry_similarity(trajectories):
        traces = []
        cell_set_name_to_trajectories = Trajectory.group_trajectories_by_cell_set(trajectories)
        cell_set_names = list(cell_set_name_to_trajectories.keys())
        for i in range(1, len(cell_set_names)):
            cell_set_name_i = cell_set_names[i]
            trajectories1 = cell_set_name_to_trajectories[cell_set_name_i]

            for j in range(i):
                cell_set_name_j = cell_set_names[j]
                trajectories2 = cell_set_name_to_trajectories[cell_set_name_j]
                x = []
                y = []
                for k in range(len(trajectories1)):
                    sim = Trajectory.ancestry_similarity_score(trajectories1[k]['p'], trajectories2[k]['p'])
                    x.append(trajectories1[k]['t'])
                    y.append(sim)

                traces.append(
                    {'x': x, 'y': y, 'name': cell_set_name_i + ' vs. ' + cell_set_name_j, 'mode': 'lines+markers',
                     'type': 'scatter'})
        return traces

    @staticmethod
    def trajectory_embedding(trajectory_results, coords):
        cell_set_name_to_trajectories = Trajectory.group_trajectories_by_cell_set(trajectory_results)
        cell_set_name_to_traces = {}
        for cell_set in cell_set_name_to_trajectories:
            highs = []
            trajectories = cell_set_name_to_trajectories[cell_set]
            traces = []
            for trajectory in trajectories:
                t = trajectory['t']
                p = trajectory['p']
                cell_ids = trajectory['cell_ids']
                joined = coords.join(pd.DataFrame(index=cell_ids, data={'p': p}), how='right')
                df_sum = joined.groupby(['x', 'y']).sum()
                traces.append({'p': p, 't': t, 'x': df_sum.index.get_level_values(0).tolist(),
                               'y': df_sum.index.get_level_values(1).tolist(),
                               'marker': {'color': df_sum['p'].values.tolist()}})

            cell_set_name_to_traces[cell_set] = traces
        return cell_set_name_to_traces

    @staticmethod
    def trajectory_for_cell_sets(transport_maps, time_to_cell_sets, cache_transport_maps=True, print_progress=False):
        """
        Args:
            transport_maps (list) A list of transport maps
            time_to_cell_sets (dict): Maps time to cell sets. A cell set is a dict containing "name" and "set"
            cache_transport_maps (bool): Whether to cache the transport maps in memory

        Returns:
            List of trajectories. A trajectory is a dict containing cell_set, p, entropy, normalized_entropy, and cell_ids
        """

        trajectory_results = []
        progress = 0
        progress_step = 1 / len(time_to_cell_sets)
        for t in time_to_cell_sets:
            cell_sets = time_to_cell_sets[t]
            trajectory_result = Trajectory.__trajectory_for_cell_sets_at_time_t(
                cell_sets=cell_sets,
                transport_maps=transport_maps,
                time=t,
                cache_transport_maps=cache_transport_maps,
                print_progress=(print_progress, progress, progress_step))
            progress += progress_step
            trajectory_results += trajectory_result

        return trajectory_results

    @staticmethod
    def __trajectory_for_cell_sets_at_time_t(cell_sets, transport_maps, time, cache_transport_maps=True,
                                             print_progress=(False, 0, 0)):
        """

        Args:
            cell_sets (list): A list of dicts containing "name" and "set"
            transport_maps (llist) A list of transport maps
            time (float): The time at which the cell sets are defined
            cache_transport_maps (bool): Whether to cache the transport maps in memory
        """

        transport_map_t1_index = -1
        transport_map_t2_index = -1
        for i in range(len(transport_maps)):
            tmap_dict = transport_maps[i]
            if tmap_dict['t1'] == time:
                transport_map_t1_index = i
            if tmap_dict['t2'] == time:
                transport_map_t2_index = i

        use_t1 = transport_map_t1_index != -1
        tmap_index = transport_map_t1_index if use_t1 else transport_map_t2_index
        tmap_dict_at_t = transport_maps[tmap_index]
        tmap_at_t = tmap_dict_at_t.get('ds')
        if tmap_at_t is None:
            tmap_at_t = wot.io.read_dataset(tmap_dict_at_t['path'])
            if cache_transport_maps:
                tmap_dict_at_t['ds'] = tmap_at_t
        results = []
        pvec_array_at_t = []
        for cell_set_index in range(len(cell_sets)):
            cell_ids_in_set = cell_sets[cell_set_index]['set']
            membership = tmap_at_t.obs.index.isin(cell_ids_in_set) if use_t1 else tmap_at_t.var.index.isin(
                cell_ids_in_set)
            membership = membership.astype(np.float)
            membership /= membership.sum()
            pvec_array_at_t.append(membership)
            entropy = np.exp(scipy.stats.entropy(membership))
            results.append(
                {'cell_set': cell_sets[cell_set_index]['name'],
                 'p': membership,
                 'entropy': entropy,
                 'normalized_entropy': entropy / len(membership), 't': time,
                 'cell_ids': tmap_at_t.obs.index.values if use_t1 else tmap_at_t.var.index.values
                 })

        trajectory_ranges = []
        if transport_map_t2_index != -1:
            trajectory_ranges.append({'is_back': True, 'range': range(transport_map_t2_index, - 1, -1)})

        if transport_map_t1_index != -1:
            trajectory_ranges.append({'is_back': False, 'range': range(transport_map_t1_index, len(transport_maps))})

        progress_count = 0
        total_progress_count = sum([len(t['range']) for t in trajectory_ranges])

        for trajectory_range in trajectory_ranges:
            is_back = trajectory_range['is_back']
            pvec_array = pvec_array_at_t  # initialize
            for transport_index in trajectory_range['range']:
                if print_progress[0]:
                    p = print_progress[1] + print_progress[2] * \
                        progress_count / total_progress_count
                    wot.io.output_progress(p)
                    progress_count += 1
                tmap_dict = transport_maps[transport_index]
                tmap_ds = tmap_dict.get('ds')
                if tmap_ds is None:
                    tmap_ds = wot.io.read_dataset(tmap_dict['path'])
                    if cache_transport_maps:
                        tmap_dict['ds'] = tmap_ds

                new_pvec_array = []
                for cell_set_index in range(len(cell_sets)):
                    p = pvec_array[cell_set_index]
                    if is_back:
                        p = tmap_ds.X.dot(p)
                    else:
                        p = p.dot(tmap_ds.X)
                    p /= p.sum()
                    entropy = np.exp(scipy.stats.entropy(p))
                    results.append(
                        {'cell_set': cell_sets[cell_set_index]['name'],
                         'p': p,
                         'entropy': entropy,
                         'normalized_entropy': entropy / len(p),
                         't': tmap_dict['t1'] if is_back else tmap_dict['t2'],
                         'cell_ids': tmap_ds.obs.index.values if is_back else tmap_ds.var.index.values})
                    # n_choose = int(np.ceil(entropy))
                    # n_choose = min(ncells, n_choose)
                    # n_choose = ncells
                    # sampled_indices = np.random.choice(len(v), n_choose, p=v, replace=True)
                    new_pvec_array.append(p)
                pvec_array = new_pvec_array

        return results

    # @staticmethod
    # def interpolate(x, xi, yi, sigma):
    #     diff = (x - xi)
    #     diff *= -diff
    #     sigma2 = 2 * np.power(sigma, 2)
    #     wi = np.exp(diff / sigma2)
    #     fx = np.sum(yi * wi) / np.sum(wi)
    #     return fx
    #
    # @staticmethod
    # def kernel_smooth(xi, yi, stop, start=0, steps=1000, sigma=0.7):
    #     xlist = np.linspace(start, stop, steps)
    #     fhat = np.zeros(len(xlist))
    #     for i in range(len(xlist)):
    #         fhat[i] = TrajectoryUtil.interpolate(xlist[i], xi, yi, sigma)
    #
    #     return xlist, fhat
