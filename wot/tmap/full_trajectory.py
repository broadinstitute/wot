# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd


def full_trajectory(transport_maps, time, ids=None, normalize=False):
    """
        Compute the trajectory of the given ids.

        Args:
            transport_maps (list): A sorted list of dictionaries
            containing 'transport_map' (anndata.AnnData), 't1', and 't2'. The
            ids (list): A list of ids to compute the trajectory for.
            time (float): The time at which the ids were measured.
            normalize (bool) Whether to normalize total mass to one at each
            timepoint.

        Returns:
            dict: A dictionary containing ancestors and descendants.
        """

    t_start_time_index = None
    for i in range(len(transport_maps)):
        if transport_maps[i]['t1'] == time:
            t_start_time_index = i
            break
    if t_start_time_index is None:
        raise RuntimeError(
            'Transport transport_map for time ' + str(time) + ' not found.')

    # transport maps have t1 on rows, t2 on columns
    transport_maps_by_start_time = list(map(lambda x: x['transport_map'], transport_maps))
    # if t=9, t_start_time_index=t9_t10

    # subset rows
    if ids is not None:
        transport_maps_by_start_time[t_start_time_index] = \
            transport_maps_by_start_time[t_start_time_index][
                transport_maps_by_start_time[t_start_time_index].index.isin(
                    ids)]


    # ancestors, go backwards in time
    # t[i−2,i] = t[i−2,i−1]t[i−1,i]
    # e.g. t=9, t7_t8, t8_t9, t9_t10 compute ancestors of t9 at t7

    if t_start_time_index >= 1:
        # subset columns to specified cell ids only
        if ids is not None:
            transport_maps_by_start_time[
                t_start_time_index - 1] = transport_maps_by_start_time[t_start_time_index - 1][ids]

        for i in range(t_start_time_index - 2, -1, -1):
            transport_maps_by_start_time[i] = transport_maps_by_start_time[i].dot(transport_maps_by_start_time[i + 1]) * \
                                              transport_maps_by_start_time[i + 1].shape[0]

    # descendants, go forwards in time
    # e.g. t=9, t_start_time_index=t9_t10, t+1=t10_t11
    for i in range(t_start_time_index + 1,
                   len(transport_maps_by_start_time)):
        transport_maps_by_start_time[i] = transport_maps_by_start_time[i - 1].dot(transport_maps_by_start_time[i]) * \
                                          transport_maps_by_start_time[i].shape[0]

    descendants = []
    descendants_summary = []
    for i in range(t_start_time_index,
                   len(transport_maps_by_start_time)):
        # sum across columns
        m = transport_maps_by_start_time[i].transpose()
        summary = m.sum(axis=1)
        if normalize:
            total = np.sum(summary.values)
            summary = summary / total
        summary = summary.to_frame(name='sum')
        descendants_summary.append(summary)
        descendants.append(m)

    descendants = pd.concat(descendants)
    descendants_summary = pd.concat(descendants_summary)

    ancestors = []
    ancestors_summary = []
    for m in transport_maps_by_start_time[t_start_time_index - 1::-1]:
        # sum across columns
        summary = m.sum(axis=1)
        if normalize:
            total = np.sum(summary.values)
            summary = summary / total
        summary = summary.to_frame(name='sum')
        ancestors_summary.append(summary)
        ancestors.append(m)

    ancestors = pd.concat(ancestors)
    ancestors_summary = pd.concat(ancestors_summary)

    return {'descendants': descendants,
            'descendants_summary': descendants_summary,
            'ancestors': ancestors,
            'ancestors_summary': ancestors_summary}
