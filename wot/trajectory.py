# -*- coding: utf-8 -*-

import pandas


def trajectory(ids, transport_maps, time):
    """
        Compute the trajectory of the given ids.

        Args:
            transport_maps (list): A sorted list of dictionaries
            containing 'transport_map' (DataFrame), 't_minus_1', and 't'. The
            ids (list): A list of ids to compute the trajectory for.
            time (string): The time at which the ids were measured.

        Returns:
            dict: A dictionary containing ancestors and descendants.
        """

    t_start_time_index = None
    for i in range(len(transport_maps)):
        if transport_maps[i]['t_minus_1'] == time:
            t_start_time_index = i
            break
    if t_start_time_index is None:
        raise RuntimeError(
            'Transport transport_map for time ' + str(time) + ' not found.')
    for i in range(len(transport_maps) - 1):
        if transport_maps[i]['t'] != transport_maps[i + 1]['t_minus_1']:
            raise RuntimeError('Matching transport transport_map not found')

    # data frames have t_minus_1 on rows, t on columns
    transport_maps_by_start_time = [None] * len(transport_maps)
    # if t=9, t_start_time_index=t9_t10
    transport_maps_by_start_time[t_start_time_index] = transport_maps[
        t_start_time_index]['transport_map']

    # subset rows
    transport_maps_by_start_time[t_start_time_index] = \
        transport_maps_by_start_time[t_start_time_index][
            transport_maps_by_start_time[t_start_time_index].index.isin(
                ids)]

    # transport maps have t-1 on rows, t on columns

    # ancestors, go backwards in time
    # t[i−2,i] = t[i−2,i−1]t[i−1,i]
    # e.g. t=9, t7_t8, t8_t9, t9_t10 compute ancestors of t9 at t7
    ancestors = None

    if t_start_time_index >= 1:
        # subset columns to specified cell ids only
        transport_maps_by_start_time[
            t_start_time_index - 1] = transport_maps[t_start_time_index - 1][
            'transport_map'][ids]
        for i in range(t_start_time_index - 2, -1, -1):
            transport_maps_by_start_time[i] = transport_maps[i]['transport_map']

            transport_maps_by_start_time[i] = \
                transport_maps_by_start_time[
                    i].dot(
                    transport_maps_by_start_time[i + 1])
        ancestors = pandas.concat(
            transport_maps_by_start_time[t_start_time_index - 1::-1])

    # descendants, go forwards in time
    # e.g. t=9, t_start_time_index=t9_t10, t+1=t10_t11

    for i in range(t_start_time_index + 1,
                   len(transport_maps)):
        transport_maps_by_start_time[i] = transport_maps[i]['transport_map']
        transport_maps_by_start_time[i] = transport_maps_by_start_time[
            i - 1].dot(
            transport_maps_by_start_time[i])
    for i in range(t_start_time_index,
                   len(transport_maps)):
        transport_maps_by_start_time[i] = transport_maps_by_start_time[
            i].transpose()

    descendants = pandas.concat(
        transport_maps_by_start_time[t_start_time_index:])
    return {'ancestors': ancestors, 'descendants': descendants}
